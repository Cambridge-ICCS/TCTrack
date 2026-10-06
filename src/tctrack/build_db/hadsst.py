"""Nearest-neighbour lookup for HadSST NetCDF datasets.

Requires the median monthly anomalies NetCDF, e.g. HadSST.4.2.0.0_median.nc via
https://www.metoffice.gov.uk/hadobs/hadsst4/data/download.html
"""

from pathlib import Path
from typing import Protocol, cast

import netCDF4
import numpy as np


class ObservationTime(Protocol):
    """The date-time fields required to find a matching HadSST time."""

    year: int
    month: int


class HadSST:
    """Read HadSST data and provide nearest-grid-cell temperature anomalies."""

    def __init__(self, filepath: str | Path) -> None:
        """
        Load and validate a HadSST dataset.

        Parameters
        ----------
        filepath
            Path to a HadSST NetCDF dataset containing a `tos` variable.

        Raises
        ------
        ValueError
            If the path cannot be opened as a compatible HadSST dataset.
        """
        path = Path(filepath)
        if not path.is_file():
            msg = f"HadSST dataset does not exist or is not a file: {path}"
            raise ValueError(msg)

        try:
            with netCDF4.Dataset(path, "r") as dataset:
                self._load_dataset(dataset, path)
        except OSError as exc:
            msg = f"Could not open HadSST dataset {path}: {exc}"
            raise ValueError(msg) from exc

    def _load_dataset(self, dataset: netCDF4.Dataset, path: Path) -> None:
        try:
            self._time = self._coordinate(dataset, "time")
            self._latitude = self._coordinate(dataset, "latitude")
            self._longitude = self._coordinate(dataset, "longitude")
            self._ssta = dataset.variables["tos"]
        except KeyError as exc:
            msg = f"HadSST dataset {path} is missing required variable {exc.args[0]!r}"
            raise ValueError(msg) from exc

        self._time_dimension = self._time.dimensions[0]
        self._latitude_dimension = self._latitude.dimensions[0]
        self._longitude_dimension = self._longitude.dimensions[0]
        coordinate_dimensions = {
            self._time_dimension,
            self._latitude_dimension,
            self._longitude_dimension,
        }
        if any(
            coordinate.ndim != 1
            for coordinate in (self._time, self._latitude, self._longitude)
        ):
            msg = f"HadSST dataset {path} coordinates must be one-dimensional"
            raise ValueError(msg)

        number_of_dimensions = 3
        if (
            self._ssta.ndim != number_of_dimensions
            or set(self._ssta.dimensions) != coordinate_dimensions
        ):
            msg = (
                f"HadSST dataset {path} variable 'tos' must have time, latitude, "
                "and longitude dimensions"
            )
            raise ValueError(msg)

        self._time_units = self._attribute(self._time, "units", path)
        self._time_calendar = getattr(self._time, "calendar", "standard")
        self._times = np.asarray(self._time[:], dtype=float)
        self._time_indices_by_month = self._monthly_time_indices(dataset, path)
        self._month_ordinals = np.asarray(
            [year * 12 + month for year, month in self._time_indices_by_month],
            dtype=int,
        )
        self._latitudes = np.asarray(self._latitude[:], dtype=float)
        self._longitudes = np.asarray(self._longitude[:], dtype=float)
        self._ssta_values = self._ssta[:]
        self._dimensions = self._ssta.dimensions
        self._fill_value = getattr(self._ssta, "_FillValue", None)

        if any(
            values.size == 0
            for values in (self._times, self._latitudes, self._longitudes)
        ):
            msg = f"HadSST dataset {path} has an empty coordinate"
            raise ValueError(msg)

    def _monthly_time_indices(
        self, dataset: netCDF4.Dataset, path: Path
    ) -> dict[tuple[int, int], int]:
        """Return HadSST time indices keyed by their calendar year and month."""
        bounds_name = getattr(self._time, "bounds", None)
        if bounds_name not in dataset.variables:
            msg = f"HadSST dataset {path} must define monthly time bounds"
            raise ValueError(msg)

        bounds = dataset.variables[bounds_name]
        number_of_bounds_dimensions = 2
        if (
            bounds.ndim != number_of_bounds_dimensions
            or bounds.dimensions[0] != self._time_dimension
        ):
            msg = f"HadSST dataset {path} has invalid monthly time bounds"
            raise ValueError(msg)

        indices = {}
        for index, time_value in enumerate(self._times):
            date = cast(
                ObservationTime,
                netCDF4.num2date(
                    np.asarray([time_value]),
                    units=self._time_units,
                    calendar=self._time_calendar,
                )[0],
            )
            indices[(date.year, date.month)] = index
        if len(indices) != len(self._times):
            msg = f"HadSST dataset {path} must contain at most one value per month"
            raise ValueError(msg)

        return indices

    @staticmethod
    def _coordinate(dataset: netCDF4.Dataset, name: str) -> netCDF4.Variable:
        """Return a coordinate variable by its CF standard name or variable name."""
        if name in dataset.variables:
            return dataset.variables[name]

        for variable in dataset.variables.values():
            if getattr(variable, "standard_name", None) == name:
                return variable

        raise KeyError(name)

    @staticmethod
    def _attribute(variable: netCDF4.Variable, name: str, path: Path) -> str:
        """Return a required NetCDF variable attribute."""
        if not hasattr(variable, name):
            msg = f"HadSST dataset {path} time coordinate is missing {name!r}"
            raise ValueError(msg)

        return str(getattr(variable, name))

    def lookup(
        self, observation_time: ObservationTime, latitude: float, longitude: float
    ) -> float | None:
        """Return the nearest HadSST anomaly, rounded to 3 decimals or `None`."""
        # HadSST values represent whole calendar months, so this also aligns
        # 360_day timestamps such as 2000-02-30 with February 2000 HadSST.
        time_index = self._nearest_month_index(observation_time)
        latitude_index = self._nearest_index(self._latitudes, latitude)
        longitude_index = self._nearest_longitude_index(
            self._normalise_longitude(longitude),
        )
        indices = {
            self._time_dimension: time_index,
            self._latitude_dimension: latitude_index,
            self._longitude_dimension: longitude_index,
        }
        ssta = self._ssta_values[
            tuple(indices[dimension] for dimension in self._dimensions)
        ]

        if np.ma.is_masked(ssta):
            return None

        ssta = float(ssta)
        if not np.isfinite(ssta) or (
            self._fill_value is not None and ssta == float(self._fill_value)
        ):
            return None

        return round(ssta, 3)

    def _nearest_month_index(self, observation_time: ObservationTime) -> int:
        """Return the matching HadSST month, or the nearest available month."""
        month = (observation_time.year, observation_time.month)
        if (time_index := self._time_indices_by_month.get(month)) is not None:
            return time_index

        month_ordinal = observation_time.year * 12 + observation_time.month
        return int(np.abs(self._month_ordinals - month_ordinal).argmin())

    def _normalise_longitude(self, longitude: float) -> float:
        """Convert a longitude to the convention used by the HadSST grid."""
        if self._longitudes.min() < 0:
            return (longitude + 180) % 360 - 180

        return longitude % 360

    def _nearest_longitude_index(self, longitude: float) -> int:
        """Return the nearest longitude index, including across the dateline."""
        distance = np.abs(self._longitudes - longitude)
        return int(np.minimum(distance, 360 - distance).argmin())

    @staticmethod
    def _nearest_index(values: np.ndarray, target: float) -> int:
        """Return the index of the nearest coordinate value."""
        return int(np.abs(values - target).argmin())
