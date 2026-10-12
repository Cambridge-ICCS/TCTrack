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
        """Load and validate a monthly HadSST dataset containing a `tos` variable."""
        path = Path(filepath)
        if not path.is_file():
            msg = f"HadSST dataset does not exist or is not a file: {path}"
            raise ValueError(msg)

        try:
            with netCDF4.Dataset(path) as dataset:
                try:
                    time = self._coordinate(dataset, "time")
                    latitude = self._coordinate(dataset, "latitude")
                    longitude = self._coordinate(dataset, "longitude")
                    ssta = dataset.variables["tos"]
                except KeyError as exc:
                    msg = (
                        f"HadSST dataset {path} is missing required variable "
                        f"{exc.args[0]!r}"
                    )
                    raise ValueError(msg) from exc

                coordinates = (time, latitude, longitude)
                if any(coordinate.ndim != 1 for coordinate in coordinates):
                    msg = f"HadSST dataset {path} coordinates must be one-dimensional"
                    raise ValueError(msg)

                dimensions = tuple(
                    coordinate.dimensions[0] for coordinate in coordinates
                )
                if ssta.ndim != len(dimensions) or set(ssta.dimensions) != set(
                    dimensions
                ):
                    msg = (
                        f"HadSST dataset {path} variable 'tos' must have time, "
                        "latitude, and longitude dimensions"
                    )
                    raise ValueError(msg)

                if not hasattr(time, "units"):
                    msg = f"HadSST dataset {path} time coordinate is missing 'units'"
                    raise ValueError(msg)

                bounds_name = getattr(time, "bounds", None)
                if bounds_name not in dataset.variables:
                    msg = f"HadSST dataset {path} must define monthly time bounds"
                    raise ValueError(msg)

                bounds = dataset.variables[bounds_name]
                number_of_bounds_dimensions = 2
                if (
                    bounds.ndim != number_of_bounds_dimensions
                    or bounds.dimensions[0] != time.dimensions[0]
                ):
                    msg = f"HadSST dataset {path} has invalid monthly time bounds"
                    raise ValueError(msg)

                if any(coordinate.size == 0 for coordinate in coordinates):
                    msg = f"HadSST dataset {path} has an empty coordinate"
                    raise ValueError(msg)

                dates = cast(
                    np.ndarray,
                    netCDF4.num2date(
                        time[:], time.units, getattr(time, "calendar", "standard")
                    ),
                )
                self._month_ordinals = np.asarray(
                    [self._month_ordinal(date) for date in dates], dtype=int
                )
                if np.unique(self._month_ordinals).size != self._month_ordinals.size:
                    msg = (
                        f"HadSST dataset {path} must contain "
                        "at most one value per month"
                    )
                    raise ValueError(msg)

                self._latitudes = np.asarray(latitude[:], dtype=float)
                self._longitudes = np.asarray(longitude[:], dtype=float)
                # Arrange axes once so lookup always uses time, latitude, longitude.
                # For source dimensions (longitude, time, latitude), the axis
                # indices below are (1, 2, 0).
                self._ssta_values = ssta[:].transpose(
                    tuple(ssta.dimensions.index(dimension) for dimension in dimensions)
                )
                self._fill_value = getattr(ssta, "_FillValue", None)
        except OSError as exc:
            msg = f"Could not open HadSST dataset {path}: {exc}"
            raise ValueError(msg) from exc

    @staticmethod
    def _coordinate(dataset: netCDF4.Dataset, name: str) -> netCDF4.Variable:
        """Return a coordinate variable by its CF standard name or variable name."""
        if name in dataset.variables:
            return dataset.variables[name]

        for variable in dataset.variables.values():
            if getattr(variable, "standard_name", None) == name:
                return variable

        raise KeyError(name)

    def lookup(
        self, observation_time: ObservationTime, latitude: float, longitude: float
    ) -> float | None:
        """Return the nearest HadSST anomaly, rounded to 3 decimals or `None`."""
        time_index = self._nearest_month_index(observation_time)
        latitude_index = self._nearest_index(self._latitudes, latitude)
        longitude_index = self._nearest_longitude_index(longitude)
        ssta = self._ssta_values[time_index, latitude_index, longitude_index]

        if np.ma.is_masked(ssta):
            return None

        ssta = float(ssta)
        if not np.isfinite(ssta) or (
            self._fill_value is not None and ssta == float(self._fill_value)
        ):
            return None

        return round(ssta, 3)

    @staticmethod
    def _month_ordinal(date: ObservationTime) -> int:
        """Express a calendar month as an integer for comparing month distances."""
        return date.year * 12 + date.month

    def _nearest_month_index(self, observation_time: ObservationTime) -> int:
        """Return the matching month, or the nearest available month."""
        # Calendar months also align 360_day dates such as February 30.
        month = self._month_ordinal(observation_time)
        return self._nearest_index(self._month_ordinals, month)

    @staticmethod
    def _nearest_index(values: np.ndarray, target: float) -> int:
        """Return the index of the nearest coordinate value."""
        distances = np.abs(values - target)
        return int(distances.argmin())

    def _nearest_longitude_index(self, longitude: float) -> int:
        """Return the nearest longitude index, including across the dateline."""
        # Wrap offsets to [-180, 180) to measure the shortest way around the globe.
        # This handles both longitude conventions without normalising the input.
        wrapped_offsets = (self._longitudes - longitude + 180) % 360 - 180
        distances = np.abs(wrapped_offsets)
        return int(distances.argmin())
