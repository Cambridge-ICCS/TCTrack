"""Monthly lookup for Relative Oceanic Niño Index datasets."""

from datetime import datetime
from pathlib import Path
from typing import cast

import cftime
import netCDF4
import numpy as np


class RONI:
    """Read monthly RONI data and provide values by calendar month."""

    def __init__(self, filepath: str | Path) -> None:
        """Load a NetCDF dataset with one-dimensional `time` and `value` variables."""
        path = Path(filepath)
        if not path.is_file():
            msg = f"RONI dataset does not exist or is not a file: {path}"
            raise ValueError(msg)

        try:
            with netCDF4.Dataset(path) as dataset:
                try:
                    time = dataset.variables["time"]
                    values = dataset.variables["value"]
                except KeyError as exc:
                    msg = (
                        f"RONI dataset {path} is missing required variable "
                        f"{exc.args[0]!r}"
                    )
                    raise ValueError(msg) from exc

                if time.ndim != 1 or values.dimensions != time.dimensions:
                    msg = (
                        f"RONI dataset {path} needs one-dimensional time and value "
                        "variables"
                    )
                    raise ValueError(msg)

                if not hasattr(time, "units"):
                    msg = f"RONI dataset {path} time coordinate is missing 'units'"
                    raise ValueError(msg)

                dates = cast(
                    np.ndarray,
                    netCDF4.num2date(
                        time[:], time.units, getattr(time, "calendar", "standard")
                    ),
                )
                self._values_by_month = {
                    (date.year, date.month): value
                    for date, value in zip(dates, values[:], strict=True)
                }
        except OSError as exc:
            msg = f"Could not open RONI dataset {path}: {exc}"
            raise ValueError(msg) from exc

        if not self._values_by_month:
            msg = f"RONI dataset {path} has an empty time coordinate"
            raise ValueError(msg)

    def lookup(self, observation_date: datetime | cftime.datetime) -> float | None:
        """Return the observation month's RONI rounded to 3 decimals, or `None`."""
        value = self._values_by_month.get(
            (observation_date.year, observation_date.month)
        )

        if value is None or np.ma.is_masked(value) or not np.isfinite(value):
            return None

        return round(float(value), 3)
