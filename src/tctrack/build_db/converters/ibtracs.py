"""Module providing the class to convert IBTrACS dataset to TCTrack format."""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from .base import Converter

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", force=True)
logger = logging.getLogger(__name__)


class IBTrACSConverter(Converter):
    """Concrete class that converts IBTrACS dataset to TCTrack format."""

    tctrack_tracker: str = "IBTrACS"

    def __init__(
        self,
        *,
        input_path: str,
        output_path: str,
        start_date: str | None,
        end_date: str | None,
        **kwargs: Any,  # noqa: ARG002
    ) -> None:
        self.output_path = output_path
        self.start_date = start_date
        self.end_date = end_date

        logger.info(f"Loading dataset from {input_path}")
        self.ds = xr.load_dataset(input_path)
        self.ds.attrs["tctrack_tracker"] = self.tctrack_tracker

    def convert(self) -> None:
        """Convert IBTrACS dataset to TCTrack format."""
        self._filter_by_date()
        self._convert_to_tctrack_format()
        self._write()

    def _filter_by_date(self) -> None:
        """Filter the dataset along the time coordinate.

        Raises
        ------
        RuntimeError
            When there is no data for the given date range.
        """
        if self._is_valid_date_range():
            logging.info(
                f"Filtering dataset for date range: {self.start_date} - {self.end_date}"
            )
            self.ds = self.ds.where(
                (self.ds["time"] >= np.datetime64(f"{self.start_date}T00:00:00"))
                & (self.ds["time"] <= np.datetime64(f"{self.end_date}T23:59:59")),
                drop=True,
            )

            if self.ds.time.size == 0:
                msg = "No data available for the given date range."
                raise RuntimeError(msg)

    def _is_valid_date_range(self) -> bool:
        """Check whether the passed date range is valid or not.

        Raises
        ------
        ValueError
            If either of start date or end date is not specified.
            If start date is greater than end date.
            If the date format is not YYYY-MM-DD.

        Returns
        -------
        Boolean value indicating the validity.
        """
        if not (self.start_date or self.end_date):
            return False

        if bool(self.start_date) != bool(self.end_date):
            msg = "Valid start date and end date are required."
            raise ValueError(msg)

        # Validating date format to be YYYY-MM-DD
        datetime.strptime(self.start_date, "%Y-%m-%d")  # type: ignore
        datetime.strptime(self.end_date, "%Y-%m-%d")  # type: ignore

        if self.start_date > self.end_date:  # type: ignore
            msg = "The start date is beyond the end date."
            raise ValueError(msg)

        return True

    def _convert_to_tctrack_format(self) -> None:
        """Workflow for converting IBTrACS to TCTrack format."""
        logging.info("Converting to TCTrack format")
        ds = self.ds.rename_dims({"storm": "trajectory", "date_time": "observation"})
        ds = ds.assign_coords(
            {
                "trajectory": range(len(ds.trajectory)),
                "observation": range(len(ds.observation)),
            }
        )
        self.ds = ds.rename({"lat": "latitude", "lon": "longitude"})
        self._filter_by_data_variables(("wmo_wind", "wmo_pres"))
        self.ds = self.ds.rename_vars(
            {"wmo_wind": "wind_speed", "wmo_pres": "air_pressure_at_sea_level"}
        )
        self._convert_units_for_data_variable("wind_speed", 0.514, "m s-1")
        self._convert_units_for_data_variable("air_pressure_at_sea_level", 100, "Pa")

    def _filter_by_data_variables(self, data_vars: tuple[str, ...]) -> None:
        """Remove all data variables except the ones specified by data_vars.

        Parameters
        ----------
        data_vars : tuple[str, ...]
            The data variables to keep.
        """
        variables_to_remove = [v for v in list(self.ds.data_vars) if v not in data_vars]
        self.ds = self.ds.drop_vars(variables_to_remove)

    def _convert_units_for_data_variable(
        self, data_var: str, scale_factor: float, target_units: str
    ) -> None:
        """Convert the values of the specified data variable to the target units.

        Parameters
        ----------
        data_var : str
            The data variable for which unit conversion is to be done.
        scale_factor : float
            The value to be multiplied with the data variable values in order to
            do the unit conversion.
        target_units : str
            The units to which the conversion is to be done.
        """
        self.ds[data_var] = self.ds[data_var] * scale_factor
        self.ds[data_var].attrs["units"] = target_units

    def _write(self) -> None:
        """Write the Xarray Dataset to NetCDF file.

        Raises
        ------
        ValueError
            If the output path is not a NetCDF file.
        """
        if Path(self.output_path).suffix != ".nc":
            msg = "Output path is not a NetCDF file."
            raise ValueError(msg)
        logging.info("Writing dataset to NetCDF file")
        self.ds.to_netcdf(self.output_path)
