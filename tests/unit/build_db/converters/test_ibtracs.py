"""Tests for IBTrACSConverter."""

# ruff: noqa: C408, ARG002, SLF001
from unittest.mock import MagicMock, call, patch

import numpy as np
import pytest
import xarray as xr

from tctrack.build_db.converters import IBTrACSConverter


@pytest.fixture
def mock_xr_ds_basic() -> xr.Dataset:
    """Mock a basic Xarray.Dataset."""
    return xr.Dataset(
        coords=dict(
            dim_1=("dim_1", ["x1", "x2"]),
            dim_2=("dim_2", ["y1", "y2", "y3"]),
            time=(
                ("dim_1", "dim_2"),
                [
                    [
                        np.datetime64("2025-03-01T12", "ns"),
                        np.datetime64("2025-03-01T13", "ns"),
                        np.datetime64("2025-03-01T14", "ns"),
                    ],
                    [
                        np.datetime64("2025-03-15T12", "ns"),
                        np.datetime64("2025-03-15T13", "ns"),
                        np.datetime64("2025-03-15T14", "ns"),
                    ],
                ],
            ),
        ),
        data_vars=dict(
            var_1=(["dim_1", "dim_2"], [[0, 1, 2], [3, 4, 5]], {"units": "a"}),
            var_2=(["dim_1", "dim_2"], [[0, 1, 2], [3, 4, 5]], {"units": "a"}),
            var_3=(["dim_1", "dim_2"], [[0, 1, 2], [3, 4, 5]], {"units": "a"}),
            var_4=(["dim_1", "dim_2"], [[0, 1, 2], [3, 4, 5]], {"units": "b"}),
            var_5=(["dim_1", "dim_2"], [[0, 1, 2], [3, 4, 5]], {"units": "b"}),
        ),
        attrs=dict(),
    )


@pytest.fixture
def mock_xr_ds_ibtracs() -> xr.Dataset:
    """Mock IBTrACS data."""
    return xr.Dataset(
        coords=dict(
            time=(
                ("storm", "date_time"),
                [
                    [
                        np.datetime64("2025-03-01T12", "ns"),
                        np.datetime64("2025-03-01T13", "ns"),
                        np.datetime64("2025-03-01T14", "ns"),
                    ],
                    [
                        np.datetime64("2025-03-15T12", "ns"),
                        np.datetime64("2025-03-15T13", "ns"),
                        np.datetime64("2025-03-15T14", "ns"),
                    ],
                ],
            ),
            lat=(("storm", "date_time"), [[10.1, 10.2, 10.3], [11.1, 11.2, 11.3]]),
            lon=(("storm", "date_time"), [[80.1, 80.2, 80.3], [90.1, 90.2, 90.3]]),
        ),
        data_vars=dict(
            wmo_wind=(["storm", "date_time"], [[0, 1, 2], [3, 4, 5]], {"units": "kts"}),
            wmo_pres=(["storm", "date_time"], [[0, 1, 2], [3, 4, 5]], {"units": "mb"}),
            var_3=(["storm", "date_time"], [[0, 1, 2], [3, 4, 5]], {"units": "a"}),
            var_4=(["storm", "date_time"], [[0, 1, 2], [3, 4, 5]], {"units": "b"}),
            var_5=(["storm", "date_time"], [[0, 1, 2], [3, 4, 5]], {"units": "c"}),
        ),
        attrs=dict(),
    )


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterInit:
    """Tests for IBTrACSConverter.__init__."""

    def test_load_dataset(self, mock_load_dataset):
        """Test that the input dataset is loaded as an Xarray Dataset."""
        mock_input_path = "/path/to/input"

        ibtracs_converter = IBTrACSConverter(
            input_path=mock_input_path,
            output_path="/path/to/output",
            start_date="2025-01-01",
            end_date="2025-12-31",
        )

        mock_load_dataset.assert_called_once_with(mock_input_path)
        assert ibtracs_converter.ds is mock_load_dataset.return_value

    def test_tctrack_tracker_attribute(self, mock_load_dataset, mock_xr_ds_basic):
        """Tests the following.

        1. IBTrACSConverter.tctrack_tracker is set to "IBTrACS".
        2. "tctrack_tracker" (with value "IBTrACS") is added to Xarray Dataset's attrs.
        """
        mock_load_dataset.return_value = mock_xr_ds_basic

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-01-01",
            end_date="2025-12-31",
        )

        assert ibtracs_converter.tctrack_tracker == "IBTrACS"
        assert ibtracs_converter.ds.attrs["tctrack_tracker"] == "IBTrACS"

    def test_tctrack_tracker_attribute_assignment(
        self, mock_load_dataset, mock_xr_ds_basic
    ):
        """Test that Xarray.Dataset.attrs["tctrack_tracker"] is set.

        It is set from the class variable `tctrack_tracker`.
        """
        mock_load_dataset.return_value = mock_xr_ds_basic
        mock_tracker = "mock_tracker"
        IBTrACSConverter.tctrack_tracker = mock_tracker

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-01-01",
            end_date="2025-12-31",
        )

        assert ibtracs_converter.tctrack_tracker == mock_tracker
        assert ibtracs_converter.ds.attrs["tctrack_tracker"] == mock_tracker


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterIsValidDateRange:
    """Tests for IBTrACSConverter._is_valid_date_range."""

    def test_no_date_range_specified(self, mock_load_dataset):
        """Test that False is returned when dates are not specified."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )

        assert not ibtracs_converter._is_valid_date_range()

    @pytest.mark.parametrize(
        "start_date,end_date",
        [
            (None, "2025-12-31"),
            ("2025-01-01", None),
        ],
    )
    def test_raise_error_on_either_date_not_specified(
        self, mock_load_dataset, start_date, end_date
    ):
        """Test for the case where either start date or end date is not specified."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=start_date,
            end_date=end_date,
        )

        with pytest.raises(ValueError) as excinfo:
            ibtracs_converter._is_valid_date_range()

        assert str(excinfo.value) == "Valid start date and end date are required."

    def test_raise_error_when_start_date_greater_than_end_date(self, mock_load_dataset):
        """Test for the case where start date is greater than end date."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-12-31",
            end_date="2025-01-01",
        )

        with pytest.raises(ValueError) as excinfo:
            ibtracs_converter._is_valid_date_range()

        assert str(excinfo.value) == "The start date is beyond the end date."

    @pytest.mark.parametrize(
        "start_date,end_date",
        [
            ("20250101", "2025-12-31"),
            ("2024-12-31", "20250101"),
        ],
    )
    def test_dates_are_in_iso_format(self, mock_load_dataset, start_date, end_date):
        """Test that the start date and end date are in ISO format (YYYY-MM-DD)."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=start_date,
            end_date=end_date,
        )

        with pytest.raises(ValueError) as excinfo:
            ibtracs_converter._is_valid_date_range()

        assert (
            str(excinfo.value)
            == "time data '20250101' does not match format '%Y-%m-%d'"
        )

    def test_valid_date_range(self, mock_load_dataset):
        """Test for valid date range case."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-01-01",
            end_date="2025-12-31",
        )

        assert ibtracs_converter._is_valid_date_range()


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterFilterByDate:
    """Tests for IBTrACSConverter._filter_by_date."""

    def test_filtering_upon_valid_date_range(self, mock_load_dataset, mock_xr_ds_basic):
        """Test that filtering is applied along the time coordinate."""
        mock_load_dataset.return_value = mock_xr_ds_basic

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-03-01",
            end_date="2025-03-02",
        )
        assert ibtracs_converter.ds["time"].size == 6
        ibtracs_converter._filter_by_date()

        assert ibtracs_converter.ds["time"].size == 3

    def test_end_date_is_exclusive_in_filtering(
        self, mock_load_dataset, mock_xr_ds_basic
    ):
        """Test that filtering is exclusive of end date."""
        mock_load_dataset.return_value = mock_xr_ds_basic

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-03-01",
            end_date="2025-03-15",
        )
        assert ibtracs_converter.ds["time"].size == 6
        ibtracs_converter._filter_by_date()

        assert ibtracs_converter.ds["time"].size == 3

    def test_nodata_for_given_valid_date_range(
        self, mock_load_dataset, mock_xr_ds_basic
    ):
        """Test when specified date range is outside of the time coordinate."""
        mock_load_dataset.return_value = mock_xr_ds_basic

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date="2025-02-01",
            end_date="2025-02-02",
        )
        assert ibtracs_converter.ds["time"].size == 6

        with pytest.raises(RuntimeError) as excinfo:
            ibtracs_converter._filter_by_date()
        assert str(excinfo.value) == "No data available for the given date range."

    def test_no_filtering_done_when_date_range_not_specified(
        self, mock_load_dataset, mocker
    ):
        """Test that filtering is not performed when dates are not specified."""
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )
        spy_where = mocker.spy(ibtracs_converter.ds, "where")

        ibtracs_converter._filter_by_date()

        assert spy_where.call_count == 0


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterFilterByDataVariables:
    """Tests for IBTrACSConverter._filter_by_data_variables."""

    def test_variables_dropped_from_dataset(self, mock_load_dataset, mock_xr_ds_basic):
        """Test that the appropriate variables are dropped from the Dataset."""
        mock_load_dataset.return_value = mock_xr_ds_basic
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )
        assert list(ibtracs_converter.ds.data_vars.keys()) == [
            "var_1",
            "var_2",
            "var_3",
            "var_4",
            "var_5",
        ]

        ibtracs_converter._filter_by_data_variables(("var_2", "var_4"))

        assert list(ibtracs_converter.ds.data_vars.keys()) == ["var_2", "var_4"]


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterConvertUnitsForDataVariable:
    """Tests for IBTrACSConverter._convert_units_for_data_variable."""

    def test_unit_conversion(self, mock_load_dataset, mock_xr_ds_basic):
        """Test the unit conversion of Dataset's data variables."""
        mock_load_dataset.return_value = mock_xr_ds_basic
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )

        assert np.all(ibtracs_converter.ds.var_1.values == [[0, 1, 2], [3, 4, 5]])
        assert ibtracs_converter.ds.var_1.attrs["units"] == "a"
        ibtracs_converter._convert_units_for_data_variable("var_1", 2, "A")

        assert np.all(ibtracs_converter.ds.var_1.values == [[0, 2, 4], [6, 8, 10]])
        assert ibtracs_converter.ds.var_1.attrs["units"] == "A"


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterConvertToTCTrackFormat:
    """Tests for IBTrACSConverter._convert_to_tctrack_format."""

    def test_conversion_workflow(self, mock_load_dataset, mock_xr_ds_ibtracs):
        """Test the workflow for converting IBTrACS to TCTrack format."""
        mock_load_dataset.return_value = mock_xr_ds_ibtracs
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )

        ibtracs_converter._convert_to_tctrack_format()

        assert list(ibtracs_converter.ds.dims) == ["trajectory", "observation"]
        assert list(ibtracs_converter.ds.coords) == [
            "time",
            "latitude",
            "longitude",
            "trajectory",
            "observation",
        ]
        assert list(ibtracs_converter.ds.data_vars) == [
            "wind_speed",
            "air_pressure_at_sea_level",
        ]
        assert np.allclose(
            ibtracs_converter.ds["wind_speed"].values,
            np.array([[0.0, 0.514444, 1.028888], [1.543332, 2.057776, 2.57222]]),
        )
        assert ibtracs_converter.ds["wind_speed"].attrs["units"] == "m s-1"
        assert np.allclose(
            ibtracs_converter.ds["air_pressure_at_sea_level"].values,
            np.array([[0, 100, 200], [300, 400, 500]]),
        )
        assert ibtracs_converter.ds["air_pressure_at_sea_level"].attrs["units"] == "Pa"


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterWrite:
    """Tests for IBTrACSConverter._write."""

    def test_raise_error_if_output_file_not_netcdf(
        self, mock_load_dataset, mock_xr_ds_basic
    ):
        """Test that the output path should be that of a NetCDF file."""
        mock_load_dataset.return_value = mock_xr_ds_basic
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output",
            start_date=None,
            end_date=None,
        )

        with pytest.raises(ValueError) as excinfo:
            ibtracs_converter._write()

        assert str(excinfo.value) == "Output path is not a NetCDF file."

    def test_writing_to_netcdf_file(self, mock_load_dataset):
        """Test writing the Dataset to NetCDF file."""
        mock_load_dataset.return_value = MagicMock()
        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output.nc",
            start_date=None,
            end_date=None,
        )

        ibtracs_converter._write()

        ibtracs_converter.ds.to_netcdf.assert_called_once_with(
            ibtracs_converter.output_path
        )


@patch("tctrack.build_db.converters.ibtracs.xr.load_dataset")
class TestIBTrACSConverterConvert:
    """Tests for IBTrACSConverter.convert."""

    def test_convert(self, mock_load_dataset):
        """Test the steps in the convert method."""
        mock_load_dataset.return_value = MagicMock()

        call_manager = MagicMock()
        call_manager._filter_by_date = IBTrACSConverter._filter_by_date = MagicMock()
        call_manager._convert_to_tctrack_format = (
            IBTrACSConverter._convert_to_tctrack_format
        ) = MagicMock()
        call_manager._write = IBTrACSConverter._write = MagicMock()

        ibtracs_converter = IBTrACSConverter(
            input_path="/path/to/input",
            output_path="/path/to/output.nc",
            start_date=None,
            end_date=None,
        )

        ibtracs_converter.convert()

        ibtracs_converter._filter_by_date.assert_called_once()
        ibtracs_converter._convert_to_tctrack_format.assert_called_once()
        ibtracs_converter._write.assert_called_once()
        # Checking the order of calls
        call_manager.assert_has_calls(
            [
                call._filter_by_date(),
                call._convert_to_tctrack_format(),
                call._write(),
            ]
        )
