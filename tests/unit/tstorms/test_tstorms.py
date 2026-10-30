"""
Unit tests for the TSTORMS Python bindings.

Note that these do not require a TSTORMS installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import importlib.metadata
import json
import os
import re
import subprocess
from dataclasses import asdict
from pathlib import Path
from typing import Tuple

import cf
import numpy as np
import pytest
from cftime import datetime
from netCDF4 import Dataset

from tctrack.tstorms import (
    TSTORMSBaseParameters,
    TSTORMSDetectParameters,
    TSTORMSStitchParameters,
    TSTORMSTracker,
    parameter_set,
)


@pytest.fixture
def tstorms_filenames() -> dict[str, str]:
    """Provide filenames for TSTORMSDetectParameters as non-optional/no default."""
    return {
        "u_in_file": "u.nc",
        "v_in_file": "v.nc",
        "vort_in_file": "vort.nc",
        "tm_in_file": "tm.nc",
        "slp_in_file": "slp.nc",
    }


@pytest.fixture
def tstorms_tracker(tmp_path: Path, tstorms_filenames) -> Tuple[TSTORMSTracker, str]:
    """Provide a TSTORMSTracker instance with basic values."""
    # Create a tempdir for tstorms_dir and tstorms_driver (assumed to exist)
    tstorms_dir = tmp_path / "tstorms"
    tstorms_driver_dir = tstorms_dir / "tstorms_driver"
    tstorms_driver_dir.mkdir(parents=True)
    tstorms_trajectory_dir = tstorms_dir / "trajectory_analysis"
    tstorms_trajectory_dir.mkdir(parents=True)
    tstorms_input_dir = tmp_path / "tstorms/input"
    tstorms_input_dir.mkdir(parents=True)

    # Create mock parameters
    tstorms_parameters = TSTORMSBaseParameters(
        tstorms_dir=str(tstorms_dir),
        output_dir=str(tmp_path / "tstorms/output"),
        input_dir=str(tstorms_input_dir),
    )
    detect_params = TSTORMSDetectParameters(
        **tstorms_filenames,
        vort_crit=3.5e-5,
        tm_crit=0.0,
        thick_crit=50.0,
        dist_crit=4.0,
        lat_bound_n=70.0,
        lat_bound_s=-70.0,
        do_spline=False,
        do_thickness=False,
        use_sfc_wind=True,
    )

    # Create a valid u_ref NetCDF file with calendar metadata
    netcdf_file = os.path.join(tstorms_parameters.input_dir, "u.nc")
    u_field = cf.Field(properties={"standard_name": "velocity", "units": "m s-1"})
    u_field.nc_set_variable("u_ref")
    domain_axis_time = cf.DomainAxis(10)
    domain_axis_time.nc_set_unlimited(True)
    _ = u_field.set_construct(domain_axis_time)
    u_data = cf.Data(np.arange(10.0))
    u_field.set_data(u_data)
    dimension_T = cf.DimensionCoordinate(
        properties={
            "standard_name": "time",
            "calendar": "360_day",
            "units": "days since 1950-01-01",
        },
        data=cf.Data(np.arange(10.0)),
    )
    u_field.set_construct(dimension_T)
    # mypy ignore cf.write raises spurious warning but used elsewhere fine.
    cf.write(u_field, netcdf_file)  # type: ignore

    return TSTORMSTracker(tstorms_parameters, detect_params), str(tstorms_dir)


class TestTSTORMSTypes:
    """Tests for the different Classes and Types defined for TSTORMS."""

    def test_base_parameters_initialization(self):
        """Check TSTORMSBaseParameters initializes correctly with valid arguments."""
        params = TSTORMSBaseParameters(
            tstorms_dir="/path/to/tstorms", output_dir="/path/to/output"
        )
        assert params.tstorms_dir == "/path/to/tstorms"
        assert params.input_dir == ""
        assert params.output_dir == "/path/to/output"

    def test_base_parameters_initialization_input_dir(self):
        """Check TSTORMSBaseParameters initializes correctly with input_dir."""
        params = TSTORMSBaseParameters(
            tstorms_dir="/path/to/tstorms",
            input_dir="/path/to/input",
            output_dir="/path/to/output",
        )
        assert params.tstorms_dir == "/path/to/tstorms"
        assert params.input_dir == "/path/to/input"
        assert params.output_dir == "/path/to/output"

    def test_base_parameters_missing_values(self):
        """Check that TypeError is raised when required arguments are missing."""
        with pytest.raises(
            TypeError,
            match=(
                "missing 2 required positional arguments: "
                "'tstorms_dir' and 'output_dir'"
            ),
        ):
            TSTORMSBaseParameters()

    def test_detect_parameters_defaults(self, tstorms_filenames: dict) -> None:
        """Check the default values for TSTORMSDetectParameters."""
        params = TSTORMSDetectParameters(**tstorms_filenames)
        # Check values of all defaults except filenames
        assert params.use_sfc_wind is True
        assert params.vort_crit == 3.5e-5
        assert params.tm_crit == 0.5
        assert params.thick_crit == 50.0
        assert params.dist_crit == 4.0
        assert params.lat_bound_n == 90.0
        assert params.lat_bound_s == -90.0
        assert params.do_spline is False
        assert params.do_thickness is False

    def test_detect_parameters_lat_bounds_error(self, tstorms_filenames: dict) -> None:
        """Check ValueError when northern lat bound is less than southern lat bound."""
        with pytest.raises(
            ValueError,
            match=r"Northern latitude bound.*is less than.*Southern latitude bound",
        ):
            TSTORMSDetectParameters(
                **tstorms_filenames, lat_bound_n=-10.0, lat_bound_s=10.0
            )

    def test_detect_parameters_do_thickness_warning(
        self, tstorms_filenames: dict
    ) -> None:
        """Check UserWarning when do_thickness is set to True."""
        with pytest.warns(
            UserWarning, match="`do_thickness` is set, but will have no effect.*"
        ):
            TSTORMSDetectParameters(**tstorms_filenames, do_thickness=True)

    def test_stitch_parameters_defaults(self) -> None:
        """Check the default values for TSTORMSStitchParameters."""
        params = TSTORMSStitchParameters()
        # Check values of all defaults
        assert params.r_crit == 900.0
        assert params.wind_crit == 17.0
        assert params.vort_crit == 3.5e-5
        assert params.tm_crit == 0.5
        assert params.thick_crit == 50.0
        assert params.n_day_crit == 2
        assert params.do_filter is True
        assert params.lat_bound_n == 40.0
        assert params.lat_bound_s == -40.0
        assert params.do_spline is False
        assert params.do_thickness is False

    def test_stitch_parameters_lat_bounds_error(self) -> None:
        """Check ValueError when northern lat bound is less than southern lat bound."""
        with pytest.raises(
            ValueError,
            match=r"Northern latitude bound.*is less than.*Southern latitude bound",
        ):
            TSTORMSStitchParameters(lat_bound_n=-10.0, lat_bound_s=10.0)

    def test_stitch_parameters_do_thickness_warning(self) -> None:
        """Check UserWarning when do_thickness is set to True."""
        with pytest.warns(
            UserWarning, match="`do_thickness` is set, but will have no effect.*"
        ):
            TSTORMSStitchParameters(do_thickness=True)


class TestTSTORMSParameterSet:
    """Tests for the parameter_set function."""

    @pytest.mark.parametrize("name", ["default", "Vitart2001"])
    @pytest.mark.parametrize("use_sfc_wind", [True, False])
    def test_parameter_set_vitart(self, name, use_sfc_wind) -> None:
        """Check the Vitart parameter set values."""
        base, detect, stitch = parameter_set(
            name, "/path/to/tstorms", "/path/to/output", use_sfc_wind
        )

        assert base.tstorms_dir == "/path/to/tstorms"
        assert base.output_dir == "/path/to/output"
        assert detect.use_sfc_wind is use_sfc_wind
        assert detect.tm_crit == 0.5
        assert detect.lat_bound_n == 70.0
        assert detect.lat_bound_s == -70.0
        assert stitch.r_crit == 800.0
        assert stitch.wind_crit == 12.0 if use_sfc_wind else 17.0
        assert stitch.tm_crit == 0.5
        assert stitch.lat_bound_n == 70.0
        assert stitch.lat_bound_s == -70.0

    def test_parameter_set_invalid_name(self) -> None:
        """Check parameter_set fails for unknown names."""
        with pytest.raises(ValueError, match="Unknown parameter set: unknown"):
            parameter_set("unknown", "/path/to/tstorms", "/path/to/output")


class TestTSTORMSTracker:
    """Tests for the TSTORMS Tracker class."""

    def test_initialisation_no_inputs(self, tmp_path):
        """Test tracker initialisation without inputs in detect parameters."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)

        assert tracker._u_in_file is None
        assert tracker._v_in_file is None
        assert tracker._vort_in_file is None
        assert tracker._tm_in_file is None
        assert tracker._slp_in_file is None

    def test_initialisation_inputs(self, tmp_path):
        """Test tracker initialisation with inputs in detect parameters."""
        input_dir = str(tmp_path / "tstorms/input")
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
            input_dir=input_dir,
        )
        detect_params = TSTORMSDetectParameters(
            u_in_file="u.nc",
            v_in_file="v.nc",
            vort_in_file="vort.nc",
            tm_in_file="tm.nc",
            slp_in_file="slp.nc",
        )
        tracker = TSTORMSTracker(base_params, detect_params)

        assert tracker._u_in_file == os.path.join(input_dir, "u.nc")
        assert tracker._v_in_file == os.path.join(input_dir, "v.nc")
        assert tracker._vort_in_file == os.path.join(input_dir, "vort.nc")
        assert tracker._tm_in_file == os.path.join(input_dir, "tm.nc")
        assert tracker._slp_in_file == os.path.join(input_dir, "slp.nc")

    @staticmethod
    def _create_valid_tstorms_files(
        tmp_path: Path,
        *,
        use_sfc_wind: bool = True,
        time_vals: np.ndarray | None = None,
        lat_vals: np.ndarray | None = None,
        unlimited_time: bool = True,
        extra_var: bool = False,
        missing_dim: str | None = None,
    ) -> list[str]:
        """Create valid NetCDF files for all 5 TSTORMS variables."""
        u_var = "u_ref" if use_sfc_wind else "u850"
        v_var = "v_ref" if use_sfc_wind else "v850"
        var_names = [u_var, v_var, "vort850", "tm", "slp"]
        n_lat, n_lon, n_time = 10, 20, 3
        lats = lat_vals if lat_vals is not None else np.linspace(-60, 60, n_lat)
        lons = np.linspace(0, 360, n_lon)
        times = time_vals if time_vals is not None else np.array([1.0, 2.0, 3.0])

        inputs = []
        for var_name in var_names:
            filename = str(tmp_path / f"{var_name}.nc")
            with Dataset(filename, "w") as nc_file:
                if missing_dim != "time":
                    nc_file.createDimension("time", None if unlimited_time else n_time)
                if missing_dim != "lat":
                    nc_file.createDimension("lat", n_lat)
                if missing_dim != "lon":
                    nc_file.createDimension("lon", n_lon)

                if missing_dim != "time":
                    t_var = nc_file.createVariable("time", "f4", ("time",))
                    t_var[:] = times
                if missing_dim != "lat":
                    la_var = nc_file.createVariable("lat", "f4", ("lat",))
                    la_var[:] = lats
                if missing_dim != "lon":
                    lo_var = nc_file.createVariable("lon", "f4", ("lon",))
                    lo_var[:] = lons

                dims = tuple(d for d in ("time", "lat", "lon") if d != missing_dim)
                shape = tuple(
                    n_time if d == "time" else (n_lat if d == "lat" else n_lon)
                    for d in dims
                )
                data_var = nc_file.createVariable(var_name, "f4", dims)
                data_var[:] = np.zeros(shape)

                if extra_var:
                    extra = nc_file.createVariable("extra_field", "f4", dims)
                    extra[:] = np.zeros(shape)
            inputs.append(filename)
        return inputs

    def test_set_input_files(self, tmp_path):
        """Test the variable input files are automatically set by set_input_files."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)

        inputs = self._create_valid_tstorms_files(tmp_path)
        tracker.set_input_files(inputs)

        # Check the attributes have been set correctly
        assert tracker._u_in_file == str(tmp_path / "u_ref.nc")
        assert tracker._v_in_file == str(tmp_path / "v_ref.nc")
        assert tracker._vort_in_file == str(tmp_path / "vort850.nc")
        assert tracker._tm_in_file == str(tmp_path / "tm.nc")
        assert tracker._slp_in_file == str(tmp_path / "slp.nc")

    def test_verify_inputs_use_sfc_wind_false(self, tmp_path):
        """Test verification passes with u850/v850 when use_sfc_wind is False."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        detect_params = TSTORMSDetectParameters(use_sfc_wind=False)
        tracker = TSTORMSTracker(base_params, detect_params=detect_params)

        inputs = self._create_valid_tstorms_files(tmp_path, use_sfc_wind=False)
        tracker.set_input_files(inputs)
        assert tracker._u_in_file == str(tmp_path / "u850.nc")
        assert tracker._v_in_file == str(tmp_path / "v850.nc")

    @pytest.mark.parametrize("missing_dim", ["time", "lat", "lon"])
    def test_verify_inputs_missing_dimension(self, tmp_path, missing_dim):
        """Test verification fails when a required dimension is missing."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(tmp_path, missing_dim=missing_dim)
        with pytest.raises(ValueError, match=f"missing dimension '{missing_dim}'"):
            tracker.set_input_files(inputs)

    def test_verify_inputs_time_not_unlimited(self, tmp_path):
        """Test verification fails when time dimension is not unlimited."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(tmp_path, unlimited_time=False)
        with pytest.raises(ValueError, match="must be an unlimited dimension"):
            tracker.set_input_files(inputs)

    def test_verify_inputs_negative_time(self, tmp_path):
        """Test verification fails when time coordinates contain negative values."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(
            tmp_path, time_vals=np.array([-1.0, 0.0, 1.0])
        )
        with pytest.raises(ValueError, match="only have positive values"):
            tracker.set_input_files(inputs)

    def test_verify_inputs_decreasing_latitude(self, tmp_path):
        """Test verification fails when latitude coordinates are decreasing."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(
            tmp_path, lat_vals=np.linspace(60, -60, 10)
        )
        with pytest.raises(ValueError, match="strictly increasing from south to north"):
            tracker.set_input_files(inputs)

    def test_verify_inputs_missing_required_variable(self, tmp_path):
        """Test verification fails when not all 5 required variables are provided."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(tmp_path)
        # Exclude slp.nc
        inputs = [f for f in inputs if not f.endswith("slp.nc")]
        with pytest.raises(ValueError, match="TSTORMS requires all 5 input variables"):
            tracker.set_input_files(inputs)

    def test_verify_inputs_grid_mismatch(self, tmp_path):
        """Test verification fails when files have different coordinate grids."""
        base_params = TSTORMSBaseParameters(
            tstorms_dir=str(tmp_path / "tstorms"),
            output_dir=str(tmp_path / "tstorms/output"),
        )
        tracker = TSTORMSTracker(base_params)
        inputs = self._create_valid_tstorms_files(tmp_path)
        # Overwrite slp.nc with different latitude coordinates
        with Dataset(str(tmp_path / "slp.nc"), "w") as nc_file:
            nc_file.createDimension("time", None)
            nc_file.createDimension("lat", 10)
            nc_file.createDimension("lon", 20)
            nc_file.createVariable("time", "f4", ("time",))[:] = [1.0, 2.0, 3.0]
            nc_file.createVariable("lat", "f4", ("lat",))[:] = np.linspace(-50, 50, 10)
            nc_file.createVariable("lon", "f4", ("lon",))[:] = np.linspace(0, 360, 20)
            nc_file.createVariable("slp", "f4", ("time", "lat", "lon"))[:] = np.zeros(
                (3, 10, 20)
            )

        with pytest.raises(
            ValueError, match=r"Latitude coordinates in .* do not match"
        ):
            tracker.set_input_files(inputs)

    def test_write_driver_namelist(self, tstorms_tracker):
        """Test the generation of the driver namelist."""
        tracker = tstorms_tracker[0]

        # Call the method
        namelist_path = tracker._write_driver_namelist()

        # Verify the file was created inside tstorms_driver/
        assert os.path.exists(namelist_path)
        assert namelist_path == os.path.join(
            tracker.tstorms_parameters.output_dir, "nml_driver"
        )

        # Verify the content of the file
        input_dir = tracker.tstorms_parameters.input_dir
        with open(namelist_path, "r") as namelist_file:
            content = namelist_file.read()
            assert "crit_vort  =  3.5000E-05" in content
            assert "crit_twc   =  0.0000" in content
            assert "crit_thick =  50.0000" in content
            assert "crit_dist  =   4.0000" in content
            assert "lat_bound_n =  70.0000" in content
            assert "lat_bound_s = -70.0000" in content
            assert "do_spline   = .false." in content
            assert "do_thickness= .false." in content
            assert f"fn_u    = '{input_dir}/u.nc'" in content
            assert f"fn_v    = '{input_dir}/v.nc'" in content
            assert f"fn_vort = '{input_dir}/vort.nc'" in content
            assert f"fn_tm   = '{input_dir}/tm.nc'" in content
            assert f"fn_slp  = '{input_dir}/slp.nc'" in content
            assert "use_sfc_wnd = .true." in content

    def test_write_trajectory_namelist(self, tstorms_tracker):
        """Test the generation of the trajectory namelist."""
        tracker = tstorms_tracker[0]
        tstorms_dir = tracker.tstorms_parameters.tstorms_dir

        # Call the method
        namelist_path = tracker._write_trajectory_analysis_namelist()

        # Verify the file was created inside tstorms_driver/
        assert os.path.exists(namelist_path)
        assert namelist_path == os.path.join(
            tracker.tstorms_parameters.output_dir, "nml_traj"
        )

        # Verify the content of the file
        with open(namelist_path, "r") as namelist_file:
            content = namelist_file.read()
            assert "rcrit      = 900.0000" in content
            assert "wcrit      = 17.0000" in content
            assert "vcrit      = 3.5000E-05" in content
            assert "twc_crit   = 0.5000" in content
            assert "thick_crit = 50.0000" in content
            assert "nwcrit     = 2.0000" in content
            assert "do_filt    = .true." in content
            assert "nlat = 40.0000" in content
            assert "slat = -40.0000" in content
            assert "do_spline    = .false." in content
            assert "do_thickness = .false." in content
            assert (
                f"landmask = '{tstorms_dir}/trajectory_analysis/landsea.map'" in content
            )
            assert f"cmask = '{tstorms_dir}/trajectory_analysis/imask_2'" in content

    def test_set_time_metadata_valid(self, tstorms_tracker):
        """Test extracting valid calendar metadata from a NetCDF file."""
        tracker = tstorms_tracker[0]

        # Assert the metadata extracted correctly
        tracker._set_time_metadata()
        assert tracker.time_metadata == {
            "calendar": "360_day",
            "units": "days since 1950-01-01",
            "start_time": datetime(1950, 1, 1, 00, calendar="360_day"),
            "end_time": datetime(1950, 1, 10, 00, calendar="360_day"),
        }

    def test_set_time_metadata_file_not_found(self, tstorms_tracker):
        """Test FileNotFoundError raised when the u NetCDF file is not found."""
        tracker = tstorms_tracker[0]

        # Remove the u_ref file to simulate file-not-found
        u_ref_path = os.path.join(
            tracker.tstorms_parameters.input_dir, tracker.detect_parameters.u_in_file
        )
        if os.path.exists(u_ref_path):
            os.remove(u_ref_path)

        with pytest.raises(FileNotFoundError):
            tracker._set_time_metadata()

    def test_set_time_metadata_no_unlimited_dim(self, tstorms_tracker):
        """Test behavior when no unlimited dimension is found in the NetCDF file."""
        tracker = tstorms_tracker[0]

        # Overwrite NetCDF file in the tracker fixture, time not unlimited
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
        u_field.nc_set_variable("u_ref")
        domain_axis_time = cf.DomainAxis(10)
        _ = u_field.set_construct(domain_axis_time)
        u_data = cf.Data(np.arange(10.0))
        u_field.set_data(u_data)
        dimension_T = cf.DimensionCoordinate(
            properties={
                "standard_name": "time",
                "calendar": "360_day",
                "units": "days since 1950-01-01",
            },
            data=cf.Data(np.arange(10.0)),
        )
        u_field.set_construct(dimension_T)
        cf.write(u_field, netcdf_file)

        # Call the method and assert the default metadata is set
        with pytest.raises(
            KeyError,
            match=(
                r"No unlimited dimension found in the u_ref NetCDF file to set "
                "calendar."
            ),
        ):
            tracker._set_time_metadata()

    def test_set_time_metadata_multiple_unlimited_dim(self, tstorms_tracker):
        """Test for error when multiple unlimited dimensions found in NetCDF file."""
        tracker = tstorms_tracker[0]

        # Overwrite NetCDF file in the tracker fixture, time not unlimited
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
        u_field.nc_set_variable("u_ref")
        domain_axis_time = cf.DomainAxis(5)
        domain_axis_time.nc_set_unlimited(True)
        domain_axis_lat = cf.DomainAxis(2)
        domain_axis_lat.nc_set_unlimited(True)
        _ = u_field.set_construct(domain_axis_time)
        _ = u_field.set_construct(domain_axis_lat)
        u_data = cf.Data(np.ones([5, 2]))
        u_field.set_data(u_data)
        dimension_T = cf.DimensionCoordinate(
            properties={
                "standard_name": "time",
                "calendar": "360_day",
                "units": "days since 1950-01-01",
            },
            data=cf.Data(np.arange(5.0)),
        )
        dimension_lat = cf.DimensionCoordinate(
            properties={
                "standard_name": "latitude",
                "units": "degree_north",
            },
            data=cf.Data(np.arange(2.0)),
        )
        u_field.set_construct(dimension_T)
        u_field.set_construct(dimension_lat)
        cf.write(u_field, netcdf_file)

        # Call the method and assert the default metadata is set
        with pytest.raises(
            ValueError,
            match=re.escape(
                "TSTORMS expects only a single unlimited variable in NetCDF files "
                "corresponding to the time dimension. "
                "Multiple found: ['time', 'latitude']."
            ),
        ):
            tracker._set_time_metadata()

    def test_set_time_metadata_no_coord_var(self, tstorms_tracker):
        """Test behavior when missing unlimited dimension coordinate data."""
        tracker = tstorms_tracker[0]

        # Overwrite temporary NetCDF file in the tracker fixture
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
        u_field.nc_set_variable("u_ref")
        domain_axis_time = cf.DomainAxis(10)
        domain_axis_time.nc_set_unlimited(True)
        _ = u_field.set_construct(domain_axis_time)
        u_data = cf.Data(np.arange(10.0))
        u_field.set_data(u_data)
        cf.write(u_field, netcdf_file)

        # Call the method and assert the default metadata is set
        with pytest.raises(
            KeyError,
            match=r"Coordinate variable for unlimited dimension 'dim' not found.",
        ):
            tracker._set_time_metadata()

    @pytest.mark.parametrize(
        "dimension_properties, expected_error, expected_warning, expected_metadata",
        [
            (
                {"standard_name": "time"},  # Missing 'units' and 'calendar'
                ValueError,
                None,
                None,
            ),
            (
                {
                    "standard_name": "time",
                    "units": "days since 2000-01-01",
                },  # Missing 'calendar'
                None,
                UserWarning,
                {
                    "calendar": "julian",
                    "units": "days since 2000-01-01",
                    "start_time": datetime(2000, 1, 1, 0, calendar="julian"),
                    "end_time": datetime(2000, 1, 10, 0, calendar="julian"),
                },
            ),
        ],
    )
    def test_set_time_metadata_missing_attributes(
        self,
        tstorms_tracker,
        dimension_properties,
        expected_error,
        expected_warning,
        expected_metadata,
    ):
        """Test behavior when missing 'units' or 'calendar' attributes from coord."""
        tracker = tstorms_tracker[0]

        # Overwrite temporary NetCDF file in the tracker fixture
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
        u_field.nc_set_variable("u_ref")
        domain_axis_time = cf.DomainAxis(10)
        domain_axis_time.nc_set_unlimited(True)
        _ = u_field.set_construct(domain_axis_time)
        u_data = cf.Data(np.arange(10.0))
        u_field.set_data(u_data)
        dimension_T = cf.DimensionCoordinate(
            properties=dimension_properties,
            data=cf.Data(np.arange(10.0)),
        )
        u_field.set_construct(dimension_T)
        cf.write(u_field, netcdf_file)

        # Ensure that code runs with errors/warnings raised as expected
        if expected_error:
            with pytest.raises(expected_error):
                tracker._set_time_metadata()
        elif expected_warning:
            with pytest.warns(
                UserWarning,
                match=(
                    "The 'calendar' attribute is missing for the coordinate variable.\n"
                    "defaulting to 'julian'"
                ),
            ):
                tracker._set_time_metadata()
                assert tracker.time_metadata == expected_metadata
        else:
            tracker._set_time_metadata()
            assert tracker.time_metadata == expected_metadata

    def test_set_metadata(self, tstorms_tracker):
        """Test setting of global and variable metadata from NetCDF files."""
        tracker = tstorms_tracker[0]
        input_dir = tracker.tstorms_parameters.input_dir

        # Create additional NetCDF files for other variables
        # Check variables with and without a long name provided
        variables = [
            (
                "slp.nc",
                "slp",
                "air_pressure_at_mean_sea_level",
                "Sea Level Pressure",
                "Pa",
            ),
            ("vort.nc", "vort850", "atmosphere_upward_relative_vorticity", None, "s-1"),
        ]

        for filename, var_name, std_name, long_name, units in variables:
            netcdf_file = os.path.join(input_dir, filename)
            properties = {"standard_name": std_name, "units": units}
            if long_name:
                properties["long_name"] = long_name
            field = cf.Field(properties=properties)
            field.nc_set_variable(var_name)
            domain_axis = cf.DomainAxis(10)
            _ = field.set_construct(domain_axis)
            data = cf.Data(np.arange(10.0))
            field.set_data(data)
            cf.write(field, netcdf_file)

        # Call the method
        tracker.set_metadata()

        # Assert the global metadata was set correctly
        assert tracker.global_metadata is not None
        assert tracker.global_metadata == {
            "tctrack_version": importlib.metadata.version("tctrack"),
            "tctrack_tracker": "TSTORMSTracker",
            "tctrack_parameters": json.dumps(
                {
                    "TSTORMSBaseParameters": asdict(tracker.tstorms_parameters),
                    "TSTORMSDetectParameters": asdict(tracker.detect_parameters),
                    "TSTORMSStitchParameters": asdict(tracker.stitch_parameters),
                }
            ),
        }

        # Assert the variable metadata was read correctly
        expected_metadata = {
            "wind_speed": {
                "standard_name": "wind_speed",
                "long_name": "Surface Wind Speed",
                "units": "m s-1",
                "cell_method": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={"comment": "lesser circle of radius 4.0 degrees"},
                ),
            },
            "air_pressure_at_mean_sea_level": {
                "standard_name": "air_pressure_at_mean_sea_level",
                "long_name": "Sea Level Pressure",
                "units": "Pa",
                "cell_method": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={"comment": "lesser circle of radius 4.0 degrees"},
                ),
            },
            "atmosphere_upward_relative_vorticity": {
                "standard_name": "atmosphere_upward_relative_vorticity",
                "long_name": "Atmosphere Upward Relative Vorticity",
                "units": "s-1",
                "cell_method": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={"comment": "lesser circle of radius 4.0 degrees"},
                ),
            },
        }

        for var_name, metadata in expected_metadata.items():
            assert var_name in tracker._variable_metadata
            assert (
                tracker._variable_metadata[var_name].properties["standard_name"]
                == metadata["standard_name"]
            )
            assert (
                tracker._variable_metadata[var_name].properties["long_name"]
                == metadata["long_name"]
            )
            assert (
                tracker._variable_metadata[var_name].properties["units"]
                == metadata["units"]
            )
            assert tracker._variable_metadata[var_name].constructs == [
                metadata["cell_method"]
            ]

    def test_set_metadata_missing_file(self, tstorms_tracker):
        """Test metadata behaviour when a required NetCDF file is missing."""
        tracker = tstorms_tracker[0]

        # Ensure removal one of the input files (slp is read first)
        missing_file = os.path.join(tracker.tstorms_parameters.input_dir, "slp.nc")
        if os.path.exists(missing_file):
            os.remove(missing_file)

        with pytest.raises(
            FileNotFoundError,
            match=(
                f"Input file '{tracker.tstorms_parameters.input_dir}/slp.nc' not found."
            ),
        ):
            tracker.set_metadata()

    def test_set_metadata_missing_variable(self, tstorms_tracker):
        """Test reading variable metadata if default TSTORMS variable missing."""
        tracker = tstorms_tracker[0]
        input_dir = tracker.tstorms_parameters.input_dir

        # Create slp file with no variable names slp
        netcdf_file = os.path.join(input_dir, "slp.nc")
        field = cf.Field(
            properties={
                "standard_name": "air_pressure_at_mean_sea_level",
                "units": "Pa",
                "long_name": "Sea Level Pressure",
            }
        )
        field.nc_set_variable("incorrectly_named_slp_for_tstorms")
        domain_axis = cf.DomainAxis(10)
        _ = field.set_construct(domain_axis)
        data = cf.Data(np.arange(10.0))
        field.set_data(data)
        cf.write(field, netcdf_file)

        with pytest.raises(
            ValueError,
            match=(rf"Variable 'slp' not found in input file {input_dir}/slp.nc."),
        ):
            tracker.set_metadata()

    def _mock_trajectories_data(self):
        """Generate mock trajectory data for testing."""
        return {
            "varnames": [
                "lon",
                "lat",
                "wind_speed",
                "air_pressure_at_mean_sea_level",
                "atmosphere_upward_relative_vorticity",
            ],
            "datenames": ["year", "month", "day", "hour"],
            "trajectories": [
                [
                    {
                        "date": ["1950", "1", "6", "0"],
                        "line": [
                            "67.32",
                            "-9.26",
                            "17.90",
                            "1002.11",
                            "0.00065708",
                        ],
                    },
                    {
                        "date": ["1950", "1", "7", "0"],
                        "line": [
                            "67.68",
                            "-9.26",
                            "15.32",
                            "1003.88",
                            "0.00044106",
                        ],
                    },
                ],
                [
                    {
                        "date": ["1950", "1", "17", "0"],
                        "line": [
                            "105.29",
                            "-17.23",
                            "17.79",
                            "997.42",
                            "0.00058948",
                        ],
                    },
                    {
                        "date": ["1950", "1", "18", "0"],
                        "line": [
                            "105.00",
                            "-19.11",
                            "17.65",
                            "999.99",
                            "0.00052948",
                        ],
                    },
                    {
                        "date": ["1950", "1", "19", "0"],
                        "line": [
                            "104.94",
                            "-16.99",
                            "18.12",
                            "998.71",
                            "0.00049719",
                        ],
                    },
                ],
            ],
        }

    @pytest.fixture
    def mock_trav_file(self, tmp_path: Path):
        """Fixture to create a mock TSTORMS 'trav' file with sample trajectory data."""
        file_path = tmp_path / "trav_filt"
        mock_data = self._mock_trajectories_data()

        content = ""
        for track in mock_data["trajectories"]:
            content += "start\t{}\t{}\n".format(len(track), "\t".join(track[0]["date"]))
            content += "\n".join(
                "\t{}\t{}".format("\t".join(obs["line"]), "\t".join(obs["date"]))
                for obs in track
            )
            content += "\n"

        file_path.write_text(content)
        return str(file_path)

    def test_read_trajectories_valid(self, tstorms_tracker, mock_trav_file):
        """Test parsing valid trajectories from a TSTORMS 'trav' file."""
        tracker = tstorms_tracker[0]

        # Update the tracker to use the temporary file
        tracker.tstorms_parameters.output_dir = str(os.path.dirname(mock_trav_file))

        trajectories = tracker.read_trajectories()

        # Assert the trajectories were parsed correctly
        assert len(trajectories) == 2

        traj1 = trajectories[0]
        assert traj1.trajectory_id == 1
        assert traj1.observations == 2
        assert traj1.start_time == datetime(
            1950,
            1,
            6,
            0,
            calendar=tracker.time_metadata["calendar"],
        )
        assert traj1.observations == 2

        # Validate the data attribute for the first trajectory
        expected_data_trajectory_1 = {
            "time": [
                datetime(
                    1950,
                    1,
                    6,
                    0,
                    calendar=tracker.time_metadata["calendar"],
                ),
                datetime(
                    1950,
                    1,
                    7,
                    0,
                    calendar=tracker.time_metadata["calendar"],
                ),
            ],
            "lon": [67.32, 67.68],
            "lat": [-9.26, -9.26],
            "wind_speed": [17.90, 15.32],
            "air_pressure_at_mean_sea_level": [1002.11, 1003.88],
            "atmosphere_upward_relative_vorticity": [0.00065708, 0.00044106],
        }

        assert traj1.data == expected_data_trajectory_1

        traj2 = trajectories[1]
        assert traj2.trajectory_id == 2
        assert traj2.observations == 3
        assert traj2.start_time == datetime(
            1950,
            1,
            17,
            0,
            calendar=tracker.time_metadata["calendar"],
        )
        assert traj2.observations == 3
        expected_data_trajectory_2 = {
            "lon": [105.29, 105.00, 104.94],
            "lat": [-17.23, -19.11, -16.99],
            "wind_speed": [17.79, 17.65, 18.12],
            "air_pressure_at_mean_sea_level": [997.42, 999.99, 998.71],
            "atmosphere_upward_relative_vorticity": [
                0.00058948,
                0.00052948,
                0.00049719,
            ],
        }
        for key, expected_values in expected_data_trajectory_2.items():
            assert traj2.data[key] == expected_values

    def test_run_tracker_success(self, mocker, tstorms_tracker, mock_trav_file) -> None:
        """Check run_tracker runs successfully."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = mocker.MagicMock(
            returncode=0, stdout="Success"
        )

        # Check run_tracker runs without error and produces an output file
        tracker = tstorms_tracker[0]
        input_dir = tracker.tstorms_parameters.input_dir

        # Create additional NetCDF files for other variables
        # Check variables with and without a long name provided
        variables = [
            (
                "slp.nc",
                "slp",
                "air_pressure_at_mean_sea_level",
                "Sea Level Pressure",
                "Pa",
            ),
            ("vort.nc", "vort850", "atmosphere_upward_relative_vorticity", None, "s-1"),
        ]

        for filename, var_name, std_name, long_name, units in variables:
            netcdf_file = os.path.join(input_dir, filename)
            properties = {"standard_name": std_name, "units": units}
            if long_name:
                properties["long_name"] = long_name
            field = cf.Field(properties=properties)
            field.nc_set_variable(var_name)
            domain_axis = cf.DomainAxis(10)
            _ = field.set_construct(domain_axis)
            data = cf.Data(np.arange(10.0))
            field.set_data(data)
            # mypy ignore cf.write raises spurious warning but used elsewhere fine.
            cf.write(field, netcdf_file)  # type: ignore

        # Update the tracker to use the temporary trav file
        # Create a dummy cyclones file
        output_dir = os.path.dirname(mock_trav_file)
        tracker.tstorms_parameters.output_dir = output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        output_file = os.path.join(output_dir, "trajectories.nc")
        tracker.run_tracker([], output_file)
        assert os.path.exists(output_file)

    def test_run_tracker_failure(self, mocker, tstorms_tracker, mock_trav_file) -> None:
        """Check run_tracker propagates RuntimeError from detect/stitch."""
        # Create tracker object
        tracker = tstorms_tracker[0]

        # Update the tracker to use the temporary trav file
        # Create a dummy cyclones file
        output_dir = str(os.path.dirname(mock_trav_file))
        tracker.tstorms_parameters.output_dir = output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        # Mock subprocess.run and define a function to mock fail for a specific command
        mock_subprocess_run = mocker.patch("subprocess.run")

        def subprocess_failure(cmd, args, **_kwargs):
            """
            Set up a mocked subprocess failure if a subprocess call contains cmd string.

            ``cmd`` is a string contained in the subprocess call for which to generate
            failure.
            Can be invoked as a mocker.patch side_effect using lambdas as follows:

            >>> patched_sp_run = mocker.patch("subprocess.run")
            >>> patched_sp_run.side_effect = lambda args, **kwargs: subprocess_failure(
            >>>     "str in failing command", args, **kwargs
            >>> )
            """
            if cmd in args[0]:
                raise subprocess.CalledProcessError(
                    returncode=1, cmd=cmd, stderr="Error occurred"
                )
            mock_output = mocker.Mock(returncode=0, stdout="Success")
            return mock_output

        # Check a RuntimeError is correctly raised for detect (tstorms_driver)
        mock_subprocess_run.side_effect = lambda args, **kwargs: subprocess_failure(
            "tstorms_driver.exe", args, **kwargs
        )
        with pytest.raises(
            RuntimeError, match="Detect failed with a non-zero exit code"
        ):
            tracker.run_tracker([], "trajectories.nc")

        # Check a RuntimeError is correctly raised for stitch (trajectory_analysis)
        mock_subprocess_run.side_effect = lambda args, **kwargs: subprocess_failure(
            "trajectory_analysis_csc.exe", args, **kwargs
        )
        with pytest.raises(
            RuntimeError, match="Stitch failed with a non-zero exit code"
        ):
            tracker.run_tracker([], "trajectories.nc")


class TestTSTORMSTrackerDetect:
    """Tests for the detect functionality of TSTORMSTracker."""

    def test_tstorms_tracker_detect(self, mocker, tstorms_tracker) -> None:
        """Checks the correct tstorms_driver call is made."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )
        # Mock warning for no cyclones file being generated.
        mocker.patch("warnings.warn")

        tracker = tstorms_tracker[0]
        tstorms_dir = tstorms_tracker[1]
        output_dir = tracker.tstorms_parameters.output_dir
        result = tracker.detect()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [f"{tstorms_dir}/tstorms_driver/tstorms_driver.exe"],
            stdin=mocker.ANY,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(output_dir),
            input=None,
        )
        stdin_file = mock_subprocess_run.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_driver"
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_detect_verbose(
        self, mocker, tstorms_tracker, capsys
    ) -> None:
        """Checks the correct tstorms_driver call is made with verbosity=2."""
        # Mock subprocess.Popen to simulate real-time output
        mock_popen = mocker.patch("subprocess.Popen")
        mock_stdout = mocker.MagicMock()
        mock_stdout.readline.side_effect = ["Line 1\n", "Line 2\n", ""]
        mock_process = mocker.MagicMock(
            stdout=mock_stdout,
            stderr="Mocked stderr output",
            returncode=0,
        )
        # Ensure communicate returns a tuple with stdout and stderr
        mock_process.communicate.return_value = ("", "Mocked stderr output")
        mock_popen.return_value = mock_process
        # Mock warning for no cyclones file being generated.
        mocker.patch("warnings.warn")

        tracker = tstorms_tracker[0]
        tstorms_dir = tstorms_tracker[1]
        output_dir = tracker.tstorms_parameters.output_dir
        result = tracker.detect(verbosity=2)

        # Check subprocess.Popen call
        mock_popen.assert_called_once_with(
            [f"{tstorms_dir}/tstorms_driver/tstorms_driver.exe"],
            stdin=mocker.ANY,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False,
            bufsize=1,
            cwd=str(output_dir),
        )

        # Verify stdin file path
        stdin_file = mock_popen.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_driver"

        # Verify returned outputs - using capsys to get the stderr feed
        captured = capsys.readouterr()
        assert captured.out == "Executing Detect...\nLine 1\nLine 2\n"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_detect_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check detect raises FileNotFoundError when executable is missing."""
        # Mock subprocess.run to simulate a FileNotFoundError
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TSTORMSTracker instance
        tracker = tstorms_tracker[0]

        # Assert that detect raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Detect failed because the executable could not be found",
        ):
            tracker.detect()

    def test_tstorms_tracker_detect_verbose_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check detect raises FileNotFoundError verbosity=2."""
        # Mock subprocess.Popen to simulate a FileNotFoundError
        mock_popen = mocker.patch("subprocess.Popen")
        mock_popen.side_effect = FileNotFoundError("Executable not found")

        tracker = tstorms_tracker[0]

        # Assert that detect raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Detect failed because the executable could not be found",
        ):
            tracker.detect(verbosity=2)

    def test_tstorms_tracker_detect_failure(self, mocker, tstorms_tracker) -> None:
        """Check detect raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="tstorms_driver", stderr="Error occurred"
        )

        # Create a TSTORMSTracker instance
        tracker = tstorms_tracker[0]

        # Assert that detect raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Detect failed with a non-zero exit code"
        ):
            tracker.detect()

    def test_tstorms_tracker_detect_verbose_failure(
        self,
        mocker,
        tstorms_tracker,
    ) -> None:
        """Check detect raises RuntimeError on with verbosity=2."""
        # Mock subprocess.Popen to simulate a failure
        mock_popen = mocker.patch("subprocess.Popen")
        mock_popen.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="tstorms_driver", stderr="Error occurred"
        )

        tracker = tstorms_tracker[0]

        # Assert that detect raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Detect failed with a non-zero exit code"
        ):
            tracker.detect(verbosity=2)


class TestTSTORMSTrackerStitch:
    """Tests for the stitch functionality of TSTORMSTracker."""

    def test_tstorms_tracker_stitch(self, mocker, tstorms_tracker) -> None:
        """Checks the correct trajectory_analysis call is made."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )
        # Mock warning for no output file being generated.
        mocker.patch("warnings.warn")

        # Set up tracker
        tracker = tstorms_tracker[0]
        tstorms_dir = tstorms_tracker[1]
        output_dir = tracker.tstorms_parameters.output_dir

        # Create a dummy cyclones file
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        result = tracker.stitch()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [f"{tstorms_dir}/trajectory_analysis/trajectory_analysis_csc.exe"],
            stdin=mocker.ANY,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(output_dir),
            input=None,
        )
        stdin_file = mock_subprocess_run.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_traj"
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_stitch_verbose(
        self, mocker, tstorms_tracker, capsys
    ) -> None:
        """Checks the correct trajectory_analysis call is made with verbosity=2."""
        # Mock subprocess.Popen to simulate real-time output
        mock_popen = mocker.patch("subprocess.Popen")
        mock_stdout = mocker.MagicMock()
        mock_stdout.readline.side_effect = ["Line 1\n", "Line 2\n", ""]
        mock_process = mocker.MagicMock(
            stdout=mock_stdout,
            stderr="Mocked stderr output",
            returncode=0,
        )
        # Ensure communicate returns a tuple with stdout and stderr
        mock_process.communicate.return_value = ("", "Mocked stderr output")
        mock_popen.return_value = mock_process
        # Mock warning for no output file being generated.
        mocker.patch("warnings.warn")

        # Set up tracker
        tracker = tstorms_tracker[0]
        tstorms_dir = tstorms_tracker[1]
        output_dir = tracker.tstorms_parameters.output_dir

        # Create a dummy cyclones file
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        result = tracker.stitch(verbosity=2)

        # Check subprocess.Popen call
        mock_popen.assert_called_once_with(
            [f"{tstorms_dir}/trajectory_analysis/trajectory_analysis_csc.exe"],
            stdin=mocker.ANY,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False,
            bufsize=1,
            cwd=str(output_dir),
        )

        # Verify stdin file path
        stdin_file = mock_popen.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_traj"

        # Verify returned outputs - using capsys to get the stderr feed
        captured = capsys.readouterr()
        assert captured.out == "Executing Stitch...\nLine 1\nLine 2\n"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_stitch_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check detect raises FileNotFoundError when executable is missing."""
        # Mock subprocess.run to simulate a FileNotFoundError
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TSTORMSTracker instance
        tracker = tstorms_tracker[0]

        # Create a dummy cyclones file
        output_dir = tracker.tstorms_parameters.output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        # Assert that stitch raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Stitch failed because the executable could not be found",
        ):
            tracker.stitch()

    def test_tstorms_tracker_stitch_verbose_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check stitch raises FileNotFoundError verbosity=2."""
        # Mock subprocess.Popen to simulate a FileNotFoundError
        mock_popen = mocker.patch("subprocess.Popen")
        mock_popen.side_effect = FileNotFoundError("Executable not found")

        tracker = tstorms_tracker[0]

        # Create a dummy cyclones file
        output_dir = tracker.tstorms_parameters.output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        # Assert that detect raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Stitch failed because the executable could not be found",
        ):
            tracker.stitch(verbosity=2)

    def test_tstorms_tracker_stitch_failure(self, mocker, tstorms_tracker) -> None:
        """Check stitch raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="trajectory_analysis", stderr="Error occurred"
        )

        # Create a TSTORMSTracker instance
        tracker = tstorms_tracker[0]

        # Create a dummy cyclones file
        output_dir = tracker.tstorms_parameters.output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        # Assert that stitch raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Stitch failed with a non-zero exit code"
        ):
            tracker.stitch()

    def test_tstorms_tracker_stitch_verbose_failure(
        self,
        mocker,
        tstorms_tracker,
    ) -> None:
        """Check stitch raises RuntimeError on with verbosity=2."""
        # Mock subprocess.Popen to simulate a failure
        mock_popen = mocker.patch("subprocess.Popen")
        mock_popen.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="trajectory_analysis", stderr="Error occurred"
        )

        tracker = tstorms_tracker[0]

        # Create a dummy cyclones file
        output_dir = tracker.tstorms_parameters.output_dir
        generated_cyclones_file = os.path.join(output_dir, "cyclones")
        with open(generated_cyclones_file, "w") as f:
            f.write("Mock cyclones content")

        # Assert that detect raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Stitch failed with a non-zero exit code"
        ):
            tracker.stitch(verbosity=2)
