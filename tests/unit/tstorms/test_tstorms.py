"""
Unit tests for the TSTORMS Python bindings.

Note that these do not require a TSTORMS installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import os

import pytest

from tctrack.tstorms import (
    TSTORMSBaseParameters,
    TSTORMSDetectParameters,
    TSTORMSStitchParameters,
    TSTORMSTracker,
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


class TestTSTORMSTypes:
    """Tests for the different Classes and Types defined for TSTORMS."""

    def test_base_parameters_initialization(self):
        """Check TSTORMSBaseParameters initializes correctly with valid arguments."""
        params = TSTORMSBaseParameters(
            tstorms_dir="/path/to/tstorms", output_dir="/path/to/output"
        )
        assert params.tstorms_dir == "/path/to/tstorms"
        assert params.input_dir is None
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


class TestTSTORMSTracker:
    """Tests for the TSTORMS Tracker class."""

    def test_write_driver_namelist(self, tmp_path, tstorms_filenames):
        """Test the generation of the driver namelist inside tstorms_driver/."""
        # Create a tempdir for tstorms_dir and tstorms_driver (assumed to exist)
        tstorms_dir = tmp_path / "tstorms"
        tstorms_driver_dir = tstorms_dir / "tstorms_driver"
        tstorms_driver_dir.mkdir(parents=True)

        # Create mock parameters
        tstorms_params = TSTORMSBaseParameters(
            tstorms_dir=str(tstorms_dir),
            output_dir=str(tmp_path / "output"),
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
        tracker = TSTORMSTracker(tstorms_params, detect_params)

        # Call the method
        namelist_path = tracker._write_driver_namelist()  # noqa: SLF001 - Private member access

        # Verify the file was created inside tstorms_driver/
        assert os.path.exists(namelist_path)
        assert namelist_path == str(tstorms_driver_dir / "nml_driver")

        # Verify the content of the file
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
            assert "fn_u    = 'u.nc'" in content
            assert "fn_v    = 'v.nc'" in content
            assert "fn_vort = 'vort.nc'" in content
            assert "fn_tm   = 'tm.nc'" in content
            assert "fn_slp  = 'slp.nc'" in content
            assert "use_sfc_wnd = .true." in content

    def test_write_driver_namelist_missing_dir(self, tmp_path, tstorms_filenames):
        """Test that FileNotFoundError is raised when tstorms_driver/ does not exist."""
        # Create a temporary directory for tstorms_dir but not tstorms_driver
        tstorms_dir = tmp_path / "tstorms"
        tstorms_dir.mkdir()

        tstorms_params = TSTORMSBaseParameters(
            tstorms_dir=str(tstorms_dir),
            output_dir=str(tmp_path / "output"),
        )
        detect_params = TSTORMSDetectParameters(**tstorms_filenames)
        tracker = TSTORMSTracker(tstorms_params, detect_params)

        with pytest.raises(
            FileNotFoundError, match=r"TSTORMS driver directory .* does not exist"
        ):
            tracker._write_driver_namelist()  # noqa: SLF001 - Private member access
