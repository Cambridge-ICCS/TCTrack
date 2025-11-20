"""
Unit tests for the TSTORMS Python bindings.

Note that these do not require a TSTORMS installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import importlib.metadata
import json
import os
import subprocess
import warnings
from dataclasses import asdict
from typing import Tuple

import cf
import numpy as np
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


@pytest.fixture
def tstorms_tracker(tmp_path, tstorms_filenames) -> Tuple[TSTORMSTracker, str]:
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
    return TSTORMSTracker(tstorms_parameters, detect_params), tstorms_dir


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


class TestTSTORMSTracker:
    """Tests for the TSTORMS Tracker class."""

    def test_tstorms_tracker_global_metadata(self, tstorms_tracker) -> None:
        """Check set_metadata correctly sets _global_metadata."""
        tracker = tstorms_tracker[0]

        tracker.set_metadata()

        assert tracker.global_metadata is not None
        assert tracker.global_metadata == {
            "tctrack_version": importlib.metadata.version("tctrack"),
            "tctrack_tracker": "TSTORMSTracker",
            "tstorms_parameters": json.dumps(asdict(tracker.tstorms_parameters)),
            "detect_parameters": json.dumps(asdict(tracker.detect_parameters)),
            "stitch_parameters": json.dumps(asdict(tracker.stitch_parameters)),
        }

    def test_write_driver_namelist(self, tstorms_tracker):
        """Test the generation of the driver namelist."""
        tracker = tstorms_tracker[0]

        # Call the method
        namelist_path = tracker._write_driver_namelist()  # noqa: SLF001 - Private member access

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
        namelist_path = tracker._write_trajectory_analysis_namelist()  # noqa: SLF001 - Private member access

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

    def test_extract_calendar_metadata_valid(self, tstorms_tracker):
        """Test extracting valid calendar metadata from a NetCDF file."""
        tracker = tstorms_tracker[0]

        # Create a temporary NetCDF file as expected by the tracker fixture
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
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
        cf.write(u_field, netcdf_file)

        # Call the method
        tracker._extract_calendar_metadata()  # noqa: SLF001 - Private member access

        # Assert the metadata was extracted correctly
        assert tracker._calendar_metadata == {  # noqa: SLF001 - Private member access
            "calendar_type": "360_day",
            "units": "days since 1950-01-01",
        }

    def test_extract_calendar_metadata_file_not_found(self, tstorms_tracker):
        """Test behavior when the NetCDF file is not found."""
        tracker = tstorms_tracker[0]

        # Do not create the file expected by the tracker fixture (u.nc)

        with pytest.warns(
            UserWarning, match="No input file for u_ref found to set calendar metadata"
        ):
            tracker._extract_calendar_metadata()  # noqa: SLF001 - Private member access
        assert tracker._calendar_metadata == {"calendar_type": "julian", "units": None}  # noqa: SLF001 - Private member access

    def test_extract_calendar_metadata_no_unlimited_dim(self, tstorms_tracker):
        """Test behavior when no unlimited dimension is found in the NetCDF file."""
        tracker = tstorms_tracker[0]

        # Create a temporary NetCDF file as expected by the tracker fixture, time not ul
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
        with pytest.warns(
            UserWarning, match="No input file for u_ref found to set calendar metadata"
        ):
            tracker._extract_calendar_metadata()  # noqa: SLF001 - Private member access
        assert tracker._calendar_metadata == {"calendar_type": "julian", "units": None}  # noqa: SLF001 - Private member access

    def test_extract_calendar_metadata_no_coord_var(self, tstorms_tracker):
        """Test behavior when missing unlimited dimension coordinate data."""
        tracker = tstorms_tracker[0]

        # Create a temporary NetCDF file as expected by the tracker fixture
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
        with pytest.warns(
            UserWarning, match="No input file for u_ref found to set calendar metadata"
        ):
            tracker._extract_calendar_metadata()  # noqa: SLF001 - Private member access
        assert tracker._calendar_metadata == {"calendar_type": "julian", "units": None}  # noqa: SLF001 - Private member access

    def test_extract_calendar_metadata_missing_attributes(self, tstorms_tracker):
        """Test behavior when missing 'units' or 'calendar' attributes from coord."""
        tracker = tstorms_tracker[0]

        # Create a temporary NetCDF file as expected by the tracker fixture
        netcdf_file = os.path.join(tracker.tstorms_parameters.input_dir, "u.nc")
        u_field = cf.Field(properties={"standard_name": "velocity"})
        u_field.nc_set_variable("u_ref")
        domain_axis_time = cf.DomainAxis(10)
        domain_axis_time.nc_set_unlimited(True)
        _ = u_field.set_construct(domain_axis_time)
        u_data = cf.Data(np.arange(10.0))
        u_field.set_data(u_data)
        dimension_T = cf.DimensionCoordinate(
            properties={
                "standard_name": "time",
            },
            data=cf.Data(np.arange(10.0)),
        )
        u_field.set_construct(dimension_T)
        cf.write(u_field, netcdf_file)

        # Ensure that code runs with no warnings for missing files raised
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            # Call the method and assert the default metadata is set
            tracker._extract_calendar_metadata()  # noqa: SLF001 - Private member access
            assert tracker._calendar_metadata == {  # noqa: SLF001 - Private member access
                "calendar_type": "julian",
                "units": None,
            }


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
        )
        stdin_file = mock_subprocess_run.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_driver"
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_detect_verbose(
        self, mocker, tstorms_tracker, capsys
    ) -> None:
        """Checks the correct tstorms_driver call is made with verbose=True."""
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
        result = tracker.detect(verbose=True)

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
        assert captured.out == "Line 1\nLine 2\n"
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

        # Assert that detect_nodes raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Detect failed because the executable could not be found",
        ):
            tracker.detect()

    def test_tstorms_tracker_detect_verbose_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check detect raises FileNotFoundError verbose=True."""
        # Mock subprocess.Popen to simulate a FileNotFoundError
        mock_popen = mocker.patch("subprocess.Popen")
        mock_popen.side_effect = FileNotFoundError("Executable not found")

        tracker = tstorms_tracker[0]

        # Assert that detect raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Detect failed because the executable could not be found",
        ):
            tracker.detect(verbose=True)

    def test_tstorms_tracker_detect_failure(self, mocker, tstorms_tracker) -> None:
        """Check detect_nodes raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="tstorms_driver", stderr="Error occurred"
        )

        # Create a TSTORMSTracker instance
        tracker = tstorms_tracker[0]

        # Assert that detect_nodes raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Detect failed with a non-zero exit code"
        ):
            tracker.detect()

    def test_tstorms_tracker_detect_verbose_failure(
        self,
        mocker,
        tstorms_tracker,
    ) -> None:
        """Check detect raises RuntimeError on with verbose=True."""
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
            tracker.detect(verbose=True)


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
        )
        stdin_file = mock_subprocess_run.call_args[1]["stdin"]
        assert stdin_file.name == f"{output_dir}/nml_traj"
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_tstorms_tracker_stitch_verbose(
        self, mocker, tstorms_tracker, capsys
    ) -> None:
        """Checks the correct trajectory_analysis call is made with verbose=True."""
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

        result = tracker.stitch(verbose=True)

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
        assert captured.out == "Line 1\nLine 2\n"
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

        # Assert that detect_nodes raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="Stitch failed because the executable could not be found",
        ):
            tracker.stitch()

    def test_tstorms_tracker_stitch_verbose_file_not_found(
        self, mocker, tstorms_tracker
    ) -> None:
        """Check stitch raises FileNotFoundError verbose=True."""
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
            tracker.stitch(verbose=True)

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

        # Assert that detect_nodes raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Stitch failed with a non-zero exit code"
        ):
            tracker.stitch()

    def test_tstorms_tracker_stitch_verbose_failure(
        self,
        mocker,
        tstorms_tracker,
    ) -> None:
        """Check stitch raises RuntimeError on with verbose=True."""
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
            tracker.stitch(verbose=True)
