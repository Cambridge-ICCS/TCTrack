"""
Unit tests for the TRACK Python bindings.

Note that these do not require a TRACK installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

from unittest import mock

import cftime
import pytest
from netCDF4 import Dataset
from numpy.testing import assert_array_equal

from tctrack.track import TRACKParameters, TRACKTracker


class TestTrackParameters:
    """Tests for the TRACKParameters dataclass."""

    def test_parameters_defaults(self):
        """Check TRACKParameters default values."""
        params = TRACKParameters(base_dir="dir", input_file="input")
        assert params.base_dir == "dir"
        assert params.input_file == "input"
        assert params.output_file is None
        assert params.filter_distance is None
        assert params.wind_var_names == ("ua", "va")
        assert params.caltype == "standard"
        assert params.binary == "bin/track.run"
        assert params.file_extension == "track_out"
        assert params.vorticity_file == "vor.dat"
        assert params.filt_vorticity_file == "vor_T63.dat"


class TestTrackTracker:
    """Tests for the TRACKTracker class and methods."""

    nx = 20
    ny = 10

    initialisation_inputs = (
        *["0", "y", "n", "y", "g", "y", "1", "0", "b", "2"],
        *["90", "0", "1", "193", "49", "96", "0", "150", "150"],
        *["-90", "0", "1", "193", "1", "48", "0", "150", "150"],
    )

    def _setup_tracker(
        self, mocker, clear_copy=True, mock_input_file=True, params_dict=None
    ) -> TRACKTracker:
        """Create the TRACKTracker object and mock the necessary functions."""
        # Mock shutil.copy, subprocess.run
        self.mock_copy = mocker.patch("shutil.copy")
        self.mock_subprocess_run = mocker.patch("subprocess.run")
        self.mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Mock the functions checking the input file
        if mock_input_file is True:
            mocker.patch(
                "tctrack.track.track.lat_lon_sizes", return_value=(self.ny, self.nx)
            )
            mocker.patch("tctrack.track.track.TRACKTracker._check_input_file")

        # Create the tracker
        if params_dict is None:
            params_dict = {}
        params_dict.setdefault("base_dir", "dir")
        params_dict.setdefault("input_file", "input")
        params = TRACKParameters(**params_dict)
        tracker = TRACKTracker(params)

        if clear_copy:
            self.mock_copy.reset_mock()
        return tracker

    def test_initialisation(self, mocker):
        """Check TRACKTracker initialises correctly."""
        tracker = self._setup_tracker(mocker, clear_copy=False)

        # Check attributes
        assert tracker._nx == self.nx  # noqa: SLF001
        assert tracker._ny == self.ny  # noqa: SLF001

        # Check calls to shutil.copy are as expected
        copy_calls = [
            mock.call("dir/data/zone.dat", "dir/data/zone.dat0"),
            mock.call("dir/data/adapt.dat", "dir/data/adapt.dat0"),
        ]
        self.mock_copy.assert_has_calls(copy_calls)

    def test_initialisation_failure(self, mocker):
        """Check TRACKTracker initialisation fails as expected."""
        with pytest.raises(FileNotFoundError, match="Input file does not exist"):
            self._setup_tracker(mocker, mock_input_file=False)

    def test_calculate_vorticity(self, mocker):
        """Check calculate_vorticity calls TRACK with the expected inputs."""
        tracker = self._setup_tracker(mocker)

        # Run the method
        tracker.calculate_vorticity()

        # Expected inputs
        input_file = tracker.parameters.input_file
        input_params = [
            *["n", "0", "4", "n", "1", "y", "ua", "n", "85000", "n", "g"],
            *["n", "1", str(self.nx), "1", str(self.ny), "y", "12", "0", "2", "y"],
            *["ua", "85000", "va", "85000", "1", "1", "10000", "1", "1"],
            *["10000", "y", "dir/indat/vor.dat", "0", "10", "0", "0", "0"],
        ]

        # Check call to subprocess.run is as expected
        self.mock_subprocess_run.assert_called_once_with(
            ["dir/bin/track.run", "-i", input_file, "-f", "track_out"],
            input="\n".join(input_params) + "\n",
            check=True,
            capture_output=True,
            text=True,
        )

        # Check call to shutil.copy is as expected
        self.mock_copy.assert_called_once_with(input_file, "dir/indat/")

    def test_spectral_filtering(self, mocker):
        """Check spectral_filtering calls TRACK with the expected inputs."""
        tracker = self._setup_tracker(mocker)

        # Run the method
        tracker.spectral_filtering()

        # Expected inputs
        input_file = tracker.parameters.vorticity_file
        input_params = [
            *["n", "0", "0", "y", "n", "n", "g", "n"],
            *["1", str(self.nx), "1", str(self.ny), "y"],
            *["4", "1", "1", "10000", "1", "n", "63", "y", "1", "63", "2", "y", "0.1"],
            *["0", "5", "63", "n"],
        ]

        # Check call to subprocess.run is as expected
        self.mock_subprocess_run.assert_called_once_with(
            ["dir/bin/track.run", "-i", input_file, "-f", "track_out"],
            input="\n".join(input_params) + "\n",
            check=True,
            capture_output=True,
            text=True,
        )

        # Check call to shutil.copy is as expected
        self.mock_copy.assert_called_once_with(
            "dir/outdat/specfil.track_out_band001",
            "dir/indat/" + tracker.parameters.filt_vorticity_file,
        )

    def test_tracking(self, mocker):
        """Check tracking calls TRACK with the expected inputs."""
        tracker = self._setup_tracker(mocker)

        # Run the method
        tracker.tracking()

        # Expected inputs
        input_file = tracker.parameters.filt_vorticity_file
        hemisphere_inputs = [
            *["y", "1.0e+5", "n", "1.0", "1", "n", "1", "1", "1", "100000", "e", "2"],
            *["7", "n", "0", "n", "n", "n", "n", "0.", "y", "1", "0", "0", "y", "d"],
            *["n", "n", "0"],
        ]
        input_params = [
            *["n", "0"],
            *self.initialisation_inputs,
            *["n", "n", "n"],
            *hemisphere_inputs,
            *hemisphere_inputs,
            *["n", "n", "0.2", "0.8", "n", "y", "y", "6.5", "1", "n", "n", "n", "y"],
            *["n", "y", "y", "0", "0", "n"],
        ]

        # Check call to subprocess.run is as expected
        self.mock_subprocess_run.assert_called_once_with(
            ["dir/bin/track.run", "-i", input_file, "-f", "track_out"],
            input="\n".join(input_params) + "\n",
            check=True,
            capture_output=True,
            text=True,
        )

    def test_filter_trajectories(self, mocker):
        """Check filter_trajectories calls TRACK with the expected inputs."""
        tracker = self._setup_tracker(mocker)

        # Run the method
        tracker.filter_trajectories()

        # Expected inputs
        input_file = tracker.parameters.filt_vorticity_file
        input_params = [
            *["n", "0"],
            *self.initialisation_inputs,
            *["y", "1", "n", "1", "1", "dir/outdat/objout.new.track_out"],
            *["dir/outdat/tdump.track_out", "1", "n", "8", "1000000", "n", "n"],
            *["n", "n", "a", "n", "n", "n", "0"],
        ]

        # Check call to subprocess.run is as expected
        self.mock_subprocess_run.assert_called_once_with(
            ["dir/bin/track.run", "-i", input_file, "-f", "track_out"],
            input="\n".join(input_params) + "\n",
            check=True,
            capture_output=True,
            text=True,
        )

    def _create_nc_output_file(self, tmp_path):
        """Create a minimal TRACK NetCDF output file with 2 tracks."""
        directory = tmp_path / "outdat"
        directory.mkdir()
        file_path = directory / "ff_trs.track_out.nc"

        with Dataset(file_path, "w", format="NETCDF4") as nc:
            # Dimensions
            nc.createDimension("tracks", 2)
            nc.createDimension("record", None)  # unlimited

            # Variables
            track_id = nc.createVariable("TRACK_ID", "i4", ("tracks",))
            first_pt = nc.createVariable("FIRST_PT", "i4", ("tracks",))
            num_pts = nc.createVariable("NUM_PTS", "i4", ("tracks",))
            index = nc.createVariable("index", "i4", ("record",))
            time = nc.createVariable("time", "i4", ("record",))
            lon = nc.createVariable("longitude", "f4", ("record",))
            lat = nc.createVariable("latitude", "f4", ("record",))
            intensity = nc.createVariable("intensity", "f4", ("record",))

            # Data
            track_id[:] = [0, 1]
            first_pt[:] = [0, 2]
            num_pts[:] = [2, 3]

            index[:] = [0, 1, 0, 1, 2]
            time[:] = [1, 2, 1, 2, 3]  # 1-indexed
            lon[:] = [0, 1, 2, 3, 4]
            lat[:] = [0, -1, -2, -3, -4]
            intensity[:] = [0, 10, 20, 30, 40]

    def _create_nc_input_file(self, tmp_path):
        """Create a minimal NetCDF input file with defined times."""
        file_path = tmp_path / "input.nc"

        with Dataset(file_path, "w", format="NETCDF4") as nc:
            # Dimensions
            nc.createDimension("lon", self.nx)
            nc.createDimension("lat", self.ny)
            nc.createDimension("time", None)  # unlimited

            # Time data
            time = nc.createVariable("time", "f4", ("time",))
            time[:] = [0, 0.5, 1, 1.5, 2]
            time.standard_name = "time"
            time.long_name = "time"
            time.units = "days since 1950-01-01"
            time.calendar = "360_day"
            time.axis = "T"

        return str(file_path)

    def test_trajectories(self, mocker, tmp_path):
        """Check trajectories correctly reads in the output netcdf file."""
        # Create the temporary input and output files
        self._create_nc_output_file(tmp_path)
        input_file = self._create_nc_input_file(tmp_path)

        # Set up the tracker pointing to the created files
        params = {"base_dir": str(tmp_path), "input_file": input_file}
        tracker = self._setup_tracker(mocker, mock_input_file=False, params_dict=params)

        # Get the trajectories
        trajectories = tracker.trajectories()

        # Check the trajectories
        t1 = trajectories[0]
        assert t1.trajectory_id == 0
        assert t1.observations == 2
        assert t1.calendar == "360_day"
        assert_array_equal(t1.data["lon"], [0, 1])
        assert_array_equal(t1.data["lat"], [0, -1])
        assert_array_equal(t1.data["intensity"], [0, 10])
        times1 = cftime.num2date(
            [0, 0.5], units="days since 1950-01-01", calendar="360_day"
        )
        assert_array_equal(t1.data["timestamp"], times1)

        t2 = trajectories[1]
        assert t2.trajectory_id == 1
        assert t2.observations == 3
        assert t2.calendar == "360_day"
        assert_array_equal(t2.data["lon"], [2, 3, 4])
        assert_array_equal(t2.data["lat"], [-2, -3, -4])
        assert_array_equal(t2.data["intensity"], [20, 30, 40])
        times2 = cftime.num2date(
            [0, 0.5, 1], units="days since 1950-01-01", calendar="360_day"
        )
        assert_array_equal(t2.data["timestamp"], times2)

    def test_trajectories_failure(self, mocker, tmp_path):
        """Check trajectories fails as expected."""
        # No output file
        tracker = self._setup_tracker(mocker)
        mocker.stopall()
        with pytest.raises(
            FileNotFoundError, match="TRACK output trajectory file does not exist"
        ):
            tracker.trajectories()

        # No input file (create the output file)
        self._create_nc_output_file(tmp_path)
        params = {"base_dir": str(tmp_path)}
        tracker = self._setup_tracker(mocker, params_dict=params)
        mocker.stopall()
        with pytest.raises(FileNotFoundError, match="Input file does not exist"):
            tracker.trajectories()
