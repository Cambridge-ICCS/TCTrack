"""
Unit tests for the TRACK Python bindings.

Note that these do not require a TRACK installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

from unittest import mock

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

    def _setup_tracker(self, mocker, clear_copy=True) -> TRACKTracker:
        """Create the TRACKTracker object and mock the necessary functions."""
        # Mock shutil.copy, subprocess.run, and lat_lon_size
        mocker.patch(
            "tctrack.track.track.lat_lon_sizes", return_value=(self.ny, self.nx)
        )
        self.mock_copy = mocker.patch("shutil.copy")
        self.mock_subprocess_run = mocker.patch("subprocess.run")
        self.mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create the tracker
        params = TRACKParameters(base_dir="dir", input_file="input")
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
