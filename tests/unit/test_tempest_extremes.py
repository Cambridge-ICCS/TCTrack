"""
Unit tests for the Tempest Extremes Python bindings.

Note that these do not require a Tempest Extremes installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import subprocess

import pytest

from tctrack.tempest_extremes import (
    DetectNodesParameters,
    StitchNodesParameters,
    TEContour,
    TEOutputCommand,
    TEThreshold,
    TETracker,
)


class TestTETypes:
    """Tests for the different Classes and Types defined for TempestExtremes."""

    def test_detect_nodes_parameters_defaults(self) -> None:
        """Check the default values for DetectNodesParameters."""
        params = DetectNodesParameters()
        # Check values of all defaults
        assert params.in_data is None
        assert params.out_header is False
        assert params.output_file is None
        assert params.search_by_min is None
        assert params.search_by_max is None
        assert params.closed_contours is None
        assert params.merge_dist == 0.0
        assert params.lat_name == "lat"
        assert params.lon_name == "lon"
        assert params.min_lat == 0.0
        assert params.max_lat == 0.0
        assert params.min_lon == 0.0
        assert params.max_lon == 0.0
        assert params.regional is False
        assert params.output_commands is None

    def test_stitch_nodes_parameters_defaults(self) -> None:
        """Check the default values for StitchNodesParameters."""
        params = StitchNodesParameters()
        # Check values of all defaults
        assert params.output_file == "tracks.txt"
        assert params.in_file is None
        assert params.in_list is None
        assert params.in_fmt is None
        assert params.allow_repeated_times is False
        assert params.caltype == "standard"
        assert params.time_begin is None
        assert params.time_end is None
        assert params.max_sep == 5.0
        assert params.max_gap == 0
        assert params.min_time == 1
        assert params.min_endpoint_dist == 0.0
        assert params.min_path_dist == 0.0
        assert params.threshold_filters is None
        assert params.prioritize is None
        assert params.add_velocity is False
        assert params.out_file_format == "gfdl"
        assert params.out_seconds is False

    def test_te_contour(self) -> None:
        """Test that TEContour creates the appropriate typed dict."""
        contour = TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0)
        # Check attributes are set properly and can be called from object
        assert contour["var"] == "psl"
        assert contour["delta"] == 200.0
        assert contour["dist"] == 5.5
        assert contour["minmaxdist"] == 0.0
        # Check item created is a dict as expected
        assert contour == {"var": "psl", "delta": 200.0, "dist": 5.5, "minmaxdist": 0.0}

    def test_te_output_command(self) -> None:
        """Test that TEOutputCommand creates the appropriate typed dict."""
        output_command = TEOutputCommand(var="psl", operator="min", dist=0.0)
        # Check attributes are set properly and can be called from object
        assert output_command["var"] == "psl"
        assert output_command["operator"] == "min"
        assert output_command["dist"] == 0.0
        # Check item created is a dict as expected
        assert output_command == {"var": "psl", "operator": "min", "dist": 0.0}

    def test_te_threshold(self) -> None:
        """Test that TEThreshold creates the appropriate typed dict."""
        threshold = TEThreshold(var="lat", op="<=", value=40, count=10)
        # Check attributes are set properly and can be called from object
        assert threshold["var"] == "lat"
        assert threshold["op"] == "<="
        assert threshold["value"] == 40
        assert threshold["count"] == 10
        # Check item created is a dict as expected
        assert threshold == {"var": "lat", "op": "<=", "value": 40, "count": 10}


class TestTETracker:
    """Tests for the TETracker class and associated methods."""

    # Test for TETracker initialization
    def test_te_tracker_initialization(self) -> None:
        """Test initialisation of TETracker with some non-default parameters."""
        params = DetectNodesParameters(output_file="output.txt")
        tracker = TETracker(detect_nodes_parameters=params)
        assert tracker.detect_nodes_parameters.output_file == "output.txt"

    def test_te_tracker_detect_nodes_defaults(self, mocker) -> None:
        """Checks the correct detect_nodes call is made for defaults."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # create a TETracker with default parameters and call detect_nodes method
        tracker = TETracker()
        result = tracker.detect_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "DetectNodes",
                "--mergedist",
                "0.0",
                "--latname",
                "lat",
                "--lonname",
                "lon",
                "--minlat",
                "0.0",
                "--maxlat",
                "0.0",
                "--minlon",
                "0.0",
                "--maxlon",
                "0.0",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_te_tracker_detect_nodes_non_defaults(self, mocker) -> None:
        """Checks the correct detect_nodes call is made for non-default values."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker with non-default parameters and call detect_nodes method
        params = DetectNodesParameters(
            merge_dist=10.0,
            lat_name="latitude",
            lon_name="longitude",
            min_lat=-90.0,
            max_lat=90.0,
            min_lon=-180.0,
            max_lon=180.0,
        )
        tracker = TETracker(detect_nodes_parameters=params)
        result = tracker.detect_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "DetectNodes",
                "--mergedist",
                "10.0",
                "--latname",
                "latitude",
                "--lonname",
                "longitude",
                "--minlat",
                "-90.0",
                "--maxlat",
                "90.0",
                "--minlon",
                "-180.0",
                "--maxlon",
                "180.0",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_te_tracker_detect_nodes_optional_parameters(self, mocker) -> None:
        """
        Checks the correct detect_nodes call is made for optional parameters.

        This will implicitly check the lod_to_te utility function.
        """
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker with optional parameters and call detect_nodes method
        params = DetectNodesParameters(
            in_data=["input1.nc", "input2.nc"],
            out_header=True,
            output_file="output.txt",
            search_by_min="psl",
            closed_contours=[
                TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0)
            ],
            output_commands=[
                TEOutputCommand(var="psl", operator="min", dist=0.0),
                TEOutputCommand(var="orog", operator="max", dist=3.0),
            ],
            regional=True,
        )
        tracker = TETracker(detect_nodes_parameters=params)
        result = tracker.detect_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "DetectNodes",
                "--in_data",
                "input1.nc;input2.nc",
                "--out_header",
                "--out",
                "output.txt",
                "--searchbymin",
                "psl",
                "--closedcontourcmd",
                "psl,200.0,5.5,0.0",
                "--mergedist",
                "0.0",
                "--latname",
                "lat",
                "--lonname",
                "lon",
                "--minlat",
                "0.0",
                "--maxlat",
                "0.0",
                "--minlon",
                "0.0",
                "--maxlon",
                "0.0",
                "--regional",
                "--outputcmd",
                "psl,min,0.0;orog,max,3.0",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_te_tracker_detect_nodes_file_not_found(self, mocker) -> None:
        """Check detect_nodes raises FileNotFoundError when executable is missing."""
        # Mock subprocess.run to simulate a FileNotFoundError
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TETracker instance
        params = DetectNodesParameters(output_file="output.txt")
        tracker = TETracker(detect_nodes_parameters=params)

        # Assert that detect_nodes raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="DetectNodes failed because the executable could not be found",
        ):
            tracker.detect_nodes()

    def test_te_tracker_detect_nodes_failure(self, mocker) -> None:
        """Check detect_nodes raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="DetectNodes", stderr="Error occurred"
        )

        # Create a TETracker instance
        params = DetectNodesParameters(output_file="output.txt")
        tracker = TETracker(detect_nodes_parameters=params)

        # Assert that detect_nodes raises RuntimeError
        with pytest.raises(
            RuntimeError, match="DetectNodes failed with a non-zero exit code"
        ):
            tracker.detect_nodes()


class TestTETrackerStitchNodes:
    """Tests for the stitch_nodes functionality of TETracker."""

    def test_stitch_nodes_defaults(self, mocker) -> None:
        """Checks the correct stitch_nodes call is made for default parameters."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with default StitchNodes parameters
        tracker = TETracker()
        result = tracker.stitch_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                "tracks.txt",
                "--caltype",
                "standard",
                "--range",
                "5",
                "--maxgap",
                "0",
                "--mintime",
                "1",
                "--min_endpoint_dist",
                "0",
                "--min_path_dist",
                "0",
                "--out_file_format",
                "gfdl",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_stitch_nodes_non_defaults(self, mocker) -> None:
        """Checks the correct stitch_nodes call is made for non-default parameters."""
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with non-default StitchNodes parameters
        params = StitchNodesParameters(
            output_file="custom_tracks.txt",
            in_file="nodes.txt",
            caltype="noleap",
            max_sep=10.0,
            max_gap=2,
            min_time="48h",
            min_endpoint_dist=5.0,
            min_path_dist=15.0,
            out_file_format="csv",
            out_seconds=True,
        )
        tracker = TETracker(stitch_nodes_parameters=params)
        result = tracker.stitch_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                "custom_tracks.txt",
                "--in",
                "nodes.txt",
                "--caltype",
                "noleap",
                "--range",
                "10.0",
                "--maxgap",
                "2",
                "--mintime",
                "48h",
                "--min_endpoint_dist",
                "5.0",
                "--min_path_dist",
                "15.0",
                "--out_file_format",
                "csv",
                "--out_seconds",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    @pytest.mark.parametrize(
        "threshold_filters,expected_cmd",
        [
            (
                [TEThreshold(var="lat", op="<=", value=40, count=10)],
                ["--threshold", "lat,<=,40,10"],
            ),
            (
                [
                    TEThreshold(var="lat", op="<=", value=40, count=10),
                    TEThreshold(var="lon", op=">=", value=-100, count="all"),
                ],
                ["--threshold", "lat,<=,40,10;lon,>=,-100,all"],
            ),
        ],
    )
    def test_stitch_nodes_threshold_filters(
        self, mocker, threshold_filters, expected_cmd
    ) -> None:
        """
        Checks the correct stitch_nodes call is made with threshold filters.

        Also checks use of string and int inputs to TEThreshold["count"].
        """
        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with threshold filters
        params = StitchNodesParameters(
            output_file="custom_tracks.txt",
            threshold_filters=threshold_filters,
        )
        tracker = TETracker(stitch_nodes_parameters=params)
        result = tracker.stitch_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                "custom_tracks.txt",
                "--caltype",
                "standard",
                "--range",
                "5",
                "--maxgap",
                "0",
                "--mintime",
                "1",
                "--min_endpoint_dist",
                "0",
                "--min_path_dist",
                "0",
                *expected_cmd,
                "--out_file_format",
                "gfdl",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result["stdout"] == "Mocked stdout output"
        assert result["stderr"] == "Mocked stderr output"
        assert result["returncode"] == 0

    def test_stitch_nodes_file_not_found(self, mocker) -> None:
        """Check stitch_nodes raises FileNotFoundError when executable is missing."""
        # Mock subprocess.run to simulate a FileNotFoundError
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TETracker instance
        params = StitchNodesParameters(output_file="custom_tracks.txt")
        tracker = TETracker(stitch_nodes_parameters=params)

        # Assert that stitch_nodes raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="StitchNodes failed because the executable could not be found",
        ):
            tracker.stitch_nodes()

    def test_stitch_nodes_failure(self, mocker) -> None:
        """Check stitch_nodes raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="StitchNodes", stderr="Error occurred"
        )

        # Create a TETracker instance
        params = StitchNodesParameters(output_file="custom_tracks.txt")
        tracker = TETracker(stitch_nodes_parameters=params)

        # Assert that stitch_nodes raises RuntimeError
        with pytest.raises(
            RuntimeError, match="StitchNodes failed with a non-zero exit code"
        ):
            tracker.stitch_nodes()
