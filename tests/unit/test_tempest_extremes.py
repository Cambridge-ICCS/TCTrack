"""
Unit tests for the Tempest Extremes Python bindings.

Note that these do not require a Tempest Extremes installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

from tctrack.tempest_extremes import (
    DetectNodesParameters,
    TEContour,
    TEOutputCommand,
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


class TestTETracker:
    """Tests for the TETracker class and associated methods."""

    # Test for TETracker initialization
    def test_te_tracker_initialization(self) -> None:
        """Test initialisation of TETracker with some non-default parameters."""
        params = DetectNodesParameters(output_file="output.txt")
        tracker = TETracker(detect_nodes_parameters=params)
        assert tracker.detect_nodes_parameters.output_file == "output.txt"

    def test_te_tracker_detect_nodes_defaults(self, mocker) -> None:
        """Checks the correct detect_nodes call is made out for defaults."""
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
