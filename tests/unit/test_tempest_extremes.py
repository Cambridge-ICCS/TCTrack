"""
Unit tests for the Tempest Extremes Python bindings.

Note that these do not require a Tempest Extremes installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import re
import subprocess

import pytest
from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap

from tctrack.tempest_extremes import (
    DetectNodesParameters,
    StitchNodesParameters,
    TEContour,
    TEOutputCommand,
    TEThreshold,
    TETracker,
    Track,
)


class TestTrack:
    """Tests for the Track class."""

    def test_track_init(self) -> None:
        """Test the Track class initialization."""
        track = Track(track_id=1, observations=0, year=1950, month=1, day=1, hour=3)
        assert track.track_id == 1
        assert track.observations == 0
        assert track.start_time.year == 1950
        assert track.start_time.month == 1
        assert track.start_time.day == 1
        assert track.start_time.hour == 3
        assert track.calendar == "gregorian"

    @pytest.mark.parametrize(
        "calendar, expected_type",
        [
            ("gregorian", DatetimeGregorian),
            ("360_day", Datetime360Day),
            ("noleap", DatetimeNoLeap),
        ],
    )
    def test_track_caltypes(self, calendar, expected_type) -> None:
        """Test the Track class initialization with different calendar types."""
        track = Track(
            track_id=1,
            observations=0,
            year=1950,
            month=1,
            day=1,
            hour=3,
            calendar=calendar,
        )
        assert track.track_id == 1
        assert track.observations == 0
        assert track.start_time.year == 1950
        assert track.start_time.month == 1
        assert track.start_time.day == 1
        assert track.start_time.hour == 3
        assert track.calendar == calendar
        assert isinstance(track.start_time, expected_type)

    def test_track_invalid_calendar(self) -> None:
        """Test that Track raises a ValueError for an invalid calendar type."""
        with pytest.raises(
            ValueError, match="Unsupported calendar type: invalid_calendar"
        ):
            Track(
                track_id=1,
                observations=0,
                year=1950,
                month=1,
                day=1,
                hour=3,
                calendar="invalid_calendar",
            )

    def test_track_add_point(self) -> None:
        """Test the Track class add_point method."""
        track = Track(track_id=1, observations=0, year=1950, month=1, day=1, hour=3)
        variables_dict = {"grid_i": 164, "grid_j": 332, "psl": 1.005377e05}
        track.add_point(1950, 1, 1, 3, variables_dict)
        assert len(track.data) == 4
        assert track.data["grid_i"] == [164]
        assert track.data["grid_j"] == [332]
        assert track.data["psl"] == [1.005377e05]

    @pytest.mark.parametrize(
        "inputs, expected_error, match",
        [
            # Missing data in the inputs
            (
                (1950, 1, 1, {"psl": 1.005377e05}),
                TypeError,
                re.escape(
                    "Track.add_point() missing 1 required positional argument: "
                    "'variables'"
                ),
            ),
            # Missing data in the inputs
            (
                (1950, 1, 1, 3, "extra", {"psl": 1.005377e05}),
                TypeError,
                re.escape(
                    "Track.add_point() takes 6 positional arguments but 7 were given"
                ),
            ),
            # Invalid variable data
            (
                (1950, 1, 1, 3, {"psl": "invalid"}),
                ValueError,
                "Invalid variable data: {'psl': 'invalid'}",
            ),
            # Invalid date/time values
            ((1950, 13, 1, 3, {"psl": 1.005377e05}), ValueError, None),
            ((1950, 1, 32, 3, {"psl": 1.005377e05}), ValueError, None),
            ((1950, 1, 1, 25, {"psl": 1.005377e05}), ValueError, None),
        ],
    )
    def test_add_point_edge_cases(self, inputs, expected_error, match):
        """Test the Track class add_point edge cases."""
        track = Track(track_id=1, observations=0, year=1950, month=1, day=1, hour=3)
        with pytest.raises(expected_error, match=match):
            track.add_point(*inputs)


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

    @pytest.mark.parametrize(
        "parameter",
        [
            {"out_file_format": "gfdl"},
            {"out_file_format": "csv"},
            {"out_file_format": "csvnohead"},
            {"caltype": "standard"},
            {"caltype": "noleap"},
            {"caltype": "360_day"},
        ],
    )
    def test_stitch_nodes_parameters_valid(self, parameter) -> None:
        """Check valid StitchNodesParameters inputs do not raise errors."""
        StitchNodesParameters(**parameter)

    @pytest.mark.parametrize(
        "parameter,msg",
        [
            (
                {"out_file_format": "invalid"},
                r"Invalid out_file_format \(invalid\). "
                + "Allowed values are 'gfdl', 'csv', or 'csvnohead'",
            ),
            (
                {"caltype": "invalid"},
                r"Invalid caltype \(invalid\). "
                + "Allowed values are 'standard', 'noleap', or '360_day'",
            ),
        ],
    )
    def test_stitch_nodes_parameters_invalid(self, parameter, msg) -> None:
        """Check invalid StitchNodesParameters inputs correctly raise errors."""
        with pytest.raises(ValueError, match=msg):
            StitchNodesParameters(**parameter)

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

    @pytest.fixture
    def mock_gfdl_file(self, tmp_path):
        """Fixture to create a mock GFDL file with two tracks."""
        file_path = tmp_path / "tracks_out_gfdl.txt"
        file_path.write_text(
            "start\t2\t1950\t1\t1\t3\n"
            "\t164\t332\t57.832031\t-12.070312\t1.005377e+05\t0.000000e+00\t1950\t1\t1\t3\n"
            "\t163\t332\t57.480469\t-12.070312\t1.005820e+05\t0.000000e+00\t1950\t1\t1\t6\n"
            "start\t2\t1950\t1\t2\t0\n"
            "\t843\t275\t296.542969\t-25.429688\t9.970388e+04\t2.633214e+02\t1950\t1\t2\t0\n"
            "\t850\t266\t299.003906\t-27.539062\t9.989988e+04\t6.951086e+01\t1950\t1\t2\t6\n"
        )
        return str(file_path)

    @pytest.fixture
    def mock_csv_file(self, tmp_path):
        """Fixture to create a mock CSV file with two tracks."""
        file_path = tmp_path / "tracks_out_csv.txt"
        file_path.write_text(
            "track_id,year,month,day,hour,i,j,lon,lat,psl,orog\n"
            "0,1950,1,1,3,164,332,57.832031,-12.070312,1.005377e+05,0.000000e+00\n"
            "0,1950,1,1,6,163,332,57.480469,-12.070312,1.005820e+05,0.000000e+00\n"
            "1,1950,1,2,0,843,275,296.542969,-25.429688,9.970388e+04,2.633214e+02\n"
            "1,1950,1,2,6,850,266,299.003906,-27.539062,9.989988e+04,6.951086e+01\n"
        )
        return str(file_path)

    @pytest.fixture
    def mock_csvnohead_file(self, tmp_path):
        """Fixture to create a mock CSV file without a header and with two tracks."""
        file_path = tmp_path / "tracks_out_csvnohead.txt"
        file_path.write_text(
            "0,1950,1,1,3,164,332,57.832031,-12.070312,1.005377e+05,0.000000e+00\n"
            "0,1950,1,1,6,163,332,57.480469,-12.070312,1.005820e+05,0.000000e+00\n"
            "1,1950,1,2,0,843,275,296.542969,-25.429688,9.970388e+04,2.633214e+02\n"
            "1,1950,1,2,6,850,266,299.003906,-27.539062,9.989988e+04,6.951086e+01\n"
        )
        return str(file_path)

    @pytest.mark.parametrize(
        "file_format, mock_file_fixture, expected_tracks",
        [
            (
                "gfdl",
                "mock_gfdl_file",
                [
                    {
                        "track_id": 1,
                        "observations": 2,
                        "data": {
                            "grid_i": [164, 163],
                            "grid_j": [332, 332],
                            "var_1": [57.832031, 57.480469],
                            "var_2": [-12.070312, -12.070312],
                            "var_3": [1.005377e05, 1.005820e05],
                            "var_4": [0.0, 0.0],
                        },
                    },
                    {
                        "track_id": 2,
                        "observations": 2,
                        "data": {
                            "grid_i": [843, 850],
                            "grid_j": [275, 266],
                            "var_1": [296.542969, 299.003906],
                            "var_2": [-25.429688, -27.539062],
                            "var_3": [9.970388e04, 9.989988e04],
                            "var_4": [2.633214e02, 6.951086e01],
                        },
                    },
                ],
            ),
            (
                "csv",
                "mock_csv_file",
                [
                    {
                        "track_id": 0,
                        "observations": 2,
                        "data": {
                            "grid_i": [164, 163],
                            "grid_j": [332, 332],
                            "lon": [57.832031, 57.480469],
                            "lat": [-12.070312, -12.070312],
                            "psl": [1.005377e05, 1.005820e05],
                            "orog": [0.0, 0.0],
                        },
                    },
                    {
                        "track_id": 1,
                        "observations": 2,
                        "data": {
                            "grid_i": [843, 850],
                            "grid_j": [275, 266],
                            "lon": [296.542969, 299.003906],
                            "lat": [-25.429688, -27.539062],
                            "psl": [9.970388e04, 9.989988e04],
                            "orog": [2.633214e02, 6.951086e01],
                        },
                    },
                ],
            ),
            (
                "csvnohead",
                "mock_csvnohead_file",
                [
                    {
                        "track_id": 0,
                        "observations": 2,
                        "data": {
                            "grid_i": [164, 163],
                            "grid_j": [332, 332],
                            "var_1": [57.832031, 57.480469],
                            "var_2": [-12.070312, -12.070312],
                            "var_3": [1.005377e05, 1.005820e05],
                            "var_4": [0.0, 0.0],
                        },
                    },
                    {
                        "track_id": 1,
                        "observations": 2,
                        "data": {
                            "grid_i": [843, 850],
                            "grid_j": [275, 266],
                            "var_1": [296.542969, 299.003906],
                            "var_2": [-25.429688, -27.539062],
                            "var_3": [9.970388e04, 9.989988e04],
                            "var_4": [2.633214e02, 6.951086e01],
                        },
                    },
                ],
            ),
        ],
    )
    def test_tracks(self, file_format, mock_file_fixture, expected_tracks, request):
        """Test the tracks() method for different file formats with multiple tracks."""
        mock_file = request.getfixturevalue(mock_file_fixture)
        tracker = TETracker()
        tracker.stitch_nodes_parameters.output_file = mock_file
        tracker.stitch_nodes_parameters.out_file_format = file_format
        tracker.stitch_nodes_parameters.in_fmt = ["var_1", "var_2", "var_3", "var_4"]
        tracks = tracker.tracks()

        # Assertions
        assert len(tracks) == len(expected_tracks)
        for track, expected in zip(tracks, expected_tracks, strict=False):
            assert track.track_id == expected["track_id"]
            assert track.observations == expected["observations"]
            for key, values in expected["data"].items():
                assert track.data[key] == values
