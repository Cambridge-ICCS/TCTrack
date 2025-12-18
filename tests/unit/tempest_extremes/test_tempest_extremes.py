"""
Unit tests for the Tempest Extremes Python bindings.

Note that these do not require a Tempest Extremes installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import importlib.metadata
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

import cf
import pytest
from cftime import datetime

from tctrack.core import TCTrackerMetadata
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
        params = DetectNodesParameters(in_data=["input_file.nc"])
        # Check values of all defaults
        assert params.out_header is False
        assert params.output_dir == ""
        assert params.output_file == "nodes.txt"
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
        assert params.output_dir == ""
        assert params.output_file == "trajectories.txt"
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
        params = DetectNodesParameters(
            in_data=["input_file.nc"],
            output_dir="outputs",
        )
        tracker = TETracker(detect_nodes_parameters=params)
        assert tracker.detect_nodes_parameters.output_dir == "outputs"

    def test_te_tracker_tempfiles(self) -> None:
        """Test the definition of temporary files when not manually defined."""
        dn_params = DetectNodesParameters(in_data=["input_file.nc"])
        tracker = TETracker(dn_params)
        tempdir = tracker._tempdir.name  # noqa: SLF001
        assert Path(tempdir).exists()
        assert tracker.detect_nodes_parameters.output_dir == tempdir
        assert tracker.stitch_nodes_parameters.in_file == tempdir + "/nodes.txt"
        assert tracker.stitch_nodes_parameters.output_dir == tempdir

        # Check the temporary directory is correctly cleaned up
        del tracker
        assert not Path(tempdir).exists()

    @pytest.fixture
    def netcdf_psl_file(self, tmp_path):
        """Create a netcdf file with metadata given by the 'properties' argument."""

        def _create_file(properties):
            file_name = str(tmp_path / "inputs.nc")
            field = cf.Field(properties=properties)
            field.nc_set_variable("psl")
            domain_axis_time = cf.DomainAxis(10)
            domain_axis_time.nc_set_unlimited(True)
            _ = field.set_construct(domain_axis_time)
            data = cf.Data([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
            field.set_data(data)
            dimension_T = cf.DimensionCoordinate(
                properties={
                    "standard_name": "time",
                    "calendar": "standard",
                    "units": "days since 1950-01-01",
                },
                data=cf.Data([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            )
            field.set_construct(dimension_T)
            cf.write(field, file_name)  # type: ignore[operator]
            return file_name

        return _create_file

    def test_te_tracker_time_metadata(self, netcdf_psl_file):
        """Test extracting valid calendar metadata from a NetCDF file."""
        properties = {
            "standard_name": "air_pressure_at_sea_level",
            "long_name": "Sea Level Pressure",
            "units": "Pa",
        }
        file_name = netcdf_psl_file(properties)
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="psl",
        )
        tracker = TETracker(dn_params)

        # Assert the metadata extracted correctly
        tracker.set_metadata()
        assert tracker.time_metadata == {
            "calendar": "standard",
            "units": "days since 1950-01-01",
            "start_time": datetime(1950, 1, 1, 00, calendar="standard"),
            "end_time": datetime(1950, 1, 10, 00, calendar="standard"),
        }

    def test_te_tracker_time_metadata_variable_not_found(self, netcdf_psl_file):
        """Test ValueError raised when the main search variable not found in inputs."""
        file_name = netcdf_psl_file({})
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="not_psl",
        )
        tracker = TETracker(dn_params)

        with pytest.raises(
            ValueError, match=r"Variable 'not_psl' not found in input files."
        ):
            tracker.set_metadata()

    def test_te_tracker_variable_metadata(self, netcdf_psl_file) -> None:
        """Test the reading in of variable metadata from input NetCDF files."""
        properties = {
            "standard_name": "air_pressure_at_sea_level",
            "long_name": "Sea Level Pressure",
            "units": "Pa",
        }
        file_name = netcdf_psl_file(properties)
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="psl",
            output_commands=[TEOutputCommand(var="psl", operator="min", dist=1)],
        )
        tracker = TETracker(dn_params)
        tracker.set_metadata()
        metadata = tracker.variable_metadata

        # Check the pressure metadata is set correctly
        assert metadata["psl"].properties == properties
        expected_cell_method = cf.CellMethod(
            "area",
            "minimum",
            qualifiers={"comment": "lesser circle of radius 1 degrees"},
        )
        assert metadata["psl"].constructs == [expected_cell_method]
        assert metadata["psl"].construct_kwargs is None

        # Check the grid index metadata is set correctly
        assert metadata["grid_i"] == TCTrackerMetadata(
            {"long_name": "longitudinal grid index"},
        )
        assert metadata["grid_j"] == TCTrackerMetadata(
            {"long_name": "latitudinal grid index"},
        )

    def test_te_tracker_variable_metadata_failure(self, netcdf_psl_file) -> None:
        """Check set_metadata raises ValueError for invalid inputs."""
        properties = {
            "standard_name": "air_pressure_at_sea_level",
            "long_name": "Sea Level Pressure",
            "units": "Pa",
        }
        file_name = netcdf_psl_file(properties)
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="psl",
            output_commands=[TEOutputCommand(var="invalid", operator="min", dist=1)],
        )
        tracker = TETracker(dn_params)
        with pytest.raises(
            ValueError,
            match=r"Variable 'invalid' not found in input files.",
        ):
            tracker.set_metadata()

    def test_te_tracker_variable_metadata_unknown(self, netcdf_psl_file) -> None:
        """Check set_metadata for missing metadata."""
        file_name = netcdf_psl_file({})
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="psl",
            output_commands=[TEOutputCommand(var="psl", operator="min", dist=1)],
        )
        tracker = TETracker(dn_params)
        tracker.set_metadata()
        metadata = tracker.variable_metadata
        assert "psl" in metadata
        assert metadata["psl"].properties["standard_name"] == "psl"
        assert metadata["psl"].properties["long_name"] == "psl"
        assert metadata["psl"].properties["units"] == "unknown"
        expected_cell_method = cf.CellMethod(
            "area",
            "minimum",
            qualifiers={"comment": "lesser circle of radius 1 degrees"},
        )
        assert metadata["psl"].constructs == [expected_cell_method]
        assert metadata["psl"].construct_kwargs is None

    def test_te_tracker_global_metadata(self, netcdf_psl_file) -> None:
        """Check set_metadata correctly sets _global_metadata."""
        file_name = netcdf_psl_file({})
        dn_params = DetectNodesParameters(
            in_data=[file_name],
            search_by_min="psl",
        )
        tracker = TETracker(dn_params)
        tracker.set_metadata()

        assert tracker.global_metadata is not None
        assert tracker.global_metadata == {
            "tctrack_version": importlib.metadata.version("tctrack"),
            "tctrack_tracker": "TETracker",
            "detect_nodes_parameters": json.dumps(
                asdict(tracker.detect_nodes_parameters)
            ),
            "stitch_nodes_parameters": json.dumps(
                asdict(tracker.stitch_nodes_parameters)
            ),
        }

    def _mock_trajectories_data(self):
        """Generate expected track data for a given track index and variable names."""
        return {
            "varnames": ["i", "j", "lon", "lat", "psl", "orog"],
            "datenames": ["year", "month", "day", "hour"],
            "trajectories": [
                [
                    {
                        "date": ["1950", "1", "1", "3"],
                        "line": [
                            "164",
                            "332",
                            "57.832031",
                            "-12.070312",
                            "1.005377e+05",
                            "0.000000e+00",
                        ],
                    },
                    {
                        "date": ["1950", "1", "1", "6"],
                        "line": [
                            "163",
                            "332",
                            "57.480469",
                            "-12.070312",
                            "1.005820e+05",
                            "0.000000e+00",
                        ],
                    },
                ],
                [
                    {
                        "date": ["1950", "1", "2", "0"],
                        "line": [
                            "843",
                            "275",
                            "296.542969",
                            "-25.429688",
                            "9.970388e+04",
                            "2.633214e+02",
                        ],
                    },
                    {
                        "date": ["1950", "1", "2", "6"],
                        "line": [
                            "850",
                            "266",
                            "299.003906",
                            "-27.539062",
                            "9.989988e+04",
                            "6.951086e+01",
                        ],
                    },
                ],
            ],
        }

    @pytest.fixture
    def mock_gfdl_file(self, tmp_path) -> Path:
        """Fixture to create a mock GFDL file with two trajectories."""
        file_path = tmp_path / "trajectories_out_gfdl.txt"
        mock_data = self._mock_trajectories_data()

        content = ""
        for track in mock_data["trajectories"]:
            # Add the start line for each track
            content += "start\t{}\t{}\n".format(len(track), "\t".join(track[0]["date"]))
            # Add the data lines for each observation in the track
            content += "\n".join(
                "\t{}\t{}".format("\t".join(obs["line"]), "\t".join(obs["date"]))
                for obs in track
            )
            content += "\n"

        file_path.write_text(content)
        return file_path

    @pytest.fixture
    def mock_csv_file(self, tmp_path) -> Path:
        """Fixture to create a mock CSV file with two trajectories."""
        file_path = tmp_path / "trajectories_out_csv.txt"
        mock_data = self._mock_trajectories_data()

        content = (
            ",".join(["track_id"] + mock_data["datenames"] + mock_data["varnames"])
            + "\n"
        )
        for i, track in enumerate(mock_data["trajectories"]):
            # Add the data lines for each observation in the track
            content += "\n".join(
                ",".join([str(i)] + obs["date"] + obs["line"]) for obs in track
            )
            content += "\n"

        file_path.write_text(content)
        return file_path

    @pytest.fixture
    def mock_csvnohead_file(self, tmp_path) -> Path:
        """Fixture to create mock CSV file with two trajectories and without header."""
        file_path = tmp_path / "trajectories_out_csvnohead.txt"
        mock_data = self._mock_trajectories_data()

        content = ""
        for i, track in enumerate(mock_data["trajectories"]):
            # Add the data lines for each observation in the track
            content += "\n".join(
                ",".join([str(i)] + obs["date"] + obs["line"]) for obs in track
            )
            content += "\n"

        file_path.write_text(content)
        return file_path

    @pytest.mark.parametrize(
        "file_format, mock_file_fixture, expected_trajectories",
        [
            (
                "gfdl",
                "mock_gfdl_file",
                [
                    {
                        "trajectory_id": 1,
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
                        "trajectory_id": 2,
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
                        "trajectory_id": 0,
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
                        "trajectory_id": 1,
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
                        "trajectory_id": 0,
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
                        "trajectory_id": 1,
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
    def test_trajectories(
        self, file_format, mock_file_fixture, expected_trajectories, request
    ):
        """Test trajectories method for different formats and multiple trajectories."""
        mock_file = request.getfixturevalue(mock_file_fixture)
        dn_params = DetectNodesParameters(in_data=["input_file.nc"])
        tracker = TETracker(dn_params)
        tracker.stitch_nodes_parameters.output_dir = str(mock_file.parent)
        tracker.stitch_nodes_parameters.output_file = mock_file.name
        tracker.stitch_nodes_parameters.out_file_format = file_format
        trajectories = tracker.trajectories()

        # Assertions
        assert len(trajectories) == len(expected_trajectories)
        for track, expected in zip(trajectories, expected_trajectories, strict=False):
            assert track.trajectory_id == expected["trajectory_id"]
            assert track.observations == expected["observations"]
            for key, values in expected["data"].items():
                assert track.data[key] == values

    @pytest.mark.parametrize(
        "file_format, mock_file_fixture",
        [
            ("gfdl", "mock_gfdl_file"),
            ("csvnohead", "mock_csvnohead_file"),
        ],
    )
    def test_trajectories_header_names(
        self, file_format, mock_file_fixture, request
    ) -> None:
        """
        Test the track header names are assigned correctly when not in file.

        Checks that trajectories() pulls correct variable names from
        DetectNodesParameters when there is no header information in the file
        (gfdl, csvnohead).
        """
        # Use DetectNodesParameters to set the variable names
        output_commands = [
            TEOutputCommand(var="v1", operator="min", dist=0.0),
            TEOutputCommand(var="v2", operator="min", dist=0.0),
        ]
        dn_params = DetectNodesParameters(
            in_data=["input_file.nc"], output_commands=output_commands
        )
        tracker = TETracker(dn_params)

        # Read in the trajectories
        mock_file = request.getfixturevalue(mock_file_fixture)
        tracker.stitch_nodes_parameters.output_dir = str(mock_file.parent)
        tracker.stitch_nodes_parameters.output_file = mock_file.name
        tracker.stitch_nodes_parameters.out_file_format = file_format
        trajectories = tracker.trajectories()

        # Check the track header names are correct
        cols = ["time", "grid_i", "grid_j", "lon", "lat", "v1", "v2"]
        for track in trajectories:
            assert [*track.data.keys()] == cols

    @pytest.mark.parametrize(
        "file_format, mock_file_fixture",
        [
            ("gfdl", "mock_gfdl_file"),
            ("csv", "mock_csv_file"),
            ("csvnohead", "mock_csvnohead_file"),
        ],
    )
    def test_to_netcdf_with_cf_read(
        self, file_format, mock_file_fixture, netcdf_psl_file, request, tmp_path
    ):
        """Test the to_netcdf method by writing out and validating with cf.read."""
        # Get the mock file and set up the TETracker
        mock_file = request.getfixturevalue(mock_file_fixture)
        file_name = netcdf_psl_file({})
        dn_params = DetectNodesParameters(in_data=[file_name], search_by_min="psl")
        sn_params = StitchNodesParameters(in_fmt=["lon", "lat", "v1", "v2"])
        tracker = TETracker(dn_params, stitch_nodes_parameters=sn_params)
        tracker.stitch_nodes_parameters.output_dir = str(mock_file.parent)
        tracker.stitch_nodes_parameters.output_file = mock_file.name
        tracker.stitch_nodes_parameters.out_file_format = file_format

        # Generate the NetCDF file
        output_file_path = tmp_path / "output.nc"
        output_file_name = str(output_file_path)
        tracker.to_netcdf(output_file_name)

        # Read the generated NetCDF file using cf-python
        fields = cf.read(output_file_name)

        # Validate the structure and content of the NetCDF file
        assert len(fields) == 4  # Ensure 4 fields (lat, lon become coordinates)
        trajectory_ids = fields[0].dimension_coordinate("trajectory").size
        observation_count = fields[0].dimension_coordinate("observation").size

        # Validate trajectory and observation dimensions
        assert trajectory_ids == 2
        assert observation_count == 2

        # Validate key variables like time, latitude, and longitude
        time_coord = fields[0].construct("time")
        lat_coord = fields[0].construct("lat")
        lon_coord = fields[0].construct("lon")

        assert time_coord.shape == (trajectory_ids, observation_count)
        assert lat_coord.shape == lon_coord.shape == (trajectory_ids, observation_count)

    def test_run_tracker_success(
        self, mocker, tmp_path, netcdf_psl_file, mock_gfdl_file
    ) -> None:
        """Check run_tracker runs successfully."""
        # Mock subprocess.run to simulate successful execution
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = mocker.MagicMock(
            returncode=0, stdout="Success"
        )

        # Set up psl file with metadata for time
        file_name = netcdf_psl_file({})
        dn_params = DetectNodesParameters(in_data=[file_name], search_by_min="psl")

        # Use the mock_gfdl_file to simulate the output of StitchNodes
        sn_in_fmt = ["lon", "lat", "psl", "orog"]
        dn_params = DetectNodesParameters(in_data=[file_name], search_by_min="psl")
        sn_params = StitchNodesParameters(
            output_dir=str(mock_gfdl_file.parent),
            output_file=mock_gfdl_file.name,
            in_fmt=sn_in_fmt,
        )
        tracker = TETracker(dn_params, stitch_nodes_parameters=sn_params)

        # Check run_tracker runs without error and produces an output file
        output_file = tmp_path / "trajectories.nc"
        tracker.run_tracker(output_file)
        assert output_file.exists()

    def test_run_tracker_failure(self, mocker) -> None:
        """Check run_tracker propagates RuntimeError from detect/stitch_nodes."""
        # Create tracker object
        dn_params = DetectNodesParameters(in_data=["input_file.nc"])
        tracker = TETracker(dn_params)

        # Mock subprocess.run and define a function to mock fail for a specific command
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")

        def subprocess_failure(cmd, args, **_kwargs):
            """
            Set up a mocked subprocess failure if a subprocess call contains cmd string.

            ``cmd`` is the first part of a subprocess call for which to generate
            failure.
            Can be invoked as a mocker.patch side_effect using lambdas as follows:

            >>> patched_sp_run = mocker.patch("subprocess.run")
            >>> patched_sp_run.side_effect = lambda args, **kwargs: subprocess_failure(
            >>>     "failing command", args, **kwargs
            >>> )
            """
            if args[0] == cmd:
                raise subprocess.CalledProcessError(
                    returncode=1, cmd=cmd, stderr="Error occurred"
                )
            mock_output = mocker.Mock(returncode=0, stdout="Success")
            return mock_output

        # Check a RuntimeError is correctly raised for detect_nodes
        mock_subprocess_run.side_effect = lambda args, **kwargs: subprocess_failure(
            "DetectNodes", args, **kwargs
        )
        with pytest.raises(
            RuntimeError, match="DetectNodes failed with a non-zero exit code"
        ):
            tracker.run_tracker("trajectories.nc")

        # Check a RuntimeError is correctly raised for stitch_nodes
        mock_subprocess_run.side_effect = lambda args, **kwargs: subprocess_failure(
            "StitchNodes", args, **kwargs
        )
        with pytest.raises(
            RuntimeError, match="StitchNodes failed with a non-zero exit code"
        ):
            tracker.run_tracker("trajectories.nc")


class TestTETrackerDetectNodes:
    """Tests for the detect_nodes functionality of TETracker."""

    def test_te_tracker_detect_nodes_defaults(self, mocker) -> None:
        """Checks the correct detect_nodes call is made for defaults."""
        # Mock the creation of the output directory
        mock_mkdir = mocker.patch("pathlib.Path.mkdir", autospec=True)

        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # create a TETracker with default parameters and call detect_nodes method
        dn_params = DetectNodesParameters(in_data=["input_file.nc"])
        tracker = TETracker(dn_params)
        result = tracker.detect_nodes()
        outdir = tracker._tempdir.name  # noqa: SLF001

        # Check the mkdir call made as expected
        mock_mkdir.assert_called_once_with(Path(outdir), parents=True, exist_ok=True)

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "DetectNodes",
                "--in_data",
                "input_file.nc",
                "--out",
                outdir + "/nodes.txt",
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
        # Mock the creation of the output directory
        mock_mkdir = mocker.patch("pathlib.Path.mkdir", autospec=True)

        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker with non-default parameters and call detect_nodes method
        params = DetectNodesParameters(
            in_data=["input_data.nc"],
            output_dir="custom_outputs",
            output_file="custom_nodes.txt",
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

        # Check the mkdir call made as expected
        mock_mkdir.assert_called_once_with(
            Path("custom_outputs"),
            parents=True,
            exist_ok=True,
        )

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "DetectNodes",
                "--in_data",
                "input_data.nc",
                "--out",
                "custom_outputs/custom_nodes.txt",
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
        """Checks the correct detect_nodes call is made for optional parameters.

        This will implicitly check the lod_to_te utility function.
        """
        # Mock subprocess.run to simulate successful execution
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker with optional parameters and call detect_nodes method
        params = DetectNodesParameters(
            in_data=["input1.nc", "input2.nc"],
            out_header=True,
            output_dir="outputs",
            search_by_min="psl",
            closed_contours=[
                TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0)
            ],
            time_filter="3hr",
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
                "--out",
                "outputs/nodes.txt",
                "--out_header",
                "--searchbymin",
                "psl",
                "--closedcontourcmd",
                "psl,200.0,5.5,0.0",
                "--mergedist",
                "0.0",
                "--timefilter",
                "3hr",
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
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TETracker instance
        params = DetectNodesParameters(
            in_data=["input_data.nc"], output_file="output.txt"
        )
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
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="DetectNodes", stderr="Error occurred"
        )

        # Create a TETracker instance
        params = DetectNodesParameters(
            in_data=["input_data.nc"], output_file="output.txt"
        )
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
        # Mock the creation of the output directory
        mock_mkdir = mocker.patch("pathlib.Path.mkdir", autospec=True)

        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with default StitchNodes parameters
        dn_params = DetectNodesParameters(in_data=["input_data.nc"])
        tracker = TETracker(dn_params)
        result = tracker.stitch_nodes()
        outdir = tracker._tempdir.name  # noqa: SLF001

        # Check the mkdir call made as expected
        mock_mkdir.assert_called_once_with(Path(outdir), parents=True, exist_ok=True)

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                outdir + "/trajectories.txt",
                "--in",
                outdir + "/nodes.txt",
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
        # Mock the creation of the output directory
        mock_mkdir = mocker.patch("pathlib.Path.mkdir", autospec=True)

        # Mock subprocess.run to simulate successful execution
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with non-default StitchNodes parameters
        dn_params = DetectNodesParameters(in_data=["input_data.nc"])
        sn_params = StitchNodesParameters(
            output_dir="custom_outputs",
            output_file="custom_trajectories.txt",
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
        tracker = TETracker(dn_params, stitch_nodes_parameters=sn_params)
        result = tracker.stitch_nodes()

        # Check the mkdir call made as expected
        mock_mkdir.assert_called_once_with(
            Path("custom_outputs"),
            parents=True,
            exist_ok=True,
        )

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                "custom_outputs/custom_trajectories.txt",
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
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.return_value = mocker.MagicMock(
            returncode=0, stdout="Mocked stdout output", stderr="Mocked stderr output"
        )

        # Create a TETracker instance with threshold filters
        dn_params = DetectNodesParameters(in_data=["input_data.nc"])
        sn_params = StitchNodesParameters(
            output_dir="outputs",
            in_file="nodes.txt",
            threshold_filters=threshold_filters,
        )
        tracker = TETracker(dn_params, stitch_nodes_parameters=sn_params)
        result = tracker.stitch_nodes()

        # Check subprocess call made as expected and returned outputs are passed back up
        mock_subprocess_run.assert_called_once_with(
            [
                "StitchNodes",
                "--out",
                "outputs/trajectories.txt",
                "--in",
                "nodes.txt",
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

    @pytest.mark.parametrize(
        "dn_params, sn_params, expected",
        [
            pytest.param(
                DetectNodesParameters(in_data=["infile.txt"]),
                None,
                [None, None],
                id="Check defaults",
            ),
            pytest.param(
                DetectNodesParameters(
                    in_data=["infile.txt"],
                    output_dir="outputs",
                    output_file="file.txt",
                    output_commands=[
                        TEOutputCommand(var="v1", operator="min", dist=0.0),
                        TEOutputCommand(var="v2", operator="max", dist=3.0),
                    ],
                ),
                None,
                ["outputs/file.txt", ["lon", "lat", "v1", "v2"]],
                id="Check the auto-assignment works",
            ),
            pytest.param(
                DetectNodesParameters(
                    in_data=["infile.txt"],
                    output_dir="outputs",
                    output_commands=[
                        TEOutputCommand(var="v1", operator="min", dist=0.0),
                        TEOutputCommand(var="v2", operator="max", dist=3.0),
                    ],
                ),
                StitchNodesParameters(in_file="file2.txt", in_fmt=["a", "b"]),
                ["file2.txt", ["a", "b"]],
                id="Check defined values aren't overriden",
            ),
        ],
    )
    def test_stitch_nodes_values_from_detect_nodes(
        self, dn_params, sn_params, expected
    ) -> None:
        """Check the values are being assigned properly from DetectNodesParameters."""
        tracker = TETracker(dn_params, sn_params)
        expected_in_file, expected_in_fmt = expected
        if expected_in_file is not None:
            assert tracker.stitch_nodes_parameters.in_file == expected_in_file
        assert tracker.stitch_nodes_parameters.in_fmt == expected_in_fmt

    def test_stitch_nodes_file_not_found(self, mocker) -> None:
        """Check stitch_nodes raises FileNotFoundError when executable is missing."""
        # Mock subprocess.run to simulate a FileNotFoundError
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = FileNotFoundError("Executable not found")

        # Create a TETracker instance
        dn_params = DetectNodesParameters(in_data=["input_data.nc"])
        tracker = TETracker(dn_params)

        # Assert that stitch_nodes raises FileNotFoundError
        with pytest.raises(
            FileNotFoundError,
            match="StitchNodes failed because the executable could not be found",
        ):
            tracker.stitch_nodes()

    def test_stitch_nodes_failure(self, mocker) -> None:
        """Check stitch_nodes raises RuntimeError on subprocess failure."""
        # Mock subprocess.run to simulate a failure
        mocker.patch("pathlib.Path.mkdir")
        mock_subprocess_run = mocker.patch("subprocess.run")
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="StitchNodes", stderr="Error occurred"
        )

        # Create a TETracker instance
        dn_params = DetectNodesParameters(in_data=["input_data.nc"])
        tracker = TETracker(dn_params)

        # Assert that stitch_nodes raises RuntimeError
        with pytest.raises(
            RuntimeError, match="StitchNodes failed with a non-zero exit code"
        ):
            tracker.stitch_nodes()
