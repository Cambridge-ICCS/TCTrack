"""Unit tests for metadata utility functions in the TCTrack utils package."""

import json
from dataclasses import asdict

import pytest
from netCDF4 import Dataset

from tctrack.tempest_extremes import TEDetectParameters, TEStitchParameters, TETracker
from tctrack.track import TRACKParameters, TRACKTracker
from tctrack.tstorms import (
    TSTORMSBaseParameters,
    TSTORMSDetectParameters,
    TSTORMSStitchParameters,
    TSTORMSTracker,
)
from tctrack.utils import load_tracker_metadata, read_tracker_metadata
from tctrack.utils.metadata import _read_metadata

METADATA_CASES = [
    pytest.param(
        TETracker,
        [TEDetectParameters(in_data=["input.nc"]), TEStitchParameters()],
        id="tempest_extremes",
    ),
    pytest.param(
        TSTORMSTracker,
        [
            TSTORMSBaseParameters(tstorms_dir="", output_dir=""),
            TSTORMSDetectParameters(
                u_in_file="u.nc",
                v_in_file="v.nc",
                vort_in_file="vort.nc",
                tm_in_file="tm.nc",
                slp_in_file="slp.nc",
            ),
            TSTORMSStitchParameters(),
        ],
        id="tstorms",
    ),
    pytest.param(
        TRACKTracker,
        [TRACKParameters(base_dir="track_dir", input_file="input.nc")],
        id="track",
    ),
]


class TestMetadata:
    """Tests for metadata utility functions."""

    def create_metadata_file(self, path, attrs):
        """Create a dummy NetCDF file with the supplied global attributes."""
        with Dataset(path, "w") as ds:
            for name, value in attrs.items():
                ds.setncattr(name, value)

    def create_file_and_dict(self, tracker_cls, parameters, ncfile):
        """Create the dummy NetCDF file and a parameter dictionary for comparison."""
        parameter_dict = {type(p).__name__: asdict(p) for p in parameters}

        # Convert tuples to lists for comparison as json.dumps doesn't preserve tuples
        for k1, v1 in parameter_dict.items():
            for k2, v2 in v1.items():
                if isinstance(v2, tuple):
                    parameter_dict[k1][k2] = list(v2)

        # Create the netcdf file
        self.create_metadata_file(
            ncfile,
            {
                "tctrack_version": "test-version",
                "tctrack_tracker": tracker_cls.__name__,
                "tctrack_parameters": json.dumps(parameter_dict),
            },
        )

        return parameter_dict

    @pytest.mark.parametrize("tracker_cls,parameters", METADATA_CASES)
    def test_read_metadata(self, tmp_path, tracker_cls, parameters):
        """Test _read_metadata returns the stored metadata unchanged."""
        ncfile = tmp_path / "metadata.nc"
        parameter_dict = self.create_file_and_dict(tracker_cls, parameters, ncfile)

        version, tracker_name, loaded_parameter_dict = _read_metadata(str(ncfile))

        assert version == "test-version"
        assert tracker_name == tracker_cls.__name__
        assert loaded_parameter_dict == parameter_dict

    @pytest.mark.parametrize("tracker_cls,parameters", METADATA_CASES)
    def test_read_tracker_metadata(self, tmp_path, capsys, tracker_cls, parameters):
        """Test read_tracker_metadata prints the metadata in the expected format."""
        ncfile = tmp_path / "metadata.nc"
        parameter_dict = self.create_file_and_dict(tracker_cls, parameters, ncfile)

        # Build up the expected standard output message
        expected_output_list = [
            "",
            "TCTrack version: test-version",
            f"TCTrack tracker: {tracker_cls.__name__}",
        ]
        for parameter_class, values in parameter_dict.items():
            expected_output_list.append("")
            expected_output_list.append(f"{parameter_class}:")
            for name, value in values.items():
                expected_output_list.append(f"\t{name}: {value}")
        expected_output = "\n".join(expected_output_list) + "\n"

        # Capture the standard output and check it is expected
        read_tracker_metadata(str(ncfile))
        assert capsys.readouterr().out == expected_output

    @pytest.mark.parametrize("tracker_cls,parameters", METADATA_CASES)
    def test_load_tracker_metadata(self, tmp_path, tracker_cls, parameters):
        """Test load_tracker_metadata reconstructs tracker and parameter objects."""
        ncfile = tmp_path / "metadata.nc"
        parameter_dict = self.create_file_and_dict(tracker_cls, parameters, ncfile)

        tracker, loaded_parameters = load_tracker_metadata(str(ncfile))

        # Check parameters are the same
        assert len(loaded_parameters) == len(parameters)
        for loaded_params, expected_params in zip(
            loaded_parameters, parameter_dict.values(), strict=True
        ):
            assert asdict(loaded_params) == expected_params

        # Check the tracker is the same
        assert tracker is tracker_cls

    def test_load_tracker_metadata_missing_attrs(self, tmp_path):
        """Test load_tracker_metadata raises a clear error for missing metadata."""
        ncfile = tmp_path / "metadata.nc"
        self.create_metadata_file(ncfile, {"tctrack_version": "test-version"})

        with pytest.raises(
            ValueError,
            match=r"is missing required TCTrack attributes\.",
        ):
            load_tracker_metadata(str(ncfile))
