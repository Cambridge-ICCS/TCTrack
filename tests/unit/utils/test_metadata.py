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
from tctrack.utils import load_tracker_metadata


class TestMetadata:
    """Tests for metadata utility functions."""

    def create_metadata_file(self, path, attrs):
        """Create a dummy NetCDF file with the supplied global attributes."""
        with Dataset(path, "w") as ds:
            for name, value in attrs.items():
                ds.setncattr(name, value)

    @pytest.mark.parametrize(
        "tracker_cls,expected_parameters",
        [
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
        ],
    )
    def test_load_tracker_metadata(self, tmp_path, tracker_cls, expected_parameters):
        """Test load_tracker_metadata reconstructs tracker and parameter objects."""
        ncfile = tmp_path / "metadata.nc"

        self.create_metadata_file(
            ncfile,
            {
                "tctrack_version": "test-version",
                "tctrack_tracker": tracker_cls.__name__,
                "tctrack_parameters": json.dumps(
                    {type(p).__name__: asdict(p) for p in expected_parameters}
                ),
            },
        )

        tracker, parameters = load_tracker_metadata(str(ncfile))

        # Check parameters are the same
        assert len(parameters) == len(expected_parameters)
        for params, expected_params in zip(
            parameters, expected_parameters, strict=True
        ):
            # Each parameter needs to be compared individually because tuples don't
            # survive the json.dumps-json.loads round-trip
            for field in fields(expected_params):
                value = getattr(params, field.name)
                expected_value = getattr(expected_params, field.name)
                if isinstance(expected_value, tuple):
                    assert tuple(value) == expected_value
                else:
                    assert value == expected_value

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
