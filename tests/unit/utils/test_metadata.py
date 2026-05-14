"""Unit tests for metadata utility functions in the TCTrack utils package."""

import json
from dataclasses import asdict

import pytest
from netCDF4 import Dataset

from tctrack.tempest_extremes import TEDetectParameters, TEStitchParameters, TETracker
from tctrack.utils import load_tracker_metadata


class TestMetadata:
    """Tests for metadata utility functions."""

    def create_metadata_file(self, path, attrs):
        """Create a dummy NetCDF file with the supplied global attributes."""
        with Dataset(path, "w") as ds:
            for name, value in attrs.items():
                ds.setncattr(name, value)

    def test_load_tracker_metadata(self, tmp_path):
        """Test load_tracker_metadata reconstructs tracker and parameter objects."""
        ncfile = tmp_path / "metadata.nc"
        expected_detect = TEDetectParameters(in_data=["input.nc"])
        expected_stitch = TEStitchParameters()
        expected_tracker = TETracker(expected_detect, expected_stitch)

        self.create_metadata_file(
            ncfile,
            {
                "tctrack_version": "test-version",
                "tctrack_tracker": "TETracker",
                "tctrack_parameters": json.dumps(
                    {
                        "TEDetectParameters": asdict(expected_detect),
                        "TEStitchParameters": asdict(expected_stitch),
                    }
                ),
            },
        )

        tracker, parameters = load_tracker_metadata(str(ncfile))

        # Check parameters are the same
        assert len(parameters) == 2
        assert parameters[0] == expected_detect
        assert parameters[1] == expected_stitch

        # Check the tracker is the same
        assert type(tracker) is type(expected_tracker)
        assert tracker._parameters == expected_tracker._parameters  # noqa:SLF001

    def test_load_tracker_metadata_missing_attrs(self, tmp_path):
        """Test load_tracker_metadata raises a clear error for missing metadata."""
        ncfile = tmp_path / "metadata.nc"
        self.create_metadata_file(ncfile, {"tctrack_version": "test-version"})

        with pytest.raises(
            ValueError,
            match=r"is missing required TCTrack attributes\.",
        ):
            load_tracker_metadata(str(ncfile))
