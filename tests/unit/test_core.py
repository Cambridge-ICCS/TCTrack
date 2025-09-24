"""Unit tests for the TCTrack Core Python package."""

import re

import pytest
from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap

from tctrack.core import Track


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
