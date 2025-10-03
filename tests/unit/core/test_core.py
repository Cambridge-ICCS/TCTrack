"""Unit tests for core.py of the TCTrack Core Python package."""

import re

import numpy as np
import pytest
from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap

from tctrack.core import Trajectory


class TestTrajectory:
    """Tests for the Trajectory class."""

    def test_trajectory_init(self) -> None:
        """Test the Trajectory class initialization."""
        trajectory = Trajectory(
            trajectory_id=1, observations=0, year=1950, month=1, day=1, hour=3
        )
        assert trajectory.trajectory_id == 1
        assert trajectory.observations == 0
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.calendar == "gregorian"

    @pytest.mark.parametrize(
        "calendar, expected_type",
        [
            ("gregorian", DatetimeGregorian),
            ("360_day", Datetime360Day),
            ("noleap", DatetimeNoLeap),
        ],
    )
    def test_trajectory_caltypes(self, calendar, expected_type) -> None:
        """Test the Trajectory class initialization with different calendar types."""
        trajectory = Trajectory(
            trajectory_id=1,
            observations=0,
            year=1950,
            month=1,
            day=1,
            hour=3,
            calendar=calendar,
        )
        assert trajectory.trajectory_id == 1
        assert trajectory.observations == 0
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.calendar == calendar
        assert isinstance(trajectory.start_time, expected_type)

    def test_trajectory_invalid_calendar(self) -> None:
        """Test that Trajectory raises a ValueError for an invalid calendar type."""
        with pytest.raises(
            ValueError, match="Unsupported calendar type: invalid_calendar"
        ):
            Trajectory(
                trajectory_id=1,
                observations=0,
                year=1950,
                month=1,
                day=1,
                hour=3,
                calendar="invalid_calendar",
            )

    def test_trajectory_add_point(self) -> None:
        """Test the Trajectory class add_point method."""
        trajectory = Trajectory(
            trajectory_id=1, observations=0, year=1950, month=1, day=1, hour=3
        )
        variables_dict = {"grid_i": 164, "grid_j": 332, "psl": 1.005377e05}
        trajectory.add_point(1950, 1, 1, 3, variables_dict)
        assert len(trajectory.data) == 4
        assert trajectory.data["grid_i"] == [164]
        assert trajectory.data["grid_j"] == [332]
        assert trajectory.data["psl"] == [1.005377e05]

    @pytest.mark.parametrize(
        "inputs, expected_error, match",
        [
            # Missing data in the inputs
            (
                (1950, 1, 1, {"psl": 1.005377e05}),
                TypeError,
                re.escape(
                    "Trajectory.add_point() missing 1 required positional argument: "
                    "'variables'"
                ),
            ),
            # Missing data in the inputs
            (
                (1950, 1, 1, 3, "extra", {"psl": 1.005377e05}),
                TypeError,
                re.escape(
                    "Trajectory.add_point() takes 6 positional arguments "
                    "but 7 were given"
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
        """Test the Trajectory class add_point edge cases."""
        trajectory = Trajectory(
            trajectory_id=1, observations=0, year=1950, month=1, day=1, hour=3
        )
        with pytest.raises(expected_error, match=match):
            trajectory.add_point(*inputs)

    def test_add_multiple_points(self):
        """Test the add_multiple_points method of the Trajectory class."""
        trajectory = Trajectory(
            trajectory_id=1,
            observations=0,
            year=2023,
            month=1,
            day=1,
            hour=0,
            calendar="gregorian",
        )

        # Add multiple points
        years = [2023, 2023]
        months = [1, 1]
        days = [1, 3]
        hours = [0, 12]
        variables = {
            "temperature": [300.0, 305.0],
            "pressure": [1000.0, 995.0],
        }

        trajectory.add_multiple_points(years, months, days, hours, variables)

        # Validate the data
        assert len(trajectory.data["timestamp"]) == 2
        assert np.allclose(trajectory.data["temperature"], [300.0, 305.0])
        assert np.allclose(trajectory.data["pressure"], [1000.0, 995.0])
