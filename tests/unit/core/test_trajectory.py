"""Unit tests for trajectory.py of the TCTrack Core Python package."""

import re

import numpy as np
import pytest
from cftime import datetime

from tctrack.core import Trajectory


class TestTrajectory:
    """Tests for the Trajectory class."""

    def test_trajectory_init_list(self) -> None:
        """Test the Trajectory class initialization with time as a list."""
        trajectory = Trajectory(trajectory_id=1, time=(1950, 1, 1, 3))
        assert trajectory.trajectory_id == 1
        assert trajectory.observations == 0
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.start_time.minute == 0
        assert trajectory.start_time.second == 0
        assert trajectory.calendar == "gregorian"

    def test_trajectory_init_list_too_long(self) -> None:
        """Test the Trajectory class throws a warning if the time list is too long."""
        with pytest.warns(UserWarning, match="The list for the time is too long."):
            trajectory = Trajectory(trajectory_id=1, time=(1950, 1, 1, 3, 15, 30, 1))
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.start_time.minute == 15
        assert trajectory.start_time.second == 30

    def test_trajectory_init_datetime(self) -> None:
        """Test the Trajectory class initialization with time as a datetime."""
        time = datetime(1950, 1, 1, 3, calendar="360_day")
        trajectory = Trajectory(trajectory_id=1, time=time, calendar="gregorian")
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.start_time.minute == 0
        assert trajectory.start_time.second == 0
        assert trajectory.calendar == "360_day"  # Check the datetime overrides this

    def test_trajectory_invalid_time(self) -> None:
        """Test the Trajectory class initialization fails if time is the wrong type."""
        with pytest.raises(TypeError, match="Invalid type for 'time'"):
            _ = Trajectory(trajectory_id=1, time="str is invalid")

    @pytest.mark.parametrize(
        "calendar",
        [
            ("gregorian"),
            ("julian"),
            ("360_day"),
            ("noleap"),
        ],
    )
    def test_trajectory_caltypes(self, calendar) -> None:
        """Test the Trajectory class initialization with different calendar types."""
        trajectory = Trajectory(
            trajectory_id=1,
            time=(1950, 1, 1, 3),
            calendar=calendar,
        )
        assert trajectory.trajectory_id == 1
        assert trajectory.observations == 0
        assert trajectory.start_time.year == 1950
        assert trajectory.start_time.month == 1
        assert trajectory.start_time.day == 1
        assert trajectory.start_time.hour == 3
        assert trajectory.calendar == calendar
        assert isinstance(trajectory.start_time, datetime)

    def test_trajectory_invalid_calendar(self) -> None:
        """Test that Trajectory raises a ValueError for an invalid calendar type."""
        with pytest.raises(
            ValueError, match="Unsupported calendar type: invalid_calendar"
        ):
            Trajectory(
                trajectory_id=1,
                time=(1950, 1, 1, 3),
                calendar="invalid_calendar",
            )

    def test_trajectory_add_point(self) -> None:
        """Test the Trajectory class add_point method."""
        trajectory = Trajectory(trajectory_id=1, time=(1950, 1, 1, 3))
        variables_dict = {"grid_i": 164, "grid_j": 332, "psl": 1.005377e05}
        trajectory.add_point((1950, 1, 1, 3), variables_dict)
        assert len(trajectory.data) == 4
        assert trajectory.data["grid_i"] == [164]
        assert trajectory.data["grid_j"] == [332]
        assert trajectory.data["psl"] == [1.005377e05]

    @pytest.mark.parametrize(
        "inputs, expected_error, match",
        [
            # Missing data in the inputs
            (
                ((1950, 1, 1, 3),),
                TypeError,
                re.escape(
                    "Trajectory.add_point() missing 1 required positional argument: "
                    "'variables'"
                ),
            ),
            # Missing data in the inputs
            (
                ((1950, 1, 1, 3), "extra", {"psl": 1.005377e05}),
                TypeError,
                re.escape(
                    "Trajectory.add_point() takes 3 positional arguments "
                    "but 4 were given"
                ),
            ),
            # Invalid variable data
            (
                ((1950, 1, 1, 3), {"psl": "invalid"}),
                ValueError,
                "Invalid variable data: {'psl': 'invalid'}",
            ),
            # Invalid date/time values
            (((1950, 13, 1, 3), {"psl": 1.005377e05}), ValueError, None),
            (((1950, 1, 32, 3), {"psl": 1.005377e05}), ValueError, None),
            (((1950, 1, 1, 25), {"psl": 1.005377e05}), ValueError, None),
            # Invalid date/time type
            (
                ("invalid time", {"psl": 1.005377e05}),
                TypeError,
                "Invalid type for 'time'.",
            ),
        ],
    )
    def test_add_point_edge_cases(self, inputs, expected_error, match):
        """Test the Trajectory class add_point edge cases."""
        trajectory = Trajectory(trajectory_id=1, time=(1950, 1, 1, 3))
        with pytest.raises(expected_error, match=match):
            trajectory.add_point(*inputs)

    def test_add_multiple_points(self):
        """Test the add_multiple_points method of the Trajectory class."""
        trajectory = Trajectory(
            trajectory_id=1,
            time=(2023, 1, 1, 0),
            calendar="gregorian",
        )

        # Add multiple points
        times = [(2023, 1, 1, 0), (2023, 1, 3, 12)]
        variables = {
            "temperature": [300.0, 305.0],
            "pressure": [1000.0, 995.0],
        }

        trajectory.add_multiple_points(times, variables)

        # Validate the data
        assert len(trajectory.data["time"]) == 2
        assert np.allclose(trajectory.data["temperature"], [300.0, 305.0])
        assert np.allclose(trajectory.data["pressure"], [1000.0, 995.0])
