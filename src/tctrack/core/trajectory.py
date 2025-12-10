"""Module providing a class for storing data for a single cyclone trajectory."""

import numbers
import warnings
from collections.abc import Sequence
from itertools import zip_longest
from typing import TypeGuard

from cftime import datetime


def _is_datetime(x: object) -> TypeGuard[datetime]:
    """Check (in a way the type-checker knows) that it is a datetime."""
    return isinstance(x, datetime)


class Trajectory:
    """
    A single Lagrangian cyclone trajectory/track with metadata and data points.

    Attributes
    ----------
    trajectory_id : int
        The unique identifier for the trajectory.
    observations : int
        Number of points in the trajectory.
    calendar : str
        The calendar type to use for datetime handling.
        Options are "gregorian", "360_day", or "noleap".
    start_time : cftime.datetime
        Start time of the trajectory as a cftime.datetime object.
    data : dict
        Dict of data for various variables along the trajectory.
        Timestamp and other variables as supplied in file.
        For compatibility elsewhere in TCTrack trajectories are assumed to contain as
        a minimum data for ``lat``, ``lon``, and ``timestep``.
    """

    def __init__(
        self,
        trajectory_id: int,
        time: Sequence[int] | datetime,
        calendar: str = "gregorian",
    ):
        """
        Initialize a Trajectory object.

        Parameters
        ----------
        trajectory_id : int
            The unique identifier for the trajectory.
        time : Sequence[int] | cftime.datetime
            The starting time of the trajectory. If using a list of integers this should
            be in the order: year, month, day, hour, minute, second. Any values not
            provided will be set to 0. This will be converted to a ``cftime.datetime``
            using the appropriate calendar.
        calendar : str, optional
            The calendar type to use for datetime handling if the time is provided as a
            list. Options are "gregorian", "julian", "360_day", or "noleap".

        Raises
        ------
        TypeError
            If time is not a ``Sequence[int]`` or ``cftime.datetime``.
        """
        self.trajectory_id = trajectory_id
        self.observations = 0
        if isinstance(time, Sequence) and all(isinstance(t, int) for t in time):
            self.calendar = calendar
            self.start_time = self._create_datetime(time)
        elif _is_datetime(time):
            self.calendar = time.calendar
            self.start_time = time
        else:
            msg = "Invalid type for 'time'. Must be a Sequence[int] or cftime.datetime."
            raise TypeError(msg)
        self.data: dict = {}

    def _create_datetime(self, time: Sequence[int]) -> datetime:
        """
        Create a cftime object based on the specified calendar attribute.

        Parameters
        ----------
        time : Sequence[int]
            Time as a list of integers in the order: year, month, day, hour, minute,
            second. Any values not provided will be set to 0.

        Returns
        -------
        datetime : datetime
            cftime.datetime object with the appropriate calendar setting

        Raises
        ------
        ValueError
            If the calendar is not one of the supported types.
        UserWarning
            If more than six values are passed as a time.
        """
        supported_types = {"360_day", "noleap", "julian", "gregorian", "standard"}
        if self.calendar in supported_types:
            time_units = ("year", "month", "day", "hour", "minute", "second")
            # Set all the time units provided in `time`, set the rest to zero
            time_dict = dict(zip_longest(time_units, time, fillvalue=0))
            if len(time) > len(time_units):
                msg = (
                    "The list for the time is too long. "
                    "Only the first six values will be used."
                )
                warnings.warn(msg, category=UserWarning, stacklevel=2)
            return datetime(
                *[time_dict[unit] for unit in time_units],
                calendar=self.calendar,
            )
        else:
            msg = (
                f"Unsupported calendar type: {self.calendar}. "
                "Supported types are "
                f"{', '.join(f'`{caltype}`' for caltype in supported_types)}."
            )
            raise ValueError(msg)

    def __str__(self) -> str:
        """Improve the representation of Trajectory class to users."""
        return (
            f"Trajectory(observations={self.observations}, "
            f"start_time={self.start_time}, "
            f"calendar={self.calendar}, "
            f"data_points={len(self.data)})"
        )

    def add_point(self, time: Sequence[int] | datetime, variables: dict):
        """
        Add a single data point to the trajectory.

        Parameters
        ----------
        time : Sequence[int] | cftime.datetime
            The time of the data point. If a list of integers this should be in the
            order: year, month, day, hour, minute, second. Any values not provided will
            be set to 0. This will be converted to a ``cftime.datetime`` using the
            appropriate calendar.
        variables : dict
            A dict containing any variables for the point as key : value pairs

        Raises
        ------
        TypeError
            If time is not a ``Sequence[int]`` or ``cftime.datetime``.
        """
        # Validate variables as int or float
        if not isinstance(variables, dict) or not all(
            isinstance(value, numbers.Real) for value in variables.values()
        ):
            msg = (
                f"Invalid variable data: {variables}."
                " Must be a dictionary with numeric values."
            )
            raise ValueError(msg)

        if _is_datetime(time):
            timestamp = time  # This assumes the time has the same calendar
        elif isinstance(time, Sequence) and all(isinstance(t, int) for t in time):
            timestamp = self._create_datetime(time)
        else:
            msg = "Invalid type for 'time'. Must be a Sequence[int] or cftime.datetime."
            raise TypeError(msg)

        # Initialize data structure if empty
        if not self.data:
            self.data = {
                "timestamp": [],
                **{key: [] for key in variables},
            }

        # Append data to the respective lists
        self.data["timestamp"].append(timestamp)
        for key, value in variables.items():
            self.data[key].append(value)

        self.observations += 1

    def add_multiple_points(
        self,
        times: Sequence[Sequence[int] | datetime],
        variables: dict,
    ):
        """
        Add multiple data points to the trajectory in one go.

        Parameters
        ----------
        times : Sequence[Sequence[int] | cftime.datetime]
            The times of the data points. If the individual times are a list of integers
            this should be in the order: year, month, day, hour, minute, second. Any
            values not provided will be set to 0. This will be converted to a
            ``cftime.datetime`` using the appropriate calendar.
        variables : dict
            A dict containing arrays of any variables for the point
            as key : list[value] pairs

        Raises
        ------
        ValueError
            If the lengths of any variable array do not match the number of times.
        """
        # Ensure all input arrays are of the same length
        if not all(len(values) == len(times) for values in variables.values()):
            err_msg = "All input arrays must have the same length."
            raise ValueError(err_msg)

        # Add each data point
        for time, variable_values in zip(
            times,
            zip(*variables.values(), strict=False),
            strict=False,
        ):
            self.add_point(
                time,
                dict(
                    zip(
                        variables.keys(),
                        variable_values,
                        strict=False,
                    )
                ),
            )
