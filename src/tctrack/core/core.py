"""Core components and utilities used throughout the TCTrack codebase."""

import numbers

from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap


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
    start_time : Datetime360Day | DatetimeNoLeap | DatetimeGregorian
        Start time of the trajectory as a datetime or cftime object.
    data : dict
        Dict of data for various variables along the trajectory.
        Timestamp and other variables as supplied in file.
    """

    def __init__(  # noqa: PLR0913 - too many arguments
        self,
        trajectory_id: int,
        observations: int,
        year: int,
        month: int,
        day: int,
        hour: int,
        calendar: str = "gregorian",
    ):
        """
        Initialize a Trajectory object.

        Parameters
        ----------
        trajectory_id : int
            The unique identifier for the trajectory.
        observations : int
            The number of points in the trajectory.
        year : int
            The starting year of the trajectory.
        month : int
            The starting month of the trajectory (1-12).
        day : int
            The starting day of the trajectory (1-31, depending on the calendar).
        hour : int
            The starting hour of the trajectory (0-23).
        calendar : str, optional
            The calendar type to use for datetime handling. Options are
            "gregorian", "360_day", or "noleap".
        """
        self.trajectory_id = trajectory_id
        self.observations = observations
        self.calendar = calendar
        self.start_time = self._create_datetime(year, month, day, hour)
        self.data: dict = {}

    def _create_datetime(self, year: int, month: int, day: int, hour: int):
        """
        Create a cftime object based on the specified calendar attribute.

        Parameters
        ----------
        year : int
            year as an integer
        month : int
            month as an integer
        day : int
            day as an integer
        hour : int
            hour as an integer

        Returns
        -------
        datetime : Datetime360Day | DatetimeNoLeap | DatetimeGregorian
        """
        if self.calendar == "360_day":
            return Datetime360Day(year, month, day, hour)
        elif self.calendar == "noleap":
            return DatetimeNoLeap(year, month, day, hour)
        elif self.calendar in {"gregorian", "standard"}:
            return DatetimeGregorian(year, month, day, hour)
        else:
            msg = (
                f"Unsupported calendar type: {self.calendar}. "
                "Supported types are '360_day', 'noleap', or 'gregorian'/'standard'."
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

    def add_point(self, year: int, month: int, day: int, hour: int, variables: dict):
        """
        Add a data point to the trajectory.

        Parameters
        ----------
        year : int
          The year of the data point.
        month : int
          The month of the data point (1-12).
        day : int
          The day of the data point (1-31, depending on the calendar).
        hour : int
          The hour of the data point (0-23).
        variables : dict
            A dict containing any variables for the point as key : value pairs
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

        timestamp = self._create_datetime(year, month, day, hour)

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

    def add_multiple_points(
        self,
        years: list[int],
        months: list[int],
        days: list[int],
        hours: list[int],
        variables: dict,
    ):
        """
        Add multiple data points to the trajectory in one go.

        Parameters
        ----------
        years : list[int]
          The year of the data point.
        months : list[int]
          The month of the data point (1-12).
        days : list[int]
          The day of the data point (1-31, depending on the calendar).
        hours : list[int]
          The hour of the data point (0-23).
        variables : dict
            A dict containing arrays of any variables for the point
            as key : list[value] pairs

        Raises
        ------
        ValueError
            If the lengths of years, months, days, hours, or any variable array
            do not match.
        """
        # Ensure all input arrays are of the same length
        num_points = len(years)
        if not all(len(lst) == num_points for lst in [months, days, hours]) or not all(
            len(values) == num_points for values in variables.values()
        ):
            err_msg = "All input arrays must have the same length."
            raise ValueError(err_msg)

        # Add each data point
        for year, month, day, hour, variable_values in zip(
            years,
            months,
            days,
            hours,
            zip(*variables.values(), strict=False),
            strict=False,
        ):
            self.add_point(
                year,
                month,
                day,
                hour,
                dict(
                    zip(
                        variables.keys(),
                        variable_values,
                        strict=False,
                    )
                ),
            )
