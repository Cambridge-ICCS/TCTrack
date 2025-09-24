"""Core components and utilities used throughout the TCTrack codebase."""

from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap


class Track:
    """
    Represents a single Lagrangian cyclone track with metadata and data points.

    Attributes
    ----------
    track_id : int
        The unique identifier for the track.
    observations : int
        Number of points in the track.
    calendar : str
        The calendar type to use for datetime handling.
        Options are "gregorian", "360_day", or "noleap".
    start_time : Datetime360Day | DatetimeNoLeap | DatetimeGregorian
        Start time of the track as a datetime or cftime object.
    data : dict
        Dict of data for various variables along the track.
        Timestamp and other variables as supplied in file.
    """

    def __init__(  # noqa: PLR0913 - too many arguments
        self,
        track_id: int,
        observations: int,
        year: int,
        month: int,
        day: int,
        hour: int,
        calendar: str = "gregorian",
    ):
        """
        Initialize a Track object.

        Parameters
        ----------
        track_id : int
            The unique identifier for the track.
        observations : int
            The number of points in the track.
        year : int
            The starting year of the track.
        month : int
            The starting month of the track (1-12).
        day : int
            The starting day of the track (1-31, depending on the calendar).
        hour : int
            The starting hour of the track (0-23).
        calendar : str, optional
            The calendar type to use for datetime handling. Options are
            "gregorian", "360_day", or "noleap".
        """
        self.track_id = track_id
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
        """Improve the representation of Track to users."""
        return (
            f"Track(observations={self.observations}, "
            f"start_time={self.start_time}, "
            f"calendar={self.calendar}, "
            f"data_points={len(self.data)})"
        )

    def add_point(self, year: int, month: int, day: int, hour: int, variables: dict):
        """
        Add a data point to the track.

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
            isinstance(value, (int, float)) for value in variables.values()
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
