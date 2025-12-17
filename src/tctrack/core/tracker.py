"""Module providing an abstract base classes for creating specific trackers."""

import importlib.metadata
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from datetime import timedelta
from typing import TypedDict

import cf
from cftime import date2num, datetime

from .trajectory import Trajectory


@dataclass
class TCTrackerParameters:
    """
    Base Data Class for containing parameters of TCTracker classes.

    Parameters for a specific algorithm should be subclassed from this base class.
    The child class should be annotated as a dataclass using ``@dataclass(repr=False)``.

    Examples
    --------
    >>> from tctrack.core import TCTrackerParameters
    >>>
    >>> @dataclass(repr=False)
    >>> class MyTrackerParameters(TCTrackerParameters):
    >>>     param_a: int
    >>>     param_b: str
    >>>
    >>> params = MyTrackerParameters(param_a=1, param_b="example")
    >>>
    >>> print(params)
    MyTrackerParameters(
        param_a     = 1
        param_b     = example
    )

    """

    def __repr__(self) -> str:
        """Provide string representation of Parameters to users."""
        attributes = "\n\t".join(
            f"{field.name} \t = {getattr(self, field.name)}" for field in fields(self)
        )
        return f"{type(self).__name__}(\n\t{attributes}\n)"

    def __str__(self) -> str:
        """Provide string representation of Parameters."""
        return self.__repr__()


class TCTrackerTimeMetadata(TypedDict):
    """Dataclass containing the time metadata for a dataset used by a Tracker."""

    calendar: str
    """The calendar type as a string."""
    units: int
    """The calendar units as a string."""
    start_time: datetime
    """The start time of the data processed as a ``cftime.datetime`` object."""
    end_time: datetime
    """The final time of the data processed as a ``cftime.datetime`` object."""


def is_typed_dict_instance(data, typed_dict_class):
    """Check if an object has all required keys of the TypedDict class."""
    return isinstance(data, dict) and all(
        key in data for key in typed_dict_class.__required_keys__
    )


@dataclass
class TCTrackerMetadata:
    """Dataclass containing the metadata for a single variable in variable_metadata."""

    properties: dict[str, str]
    """The basic metadata properties for the variable."""

    constructs: list | None = None
    """A list of any CF constructs to add, such as :class:`cf.CellMethod` constructs."""

    construct_kwargs: list[dict] | None = None
    """
    A list of kwargs (as dicts) to use when the :meth:`cf.Field.set_construct` method
    is called in :meth:`TCTracker.to_netcdf`. This can be left as ``None``, otherwise it
    must be the same length as :attr:`constructs`.
    """

    def __post_init__(self):
        """Ensure both constructs and construct_kwargs are the same length."""
        if (
            self.constructs
            and self.construct_kwargs
            and len(self.constructs) != len(self.construct_kwargs)
        ):
            msg = (
                "'constructs' and 'construct_kwargs' have mismatched lengths "
                f"(got {len(self.constructs)} and {len(self.construct_kwargs)})"
            )
            raise ValueError(msg)


class TCTracker(ABC):
    """
    Abstract Base Class representing a generic TCTracker class.

    Attributes
    ----------
    _variable_metadata : dict[str, TCTrackerMetadata] | None
        A dictionary containing metadata for variables.
        This attribute must be initialized by the subclass through the
        :meth:`set_metadata` method, prior to which it is initialised as ``None``.
    _time_metadata : dict[str, ] | None
        A dictionary containing metadata for times including calendar, units, and start
        and end times of the dataset.
        This attribute must be initialized by the subclass through the
        :meth:`set_metadata` method, prior to which it is initialised as ``None``.
    _global_metadata : dict[str, str]
        A dictionary containing global metadata about the data and TCTrack parameters.
        This attribute should be initialized by the subclass through the
        :meth:`set_metadata` method. The parameters should be output in a json format
        using ``json.dumps(asdict(parameters))``
    """

    # Private attributes
    _variable_metadata: dict[str, TCTrackerMetadata] | None = None
    _time_metadata: TCTrackerTimeMetadata | None = None
    _global_metadata: dict[str, str]

    @property
    def variable_metadata(self) -> dict:
        """
        dict: Read-only property containing NetCDF metadata for variables.

        Raises
        ------
        AttributeError
            If `_variable_metadata` has not been initialized.
        """
        if not self._variable_metadata:
            err_msg = "_variable_metadata has not been initialized."
            raise AttributeError(err_msg)
        return self._variable_metadata

    @property
    def time_metadata(self) -> TCTrackerTimeMetadata:
        """
        dict: Read-only property containing time metadata for this Tracker run.

        Raises
        ------
        AttributeError
            If `_time_metadata` has not been initialized.
        """
        if not self._time_metadata:
            err_msg = "_time_metadata has not been initialized."
            raise AttributeError(err_msg)
        if self._time_metadata and not is_typed_dict_instance(
            self._time_metadata, TCTrackerTimeMetadata
        ):
            err_msg = (
                "_time_metadata does not conform to the expected format of "
                "`TCTrackerTimeMetadata`."
            )
            raise TypeError(err_msg)

        return self._time_metadata

    @property
    def global_metadata(self) -> dict:
        """
        dict: Read-only property containing metadata for tctrack / tracker instance.

        Raises
        ------
        AttributeError
            If `_global_metadata` has not been initialized.
        """
        if not hasattr(self, "_global_metadata"):
            err_msg = "_global_metadata has not been initialized."
            raise AttributeError(err_msg)
        return self._global_metadata

    @abstractmethod
    def set_metadata(self) -> None:
        """
        Abstract method to initialize the metadata attributes.

        This method must be implemented by subclasses to populate the
        :attr:`_variable_metadata` attribute with relevant metadata for variables,
        `:attr:`_time_metadata` with metadata about the time for the dataset, and
        the :attr:`_global_metadata` attribute for metadata about the TCTrack
        parameters and name of the tracker.

        This will be called from the :meth:`to_netcdf` method.

        Notes
        -----
        The :attr:`_variable_metadata` attribute is expected to be a dictionary where
        keys are variable names and values are instances of :class:`TCTrackerMetadata`
        containing the metadata for that variable (e.g., `standard_name`, `long_name`,
        `units`).

        The :attr:`_time_metadata` attribute is expected to be a typed dictionary of
        the :class:`TCTrackerTimeMetadata` form containing the calendar type and units
        of the dataset, as well as the start and end times.

        The :attr:`_global_metadata` should contain each parameter object in a json
        format using ``json.dumps(asdict(parameters))``.

        Examples
        --------
        >>> class MyTracker(TCTracker):
        ...     def set_metadata(self):
        ...         super().set_metadata()
        ...         self._variable_metadata = {
        ...             "example_variable": TCTrackerMetadata(
        ...                 properties={
        ...                     "standard_name": "example_standard_name",
        ...                     "long_name": "Example Long Name",
        ...                     "units": "example_units",
        ...                 },
        ...                 constructs=[<CF CellMethod>],
        ...             }
        ...         }
        ...         self._time_metadata = {
        ...             "calendar": "example_calendar",
        ...             "units": "days since yyyy-mm-dd",
        ...             "start_time": cftime.datetime(
        ...                 yyyy, mm, dd, hh, calendar=self.example_calendar
        ...             ),
        ...             "end_time": cftime.datetime(
        ...                 yyyy, mm, dd, hh, calendar=self.example_calendar
        ...             ),
        ...         }
        ...         self.global_metadata["mytracker_parameters"] = json.dumps(
        ...             asdict(MyTrackerParameters)
        ...         ),
        """
        self._global_metadata = {
            "tctrack_version": importlib.metadata.version("tctrack"),
            "tctrack_tracker": type(self).__name__,
        }

    @abstractmethod
    def trajectories(self) -> list[Trajectory]:
        """
        Parse tracking algorithm outputs into list of :class:`tctrack.core.Trajectory`.

        Implementation is deferred to the specific tracking algorithm.
        For compatibility elsewhere in TCTrack trajectories are assumed to contain as
        a minimum data for ``lat``, ``lon``, and ``timestep``.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """

    def to_netcdf(self, output_file: str) -> None:  # noqa: PLR0915, PLR0912
        """
        Write track trajectories to CF-compliant NetCDF trajectory file.

        Reads in trajectories based on the parameters set for a specific implementation
        of the class and writes them to a CF-Conventions compliant NetCDF trajectory
        file using cf-python.
        Trajectories are assumed to contain as a minimum data for ``lat``, ``lon``,
        and ``timestep``.

        An ancillary field variable is added to the output file indicating any tracks
        that start/end within 1 day of the input dataset boundaries.

        Parameters
        ----------
        output_file: str
            filename for the output netCDF file
            Note: This will be placed in the local directory unless a full path is given

        Warnings
        --------
        UserWarning
            If there are no trajectories read in from the tracker outputs by
            :meth:`trajectories()` meaning no NetCDF output file can be written.

        References
        ----------
        `CF-Conventions v1.1 - H.4. Trajectory Data <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#trajectory-data>`_
        `cf-python documentation <https://ncas-cms.github.io/cf-python/index.html>`_

        Examples
        --------
        Instantiate a :class:`TCTracker` subclass instance with appropriate parameters,
        run the relevant methods to generate cyclone track trajectories, and then save
        the results to a CF-compliant trajectory file:

        >>> my_tracker = TCTracker(...)
        >>> ...
        >>> my_tracker.to_netcdf("my_netcdf_file.nc")
        """
        # Set the metadata if not done already
        if (
            not self._variable_metadata
            or not self._global_metadata
            or not self._time_metadata
        ):
            self.set_metadata()

        # Read in the trajectories generated by the tracker implementation
        # Ensure they contain assumed variables lon, lat, time
        trajectories = self.trajectories()

        if len(trajectories) == 0:
            msg = (
                "There are no trajectories in this period so no output file will be "
                "written."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)
            return

        # Validate that each trajectory contains the required keys and
        # check for trajectories starting and ending on file boundaries
        required_keys = {"time", "lat", "lon"}
        starting_trajectory = [False] * len(trajectories)
        ending_trajectory = [False] * len(trajectories)

        for i, trajectory in enumerate(trajectories):
            missing_keys = required_keys - trajectory.data.keys()
            if missing_keys:
                errmsg = (
                    f"Trajectory {i} is missing required keys: "
                    f"{', '.join(missing_keys)}"
                )
                raise ValueError(errmsg)

            # Check for trajectories starting and ending within a day of file boundaries
            if (
                trajectory.data["time"][0] - self.time_metadata["start_time"]
            ) <= timedelta(days=1):
                starting_trajectory[i] = True
            if (
                self.time_metadata["end_time"] - trajectory.data["time"][-1]
            ) <= timedelta(days=1):
                ending_trajectory[i] = True

        start_field = cf.FieldAncillary(
            data=starting_trajectory,
            properties={
                "standard_name": "status_flag",
                "long_name": "Trajectory starting at start of dataset flag.",
            },
        )
        start_field.nc_set_variable("start_flag")
        end_field = cf.FieldAncillary(
            data=ending_trajectory,
            properties={
                "standard_name": "status_flag",
                "long_name": "Trajectory finishing at enf of dataset flag.",
            },
        )
        end_field.nc_set_variable("end_flag")

        # Determine dimensions
        num_trajectories = len(trajectories)
        max_obs = max(trajectory.observations for trajectory in trajectories)

        # Define domain axes and coords based on number of trajectories and max lengths
        domain_axis_traj = cf.DomainAxis(size=num_trajectories)
        domain_axis_obs = cf.DomainAxis(size=max_obs)
        domain_axis_traj.nc_set_dimension("trajectory")
        domain_axis_obs.nc_set_dimension("observation")

        dim_traj = cf.DimensionCoordinate(
            data=cf.Data(range(num_trajectories)),
            properties={
                "standard_name": "trajectory",
                "cf_role": "trajectory_id",
                "long_name": "trajectory index",
            },
        )

        dim_obs = cf.DimensionCoordinate(
            data=cf.Data(range(max_obs)),
            properties={
                "standard_name": "observation",
                "long_name": "observation index",
            },
        )

        # Create auxiliary coordinates for time, latitude, longitude
        # Convert time from cftime to num format to write out via cf
        time_fill = -1e8
        time_data = cf.Data(
            [
                date2num(
                    trajectory.data["time"],
                    units=self.time_metadata["units"],
                    calendar=self.time_metadata["calendar"],
                ).tolist()
                + [time_fill] * (max_obs - trajectory.observations)
                for trajectory in trajectories
            ],
            fill_value=time_fill,
        )
        time_coord = cf.AuxiliaryCoordinate(
            data=time_data,
            properties={
                "standard_name": "time",
                "long_name": "time",
                "units": cf.Units(
                    self.time_metadata["units"], calendar=self.time_metadata["calendar"]
                ),
                "missing_value": time_fill,
            },
        )

        lat_lon_fill = -999.9
        lat_data = cf.Data(
            [
                trajectory.data["lat"] + [None] * (max_obs - trajectory.observations)
                for trajectory in trajectories
            ],
            fill_value=lat_lon_fill,
        )
        lat_coord = cf.AuxiliaryCoordinate(
            data=lat_data,
            properties={
                "standard_name": "lat",
                "long_name": "latitude",
                "units": "degrees_north",
                "missing_value": lat_lon_fill,
            },
        )

        lon_data = cf.Data(
            [
                trajectory.data["lon"] + [None] * (max_obs - trajectory.observations)
                for trajectory in trajectories
            ],
            fill_value=lat_lon_fill,
        )
        lon_coord = cf.AuxiliaryCoordinate(
            data=lon_data,
            properties={
                "standard_name": "lon",
                "long_name": "longitude",
                "units": "degrees_east",
                "missing_value": lat_lon_fill,
            },
        )

        # Create a cf.Field for each non-coordinate variable in Track.data
        # Assumes all trajectories contain the same variables
        fields = []
        for variable in trajectories[0].data:
            if variable in {"time", "lat", "lon"}:
                continue

            # Define the variable metadata
            metadata = self.variable_metadata.get(variable, TCTrackerMetadata({}))
            metadata.properties["featureType"] = "trajectory"

            field = cf.Field(properties=metadata.properties)

            # By default cf-python sets variable name as the standard name.
            # If there is no standard name in the metadata for this variable we set
            # it manually to something meaningful (cf-python default is `data_n`).
            if "standard_name" not in metadata.properties:
                field.nc_set_variable(variable)
            # Add any metadata constructs
            if metadata.constructs:
                for ic, construct in enumerate(metadata.constructs):
                    kwargs = {}
                    if metadata.construct_kwargs:
                        kwargs = metadata.construct_kwargs[ic]
                    field.set_construct(construct, **kwargs)

            # Add the axes / coordinates
            axis_traj = field.set_construct(domain_axis_traj)
            axis_obs = field.set_construct(domain_axis_obs)
            field.set_construct(dim_traj, axes=(axis_traj,))
            field.set_construct(dim_obs, axes=(axis_obs,))
            field.set_construct(time_coord, axes=(axis_traj, axis_obs))
            field.set_construct(lat_coord, axes=(axis_traj, axis_obs))
            field.set_construct(lon_coord, axes=(axis_traj, axis_obs))
            field.set_construct(start_field, axes=(axis_traj))
            field.set_construct(end_field, axes=(axis_traj))

            field_fill = -1e10
            variable_data = cf.Data(
                [
                    trajectory.data[variable]
                    + [None] * (max_obs - trajectory.observations)
                    for trajectory in trajectories
                ],
                fill_value=field_fill,
            )

            # Add the variable coordinate to the field
            field.set_data(variable_data, axes=(axis_traj, axis_obs))
            field.set_property("missing_value", field_fill)

            # Add the global metadata
            if self.global_metadata:
                field.nc_set_global_attributes(self.global_metadata)

            fields.append(field)

        # Write to file
        cf.write(fields, output_file)  # type: ignore[operator]

    @abstractmethod
    def run_tracker(self, output_file: str) -> None:
        """
        Run the tracker to obtain tropical cyclone track trajectories as NetCDF file.

        Implementation is deferred to the specific tracking algorithm.

        This should first run any relevant methods for generating cyclone trajectories
        from the input data files.
        The trajectories output is then saved as a CF-compliant trajectory netCDF file
        by calling the :meth:`to_netcdf()` method of this class.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone trajectories are saved.
        """
