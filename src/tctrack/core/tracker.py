"""Module providing an abstract base classes for creating specific trackers."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, fields

import cf
from cftime import date2num

from .core import Trajectory


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


class TCTracker(ABC):
    """
    Abstract Base Class representing a generic TCTracker class.

    Attributes
    ----------
    _variable_metadata : dict
        A dictionary containing metadata for variables.
        This attribute must be initialized by the subclass through the
        :meth:`read_variable_metadata` method.
    """

    # Private attributes
    _variable_metadata: dict

    @property
    def variable_metadata(self) -> dict:
        """
        dict: Read-only property containing NetCDF metadata for variables.

        Raises
        ------
        AttributeError
            If `_variable_metadata` has not been initialized.
        """
        if not hasattr(self, "_variable_metadata"):
            err_msg = "_variable_metadata has not been initialized."
            raise AttributeError(err_msg)
        return self._variable_metadata

    @abstractmethod
    def read_variable_metadata(self) -> None:
        """
        Abstract method to initialize the `_variable_metadata` attribute.

        Reads and sets metadata for each variable used in tracking to be written out.
        These will be stored in the :attr:`variable_metadata` attribute as a
        dictionary of dictionaries. This will be called from the :meth:`to_netcdf`
        method.

        This method must be implemented by subclasses to populate the
        :attr:`_variable_metadata` attribute with relevant metadata for variables.

        Notes
        -----
        The :attr:`_variable_metadata` attribute is expected to be a dictionary
        where keys are variable names and values are dictionaries containing
        metadata (e.g., `standard_name`, `long_name`, `units`, `cell_method`).

        Examples
        --------
        >>> class MyTracker(TCTracker):
        ...     def read_variable_metadata(self):
        ...         self._variable_metadata = {
        ...             "example_variable": {
        ...                 "standard_name": "example_standard_name",
        ...                 "long_name": "Example Long Name",
        ...                 "units": "example_units",
        ...                 "cell_method": <CF CellMethod>,
        ...             }
        ...         }
        """

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

    def to_netcdf(self, output_file: str) -> None:
        """
        Write track trajectories to CF-compliant NetCDF trajectory file.

        Reads in trajectories based on the parameters set for a specific implementation
        of the class and writes them to a CF-Conventions compliant NetCDF trajectory
        file using cf-python.
        Trajectories are assumed to contain as a minimum data for ``lat``, ``lon``,
        and ``timestep``.

        Parameters
        ----------
        output_file: str
            filename for the output netCDF file

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
        # Read in the variable metadata from input files if not done already
        if not self.variable_metadata:
            self.read_variable_metadata()

        # Read in the trajectories generated by the tracker implementation
        # Ensure they contain assumed variables lon, lat, timestamp
        trajectories = self.trajectories()
        # Validate that each trajectory contains the required keys
        required_keys = {"timestamp", "lat", "lon"}
        for i, trajectory in enumerate(trajectories):
            missing_keys = required_keys - trajectory.data.keys()
            if missing_keys:
                errmsg = (
                    f"Trajectory {i} is missing required keys: "
                    f"{', '.join(missing_keys)}"
                )
                raise ValueError(errmsg)

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
        time_fill = -1e10
        time_data = cf.Data(
            [
                date2num(
                    trajectory.data["timestamp"],
                    units="days since 1970-01-01",
                    calendar=trajectories[
                        0
                    ].calendar,  # Assuming all trajectories use the same calendar
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
                    "days since 1970-01-01", calendar=trajectories[0].calendar
                ),
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
            },
        )

        # Create a cf.Field for each non-coordinate variable in Track.data
        # Assumes all trajectories contain the same variables
        fields = []
        for variable in trajectories[0].data:
            if variable in {"timestamp", "lat", "lon"}:
                continue

            # Define the variable metadata
            properties = self.variable_metadata.get(variable, {})
            properties["featureType"] = "trajectory"
            # Remove the cell method (if there is one) to add separately
            cell_method = properties.pop("cell_method", None)

            field = cf.Field(properties=properties)

            if cell_method is not None:
                field.set_construct(cell_method)

            axis_traj = field.set_construct(domain_axis_traj)
            axis_obs = field.set_construct(domain_axis_obs)
            field.set_construct(dim_traj, axes=(axis_traj,))
            field.set_construct(dim_obs, axes=(axis_obs,))
            field.set_construct(time_coord, axes=(axis_traj, axis_obs))
            field.set_construct(lat_coord, axes=(axis_traj, axis_obs))
            field.set_construct(lon_coord, axes=(axis_traj, axis_obs))

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
