"""Module providing a class to run a tracking algorithm.

Replace all instances of NAME with the name of the tracking algorithm.

Refer to the detailed instructions linked below and remember to add documentation and
tests for the tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

from dataclasses import dataclass

import cf

from tctrack.core import (
    TCTracker,
    TCTrackerMetadata,
    TCTrackerParameters,
    Trajectory,
)


@dataclass(repr=False)
class NAMEParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by NAMETracker."""

    input_file: str
    """CF-NetCDF file containing the input variables (as an example - not required)."""


class NAMETracker(TCTracker):
    """Class used to run the NAME tracking algorithm.

    Attributes
    ----------
    parameters : NAMEParameters
        Class containing the parameters for the tracking algorithm.
    """

    def __init__(self, parameters: NAMEParameters):
        """Construct the tracker.

        Multiple classes can be used to group the parameters.

        Parameters
        ----------
        parameters : NAMEParameters
            Class containing the parameters for the tracking algorithm(s).
            Defaults to the default values in NAMEParameters class.
        """
        self.parameters: NAMEParameters = parameters

        # Set any private attributes, make changes to parameters, setup, etc

    @property
    def _parameters(self) -> list[TCTrackerParameters]:
        """A list of the parameter objects that is accessible from the base class."""
        return [self.parameters]

    def read_trajectories(self) -> list[Trajectory]:
        """Parse tracker outputs into a list of :class:`tctrack.core.Trajectory`.

        This can be left empty if the tracking algorithm creates the list rather than
        writing an intermediate output file. In this case, the list should be passed to
        the `trajectories` argument of :meth:`to_netcdf`.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        trajectories: list[Trajectory] = []
        return trajectories

    def _set_metadata(self) -> None:
        """Set the time and variable metadata attributes.

        Raises
        ------
        ValueError
            If a variable with time coordinate is not found in the input file.
        """
        # Set the _time_metadata dictionary using information from the input files
        # For example, if the input is a CF-netcdf file with a 'time' coordinate:
        vars_with_time = cf.read(self.parameters.input_file).select_by_construct("time")  # type: ignore[operator]
        if len(vars_with_time) == 0:
            msg = r"No variable with 'time' coordinate found in the input file."
            raise ValueError(msg)
        time_coord = vars_with_time[0].coordinate("time")
        self._time_metadata = {
            "calendar": time_coord.calendar,
            "units": time_coord.units,
            "start_time": time_coord.data[0],
            "end_time": time_coord.data[-1],
        }

        # Set the _variable_metadata dictionary for each output variable
        # The variable names should match those in the Trajectory objects
        self._variable_metadata = {}
        self._variable_metadata["variable_name"] = TCTrackerMetadata(
            properties={
                "standard_name": "...",
                "long_name": "...",
                "units": "...",
            }
        )

    def run_tracker(self, output_file: str) -> None:
        """Run the tracker to obtain the tropical cyclone trajectories.

        This runs each of the individual steps of the tracking algorithm.

        Finally, the output is saved as a CF-netCDF file by calling :meth:`to_netcdf`.
        If the tracking algorithm doesn't produce an intermediate output file then the
        trajectories should be passed to it as an argument.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone track trajectories are saved.

        Examples
        --------
        To set the parameters, instantiate a :class:`NAMETracker` instance and run the
        algorithm to generate trajectories:

        >>> params = NAMEParameters(...)
        >>> tracker = NAMETracker(params)
        >>> tracker.run_tracker("trajectories.nc")
        """
        # Run the steps for the tracking algorithm, eg:
        # self.preprocess()
        # self.detect()
        # self.stitch()

        # Output the trajectories as a CF-netcdf file.
        # self.to_netcdf(output_file)
