"""Module providing a class to run a machine-learning tracking algorithm.

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
class MLParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by MLTracker."""

    input_file: str
    """CF-NetCDF file containing the input variables (as an example - not required)."""


class MLTracker(TCTracker):
    """Class used to run the machine-learning tracking algorithm.

    Attributes
    ----------
    parameters : MLParameters
        Class containing the parameters for the tracking algorithm.
    """

    def __init__(self, parameters: MLParameters):
        """Construct the tracker.

        Parameters
        ----------
        parameters : MLParameters
            Class containing the parameters for the tracking algorithm.
            Defaults to the default values in MLParameters class.
        """
        self.parameters: MLParameters = parameters

        # Set any private attributes, make changes to parameters, setup, etc

    @property
    def _parameters(self) -> list[TCTrackerParameters]:
        """A list of the parameter objects that is accessible from the base class."""
        return [self.parameters]

    def read_trajectories(self) -> list[Trajectory]:
        """Parse tracker outputs into a list of :class:`tctrack.core.Trajectory`.

        This is left empty as the tracking algorithm creates the list directly.

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

    def detect(self) -> None:
        """Run the ML detector and output the TC candidates."""

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
        To set the parameters, instantiate a :class:`MLTracker` instance and run the
        algorithm to generate trajectories:

        >>> params = MLParameters(...)
        >>> tracker = MLTracker(params)
        >>> tracker.run_tracker("trajectories.nc")
        """
        # Run the steps for the tracking algorithm, eg:
        # self.preprocess()
        # self.detect()
        # trajectories = self.stitch()

        # Output the trajectories as a CF-netcdf file.
        # self.to_netcdf(output_file, trajectories)
