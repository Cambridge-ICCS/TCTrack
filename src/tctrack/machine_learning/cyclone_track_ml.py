"""Module providing a class to run a machine-learning tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

from dataclasses import dataclass

import cf

from tctrack.core import (
    TCTrackerMetadata,
    Trajectory,
)
from tctrack.core.ml_tracker import TCMLParameters, TCMLTracker


@dataclass(repr=False)
class MLParameters(TCMLParameters):
    """Dataclass containing values for parameters used by MLTracker.

    Extends :class:`~tctrack.core.ml_tracker.TCMLParameters` with the input
    file path and sets the default HuggingFace repository for this model.

    See Also
    --------
    MLTracker : The tracker class that uses these parameters.
    """

    input_file: str = ""
    """CF-NetCDF file containing the input variables for the model."""

    hf_repo_id: str = "surbhigoel456/cyclone-TC-ML"
    """HuggingFace Hub repository ID for the cyclone detection model."""


class MLTracker(TCMLTracker):
    """Class used to run the machine-learning cyclone tracking algorithm.

    The model weights are loaded from a HuggingFace Hub repository on first
    use and cached locally by ``huggingface_hub``.

    Attributes
    ----------
    parameters : MLParameters
        Class containing the parameters for the tracking algorithm.
    model : torch.jit.ScriptModule
        The loaded TorchScript model in evaluation mode.

    Examples
    --------
    >>> params = MLParameters(input_file="era5_2020.nc")
    >>> tracker = MLTracker(params)
    >>> tracker.run_tracker("trajectories.nc")
    """

    @property
    def _hf_filename(self) -> str:
        """Filename of the model weights in the HuggingFace repository."""
        return "cyclone-detect-ml-scripted.pt"

    def __init__(self, parameters: MLParameters):
        """Construct the MLTracker and load the model.

        Parameters
        ----------
        parameters : MLParameters
            Class containing the parameters for the tracking algorithm.

        Raises
        ------
        OSError
            If ``parameters.model_path`` is given but the file does not exist.
        huggingface_hub.errors.RepositoryNotFoundError
            If the HuggingFace repository cannot be found or accessed.
        """
        self.parameters: MLParameters = parameters
        self._trajectories: list[Trajectory] = []
        self._load_model(parameters)

    @property
    def _parameters(self) -> list[MLParameters]:
        """A list of the parameter objects accessible from the base class."""
        return [self.parameters]

    def read_trajectories(self) -> list[Trajectory]:
        """Parse tracker outputs into a list of :class:`tctrack.core.Trajectory`.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        return self._trajectories

    def _set_metadata(self) -> None:
        """Set the time and variable metadata attributes.

        Raises
        ------
        ValueError
            If a variable with time coordinate is not found in the input file.
        """
        vars_with_time = cf.read(self.parameters.input_file).select_by_construct("time")  # type: ignore[operator]
        if len(vars_with_time) == 0:
            msg = "No variable with 'time' coordinate found in the input file."
            raise ValueError(msg)
        time_coord = vars_with_time[0].coordinate("time")
        self._time_metadata = {
            "calendar": time_coord.calendar,
            "units": time_coord.units,
            "start_time": time_coord.data[0],
            "end_time": time_coord.data[-1],
        }

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
