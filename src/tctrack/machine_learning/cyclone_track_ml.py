"""Module providing a class to run a machine-learning tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

import os
from dataclasses import dataclass

import cf
import torch
from huggingface_hub import hf_hub_download

from tctrack.core import (
    TCTracker,
    TCTrackerMetadata,
    TCTrackerParameters,
    Trajectory,
)


@dataclass(repr=False)
class MLParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by MLTracker.

    See Also
    --------
    MLTracker : The tracker class that uses these parameters.
    """

    input_file: str
    """CF-NetCDF file containing the input variables for the model."""

    model_path: str | None = None
    """
    Local path to the ``.pt`` model file. If ``None``, the model is downloaded
    from the HuggingFace Hub repository given by :attr:`hf_repo_id` and cached
    locally by ``huggingface_hub``.
    """

    hf_repo_id: str = "surbhigoel456/cyclone-TC-ML"
    """HuggingFace Hub repository ID from which the model is downloaded."""

    hf_token: str | None = None
    """
    HuggingFace access token for the private model repository. If ``None``,
    the value of the ``HF_TOKEN`` environment variable is used instead.
    """

    device: str = "cpu"
    """PyTorch device on which to run inference (e.g. ``"cpu"`` or ``"cuda"``)."""

    threshold: float = 0.5
    """
    Confidence threshold for retaining detections. Candidates with a model
    score below this value are discarded.
    """

    def __post_init__(self):
        """Validate parameters."""
        if not (0.0 <= self.threshold <= 1.0):
            msg = f"threshold must be in [0, 1], got {self.threshold}"
            raise ValueError(msg)


class MLTracker(TCTracker):
    """Class used to run the machine-learning cyclone tracking algorithm.

    The model weights are loaded from a private HuggingFace Hub repository on
    first use and cached locally by ``huggingface_hub``.

    Attributes
    ----------
    parameters : MLParameters
        Class containing the parameters for the tracking algorithm.
    model : torch.nn.Module
        The loaded PyTorch model in evaluation mode.

    Examples
    --------
    >>> params = MLParameters(input_file="era5_2020.nc")
    >>> tracker = MLTracker(params)
    >>> tracker.run_tracker("trajectories.nc")
    """

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

        if parameters.model_path is not None:
            model_file = parameters.model_path
        else:
            token = parameters.hf_token or os.environ.get("HF_TOKEN")
            model_file = hf_hub_download(
                repo_id=parameters.hf_repo_id,
                filename="cyclone-detect-ml-scripted.pt",
                token=token,
            )

        self.model: torch.jit.ScriptModule = torch.jit.load(
            model_file, map_location=parameters.device
        )
        self.model.eval()

    @property
    def _parameters(self) -> list[TCTrackerParameters]:
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
