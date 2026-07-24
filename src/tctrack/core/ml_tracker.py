"""Module providing abstract base classes for ML-based tracking algorithms."""

import os
from abc import abstractmethod
from dataclasses import dataclass

import torch
from huggingface_hub import hf_hub_download

from .tracker import TCTracker, TCTrackerParameters


@dataclass(repr=False)
class TCMLParameters(TCTrackerParameters):
    """Base dataclass for parameters common to all ML-based trackers."""

    model_path: str | None = None
    """
    Local path to the ``.pt`` model file. If ``None``, the model is downloaded
    from the HuggingFace Hub repository given by :attr:`hf_repo_id`.
    """

    hf_repo_id: str = ""
    """HuggingFace Hub repository ID from which the model is downloaded."""

    hf_token: str | None = None
    """
    HuggingFace access token for the model repository. If ``None``, the value
    of the ``HF_TOKEN`` environment variable is used instead.
    """

    device: str = "cpu"
    """PyTorch device to run inference (e.g. ``"cpu"`` or ``"cuda"``)."""

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


class TCMLTracker(TCTracker):
    """Abstract base class for ML-based tropical cyclone trackers.

    Extends `TCTracker` functionality common to all
    ML-based trackers e.g. loading a TorchScript model from
    HuggingFace Hub repository, and the `detect` step.

    Attributes
    ----------
    model : torch.jit.ScriptModule
        The loaded TorchScript model in evaluation mode.
    """

    model: torch.jit.ScriptModule

    @property
    @abstractmethod
    def _hf_filename(self) -> str:
        """Filename of the model weights in the HuggingFace repository."""

    def _load_model(self, parameters: TCMLParameters) -> None:
        """Load the TorchScript model and set it to evaluation mode.

        Loads model from TCMLParameters.model_path or downloads from HuggingFace Hub

        Parameters
        ----------
        parameters : TCMLParameters
            The parameter object containing model location and device type.

        Raises
        ------
        OSError
            If ``parameters.model_path`` is given but the file does not exist.
        huggingface_hub.errors.RepositoryNotFoundError
            If the HuggingFace repository cannot be found or accessed.
        """
        if parameters.model_path is not None:
            if not os.path.exists(parameters.model_path):
                msg = f"Model file not found: {parameters.model_path}"
                raise OSError(msg)
            model_file = parameters.model_path
        else:
            token = parameters.hf_token or os.environ.get("HF_TOKEN")
            model_file = hf_hub_download(
                repo_id=parameters.hf_repo_id,
                filename=self._hf_filename,
                token=token,
            )

        self.model: torch.jit.ScriptModule = torch.jit.load(
            model_file, map_location=parameters.device
        )
        self.model.eval()

    @abstractmethod
    def detect(self) -> None:
        """Run the ML model on the preprocessed input to obtain TC candidates.

        This method should populate the internal state (e.g. a list of
        candidate detections) that can be used for reading trajectories or the stitching
        step.
        """
