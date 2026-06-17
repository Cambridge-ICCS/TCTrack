"""Module providing a class to run a machine-learning tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

from dataclasses import dataclass

from tctrack.core import (
    TCTrackerParameters,
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
