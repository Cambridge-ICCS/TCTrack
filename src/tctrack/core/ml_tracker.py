"""Module providing abstract base classes for ML-based tracking algorithms."""

from dataclasses import dataclass

from .tracker import TCTrackerParameters


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
