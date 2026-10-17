"""Module providing the abstract base class for the converters.

The converters take in a dataset, such as IBTrACS, and convert
it to TCTrack format.
"""

from abc import ABC, abstractmethod
from typing import Any


class Converter(ABC):
    """Abstract base class for the converters."""

    tctrack_tracker: str = ""

    @abstractmethod
    def __init__(self, *, input_path: str, output_path: str, **kwargs: Any) -> None: ...

    @abstractmethod
    def convert(self) -> None:
        """Abstract method for the conversion logic."""
