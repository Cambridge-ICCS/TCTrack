"""Module providing an abstract base classes for creating specific trackers."""

from dataclasses import dataclass, fields


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
    >>> @dataclass
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
