"""Package providing tropical cyclone tracking utilities."""

from contextlib import suppress

with suppress(FileNotFoundError):
    from tctrack import (
        core,
        geographic,
        preprocessing,
        tempest_extremes,
        track,
        tstorms,
        utils,
    )

__all__ = [
    "core",
    "geographic",
    "preprocessing",
    "tempest_extremes",
    "track",
    "tstorms",
    "utils",
]
