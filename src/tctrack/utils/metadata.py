"""Module providing utility functions for reading metadata from output files."""

import json

from netCDF4 import Dataset

from tctrack.core import TCTracker, TCTrackerParameters
from tctrack.tempest_extremes import TEDetectParameters, TEStitchParameters, TETracker
from tctrack.track import TRACKParameters, TRACKTracker
from tctrack.tstorms import (
    TSTORMSBaseParameters,
    TSTORMSDetectParameters,
    TSTORMSStitchParameters,
    TSTORMSTracker,
)

_TRACKER_CLASSES: dict[str, type] = {
    "TETracker": TETracker,
    "TRACKTracker": TRACKTracker,
    "TSTORMSTracker": TSTORMSTracker,
}

_PARAMETER_CLASSES: dict[str, type] = {
    "TEDetectParameters": TEDetectParameters,
    "TEStitchParameters": TEStitchParameters,
    "TRACKParameters": TRACKParameters,
    "TSTORMSBaseParameters": TSTORMSBaseParameters,
    "TSTORMSDetectParameters": TSTORMSDetectParameters,
    "TSTORMSStitchParameters": TSTORMSStitchParameters,
}


def _read_metadata(filename: str) -> tuple[str, str, dict[str, dict]]:
    """Read the TCTrack metadata from the output file.

    Returns the TCTrack version, the tracker used, and the parameters (as a dictionary
    of parameter classes, each containing a dictionary of parameters).
    """
    with Dataset(filename) as ds:
        try:
            version = ds.getncattr("tctrack_version")
            tracker_name = ds.getncattr("tctrack_tracker")
            parameter_dicts = json.loads(ds.getncattr("tctrack_parameters"))
        except AttributeError as e:
            msg = f"File {filename!r} is missing required TCTrack attributes."
            raise ValueError(msg) from e

    return version, tracker_name, parameter_dicts


def read_tracker_metadata(filename: str) -> None:
    """Print the TCTrack metadata from an output file in a readable format.

    Parameters
    ----------
    filename : str
        The filename of the output netcdf file containing TC tracks.
    """
    version, tracker_name, parameter_dicts = _read_metadata(filename)

    print(f"\nTCTrack version: {version}")
    print(f"TCTrack tracker: {tracker_name}")
    for parameter_class, parameters in parameter_dicts.items():
        print(f"\n{parameter_class}:")
        for name, parameter in parameters.items():
            print(f"\t{name}: {parameter}")


def load_tracker_metadata(filename: str) -> tuple[TCTracker, list[TCTrackerParameters]]:
    """Read the TCTrack metadata from an output file, returns the parameters.

    Parameters
    ----------
    filename : str
        The filename of the output netcdf file containing TC tracks.

    Returns
    -------
    tracker : TCTracker
        A tracker object that uses the parameters from the metadata.
    parameters : list[TCTrackerParameters]
        A list of parameter objects that duplicate the information in the metadata.
    """
    _, tracker_name, parameter_dicts = _read_metadata(filename)

    # Instantiate parameter classes from each parameter dictionary
    parameters = []
    for name, parameter_dict in parameter_dicts.items():
        parameter_cls = _PARAMETER_CLASSES[name]
        parameters.append(parameter_cls(**parameter_dict))

    tracker = _TRACKER_CLASSES[tracker_name](*parameters)

    return tracker, parameters
