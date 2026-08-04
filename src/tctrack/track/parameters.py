"""Module containing parameter classes for the TRACK tracker."""

from dataclasses import dataclass

from tctrack.core import TCTrackerParameters


@dataclass(repr=False)
class TRACKParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by TRACK."""

    base_dir: str
    """The filepath to the directory where TRACK is installed."""

    filter_distance: float | None = None
    """The minimum start-to-end distance which trajectories must travel, in degrees."""

    wind_var_names: tuple[str, str] = ("ua", "va")
    """The variable names for the Eastward and Northward Wind in the input file."""

    pressure_level: int = 85000
    """The pressure level at which to calculate the vorticity. In Pa."""

    binary: str = "bin/track.run"
    """The filepath of the main TRACK compiled binary relative to :attr:`base_dir`."""

    file_extension: str = "track_out"
    """
    The file extension to use for intermediate TRACK output files. This cannot be the
    same as part of the /path/to/TRACK/outdat.
    """

    vorticity_file: str = "vor.dat"
    """The filename for the vorticity intermediate output file."""

    filt_vorticity_file: str = "vor_T63.dat"
    """The filename for the spectral filtered vorticity intermediate output file."""

    export_inputs: bool = False
    """Flag to export the TRACK command line inputs to files."""

    read_inputs: bool = False
    """
    Flag to read the command line inputs from files instead of TRACKParameters. There
    should be a file for each step, named as ``calculate_vorticity.in``,
    ``spectral_filtering.in``, ``tracking.in``, and ``filter_trajectories.in``. If a
    file dose not exist it will generate inputs from TRACKParameters as normal.
    """

    inputs_directory: str | None = None
    """
    Directory containing files for exporting / reading command line inputs. This will be
    created if it does not exist. By default, the current directory is used.
    """
