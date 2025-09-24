"""Module providing a class that binds the TRACK code.

References
----------
- TRACK code on the Reading University GitLab: https://gitlab.act.reading.ac.uk/track/track
"""

import subprocess
import tempfile
from dataclasses import dataclass

from tctrack.core import TCTracker, TCTrackerParameters, Trajectory


@dataclass(repr=False)
class TRACKParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by TRACK."""

    in_data: list[str] | None = None
    """List of strings of NetCDF input files."""

    output_file: str | None = None
    """
    Output file to write to. If ``None``, a temporary file will be created for the
    lifetime of the :class:`TRACKTracker` instance.
    """


class TRACKTracker(TCTracker):
    """Class containing bindings to the TRACK code.

    Attributes
    ----------
    track_parameters : TRACKParameters
        Class containing the parameters for the TRACK algorithm(s)
    """

    # Private attributes
    _tempdir: tempfile.TemporaryDirectory

    def __init__(
        self,
        track_parameters: TRACKParameters | None = None,
    ):
        """
        Construct the TRACK class.

        Parameters
        ----------
        track_parameters : TRACKParameters
            Class containing the parameters for the TRACK algorithm(s)
            Defaults to the default values in TRACKParameters Class
        """
        if track_parameters is not None:
            self.track_parameters: TRACKParameters = track_parameters
        else:
            self.track_parameters = TRACKParameters()

        self._variable_metadata = {}

    def _run_track_process(self, command_name: str, command_list: list[str]):
        """Run a TRACK command.

        Parameters
        ----------
        command_name : str
            The name of the command to be used in the log and error reporting.
        command_list : list[str]
            The list of strings that produce the command as given by
            other private class methods.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the TRACK executables cannot be found on the system.
        RuntimeError
            If TRACK executable returns a non-zero exit code.
        """
        try:
            result = subprocess.run(  # noqa: S603 - no shell
                command_list,
                check=True,
                capture_output=True,
                text=True,
            )
            print(f"{command_name} completed successfully.")
            print(
                f"First 12 lines of output:\n"
                f"{''.join(result.stdout.splitlines(True)[:12])}"
                f"\n...\n\n"
                f"Last 12 lines of output:\n"
                f"{''.join(result.stdout.splitlines(True)[-12:])}"
            )
            return {
                "stdout": result.stdout,
                "stderr": result.stderr,
                "returncode": result.returncode,
            }
        except FileNotFoundError as exc:
            msg = (
                f"{command_name} failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc
        except subprocess.CalledProcessError as exc:
            msg = (
                f"{command_name} failed with a non-zero exit code:"
                f"({exc.returncode}):\n"
                f"{exc.stderr}"
            )
            raise RuntimeError(msg) from exc

    def track(self):
        """
        Call the TRACK utility of TRACK.

        This will make a system call out to TRACK (provided it has been installed as an
        external dependency).
        It will be run according to the parameters in the :attr:`track_parameters`
        attribute that were set when the :class:`TRACKTracker` instance was created.

        TODO: Detail the output format

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the TRACK executeables cannot be found.
        RuntimeError
            If a TRACK executable returns a non-zero exit code.

        References
        ----------
        TODO: Add reference to Source Code?

        Examples
        --------
        To set the parameters, instantiate a :class:`TRACKTracker` instance and run
        DetectNodes:

        >>> my_params = TRACKParameters(...)
        >>> my_tracker = TRACKTracker(track_params=my_params)
        >>> result = my_tracker.track()
        """
        return self._run_track_process("Track")

    def trajectories(self) -> list[Trajectory]:
        """
        Parse outputs from TRACK to list of :class:`tctrack.core.Trajectory`.

        The file to be read and its properties are based on the values in the
        TODO

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        # TODO: Define trajectories based on TRACK output data
        trajectories: list[Trajectory] = []

        return trajectories

    def _read_variable_metadata(self) -> None:
        """Read in the metadata from the input files for each variable."""
        self._variable_metadata = {}

        # TODO: Complete this method based on input data files and TRACK code

    def run_tracker(self, output_file: str):
        """Run the TRACK tracker to obtain the tropical cyclone track trajectories.

        This runs the relevant methods in order :meth:`TODO`
        The TRACK output is then saved as a CF-compliant trajectory netCDF file.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone track trajectories are saved.

        Raises
        ------
        FileNotFoundError
            - If the TRACK executables cannot be found.
        RuntimeError
            If the TRACK commands return a non-zero exit code.

        Examples
        --------
        To set the parameters, instantiate a :class:`TRACKTracker` instance and run the
        algorithms to generate trajectories:

        >>> track_params = TRACKParameters(...)
        >>> my_tracker = TRACKTracker(track_params)
        >>> my_tracker.run_tracker()
        """
        # TODO: Add calls to relevant methods to perform algorithm

        self.to_netcdf(output_file)
