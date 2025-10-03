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

    filter_distance: float | None = None
    """The minimum start-to-end distance which trajectories must travel."""


class TRACKTracker(TCTracker):
    """Class containing bindings to the TRACK code.

    Attributes
    ----------
    parameters : TRACKParameters
        Class containing the parameters for the TRACK algorithm(s)
    """

    # Private attributes
    _tempdir: tempfile.TemporaryDirectory
    _track_exe: str
    _track_fext: str

    def __init__(
        self,
        parameters: TRACKParameters | None = None,
    ):
        """
        Construct the TRACK class.

        Parameters
        ----------
        parameters : TRACKParameters
            Class containing the parameters for the TRACK algorithm(s)
            Defaults to the default values in TRACKParameters Class
        """
        if parameters is not None:
            self.parameters: TRACKParameters = parameters
        else:
            self.parameters = TRACKParameters()

        self._variable_metadata = {}

    def _get_initialisation_inputs(self, inputs: list[str]):
        """Adds "initialisation" inputs in common between tracking and filter_tracks.

        Parameters
        ----------
        inputs : list[str]
            The list of inputs which will be appended to in-place.
        """

        inputs.append("0") # Binary input
        inputs.append("y") # Frames separated by newline
        inputs.append("n") # Translate the grid?
        inputs.append("y") # Make periodic
        inputs.append("g") # Geodesic distance metric
        inputs.append("y") # Change projection
        inputs.append("1") # Azimuthal projection
        inputs.append("0") # Not already azimuthal
        inputs.append("b") # Both hemispheres
        inputs.append("2") # Stereogrphic azimuthal projection

        # Hemisphere 1
        inputs.append("90") # Origin lat
        inputs.append("0") # Origin lon
        inputs.append("1") # Start x
        inputs.append("193") # End x
        inputs.append("49") # Start y
        inputs.append("96") # End y
        inputs.append("0") # Define by number of points
        inputs.append("150") # New grid points in x
        inputs.append("150") # New grid points in y

        # Hemisphere 2
        inputs.append("-90") # Origin lat
        inputs.append("0") # Origin lon
        inputs.append("1") # Start x
        inputs.append("193") # End x
        inputs.append("1") # Start y
        inputs.append("48") # End y
        inputs.append("0") # Define by number of points
        inputs.append("150") # New grid points in x
        inputs.append("150") # New grid points in y

    def _get_calculate_vorticity_inputs(self) -> list[str]:
        """Builds the list of TRACK input parameters to calculate the vorticity.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        # Initialisation
        inputs.append("n") # No country map
        inputs.append("0") # No initialisation file
        inputs.append("4") # NetCDF input
        inputs.append("n") # No file summary
        inputs.append("1") # Use netcdf names
        inputs.append("y") # Uses COARDS convention
        inputs.append("ua") # Variable to use
        inputs.append("n") # Translate the grid?
        inputs.append("85000") # Pressure level to use (in Pa)
        inputs.append("n") # Make periodic? (Only needed for tracking)
        inputs.append("g") # Geodesic distance metric
        inputs.append("n") # Don't change projection
        inputs.append("1") # Start x domain at 1
        inputs.append("192") # End x domain at final lon id
        inputs.append("1") # Start y domain at 1
        inputs.append("96") # End y domain at final lat id
        inputs.append("y") # Perform analysis

        # Vorticity calculation
        inputs.append("12") # Compute vorticity
        inputs.append("0") # Use B-splines
        inputs.append("2") # Vorticity from winds
        inputs.append("y") # File has both u and v
        inputs.append("ua") # Field for u
        inputs.append("85000") # Pressure level
        inputs.append("va") # Field for v
        inputs.append("85000") # Pressure level
        inputs.append("1") # U start frame id
        inputs.append("1") # U frame step(?)
        inputs.append("10000") # U end frame id
        inputs.append("1") # V start frame id
        inputs.append("1") # V frame step (?)
        inputs.append("10000") # V end frame id
        inputs.append("y") # Continue despite missing data
        inputs.append("indat/vor850.dat") # Output file name
        inputs.append("0") # Smoopy interpolation - no spherical continuity
        inputs.append("10") # 10 wrap around (x?)
        inputs.append("0") # No smoothing - interpolation only
        inputs.append("0") # No smoothing - interpolation only
        inputs.append("0") # Exit

        return inputs

    def _get_spectral_filtering_inputs(self) -> list[str]:
        """Builds the list of TRACK input parameters to perform the spectral filtering.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        # Initialisation
        inputs.append("n") # No country map
        inputs.append("0") # No initialisation file
        inputs.append("0") # Binary input
        inputs.append("y") # Frames separated by newline
        inputs.append("n") # Translate the grid?
        inputs.append("n") # Make periodic? (Only needed for tracking)
        inputs.append("g") # Geodesic distance metric
        inputs.append("n") # Don't change projection
        inputs.append("1") # Start x domain at 1
        inputs.append("192") # End x domain at final lon id
        inputs.append("1") # Start y domain at 1
        inputs.append("96") # End y domain at final lat id
        inputs.append("y") # Perform analysis

        # Spectral Filtering
        inputs.append("4") # Spectral Filtering
        inputs.append("1") # Start frame ID
        inputs.append("1") # Frame step
        inputs.append("10000") # End frame ID
        inputs.append("1") # Fast spectral transform
        inputs.append("n") # Derived fields are not required
        inputs.append("63") # If ny>=96 use T63 truncation, otherwise T42 truncation
        inputs.append("y") # Output on new grid
        inputs.append("1") # Use a Gaussian grid
        inputs.append("63") # Truncation for output grid
        inputs.append("2") # Number of filter bands
        inputs.append("y") # Use Hoskins filter
        inputs.append("0.1") # Cutoff constant value
        inputs.append("0") # Band boundary wavenumber 1
        inputs.append("5") # Band boundary wavenumber 2
        inputs.append("63") # Band boundary wavenumber 3
        inputs.append("n") # Do not restrict values

        return inputs

    def _get_tracking_inputs(self) -> list[str]:
        """Builds the list of TRACK input parameters to perform the tracking.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        inputs.append("n") # No country map
        inputs.append("0") # No initialisation file
        self._get_initialisation_inputs(inputs)
        inputs.append("n") # No analysis
        inputs.append("n") # No existing set of object / feature data
        inputs.append("n") # Don't compute tendency

        # Repeat twice for each hemisphere
        for _ in range(2):
            inputs.append("y") # Scale the field
            inputs.append("1.0e+5") # Scaling
            inputs.append("n") # No offset subtraction
            inputs.append("1.0") # Required threshold
            inputs.append("1") # Sphery smoothing (spherical continuity)
            inputs.append("n") # Don't add another field
            inputs.append("1") # MAX thresholding
            inputs.append("1") # Start at frame 1
            inputs.append("1") # Frame interval
            inputs.append("100000") # End frame
            inputs.append("e") # Edge connectivity only (not vertex)
            inputs.append("2") # Lower limiting size of objects
            inputs.append("7") # Surface fitting of region of interest and local optimisation [0,3,4,7,8,9]
            inputs.append("n") # Don't compute anisotropy / orientation / area
            inputs.append("0") # Smoopy interpolation
            inputs.append("n") # Don't filter object feature points for unphysical values
            inputs.append("n") # Don't filter to retain the point with the largest value
            inputs.append("n") # Don't filter for points too close together
            inputs.append("n") # Don't use time average subtraction / thresholding
            inputs.append("0.") # Sphery smoothing factor (0 = spline interpolation)
            inputs.append("y") # Specify continuity properties at the poles
            inputs.append("1") # C1 continuity
            inputs.append("0") # Don't exclude boundary maxima (close to threshold boundary)
            inputs.append("0") # Smoopy smoothing factor
            inputs.append("y") # Constrained optimisation
            inputs.append("d") # Use default constraints (see data/constraints.dat)
            inputs.append("n") # Don't write object data to file as well as feature data
            inputs.append("n") # Don't plot first chosen frame
            inputs.append("0") # No other plotting

        inputs.append("n") # Don't rotate feature locations (to account for a rotated pole)
        inputs.append("n") # Don't load an existing file of track data
        inputs.append("0.2") # Input weight 1 for the modified greedy exchange algorithm
        inputs.append("0.8") # Input weight 2
        inputs.append("n") # Don't make tracking aware of any missing frames
        inputs.append("y") # Use regional upper bound displacements (uses data/zone.dat0)
        inputs.append("y") # Use adaptive tracking
        inputs.append("6.5") # Max upperbound displacement for a phantom point (dmax)
        inputs.append("1") # Penalty value of the cost function for a phantom point (phimax)
        inputs.append("n") # Don't plot initial track data
        inputs.append("n") # Don't use a different initialisation
        inputs.append("n") # Don't do a missing frame search
        inputs.append("y") # Apply zonal upper bound displacements / adaptive track smoothness
        inputs.append("n") # Don't make tracking aware of the missing frames
        inputs.append("y") # Use regional upper bound displacments
        inputs.append("y") # Use adaptive tracking
        inputs.append("0") # No further plotting (min points = 0)
        inputs.append("0") # No track plotted
        inputs.append("n") # Don't repeat with a different parametisation

        return inputs

    def _get_filter_tracks_inputs(self) -> list[str]:
        """Builds the list of TRACK input parameters to perform the track filtering.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        inputs.append("n") # No country map
        inputs.append("0") # No initialisation file
        self._get_initialisation_inputs(inputs)
        inputs.append("y") # Perform combination / analysis / etc
        inputs.append("1") # Combine track data
        inputs.append("n") # Don't use an existing combined file
        inputs.append("1") # Use matching
        inputs.append("1") # Number of batches
        inputs.append("<objout.new file>") # object file to use
        inputs.append("<tdump file>") # track file to use
        inputs.append("1") # First frame number
        inputs.append("n") # Don't check for and merge split tracks
        inputs.append("8") # Min number of points in a track
        inputs.append("1000000") # Max number of points in a track
        if (self.parameters.filter_distance is None):
            inputs.append("n") # Don't filter tracks according to distance
        else:
            inputs.append("y") # Filter tracks according to distance
            inputs.append(self.parameters.filter_distance)
            inputs.append("s") # Distance between start and finish (or along track 't')
        inputs.append("n") # Don't filter by system strength
        inputs.append("n") # Don't filter by propogation direction
        inputs.append("n") # Don't restrict the frame range
        inputs.append("a") # Use all tracks
        inputs.append("n") # We don't want the nearest grid point positions
        inputs.append("n") # Don't plot the combined track data
        inputs.append("n") # Don't redo with different minimum lifetime / data set / area of interest
        inputs.append("0") # Exit

        return inputs

    def _run_track_process(self, command_name: str, command_list: list[str], inputs:
                           list[str]):
        """Run a TRACK command.

        Parameters
        ----------
        command_name : str
            The name of the command to be used in the log and error reporting.
        command_list : list[str]
            The list of strings that produce the command as given by
            other private class methods.
        inputs : list[str]
            The list of input parameters to pass to the TRACK program.

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
                input="\n".join(inputs) + "\n",
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

    def calculate_vorticity(self):
        command = [self._track_exe, "-i", uv_file, "f", self._track_fext]
        inputs = self._get_calculate_vorticity_inputs()
        return self._run_track_process("calculate_vorticity", command, inputs)

    def spectral_filtering(self):
        command = [self._track_exe, "-i", vorticity_file, "f", self._track_fext]
        inputs = self._get_spectral_filtering_inputs()
        return self._run_track_process("spectral_filtering", command, inputs)

    def tracking(self):
        """
        Call the tracking utility of TRACK.

        This will make a system call out to TRACK (provided it has been installed as an
        external dependency) to perform the detection and stitching of tropical cyclone
        trajectories according to the :attr:`parameters` attribute that was set when the
        :class:`TRACKTracker` instance was created.

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
        >>> my_tracker = TRACKTracker(parameters=my_params)
        >>> result = my_tracker.tracking()
        """

        command = [self._track_exe, '-i', filt_vorticity_file, 'f', self._track_fext]
        inputs = self._get_tracking_inputs()
        return self._run_track_process("tracking", command, inputs)

    def filter_tracks(self):
        command = [self._track_exe, "-i", "N/A", "f", self._track_fext]
        inputs = self._get_filter_tracks_inputs()
        return self._run_track_process("filter_tracks", command, inputs)

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
        self.calculate_vorticity()
        self.spectral_filtering()
        self.tracking()
        self.filter_tracks()

        self.to_netcdf(output_file)
