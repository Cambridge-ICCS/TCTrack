"""Module providing a class that binds the TRACK code.

References
----------
- `TRACK code on the University of Reading GitLab
  <https://gitlab.act.reading.ac.uk/track/track>`__
- Publication describing the algorithm: `Hodges et al., 2017
  <https://doi.org/10.1175/JCLI-D-16-0557.1>`__
"""

import json
import shutil
import subprocess
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import cf

from tctrack.core import TCTracker, TCTrackerMetadata, TCTrackerParameters, Trajectory
from tctrack.utils import lat_lon_sizes


@dataclass(repr=False)
class TRACKParameters(TCTrackerParameters):
    """Dataclass containing values for parameters used by TRACK."""

    base_dir: str
    """The filepath to the directory where TRACK is installed."""

    input_file: str
    """
    NetCDF input file containing the north and easterly wind speeds on a Gaussian grid
    in m/s.
    """

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


class TRACKTracker(TCTracker):
    """Class containing bindings to the TRACK code.

    Attributes
    ----------
    parameters : TRACKParameters
        Class containing the parameters for the TRACK algorithm(s)

    References
    ----------
    `TRACK code on the University of Reading GitLab <https://gitlab.act.reading.ac.uk/track/track>`__
    """

    # Private attributes
    _nx: int  # Length of the longitude dimension of the data
    _ny: int  # Length of the latitude dimension of the data

    def __init__(self, parameters: TRACKParameters):
        """Construct the TRACK class.

        Parameters
        ----------
        parameters : TRACKParameters
            Class containing the parameters for the TRACK algorithm(s)
            Defaults to the default values in TRACKParameters Class

        Raises
        ------
        FileNotFoundError
            If the input file does not exist.
        """
        self.parameters: TRACKParameters = parameters

        self._variable_metadata = {}

        # Get sizes from input file
        input_file = self.parameters.input_file
        if not Path(input_file).exists():
            msg = f"Input file does not exist ({input_file})."
            raise FileNotFoundError(msg)
        self._ny, self._nx = lat_lon_sizes(input_file)

        # Set up files in the TRACK directory (if not already)
        base_dir = self.parameters.base_dir
        shutil.copy(base_dir + "/data/zone.dat", base_dir + "/data/zone.dat0")
        shutil.copy(base_dir + "/data/adapt.dat", base_dir + "/data/adapt.dat0")

    def _get_initialisation_inputs(self, inputs: list[str]):
        """Add "initialisation" inputs common to both tracking and filter_tracks calls.

        Parameters
        ----------
        inputs : list[str]
            The list of inputs which will be appended to in-place.
        """
        inputs.append("0")  # Binary input
        inputs.append("y")  # Frames separated by newline
        inputs.append("n")  # Translate the grid?
        inputs.append("y")  # Make periodic
        inputs.append("g")  # Geodesic distance metric
        inputs.append("y")  # Change projection
        inputs.append("1")  # Azimuthal projection
        inputs.append("0")  # Not already azimuthal
        inputs.append("b")  # Both hemispheres
        inputs.append("2")  # Stereogrphic azimuthal projection

        # Hemisphere 1
        inputs.append("90")  # Origin lat
        inputs.append("0")  # Origin lon
        inputs.append("1")  # Start x
        inputs.append("193")  # End x
        inputs.append("49")  # Start y
        inputs.append("96")  # End y
        inputs.append("0")  # Define by number of points
        inputs.append("150")  # New grid points in x
        inputs.append("150")  # New grid points in y

        # Hemisphere 2
        inputs.append("-90")  # Origin lat
        inputs.append("0")  # Origin lon
        inputs.append("1")  # Start x
        inputs.append("193")  # End x
        inputs.append("1")  # Start y
        inputs.append("48")  # End y
        inputs.append("0")  # Define by number of points
        inputs.append("150")  # New grid points in x
        inputs.append("150")  # New grid points in y

    def _get_calculate_vorticity_inputs(self) -> list[str]:
        """Build the list of TRACK input parameters to calculate the vorticity.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        params = self.parameters

        inputs = []

        # Initialisation
        inputs.append("n")  # No country map
        inputs.append("0")  # No initialisation file
        inputs.append("4")  # NetCDF input
        inputs.append("n")  # No file summary
        inputs.append("1")  # Use netcdf names
        inputs.append("y")  # Uses COARDS convention
        inputs.append(params.wind_var_names[0])  # Variable to use
        inputs.append("n")  # Translate the grid?
        inputs.append(str(params.pressure_level))  # Pressure level to use (in Pa)
        inputs.append("n")  # Make periodic? (Only needed for tracking)
        inputs.append("g")  # Geodesic distance metric
        inputs.append("n")  # Don't change projection
        inputs.append("1")  # Start x domain at 1
        inputs.append(str(self._nx))  # End x domain at final lon id
        inputs.append("1")  # Start y domain at 1
        inputs.append(str(self._ny))  # End y domain at final lat id
        inputs.append("y")  # Perform analysis

        # Vorticity calculation
        inputs.append("12")  # Compute vorticity
        inputs.append("0")  # Use B-splines
        inputs.append("2")  # Vorticity from winds
        inputs.append("y")  # File has both u and v
        inputs.append(params.wind_var_names[0])  # Field for u
        inputs.append(str(params.pressure_level))  # Pressure level
        inputs.append(params.wind_var_names[1])  # Field for v
        inputs.append(str(params.pressure_level))  # Pressure level
        inputs.append("1")  # U start frame id
        inputs.append("1")  # U frame step(?)
        inputs.append("10000")  # U end frame id
        inputs.append("1")  # V start frame id
        inputs.append("1")  # V frame step (?)
        inputs.append("10000")  # V end frame id
        inputs.append("y")  # Continue despite missing data
        output_file = params.base_dir + "/indat/" + params.vorticity_file
        inputs.append(output_file)  # Output file name
        inputs.append("0")  # Smoopy interpolation - no spherical continuity
        inputs.append("10")  # 10 wrap around (x?)
        inputs.append("0")  # No smoothing - interpolation only
        inputs.append("0")  # No smoothing - interpolation only
        inputs.append("0")  # Exit

        return inputs

    def _get_spectral_filtering_inputs(self) -> list[str]:
        """Build the list of TRACK input parameters to perform the spectral filtering.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        # Initialisation
        inputs.append("n")  # No country map
        inputs.append("0")  # No initialisation file
        inputs.append("0")  # Binary input
        inputs.append("y")  # Frames separated by newline
        inputs.append("n")  # Translate the grid?
        inputs.append("n")  # Make periodic? (Only needed for tracking)
        inputs.append("g")  # Geodesic distance metric
        inputs.append("n")  # Don't change projection
        inputs.append("1")  # Start x domain at 1
        inputs.append(str(self._nx))  # End x domain at final lon id
        inputs.append("1")  # Start y domain at 1
        inputs.append(str(self._ny))  # End y domain at final lat id
        inputs.append("y")  # Perform analysis

        # Spectral Filtering
        inputs.append("4")  # Spectral Filtering
        inputs.append("1")  # Start frame ID
        inputs.append("1")  # Frame step
        inputs.append("10000")  # End frame ID
        inputs.append("1")  # Fast spectral transform
        inputs.append("n")  # Derived fields are not required
        inputs.append("63")  # If ny>=96 use T63 truncation, otherwise T42 truncation
        inputs.append("y")  # Output on new grid
        inputs.append("1")  # Use a Gaussian grid
        inputs.append("63")  # Truncation for output grid
        inputs.append("2")  # Number of filter bands
        inputs.append("y")  # Use Hoskins filter
        inputs.append("0.1")  # Cutoff constant value
        inputs.append("0")  # Band boundary wavenumber 1
        inputs.append("5")  # Band boundary wavenumber 2
        inputs.append("63")  # Band boundary wavenumber 3
        inputs.append("n")  # Do not restrict values

        return inputs

    def _get_tracking_inputs(self) -> list[str]:  # noqa: PLR0915
        """Build the list of TRACK input parameters to perform the tracking.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        inputs = []

        inputs.append("n")  # No country map
        inputs.append("0")  # No initialisation file
        self._get_initialisation_inputs(inputs)
        inputs.append("n")  # No analysis
        inputs.append("n")  # No existing set of object / feature data
        inputs.append("n")  # Don't compute tendency

        # Repeat twice for each hemisphere
        for _ in range(2):
            inputs.append("y")  # Scale the field
            inputs.append("1.0e+5")  # Scaling
            inputs.append("n")  # No offset subtraction
            inputs.append("1.0")  # Required threshold
            inputs.append("1")  # Sphery smoothing (spherical continuity)
            inputs.append("n")  # Don't add another field
            inputs.append("1")  # MAX thresholding
            inputs.append("1")  # Start at frame 1
            inputs.append("1")  # Frame interval
            inputs.append("100000")  # End frame
            inputs.append("e")  # Edge connectivity only (not vertex)
            inputs.append("2")  # Lower limiting size of objects
            # Surface fitting of region of interest and local optimisation [0,3,4,7,8,9]
            inputs.append("7")
            inputs.append("n")  # Don't compute anisotropy / orientation / area
            inputs.append("0")  # Smoopy interpolation
            # Don't filter object feature points for unphysical values
            inputs.append("n")
            # Don't filter to retain the point with the largest value
            inputs.append("n")
            inputs.append("n")  # Don't filter for points too close together
            inputs.append("n")  # Don't use time average subtraction / thresholding
            inputs.append("0.")  # Sphery smoothing factor (0 = spline interpolation)
            inputs.append("y")  # Specify continuity properties at the poles
            inputs.append("1")  # C1 continuity
            # Don't exclude boundary maxima (close to threshold boundary)
            inputs.append("0")
            inputs.append("0")  # Smoopy smoothing factor
            inputs.append("y")  # Constrained optimisation
            inputs.append("d")  # Use default constraints (see data/constraints.dat)
            inputs.append("n")  # Don't save object data to file as well as feature data
            inputs.append("n")  # Don't plot first chosen frame
            inputs.append("0")  # No other plotting

        # Don't rotate feature locations (to account for a rotated pole)
        inputs.append("n")
        inputs.append("n")  # Don't load an existing file of track data
        inputs.append("0.2")  # Input weight 1 for modified greedy exchange algorithm
        inputs.append("0.8")  # Input weight 2
        inputs.append("n")  # Don't make tracking aware of any missing frames
        inputs.append("y")  # Use regional upper bound displacements (data/zone.dat0)
        inputs.append("y")  # Use adaptive tracking
        inputs.append("6.5")  # Max upperbound displacement for a phantom point (dmax)
        inputs.append("1")  # Cost function penalty for a phantom point (phimax)
        inputs.append("n")  # Don't plot initial track data
        inputs.append("n")  # Don't use a different initialisation
        inputs.append("n")  # Don't do a missing frame search
        # Apply zonal upper bound displacements / adaptive track smoothness
        inputs.append("y")
        inputs.append("n")  # Don't make tracking aware of the missing frames
        inputs.append("y")  # Use regional upper bound displacments
        inputs.append("y")  # Use adaptive tracking
        inputs.append("0")  # No further plotting (min points = 0)
        inputs.append("0")  # No track plotted
        inputs.append("n")  # Don't repeat with a different parametisation

        return inputs

    def _get_filter_trajectories_inputs(self) -> list[str]:
        """Build the list of TRACK input parameters to perform the track filtering.

        Returns
        -------
        list[str]
            The list of input parameters.
        """
        params = self.parameters

        inputs = []

        inputs.append("n")  # No country map
        inputs.append("0")  # No initialisation file
        self._get_initialisation_inputs(inputs)
        inputs.append("y")  # Perform combination / analysis / etc
        inputs.append("1")  # Combine track data
        inputs.append("n")  # Don't use an existing combined file
        inputs.append("1")  # Use matching
        inputs.append("1")  # Number of batches
        objout_file = params.base_dir + "/outdat/objout.new." + params.file_extension
        tdump_file = params.base_dir + "/outdat/tdump." + params.file_extension
        inputs.append(objout_file)  # object file to use
        inputs.append(tdump_file)  # track file to use
        inputs.append("1")  # First frame number
        inputs.append("n")  # Don't check for and merge split tracks
        inputs.append("8")  # Min number of points in a track
        inputs.append("1000000")  # Max number of points in a track
        if params.filter_distance is None:
            inputs.append("n")  # Don't filter tracks according to distance
        else:
            inputs.append("y")  # Filter tracks according to distance
            inputs.append(str(params.filter_distance))
            inputs.append("s")  # Distance between start and finish (or along track 't')
        inputs.append("n")  # Don't filter by system strength
        inputs.append("n")  # Don't filter by propogation direction
        inputs.append("n")  # Don't restrict the frame range
        inputs.append("a")  # Use all tracks
        inputs.append("n")  # We don't want the nearest grid point positions
        inputs.append("n")  # Don't plot the combined track data
        # Don't redo with different minimum lifetime / data set / area of interest
        inputs.append("n")
        inputs.append("0")  # Exit

        return inputs

    def _prepare_inputs(self, command_name: str, inputs: list[str]) -> str:
        """Get the inputs as a string. Optionally reading from / exporting to a file.

        The file is a text file containing the raw, unlabelled inputs in a specific
        order. Refer to the _get_<step>_inputs functions for descriptions of the
        variables.

        If :attr:`~TRACKParameters.read_inputs` is ``True`` but the file doesn't exist a
        warning will be thrown and it will continue with the generated inputs.

        Parameters
        ----------
        command_name : str
            The name of the TRACK step. Used in the file name:
            ``<inputs_directory>/<command_name>.in``.
        inputs : list[str]
            The list of input parameters to pass to the TRACK program.

        Returns
        -------
        str
            A string containing the (unlabelled) input parameters separated by newline
            characters.
        """
        input_commands = "\n".join(inputs) + "\n"

        inputs_file = Path(f"{command_name}.in")
        inputs_directory = self.parameters.inputs_directory
        if inputs_directory is not None:
            inputs_file = Path(inputs_directory) / inputs_file

        # Read from file
        if self.parameters.read_inputs:
            try:
                with open(inputs_file, "r") as f:
                    input_commands = f.read()
                warnings.warn(
                    f"TRACK inputs are being read from file for '{command_name}'. "
                    "TRACKParameter values will be ignored.",
                    stacklevel=4,
                )
            except FileNotFoundError:
                warnings.warn(
                    f"Exported TRACK inputs file for '{command_name}' does not exist "
                    "despite read_inputs=True. Continuing with generated inputs.",
                    stacklevel=4,
                )

        # Export to file
        if self.parameters.export_inputs:
            if inputs_directory is not None:
                Path(inputs_directory).mkdir(parents=True, exist_ok=True)
            with open(inputs_file, "w") as f:
                f.write(input_commands)

        return input_commands

    def _run_track_process(self, command_name: str, input_file: str, inputs: list[str]):
        """Run a TRACK command.

        Parameters
        ----------
        command_name : str
            The name of the command to be used in the log and error reporting.
        input_file : str
            The filename containing the input data for the command. The file should
            first be copied into the 'indat' directory in the TRACK root directory.
        inputs : list[str]
            The list of input parameters to pass to the TRACK program.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found on the system.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.
        """
        try:
            params = self.parameters
            command = [
                params.base_dir + "/" + params.binary,
                "-i",
                input_file,
                "-f",
                params.file_extension,
            ]

            input_commands = self._prepare_inputs(command_name, inputs)

            print(f"{command_name} running...")
            result = subprocess.run(  # noqa: S603 - no shell
                command,
                input=input_commands,
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
        """Use TRACK to calculate the vorticity from the wind components.

        This method requires the wind speed components to be on a Gaussian grid in a
        netcdf file, given by :attr:`~TRACKParameters.input_file`. This is copied to the
        TRACK 'indat' folder and used as an input for a system call to TRACK to
        calculate the vorticity.

        The vorticity is written to a new file in the 'indat' folder given by
        :attr:`~TRACKParameters.vorticity_file`. This uses the following layout::

            <n_lon> <n_lat> <n_frames>
            <List of longitudes spaced 10 to a line>
            <List of latitudes spaced 10 to a line>
            FRAME 1
            <Binary block of floats for the vorticities in frame 1>
            FRAME 2
            <Binary block of floats for the vorticities in frame 2>
            ...

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found on the system.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.
        """
        params = self.parameters
        # Copy the input file to TRACK
        shutil.copy(params.input_file, params.base_dir + "/indat/")
        input_filename = Path(params.input_file).name
        # Run TRACK to perform the calculation
        inputs = self._get_calculate_vorticity_inputs()
        self._run_track_process("calculate_vorticity", input_filename, inputs)

    def spectral_filtering(self):
        """Use TRACK to perform the spectral filtering of the vorticity.

        This makes a system call to TRACK to spectrally filter the vorticity (in the
        file :attr:`~TRACKParameters.vorticity_file`) to keep only wavenumbers 6-63.

        The output file is copied to :attr:`~TRACKParameters.filt_vorticity_file` in the
        TRACK 'indat' folder, so that it can be used in :meth:`tracking`. This file
        has the same format as described in :meth:`calculate_vorticity`.

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found on the system.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.
        """
        params = self.parameters
        inputs = self._get_spectral_filtering_inputs()
        self._run_track_process("spectral_filtering", params.vorticity_file, inputs)
        shutil.copy(
            f"{params.base_dir}/outdat/specfil.{params.file_extension}_band001",
            f"{params.base_dir}/indat/{params.filt_vorticity_file}",
        )

    def tracking(self):
        """Call the tracking utility of TRACK.

        This will make a system call out to TRACK to perform the detection of tropical
        cyclone candidates using the spectrally filtered vorticity (stored in the file
        :attr:`~TRACKParameters.filt_vorticity_file`). These candidates are then
        combined into trajectories by minimising a cost function.

        The trajectories are output in the ``objout.new.<ext>`` and ``tdump.<ext>`` text
        files in the TRACK 'outdat' folder.

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found on the system.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.
        """
        params = self.parameters
        inputs = self._get_tracking_inputs()
        self._run_track_process("tracking", params.filt_vorticity_file, inputs)

    def filter_trajectories(self):
        """Use TRACK to filter the identified trajectories.

        This step takes the identified trajectories from the :meth:`tracking` step
        stored in the ``objout.new.<ext>`` and ``tdump.<ext>`` files. It then filters
        out trajectories with fewer than 8 points (i.e. based on the duration) or a
        total distance less than :attr:`~TRACKParameters.filter_distance` (if set).

        This step could also be used to combine multiple outputs if the previous
        tracking step was batched. However, this is not currently used.

        The output of this step is a ``ff_trs.<ext>.nc`` file in the TRACK 'outdat'
        folder containing the data along the filtered trajectories. To convert this into
        a CF-compliant trajectory format the :meth:`to_netcdf` function should be used.

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found on the system.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.
        """
        input_file = self.parameters.filt_vorticity_file
        inputs = self._get_filter_trajectories_inputs()
        self._run_track_process("filter_trajectories", input_file, inputs)

    def trajectories(self) -> list[Trajectory]:
        """Parse outputs from TRACK to list of :class:`tctrack.core.Trajectory`.

        This reads the output file from the filter_trajectories step, i.e.
        ``ff_trs.<ext>.nc`` in the TRACK 'outdat' folder. It also takes the times from
        the data in :attr:`TRACKParameters.input_file`.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.

        Raises
        ------
        FileNotFoundError
            If the TRACK output file does not exist.
        """
        params = self.parameters
        trajectory_file = f"{params.base_dir}/outdat/ff_trs.{params.file_extension}.nc"

        try:
            fields = cf.read(trajectory_file)  # type: ignore[operator]
        except FileNotFoundError as e:
            msg = (
                f"TRACK output trajectory file does not exist ({trajectory_file}).\n"
                "Check if TRACK completed successfully."
            )
            raise FileNotFoundError(msg) from e

        # Get the number of points for each track
        num_pts_all = fields.select_field("ncvar%NUM_PTS").array
        (ids,) = num_pts_all.nonzero()
        num_pts = num_pts_all[ids]
        i_start = fields.select_field("ncvar%FIRST_PT").array[ids]

        # Get the data for each track
        time_idx = fields.select_field("ncvar%time").array - 1  # This is 1-indexed
        lon = fields.select_field("ncvar%longitude").array
        lat = fields.select_field("ncvar%latitude").array
        intensity = fields.select_field("ncvar%intensity").array

        # Convert the time indicies to datetimes using the input file data
        all_times = cf.read(params.input_file)[0].coordinate("time")  # type: ignore[operator]
        datetimes = all_times.datetime_array[time_idx]
        years = [dt.year for dt in datetimes]
        months = [dt.month for dt in datetimes]
        days = [dt.day for dt in datetimes]
        hours = [dt.hour for dt in datetimes]

        trajectories: list[Trajectory] = []
        for it in range(len(ids)):
            i1 = i_start[it]
            i2 = i1 + num_pts[it]
            trajectory = Trajectory(
                it,
                num_pts[it],
                years[i1],
                months[i1],
                days[i1],
                hours[i1],
                calendar=all_times.calendar,
            )
            variables = {
                "lon": lon[i1:i2],
                "lat": lat[i1:i2],
                "intensity": intensity[i1:i2],
            }
            trajectory.add_multiple_points(
                years[i1:i2],
                months[i1:i2],
                days[i1:i2],
                hours[i1:i2],
                variables,
            )
            trajectories.append(trajectory)

        return trajectories

    def set_metadata(self) -> None:
        """Set the global and variable metadata attributes."""
        self._global_metadata = {
            "tctrack_tracker": type(self).__name__,
            "track_parameters": json.dumps(asdict(self.parameters)),
        }

        self._variable_metadata = {}

        plev_domain = cf.DomainAxis(size=1)
        plev = cf.Data(self.parameters.pressure_level, units="Pa")
        plev_coord = cf.AuxiliaryCoordinate(
            data=plev,
            properties={
                "standard_name": "air_pressure",
                "long_name": "pressure",
            },
        )

        self._variable_metadata["intensity"] = TCTrackerMetadata(
            properties={
                "standard_name": "atmosphere_upward_relative_vorticity",
                "long_name": f"Relative vorticity at {plev}",
                "units": "s-1",
            },
            constructs=[plev_domain, plev_coord],
            construct_kwargs=[{"key": "plev"}, {"axes": "plev"}],
        )

    def run_tracker(self, output_file: str):
        """Run the TRACK tracker to obtain the tropical cyclone track trajectories.

        This runs the relevant methods in order:

        * :meth:`calculate_vorticity`
        * :meth:`spectral_filtering`
        * :meth:`tracking`
        * :meth:`filter_trajectories`

        Finally the output is then saved as a CF-compliant trajectory netCDF file by
        calling :meth:`to_netcdf`.

        This method will overwrite the files in the 'outdat' folder in the TRACK source
        location. Therefore, any outputs from running TRACK manually should be moved
        to prevent loss of data.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone track trajectories are saved.

        Raises
        ------
        FileNotFoundError
            If the TRACK executable cannot be found.
        RuntimeError
            If the TRACK executable returns a non-zero exit code.

        Examples
        --------
        To set the parameters, instantiate a :class:`TRACKTracker` instance and run the
        algorithms to generate trajectories:

        >>> track_params = TRACKParameters(...)
        >>> my_tracker = TRACKTracker(track_params)
        >>> my_tracker.run_tracker("trajectories.nc")
        """
        self.calculate_vorticity()
        self.spectral_filtering()
        self.tracking()
        self.filter_trajectories()

        self.to_netcdf(output_file)
