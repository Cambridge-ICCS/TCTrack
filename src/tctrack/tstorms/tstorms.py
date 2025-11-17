"""Module providing classes that bind the TSTORMS code.

References
----------
- TSTORMS NOAA page: https://www.gfdl.noaa.gov/tstorms/
"""

import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import warnings
from dataclasses import asdict, dataclass

from tctrack.core import TCTracker, TCTrackerParameters


@dataclass(repr=False)
class TSTORMSBaseParameters(TCTrackerParameters):
    """
    Dataclass containing values used in configuring the TSTORMS install.

    This class is intended to configure the install so that executables can be run
    correctly and data stored in a desired location.
    Configuration of the actual cyclone detection algorithm is done using
    :class:`TSTORMSDetectParameters` and :class:`TSTORMSStitchParameters`.
    """

    tstorms_dir: str
    """
    Full path to the TSTORMS installation directory. This will be likely be a directory
    named `tropical_storms_pub/` and should contain `tstorms_driver/` and
    `trajectory_analysis/` subdirectories.
    """

    output_dir: str
    """
    Full path to the directory where TSTORMS outputs should be deposited.
    """

    input_dir: str | None = None
    """
    Full path to the directory where TSTORMS input files can be found.
    Defaults to None.
    """


@dataclass(repr=False)
class TSTORMSDetectParameters(TCTrackerParameters):
    """
    Dataclass of values used by the Driver operation of TSTORMS for candidate detection.

    Default values are set to match those in the TSTORMS tstorms_driver source-code.

    Raises
    ------
    ValueError
        - If northern latitude bound is less than the southern latitude bound.
        - If any of the input file filenames are left as ``None``.
    UserWarning
        - If do_thickness is set to True as this has no effect.
    """

    u_in_file: str
    """
    Filename of the u (zonal) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the u field named
    ``u_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``u850``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    """

    v_in_file: str
    """
    Filename of the v (meridional) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the v field named
    ``v_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``v850``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    """

    vort_in_file: str
    """
    Filename of the vorticity input file.
    This should be a NetCDF file containing a single lat-lon slice of the vorticity
    field named ``vort850`` at the 850 hPa level.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    """

    tm_in_file: str
    """
    Filename of the temperature input file.
    This should be a NetCDF file containing a single lat-lon slice of the mean
    temperature of the warm-core layer named ``tm``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    """

    slp_in_file: str
    """
    Filename of the sea-level pressure input file.
    This should be a NetCDF file containing a single lat-lon slice of sea-level pressure
    named ``slp``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    """

    use_sfc_wind: bool = True
    """
    Whether to use surface winds (``True``), or winds at 850 hPa level.
    """

    vort_crit: float = 3.5e-5
    """
    Critical vorticity threshold [s-1] to be exceeded to qualify as a candidate storm.
    """

    tm_crit: float = 0.5
    """
    Critical warm core threshold to be exceeded to qualify as a candidate storm.
    """

    thick_crit: float = 50.0
    """
    Critical thickness threshold to be exceeded to qualify as a candidate storm.
    Note that thickness calculations are not yet implemented in TSTORMS.
    """

    dist_crit: float = 4.0
    """
    Critical radius [degrees] within which vorticity, sea-level pressure, and other
    maxima/minima must lie within each other to qualify as a candidate storm.
    """

    lat_bound_n: float = 90.0
    """
    Northern latitude bound [degrees] below which storm detection should occur.
    """

    lat_bound_s: float = -90.0
    """
    Southern latitude bound [degrees] above which storm detection should occur.
    """

    do_spline: bool = False
    """
    Whether to use splines for detecting minma instead of gridpointwise search.
    """

    do_thickness: bool = False
    """
    Whether to use thickness of the 200-1000 hPa layer as a variable for detecting
    candidate storms. Note that this functionality is not yet implemented in TSTORMS.
    """

    def __post_init__(self):
        """Validate parameters."""
        if self.lat_bound_n < self.lat_bound_s:
            msg = (
                f"Northern latitude bound ({self.lat_bound_n}) is less than "
                f"Southern latitude bound ({self.lat_bound_s})."
            )
            raise ValueError(msg)
        if self.do_thickness:
            msg = (
                "`do_thickness` is set, but will have no effect as this feature is not "
                "implemented in TSTORMS."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)


@dataclass(repr=False)
class TSTORMSStitchParameters(TCTrackerParameters):
    """
    Dataclass containing values used by the stitching trajectory operation of TSTORMS.

    Default values are set to match those in the TSTORMS trajectory_analysis
    source-code.

    Raises
    ------
    ValueError
        If northern latitude bound is less than the southern latitude bound.
    UserWarning
        If do_thickness is set to True as this has no effect.
    """

    r_crit: float = 900.0
    """Maximum daily track length [km] between succcessive points in a trajectory."""

    wind_crit: float = 17.0
    """Critical wind speed [m/s] for trajectory calculations."""

    vort_crit: float = 3.5e-5
    """Critical vorticity threshold [s-1] for trajectory calculations."""

    tm_crit: float = 0.5
    """Critical warm core threshold for trajectory calculations."""

    thick_crit: float = 50.0
    """Critical thickness threshold for trajectory calculations."""

    n_day_crit: int = 2
    """Minimum number of days a trajectory must last to be valid."""

    do_filter: bool = True
    """
    Whether to apply filtering of trajectories.
    Filtering is based on landmask and lat-lon bounds to generate *_filt output files.
    """

    lat_bound_n: float = 40.0
    """Northern latitude bound [degrees] for trajectory filtering."""

    lat_bound_s: float = -40.0
    """Southern latitude bound [degrees] for trajectory filtering."""

    do_spline: bool = False
    """
    Whether to use splines for trajectory calculations.
    Should match the value used in the Driver routine.
    If splines used then :attr:`twc_crit` and :attr:`thick_crit` will not be used in
    comparisons.
    """

    do_thickness: bool = False
    """
    Whether to use thickness of the 200-1000 hPa layer as a variable for detecting
    candidate storms. Note that this functionality is not yet implemented in TSTORMS.
    """

    def __post_init__(self):
        """Validate parameters."""
        if self.lat_bound_n < self.lat_bound_s:
            msg = (
                f"Northern latitude bound ({self.lat_bound_n}) is less than "
                f"Southern latitude bound ({self.lat_bound_s})."
            )
            raise ValueError(msg)
        if self.do_thickness:
            msg = (
                "`do_thickness` is set, but will have no effect as this feature is not "
                "implemented in TSTORMS."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)


class TSTORMSTracker(TCTracker):
    """Class containing bindings to the TSTORMS code.

    Attributes
    ----------
    driver_parameters : TSTORMSDetectParameters | None
        Class containing the parameters for the driver detection algorithm
    trajectory_parameters : TSTORMSStitchParameters | None
        Class containing the parameters for the trajectory stitching algorithm
    """

    # Private attributes
    _tempdir: tempfile.TemporaryDirectory

    def __init__(
        self,
        tstorms_parameters: TSTORMSBaseParameters,
        detect_parameters: TSTORMSDetectParameters,
        stitch_parameters: TSTORMSStitchParameters | None = None,
    ):
        """
        Construct the TSTORMSTracker class.

        Parameters
        ----------
        tstorms_parameters : TSTORMSBaseParameters
            Class containing the parameters for setting up TSTORMS usage.
        detect_parameters : TSTORMSDetectParameters
            Class containing the parameters for the node detection in TSTORMS.
            Used to create the namelists for `tstorms_driver`.
        stitch_parameters : TSTORMSStitchParameters | None
            Class containing the parameters for the stitching algorithm of TSTORMS
            Defaults to the default values in StitchNodesParameters Class
        """
        self.tstorms_parameters: TSTORMSBaseParameters = tstorms_parameters
        self.detect_parameters: TSTORMSDetectParameters = detect_parameters

        if stitch_parameters is not None:
            self.stitch_parameters: TSTORMSStitchParameters = stitch_parameters
        else:
            self.stitch_parameters = TSTORMSStitchParameters()

        # Ensure the output directory exists, create if not.
        output_dir = self.tstorms_parameters.output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Check TSTORMSStitchParameters input arguments according to
        # TSTORMSDetectParameters
        # TODO: Decide which parameters we want to synchronise

        self._variable_metadata = {}

    def _run_tstorms_process(
        self,
        command_name: str,
        command_list: list[str],
        input_file: str,
        verbose: bool = False,
    ):
        """Run a TSTORMS command.

        Parameters
        ----------
        command_name : str
            The name of the command to be used in the log and error reporting.
        command_list : list[str]
            The list of strings that produce the command as given by
            _make_driver_call and _make_trajectories_call.
        input_file : str
            Path to an input file to be passed to the command's stdin (e.g. namelist).
            Defaults to None.
        verbose : bool
            Whether to print the entire TSTORMS output to screen in real-time or just
            the start/end summary. Defaults to False.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the tstorms executable from cannot be found.
        RuntimeError
            If tstorms executable returns a non-zero exit code.
        """
        try:
            if verbose:
                with open(input_file, "r") as stdin:
                    # Real-time output with Popen
                    process = subprocess.Popen(  # noqa: S603 - no shell
                        command_list,
                        stdin=stdin,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        shell=False,
                        bufsize=1,  # Line-buffered output
                    )

                    # Print stdout in real-time
                    for line in iter(process.stdout.readline, ""):  # type: ignore[union-attr]
                        print(line, end="")

                    stdout, stderr = process.communicate()
                    returncode = process.returncode

                    if returncode != 0:
                        err_msg = (
                            f"{command_name} failed with a non-zero exit code: "
                            f"{returncode}:\n{stderr}"
                        )
                        raise RuntimeError(err_msg)
            else:
                # verbose=False: Concise output print start/end summary
                with open(input_file, "r") as stdin:
                    print(f"Executing {command_name}...")
                    result = subprocess.run(  # noqa: S603 - no shell
                        command_list,
                        stdin=stdin,
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    stdout, stderr, returncode = (
                        result.stdout,
                        result.stderr,
                        result.returncode,
                    )

                    print(f"{command_name} completed successfully.")
                    print(
                        f"First 12 lines of output:\n"
                        f"{''.join(stdout.splitlines(True)[:12])}"
                        f"\n...\n\n"
                        f"Last 12 lines of output:\n"
                        f"{''.join(stdout.splitlines(True)[-12:])}"
                    )

            return {
                "stdout": stdout,
                "stderr": stderr,
                "returncode": returncode,
            }
        except FileNotFoundError as exc:
            msg = f"{command_name} failed because the executable could not be found."
            raise FileNotFoundError(msg) from exc
        except subprocess.CalledProcessError as exc:
            msg = (
                f"{command_name} failed with a non-zero exit code:"
                f"({exc.returncode}):\n"
                f"{exc.stderr}"
            )
            raise RuntimeError(msg) from exc

    def _write_driver_namelist(self) -> str:
        """Generate the namelist file for the tstorms_driver routine.

        Returns
        -------
        str
            The full path to the generated namelist file.

        Raises
        ------
        FileNotFoundError
            If the `tstorms_dir` does not exist.
        """
        # Ensure the tstorms_dir exists
        tstorms_driver_dir = os.path.join(
            self.tstorms_parameters.tstorms_dir, "tstorms_driver"
        )
        if not os.path.exists(tstorms_driver_dir):
            err_msg = f"TSTORMS driver directory '{tstorms_driver_dir}' does not exist."
            raise FileNotFoundError(err_msg)

        # Define the namelist file path
        namelist_path = os.path.join(tstorms_driver_dir, "nml_driver")

        # Format the namelist content
        if self.tstorms_parameters.input_dir:
            input_dir = self.tstorms_parameters.input_dir
        else:
            input_dir = ""
        detect_params = self.detect_parameters
        namelist_content = textwrap.dedent(f"""
         &nml_tstorms
           crit_vort  =  {detect_params.vort_crit:.4E}
           crit_twc   =  {detect_params.tm_crit:.4f}
           crit_thick =  {detect_params.thick_crit:.4f}
           crit_dist  =   {detect_params.dist_crit:.4f}
          lat_bound_n =  {detect_params.lat_bound_n:.4f}
          lat_bound_s = {detect_params.lat_bound_s:.4f}
          do_spline   = {".true." if detect_params.do_spline else ".false."}
          do_thickness= {".true." if detect_params.do_thickness else ".false."}
         &end
         &input
           fn_u    = '{os.path.join(input_dir, detect_params.u_in_file)}'
           fn_v    = '{os.path.join(input_dir, detect_params.v_in_file)}'
           fn_vort = '{os.path.join(input_dir, detect_params.vort_in_file)}'
           fn_tm   = '{os.path.join(input_dir, detect_params.tm_in_file)}'
           fn_slp  = '{os.path.join(input_dir, detect_params.slp_in_file)}'
           use_sfc_wnd = {".true." if detect_params.use_sfc_wind else ".false."}
         &end
        """)

        # Write the namelist content to the file
        with open(namelist_path, "w") as namelist_file:
            namelist_file.write(namelist_content)

        return namelist_path

    def _make_driver_call(self):
        """
        Construct a driver call based on options set in driver_parameters.

        Makes a call to write the namelist and then construct and call the driver
        command.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a driver command
            based on the parameters set in self.detect_parameters
        """
        tstorms_driver_dir = os.path.join(
            self.tstorms_parameters.tstorms_dir, "tstorms_driver"
        )
        dn_argslist = [os.path.join(tstorms_driver_dir, "tstorms_driver.exe")]

        return dn_argslist

    def detect(self, verbose=False):
        """
        Call the driver utility of TSTORMS.

        This will make a system call out to the tstorms_driver code from TSTORMS
        (provided it has been installed as an external dependency). This will be
        run according to the parameters in the :attr:`detect_parameters` attribute
        that were set when the :class:`TSTORMSTracker` instance was created.

        The output file is a plain text file named ``cyclones`` containing each of the
        TC candidates at each time from the input files. This will be generated
        in the current directory before being copied over to the desired output location
        provided in the base parameters attribute :attr:`tstorms_parameters`.
        Cyclones in the file are listed for each time in the format:

        .. code-block:: text

                <day>       <month>       <year>       <number>       <hour>
            <i_index>     <j_index>    <lon_slpmin>      <lat_slpmin>
           <max_wind>    <max_vort>    <min_slp>    <exist_warm_core> <exist_thickness>
           <max_warm_core>    <max_thickness>
                  ...

        - ``number`` is the number of nodes at that time.
        - ``i_index``, ``j_index`` are the grid indices of the node.
        - ``lon_slpmin``, ``lat_slpmin`` are the coordinates of the slp minimum.
        - Other values are the variable values for the cyclone output from the code.

        Parameters
        ----------
        verbose : bool
            Whether to print the entire TSTORMS output to screen in real-time or just
            the start/end summary.
            Defaults to False.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the tstorms_driver executeable from TSTORMS cannot be found.
        RuntimeError
            If the tstorms_driver executeable from TSTORMS returns a non-zero exit code.

        Examples
        --------
        To set the parameters, instantiate a :class:`TSTORMSTracker` instance and run
        detect:

        >>> my_params = TSTORMSDetectParameters(...)
        >>> my_tracker = TSTORMSTracker(driver_parameters=my_params)
        >>> result = my_tracker.detect()
        """
        namelist_filepath = self._write_driver_namelist()
        driver_call_list = self._make_driver_call()
        process_output = self._run_tstorms_process(
            "Detect", driver_call_list, namelist_filepath, verbose=verbose
        )

        # Copy the cyclones file to the output directory
        output_dir = self.tstorms_parameters.output_dir
        cyclones_file = "cyclones"
        destination_file = os.path.join(output_dir, cyclones_file)
        if os.path.exists(cyclones_file):
            # Copy to output_dir if not cwd
            if os.path.abspath(cyclones_file) != os.path.abspath(destination_file):
                shutil.copy(cyclones_file, destination_file)
        else:
            warnings.warn("cyclones file not found.", stacklevel=2)

        return process_output

    def _write_trajectory_analysis_namelist(self) -> str:
        """Generate the namelist file for the trajectory_analysis routine.

        Returns
        -------
        str
            The full path to the generated namelist file.

        Raises
        ------
        FileNotFoundError
            If the `trajectory_analysis` directory does not exist.
        """
        # Ensure the tstorms_dir exists
        tstorms_trajectory_dir = os.path.join(
            self.tstorms_parameters.tstorms_dir, "trajectory_analysis"
        )
        if not os.path.exists(tstorms_trajectory_dir):
            err_msg = (
                f"TSTORMS trajectory directory '{tstorms_trajectory_dir}' "
                "does not exist."
            )
            raise FileNotFoundError(err_msg)

        # Define the namelist file path
        namelist_path = os.path.join(tstorms_trajectory_dir, "nml_traj")

        # Format the namelist content
        stitch_params = self.stitch_parameters
        namelist_content = textwrap.dedent(f"""
         &input
            rcrit      = {stitch_params.r_crit:.4E}
            wcrit      = {stitch_params.wind_crit:.4f}
            vcrit      = {stitch_params.vort_crit:.4f}
            twc_crit   = {stitch_params.tm_crit:.4f}
            thick_crit =   {stitch_params.thick_crit:.4f}
            nwcrit     = {stitch_params.n_day_crit:.4f}
            do_filt    = {".true." if stitch_params.do_filter else ".false."}
            nlat = {stitch_params.lat_bound_n:.4f}
            slat = {stitch_params.lat_bound_s:.4f}
            do_spline    = {".true." if stitch_params.do_spline else ".false."}
            do_thickness = {".true." if stitch_params.do_thickness else ".false."}
         &end
        """)

        # Write the namelist content to the file
        with open(namelist_path, "w") as namelist_file:
            namelist_file.write(namelist_content)

        return namelist_path

    def _make_trajectory_analysis_call(self):
        """
        Construct a trajectory_analysis call based on options set in stitch_parameters.

        Makes a call to write the namelist and then construct and call the
        trajectory_analysis command.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a trajectory_analysis command
            based on the parameters set in self.stitch_parameters
        """
        tstorms_trajectory_analysis_dir = os.path.join(
            self.tstorms_parameters.tstorms_dir, "trajectory_analysis"
        )
        stitch_argslist = [
            os.path.join(tstorms_trajectory_analysis_dir, "trajectory_analysis_csc.exe")
        ]

        return stitch_argslist

    def stitch(self, verbose=False):
        """
        Call the trajectory analysis utility of TSTORMS to stitch candidate storms.

        This will make a system call out to the trajectory_analysis code from TSTORMS
        (provided it has been installed as an external dependency). This will be
        run according to the parameters in the :attr:`stitch_parameters` attribute
        that were set when the :class:`TSTORMSTracker` instance was created.
        It assumes that the candidate storms are contained in the current working
        directory in a file ``cyclones`` text file.

        The outputs are plain text files named:
            - ``ori`` containing the origins of each storm in the form
              ``lon, lat, YY, MM, DD, HH``,
            - ``traj`` containing trajectory point data in the form
              ``lon, lat, wind, psl, YY, MM, DD, HH`` for each trajectory,
            - ``trav`` containing trajectory point data and vorticity in the form
              ``lon, lat, wind, psl, vort_max, YY, MM, DD, HH`` for each trajectory,
            - ``stats`` containing information about the number of storms per month/year
              in each basin.
        If filtering is applied (:attr:`stitch_parameters` ``do_filter``) there will be
        additional files ``ori_filt``, ``traj_filt``, and ``trav_filt`` from after this
        takes place.
        These will be generated in the current directory before being copied over to
        the desired output location provided in the base parameters attribute
        :attr:`tstorms_parameters`.

        Parameters
        ----------
        verbose : bool
            Whether to print the entire TSTORMS output to screen in real-time or just
            the start/end summary.
            Defaults to False.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the trajectory_analysis executeable from TSTORMS cannot be found.
        RuntimeError
            If the trajectory_analysis executeable from TSTORMS returns a non-zero
            exit code.

        Examples
        --------
        To set the parameters, instantiate a :class:`TSTORMSTracker` instance and run
        stitch:

        >>> my_params = TSTORMSStitchParameters(...)
        >>> my_tracker = TSTORMSTracker(stitch_parameters=my_params)
        >>> result = my_tracker.stitch()
        """
        # Check if the cyclones file exists before proceeding
        cyclones_file = "cyclones"
        if not os.path.exists(cyclones_file):
            err_msg = (
                "No cyclones file found in the current directory. "
                "Did you run `detect` or remember to place the file locally?"
            )
            raise FileNotFoundError(err_msg)

        namelist_filepath = self._write_trajectory_analysis_namelist()
        trajectory_call_list = self._make_trajectory_analysis_call()
        process_output = self._run_tstorms_process(
            "Stitch", trajectory_call_list, namelist_filepath, verbose=verbose
        )

        # Copy the cyclones file to the output directory
        output_dir = self.tstorms_parameters.output_dir
        output_files = ["ori", "trav", "traj", "stats"]
        if self.stitch_parameters.do_filter:
            output_files.extend(["ori_filt", "trav_filt", "traj_filt"])

        for output_file in output_files:
            destination_file = os.path.join(output_dir, output_file)
            if os.path.exists(output_file):
                # Copy to output_dir if not cwd
                if os.path.abspath(output_file) != os.path.abspath(destination_file):
                    shutil.copy(output_file, destination_file)
            else:
                warnings.warn(f"output file `{output_file}` not found.", stacklevel=2)

        return process_output

    def set_metadata(self):
        """Create placeholder for abstract method."""
        super().set_metadata()

        tstorms_params_json = json.dumps(asdict(self.tstorms_parameters))
        detect_params_json = json.dumps(asdict(self.detect_parameters))
        stitch_params_json = json.dumps(asdict(self.stitch_parameters))
        self.global_metadata["tstorms_parameters"] = tstorms_params_json
        self.global_metadata["detect_parameters"] = detect_params_json
        self.global_metadata["stitch_parameters"] = stitch_params_json

    def run_tracker(self):
        """Create placeholder for abstract method."""

    def trajectories(self):
        """Create placeholder for abstract method."""
