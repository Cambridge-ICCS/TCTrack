"""Module providing classes that bind the TSTORMS code.

References
----------
- TSTORMS NOAA page: https://www.gfdl.noaa.gov/tstorms/
"""

import json
import os
import subprocess
import tempfile
import textwrap
import warnings
from dataclasses import asdict, dataclass

import cf
from cftime import num2date
from netCDF4 import Dataset

from tctrack.core import TCTracker, TCTrackerMetadata, TCTrackerParameters, Trajectory


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

    input_dir: str = ""
    """
    Full path to the directory where TSTORMS input files can be found.
    Defaults to empty string which will load from current directory.
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
        # Execute subproces commands from output_dir so outputs appear there.
        output_dir = self.tstorms_parameters.output_dir
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
                        cwd=output_dir,
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
                        cwd=output_dir,
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
        """
        Generate the namelist file for the tstorms_driver routine.

        The namelist file will be written to the ``output_dir`` specified in
        :attr:`tstorms_parameters`.

        Returns
        -------
        str
            The full path to the generated namelist file.

        Raises
        ------
        FileNotFoundError
            If the `tstorms_dir` does not exist.
        """
        # Ensure the output_dir exists to place namelist in
        output_dir = self.tstorms_parameters.output_dir
        if not os.path.exists(output_dir):
            err_msg = f"TSTORMS output directory '{output_dir}' does not exist."
            raise FileNotFoundError(err_msg)

        namelist_path = os.path.join(output_dir, "nml_driver")

        # Format the namelist content
        input_dir = self.tstorms_parameters.input_dir
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
        in the desired output location provided in the base parameters attribute
        :attr:`tstorms_parameters`.
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

        >>> base_params = TSTORMSBaseParameters(...)
        >>> detect_params = TSTORMSDetectParameters(...)
        >>> my_tracker = TSTORMSTracker(
                tstorms_parameters=base_params,
                detect_parameters=detect_params
            )
        >>> result = my_tracker.detect()
        """
        namelist_filepath = self._write_driver_namelist()
        driver_call_list = self._make_driver_call()
        process_output = self._run_tstorms_process(
            "Detect", driver_call_list, namelist_filepath, verbose=verbose
        )

        return process_output

    def _write_trajectory_analysis_namelist(self) -> str:
        """Generate the namelist file for the trajectory_analysis routine.

        The namelist file will be written to the ``output_dir`` specified in
        :attr:`tstorms_parameters`.

        Returns
        -------
        str
            The full path to the generated namelist file.

        Raises
        ------
        FileNotFoundError
            If the `trajectory_analysis` directory does not exist.
        """
        tstorms_dir = self.tstorms_parameters.tstorms_dir
        # Ensure the output_dir exists to place namelist in
        output_dir = self.tstorms_parameters.output_dir
        if not os.path.exists(output_dir):
            err_msg = f"TSTORMS output directory '{output_dir}' does not exist."
            raise FileNotFoundError(err_msg)

        # Define the namelist file path
        namelist_path = os.path.join(output_dir, "nml_traj")

        # Format the namelist content
        stitch_params = self.stitch_parameters
        namelist_content = textwrap.dedent(f"""
         &input
            rcrit      = {stitch_params.r_crit:.4f}
            wcrit      = {stitch_params.wind_crit:.4f}
            vcrit      = {stitch_params.vort_crit:.4E}
            twc_crit   = {stitch_params.tm_crit:.4f}
            thick_crit = {stitch_params.thick_crit:.4f}
            nwcrit     = {stitch_params.n_day_crit:.4f}
            do_filt    = {".true." if stitch_params.do_filter else ".false."}
            nlat = {stitch_params.lat_bound_n:.4f}
            slat = {stitch_params.lat_bound_s:.4f}
            do_spline    = {".true." if stitch_params.do_spline else ".false."}
            do_thickness = {".true." if stitch_params.do_thickness else ".false."}
            landmask = '{os.path.join(tstorms_dir, "trajectory_analysis/landsea.map")}'
            cmask = '{os.path.join(tstorms_dir, "trajectory_analysis/imask_2")}'
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
        It assumes that the candidate storms are contained in the output_dir
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
        These files will be generated in the desired output location provided in the
        base parameters attribute :attr:`tstorms_parameters`.

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

        >>> base_params = TSTORMSBaseParameters(...)
        >>> stitch_params = TSTORMSStitchParameters(...)
        >>> stitch_params = TSTORMSTracker(
                tstorms_parameters=base_params,
                stitch_parameters=stitch_params
            )
        >>> result = my_tracker.stitch()
        """
        # Check if the cyclones file exists before proceeding
        cyclones_file = os.path.join(self.tstorms_parameters.output_dir, "cyclones")
        if not os.path.exists(cyclones_file):
            err_msg = (
                "No cyclones file found in the output directory. "
                "Did you forget to run `detect` or delete the file?"
            )
            raise FileNotFoundError(err_msg)

        namelist_filepath = self._write_trajectory_analysis_namelist()
        trajectory_call_list = self._make_trajectory_analysis_call()
        process_output = self._run_tstorms_process(
            "Stitch", trajectory_call_list, namelist_filepath, verbose=verbose
        )

        return process_output

    def set_metadata(self) -> None:
        """
        Set the global and variable (reading from input files) metadata attributes.

        Reads metadata for each variable from the input NetCDF files
        defined in :attr:`detect_parameters`.
        These will be stored in the :attr:`variable_metadata` attribute as a
        dictionary of :class:`TCTrackerMetadata` objects.
        This will be called from the :meth:`to_netcdf` method.

        Raises
        ------
        ValueError
            If a variable is not found in the input files.

        Examples
        --------
        To generate in the metadata for variables from parameters and inputs:

        >>> tracker = TSTORMSTracker(tstorms_params, detect_params, stitch_params)
        >>> tracker.set_metadata()
        >>> tracker.variable_metadata
        {
            "wind_speed": TCTrackerMetadata(
                properties={
                "standard_name": "wind_speed",
                "long_name": "Wind Speed",
                ...
                },
            )
            ...
        }
        """
        super().set_metadata()

        tstorms_params_json = json.dumps(asdict(self.tstorms_parameters))
        detect_params_json = json.dumps(asdict(self.detect_parameters))
        stitch_params_json = json.dumps(asdict(self.stitch_parameters))
        self.global_metadata["tstorms_parameters"] = tstorms_params_json
        self.global_metadata["detect_parameters"] = detect_params_json
        self.global_metadata["stitch_parameters"] = stitch_params_json

        if not self._time_metadata:
            self._set_time_metadata()

        if not self._variable_metadata:
            self._set_variable_metadata()

    def _set_variable_metadata(self) -> None:
        """
        Extract variable metadata from inputs and return dict for setting.

        TSTORMS calculates windspeed as Euclidian norm of u and v, so read
        key units and attributes from u and assume they apply.
        Vorticity, Temperature, and Sea-level-pressure all read in directly from file.

        Note that there is some slight misdirection here in that wind, slp, and tm are
        detected relative to the vorticity location, but the slp location is used to
        provide coordinates. This means a situation could arise in which points are
        in fact from an area more than quoted (up to twice the radius).

        Raises
        ------
        FileNotFoundError
            If an input file cannot be found.
        ValueError
            If a variable is not found in the input files.
        """
        input_dir = self.tstorms_parameters.input_dir
        detect_params = self.detect_parameters
        var_outputs: list[dict] = [
            {
                "std_name": "wind_speed",
                "long_name": (
                    "Surface Wind Speed" if detect_params.use_sfc_wind else "Wind Speed"
                ),
                "tstorms_name": ("u_ref" if detect_params.use_sfc_wind else "u850"),
                "filename": os.path.join(input_dir, detect_params.u_in_file),
                "cellmethod": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={
                        "comment": (
                            f"lesser circle of radius {detect_params.dist_crit} degrees"
                        ),
                    },
                ),
            },
            {
                "std_name": "air_pressure_at_mean_sea_level",
                "long_name": "Air Pressure at Mean Sea Level",
                "tstorms_name": "slp",
                "filename": os.path.join(input_dir, detect_params.slp_in_file),
                "cellmethod": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={
                        "comment": (
                            f"lesser circle of radius {detect_params.dist_crit} degrees"
                        ),
                    },
                ),
            },
            {
                "std_name": "atmosphere_upward_relative_vorticity",
                "long_name": "Atmosphere Upward Relative Vorticity",
                "tstorms_name": "vort850",
                "filename": os.path.join(input_dir, detect_params.vort_in_file),
                "cellmethod": cf.CellMethod(
                    "area",
                    "maximum",
                    qualifiers={
                        "comment": (
                            f"lesser circle of radius {detect_params.dist_crit} degrees"
                        ),
                    },
                ),
            },
        ]

        # Initialise variable metadata as empty dict to populate
        self._variable_metadata = {}

        for output in var_outputs:
            # Check input file exists to read from
            if not os.path.exists(output["filename"]):
                errmsg = f"Input file '{output['filename']}' not found."
                raise FileNotFoundError(errmsg)

            var_name = output["std_name"]
            tstorms_name = output["tstorms_name"]

            # Get the variable field from the netcdf file
            fields = cf.read(output["filename"], select=f"ncvar%{tstorms_name}")  # type: ignore[operator]
            if not fields:
                msg = (
                    f"Variable '{tstorms_name}' not found "
                    f"in input file {output['filename']}."
                )
                raise ValueError(msg)
            field = fields[0]

            # Read and store the relevant metadata
            # Override wind speed as reading from u file but this is zonal velocity
            self._variable_metadata[var_name] = TCTrackerMetadata(
                {
                    "standard_name": (
                        output["std_name"]
                        if var_name == "wind_speed"
                        else field.get_property("standard_name", output["std_name"])
                    ),
                    "long_name": (
                        output["long_name"]
                        if var_name == "wind_speed"
                        else field.get_property("long_name", output["long_name"])
                    ),
                    "units": field.get_property("units", "unknown"),
                }
            )
            self._variable_metadata[var_name].constructs = [output["cellmethod"]]

    def _set_time_metadata(self):
        """
        Extract time metadata from u_ref input and set attribute.

        Calendar type, units, start time, and end time extracted from u_ref as for
        TSTORMS by identifying the unlimited dimension and fetching coordinate
        attributes. Assumes other files match.
        If file not found or variable cannot be identified errors are raised.
        The calendar defaults to Julian as for TSTORMS.

        Raises
        ------
        KeyError
            If no unlimited dimension or corresponding coordinate is found in the u
            input file.
        ValueError
            If multiple unlimited dimensions are found in the u input file.
        ValueError
            If the 'units' attribute is missing for the time coordinate.

        Warnings
        --------
        UserWarning
            If the 'calendar' attribute is missing for the coordinate variable.
            Defaults to 'julian'.
        """
        input_u_file = os.path.join(
            self.tstorms_parameters.input_dir, self.detect_parameters.u_in_file
        )

        with Dataset(input_u_file, "r") as nc_file:
            unlimited_dims = [
                dim
                for dim in nc_file.dimensions
                if nc_file.dimensions[dim].isunlimited()
            ]
            if not unlimited_dims:
                errmsg = (
                    "No unlimited dimension found in the u_ref NetCDF file "
                    "to set calendar."
                )
                raise KeyError(errmsg)
            if len(unlimited_dims) > 1:
                errmsg = (
                    "TSTORMS expects only a single unlimited variable in NetCDF "
                    "files corresponding to the time dimension. "
                    f"Multiple found: {unlimited_dims}."
                )
                raise ValueError(errmsg)

            unlimited_dim = unlimited_dims[0]
            if unlimited_dim not in nc_file.variables:
                errmsg = (
                    "Coordinate variable for unlimited dimension "
                    f"'{unlimited_dim}' not found."
                )
                raise KeyError(errmsg)

            coord_var = nc_file.variables[unlimited_dim]
            units = getattr(coord_var, "units", None)
            if units is None:
                msg = "The 'units' attribute is required for the time coordinate."
                raise ValueError(msg)

            calendar = getattr(coord_var, "calendar", None)
            if calendar is None:
                msg = (
                    "The 'calendar' attribute is missing for the coordinate variable.\n"
                    "defaulting to 'julian'"
                )
                warnings.warn(msg, category=UserWarning, stacklevel=2)
                calendar = "julian"

            start_time_num = coord_var[0]
            end_time_num = coord_var[-1]

            self._time_metadata = {
                "calendar": calendar,
                "units": units,
                "start_time": num2date(start_time_num, units=units, calendar=calendar),
                "end_time": num2date(end_time_num, units=units, calendar=calendar),
            }

    def trajectories(self):
        """
        Parse outputs from TSTORMS to list of :class:`tctrack.core.Trajectory`.

        This will read from the ``trav`` file output from trajectory stitching, or
        ``trav_filt`` if filtering was applied (see ``do_filter`` in the
        :attr:`tstorms_parameters` attribute).

        Returns
        -------
        trajectories : list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        trajectories = []
        current_trajectory_id = 0  # Initialize trajectory ID

        # We need time metadata, so set the metadata if not done already
        if not self._time_metadata:
            self._set_time_metadata()

        # Use the trav file to get full data (including vorticity)
        # If filtering was applied use filtered file
        output_dir = self.tstorms_parameters.output_dir
        trav_file = "trav_filt" if self.stitch_parameters.do_filter else "trav"
        output_trav_filepath = os.path.join(output_dir, trav_file)

        with open(output_trav_filepath, "r") as file:
            for line in file:
                items = line.split()
                if items[0] == "start":
                    # Start of new trajectory.
                    # Extract metadata and add Trajectory to dict
                    current_trajectory_id += 1
                    time = list(map(int, items[2:6]))

                    trajectories.append(
                        Trajectory(
                            current_trajectory_id,
                            time,
                            calendar=self.time_metadata["calendar"],
                        )
                    )

                # Continue processing ongoing trajectory
                # Trajectories indexed from 1 so subtract when assigning
                else:
                    trajectories[current_trajectory_id - 1].add_point(
                        *self._parse_tstorms_trav_line_to_point(items)
                    )
        return trajectories

    @staticmethod
    def _parse_tstorms_trav_line_to_point(
        line: list[str], variable_names: list[str] | None = None
    ) -> tuple[list[int], dict[str, int | float]]:
        """
        Parse line from TSTORMS trav output into a trajectory data point.

        Data point format is that expected by a :class:`tctrack.core.Trajectory`.

        Data from a traj file can be parsed by providing a different set of
        ``variable_names``

        Parameters
        ----------
        line : list[str]
            A list of strings representing the line split into parts.
        variable_names : list[str] | None
            List of variable names for the data columns.
            Defaults to None and then set to those for a trav file.

        Returns
        -------
        tuple
            A tuple containing the time as an integer list of [year, day, month, hour]
            and a dict of variables.
        """
        if not variable_names:
            variable_names = [
                "lon",
                "lat",
                "wind_speed",
                "air_pressure_at_mean_sea_level",
                "atmosphere_upward_relative_vorticity",
            ]

        return_vars: dict[str, int | float] = {}
        return_vars.update(
            {
                name: float(value)
                for name, value in zip(variable_names, line[:-4], strict=True)
            }
        )
        time = list(map(int, line[-4:]))
        return time, return_vars

    def run_tracker(self, output_file: str):
        """Run TSTORMS tracker to obtain tropical cyclone track trajectories.

        This first runs :meth:`detect` to get TC candidates at each time. Then
        these are combined into trajectories using :meth:`stitch`.
        The output is then saved as a CF-compliant NetCDF trajectory file.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone trajectories are saved.

        Raises
        ------
        FileNotFoundError
            - If the TSTORMS executables cannot be found.
        RuntimeError
            If the TSTORMS commands return a non-zero exit code.

        Examples
        --------
        To set the parameters, instantiate a :class:`TSTORMSTracker` instance and run
        `run_tracker()`:

        >>> tstorms_params = TSTORMSBaseParameters(...)
        >>> detect_params = TSTORMSDetectParameters(...)
        >>> stitch_params = TSTORMSStitchParameters(...)
        >>> my_tracker = TSTORMSTracker(tstorms_params, detect_params, stitch_params)
        >>> my_tracker.run_tracker()
        """
        self.detect()
        self.stitch()
        self.to_netcdf(output_file)
