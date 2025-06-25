"""Module providing a class that binds the Tempest Extremes code.

References
----------
- Tempest Extremes code on GitHub: https://github.com/ClimateGlobalChange/tempestextremes
- Tempest Extremes User Guide: https://climate.ucdavis.edu/tempestextremes.php
- GMD paper on Tempest Extremes v2.1: https://doi.org/10.5194/gmd-14-5023-2021
"""

import csv
import subprocess
from dataclasses import asdict, dataclass
from typing import TypedDict

from cftime import Datetime360Day, DatetimeGregorian, DatetimeNoLeap


class Track:
    """
    Represents a single Lagrangian cyclone track with metadata and data points.

    Attributes
    ----------
    track_id : int
        The unique identifier for the track.
    num_points : int
        Number of points in the track.
    calendar : str
        The calendar type to use for datetime handling.
        Options are "gregorian", "360_day", or "noleap".
    start_time : Datetime360Day | DatetimeNoLeap | DatetimeGregorian
        Start time of the track as a datetime or cftime object.
    data : list
        List of data points, where each point is a dictionary containing:
            grid_i, grid_j, variables, timestamp
    """

    def __init__(
        self,
        track_id: int,
        num_points: int,
        year: int,
        month: int,
        day: int,
        hour: int,
        calendar: str = "gregorian",
    ):
        """
        Initialize a Track object.

        Parameters
        ----------
        track_id : int
            The unique identifier for the track.
        num_points : int
            The number of points in the track.
        year : int
            The starting year of the track.
        month : int
            The starting month of the track (1-12).
        day : int
            The starting day of the track (1-31, depending on the calendar).
        hour : int
            The starting hour of the track (0-23).
        calendar : str, optional
            The calendar type to use for datetime handling. Options are
            "gregorian", "360_day", or "noleap".
        """
        self.track_id = track_id
        self.num_points = num_points
        self.calendar = calendar
        self.start_time = self._create_datetime(year, month, day, hour)
        self.data: list = []

    def _create_datetime(self, year: int, month: int, day: int, hour: int):
        """
        Create a cftime object based on the specified calendar attribute.

        Parameters
        ----------
        year : int
            year as an integer
        month : int
            month as an integer
        day : int
            day as an integer
        hour : int
            hour as an integer

        Returns
        -------
        datetime : Datetime360Day | DatetimeNoLeap | DatetimeGregorian
        """
        if self.calendar == "360_day":
            return Datetime360Day(year, month, day, hour)
        elif self.calendar == "noleap":
            return DatetimeNoLeap(year, month, day, hour)
        elif self.calendar in {"gregorian", "standard"}:
            return DatetimeGregorian(year, month, day, hour)
        else:
            msg = (
                f"Unsupported calendar type: {self.calendar}. "
                "Supported types are '360_day', 'noleap', or 'gregorian'/'standard'."
            )
            raise ValueError(msg)

    def __str__(self) -> str:
        """Improve the representation of Track to users."""
        return (
            f"Track(num_points={self.num_points}, "
            f"start_time={self.start_time}, "
            f"calendar={self.calendar}, "
            f"data_points={len(self.data)})"
        )

    def add_point(self, point: tuple):
        """
        Add a data point to the track.

        Parameters
        ----------
        point : tuple
            A tuple containing the following elements:
            - grid_i : int
                The grid index i.
            - grid_j : int
                The grid index j.
            - variables : dict
                A dictionary of variable names and their corresponding values.
            - year : int
                The year of the data point.
            - month : int
                The month of the data point (1-12).
            - day : int
                The day of the data point (1-31, depending on the calendar).
            - hour : int
                The hour of the data point (0-23).
        """
        grid_i, grid_j, variables, year, month, day, hour = point
        timestamp = self._create_datetime(year, month, day, hour)
        self.data.append(
            {
                "grid_i": grid_i,
                "grid_j": grid_j,
                "variables": variables,
                "timestamp": timestamp,
            }
        )


def lod_to_te(inputs: list[dict]) -> str:
    """
    Convert a sane list of dicts to input format used by Tempest Extremes commands.

    The input list of dictionaries will be concatenated to a single string of values
    separated by `,` and dicts separated by `;`.

    Parameters
    ----------
    inputs : list[dict]
        List of the dictionaries to be concatenated.

    Returns
    -------
    str
        Single string of values separated by `,` and dicts separated by `;`

    Examples
    --------
    >>> lod_to_te([{"a": 1, "b": 2, "c": 3}, {"d": 4, "e": 5, "f": 6}])
    '1,2,3;4,5,6'
    """
    return ";".join(",".join(str(value) for value in d.values()) for d in inputs)


class TEContour(TypedDict):
    """
    Data required for checking a contour of a single variable in DetectNodes.

    Points will be eliminated in a DetectNodes search if they fail this criterion.
    The closed contour is determined by breadth first search: if any paths exist from
    the candidate point (or nearby minima/maxima if minmaxdist is specified) that
    reach the specified distance before achieving the specified delta then we say no
    closed contour is present.
    Each contour takes the form of a ``dict`` with keys ``"var"``, ``"delta"``,
    ``"dist"``, and ``"minmaxdist"``.

    See Also
    --------
    TETracker.detect_nodes : The Detect Nodes call from the TETracker object
    DetectNodesParameters : The DetectNodes parameter class

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`_
    and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_

    Examples
    --------
    To add a contour requirement on ``"psl"`` with a change of ``200.0`` within
    ``5.5`` degrees of the candidate we create a TEContour as follows:

    >>> TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0)
    {'var': 'psl', 'delta': 200.0, 'dist': 5.5, 'minmaxdist': 0.0}
    """

    var: str
    """Name of the variable to contour in NetCDF files."""

    delta: float
    """
    Amount by which the field must change from the pivot value.
    If positive (negative) the field must increase (decrease) by this value along
    the contour.
    """

    dist: float
    """
    Great-circle distance (degrees) from the pivot within which the
    criteria must be satisfied.
    """

    minmaxdist: float
    """
    Great-circle distance away from the candidate to search for the minima/maxima.
    If delta is positive (negative), the pivot is a local minimum (maximum).
    """


class TEOutputCommand(TypedDict):
    """
    Data required to specify an additional column in the DetectNodes output.

    Each output command takes the form of a ``dict`` with keys ``"var"``, ``"op"``, and
    ``"dist"``.

    See Also
    --------
    TETracker.detect_nodes : The Detect Nodes call from the TETracker object
    DetectNodesParameters : The DetectNodes parameter class

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`_
    and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_

    Examples
    --------
    To add an output column with the minimum ``"psl"`` at the candidate point:

    >>> TEOutputCommand(var="psl", operator="min", dist=0.0)
    {'var': 'psl', 'operator': 'min', 'dist': 0.0}
    """

    var: str
    """Name of the variable to write in output files."""

    operator: str
    """
    Operator that is applied over all points within the specified distance of the
    candidate (options include ``"max"``, ``"min"``, ``"avg"``, ``"maxdist"``, and
    ``"mindist"``).
    """

    dist: float
    """
    Great-circle distance (degrees) from the candidate within which the
    operator is applied.
    """


class TEThreshold(TypedDict):
    """Data required for a threshold filter for a track in StitchNodes.

    Any tracks that do not satisfy the threshold value for a given number of points will
    be filtered out. Each condition is of the form "var,op,value,count" and conditions
    are separated by ";".

    See Also
    --------
    TETracker.stitch_nodes : The StitchNodes call from the TETracker object.
    StitchNodesParameters : The StitchNodes parameter class.

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`_
    and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_

    Examples
    --------
    To add a filter requiring latitude ``"lat"`` to be less than 40 degrees for
    10 or more points in each track:

    >>> TEThreshold(var="lat", op="<=", value=40, count=10)
    {"var": "lat", "op": "<=", "value": 40, "count": 10},
    """

    var: str
    """Name of the variable being tested. Called "col" in TempestExtremes."""

    op: str
    """Operator used for the comparison (options include >,>=,<,<=,=,!=,|>=,|<=)."""

    value: float
    """Value on the right-hand-side of the comparison."""

    count: int | str
    """
    Either the minimum number of points where the threshold must be satisfied or the
    instruction ``"all"``, ``"first"``, or ``"last"``. ``"all"`` for all points along
    the path, ``"first"`` for just the first point, and ``"last"`` for only the last
    point.
    """


@dataclass
class DetectNodesParameters:
    """
    Dataclass containing values used by the DetectNodes operation of TE.

    See Also
    --------
    TEContour : The class used to define contour criteria
    TEOutputCommand : The class used to define additional outputs

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`_
    and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_
    """

    in_data: list[str] | None = None
    """List of strings of NetCDF input files."""

    out_header: bool = False
    """Include header at the top of the output file?"""

    output_file: str | None = None
    """Output nodefile to write to."""

    search_by_min: str | None = None
    """
    Input variable in NetCDF files for selecting candidate points (defined as local
    minima). If ``None``, then uses ``"PSL"`` in Tempest Extremes.
    """

    search_by_max: str | None = None
    """
    Input variable in NetCDF files for selecting candidate points (defined as local
    maxima).
    """

    closed_contours: list[TEContour] | None = None
    """
    Criteria for candidates to be eliminated if they do not have a closed contour
    as a list of separate :class:`TEContour` criteria.
    """

    merge_dist: float = 0.0
    """
    DetectNodes merges candidate points with a distance (in degrees
    great-circle-distance) shorter than the specified value. Among two candidates
    within the merge distance, only the candidate with the lowest value of the
    search_by_min field or highest value of the search_by_max field are retained.
    """

    lat_name: str = "lat"
    """String for the latitude dimension in the NetCDF files."""

    lon_name: str = "lon"
    """String for the longitude dimension in the NetCDF files."""

    min_lat: float = 0.0
    """Minimum latitude for candidate points."""

    max_lat: float = 0.0
    """
    Maximum latitude for candidate points.
    If max_lat and min_lat are equal then these arguments are ignored.
    """

    min_lon: float = 0.0
    """Minimum longitude for candidate points."""

    max_lon: float = 0.0
    """
    Maximum longitude for candidate points.
    If ``max_lon`` and ``min_lon`` are equal then these arguments are ignored.
    """

    regional: bool = False
    """Should lat-lon grid be periodic in longitude."""

    output_commands: list[TEOutputCommand] | None = None
    """
    Criteria for any additional columns to be added to the output. Criteria are provided
    as a list of separate :class:`TEOutputCommand` criteria.
    """

    def __str__(self) -> str:
        """Improve the representation of DetectNodesParameters to users."""
        attributes = "\n\t".join(
            f"{key} \t = {value}" for key, value in asdict(self).items()
        )
        return f"DetectNodesParameters(\n\t{attributes}\n)"


@dataclass
class StitchNodesParameters:
    """Dataclass containing values used by the StitchNodes operation of TE.

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`_
    and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_
    """

    output_file: str = "tracks.txt"
    """The output filename to save the tracks to. Called "out" in TempestExtremes."""

    in_file: str | None = None
    """
    Filename of the DetectNodes output file. If this and ``in_list`` are ``None``, it
    will be taken from :attr:`DetectNodesParameters.output_file`. Called "in" in
    TempestExtremes.
    """

    in_list: str | None = None
    """
    File containing a list of input files to be processed together. This is unadvised to
    use at present as it is likely to be changed.
    """

    in_fmt: str | None = None
    """
    Comma-separated list of the variables in the input file. If ``None``, it will be
    ``"lon,lat"`` and any others defined in
    :attr:`DetectNodesParameters.output_commands`.
    """

    allow_repeated_times: bool = False
    """
    If ``False``, an error is thrown if there are multiple sections in the input
    nodefile with the same time.
    """

    caltype: str = "standard"
    """
    The type of calendar to use. Options are ``"standard"`` (365 days with leap years),
    ``"noleap"``, ``"360_day"``.
    """

    time_begin: str | None = None
    """Starting date / time for stitching tracks. Earlier times will be ignored."""

    time_end: str | None = None
    """Ending date / time for stitching tracks. Later times will be ignored."""

    max_sep: float = 5
    """
    The maximum distance allowed between candidates (degrees). Called "range" in
    TempestExtremes.
    """

    max_gap: int = 0
    """The number of missing points allowed between candidates."""

    min_time: int | str = 1
    """
    The minimum required length of a path. Either as an integer for the number of
    candidates, or a string for total duration, e.g. ``"24h"``.
    """

    min_endpoint_dist: float = 0
    """The minimum required distance between the first and last candidates (degrees)."""

    min_path_dist: float = 0
    """The minimum required acumulated distance along the path (degrees)."""

    threshold_filters: list[TEThreshold] | None = None
    """
    Filters for paths based on the number of nodes that satisfy a threshold. Uses a list
    of :class:`TEThreshold` objects.  Called "thresholdcmd" in TempestExtremes.
    """

    prioritize: str | None = None
    """
    The variable to use to determine the precedence (lowest to highest) of nodes for
    matching to the next position.
    """

    add_velocity: bool = False
    """
    Whether to include the velocity components (m/s) of the movement of the TC to the
    output file.
    """

    out_file_format: str = "gfdl"
    """
    Format of the output file. ``"gfdl"``, ``"csv"``, or ``"csvnoheader"``.
    See :meth:`TETracker.stitch_nodes` for details.
    """

    out_seconds: bool = False
    """
    For GFDL output file types, determines whether to report the sub-daily time in
    seconds (``True``) or hours (``False``).
    """

    def __str__(self) -> str:
        """Improve the representation to users."""
        attributes = "\n\t".join(
            f"{key} \t = {value}" for key, value in asdict(self).items()
        )
        return f"StitchNodesParameters(\n\t{attributes}\n)"


class TETracker:
    """Class containing bindings to the Tempest Extremes code.

    Attributes
    ----------
    detect_nodes_parameters : DetectNodesParameters
        Class containing the parameters for the DetectNodes algorithm
    stitch_nodes_parameters : StitchNodesParameters | None
        Class containing the parameters for the StitchNodes algorithm
    """

    def __init__(
        self,
        detect_nodes_parameters: DetectNodesParameters | None = None,
        stitch_nodes_parameters: StitchNodesParameters | None = None,
    ):
        """
        Construct the TempestExtremes class.

        Parameters
        ----------
        detect_nodes_parameters : DetectNodesParameters
            Class containing the parameters for the DetectNodes algorithm
            Defaults to the default values in DetectNodesParameters Class
        stitch_nodes_parameters : StitchNodesParameters | None
            Class containing the parameters for the StitchNodes algorithm
            Defaults to the default values in StitchNodesParameters Class
        """
        if detect_nodes_parameters is not None:
            self.detect_nodes_parameters: DetectNodesParameters = (
                detect_nodes_parameters
            )
        else:
            self.detect_nodes_parameters = DetectNodesParameters()

        if stitch_nodes_parameters is not None:
            self.stitch_nodes_parameters: StitchNodesParameters = (
                stitch_nodes_parameters
            )
        else:
            self.stitch_nodes_parameters = StitchNodesParameters("tracks.txt")

        # Set StitchNodes input arguments according to DetectNodes parameters,
        # if not provided
        sn_params = self.stitch_nodes_parameters
        dn_params = self.detect_nodes_parameters
        sn_input_none = sn_params.in_file is None and sn_params.in_list is None
        if sn_input_none and dn_params.output_file is not None:
            sn_params.in_file = dn_params.output_file
        if sn_params.in_fmt is None and dn_params.output_commands is not None:
            variables = [output["var"] for output in dn_params.output_commands]
            sn_params.in_fmt = ",".join(["lon", "lat", *variables])

    def _run_te_process(self, command_name: str, command_list: list[str]):
        """Run a TempestExtremes command (DetectNodes or StitchNodes).

        Parameters
        ----------
        command_name : str
            The name of the command to be used in the log and error reporting.
        command_list : list[str]
            The list of strings that produce the command as given by
            _make_detect_nodes_call and _make_stitch_nodes_call.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the DetectNodes executeable from TempestExtremes cannot be found.
        RuntimeError
            If Tempest Extremes DetectNodes returns a non-zero exit code.
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

    def _make_detect_nodes_call(self):  # noqa: PLR0912 - all branches same logic
        """
        Construct a DetectNodes call based on options set in parameters.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a DetectNodes command
            based on the parameters set in self.detect_nodes_parameters
        """
        dn_argslist = ["DetectNodes"]

        if self.detect_nodes_parameters.in_data is not None:
            dn_argslist.extend(
                [
                    "--in_data",
                    ";".join(self.detect_nodes_parameters.in_data),
                ]
            )
        if self.detect_nodes_parameters.out_header:
            dn_argslist.extend(["--out_header"])
        if self.detect_nodes_parameters.output_file is not None:
            dn_argslist.extend(["--out", self.detect_nodes_parameters.output_file])
        if self.detect_nodes_parameters.search_by_min is not None:
            dn_argslist.extend(
                [
                    "--searchbymin",
                    self.detect_nodes_parameters.search_by_min,
                ]
            )
        if self.detect_nodes_parameters.search_by_max is not None:
            dn_argslist.extend(
                [
                    "--searchbymax",
                    self.detect_nodes_parameters.search_by_max,
                ]
            )
        if self.detect_nodes_parameters.closed_contours is not None:
            dn_argslist.extend(
                [
                    "--closedcontourcmd",
                    lod_to_te(self.detect_nodes_parameters.closed_contours),
                ]
            )
        dn_argslist.extend(
            [
                "--mergedist",
                str(self.detect_nodes_parameters.merge_dist),
            ]
        )
        if self.detect_nodes_parameters.lat_name is not None:
            dn_argslist.extend(
                [
                    "--latname",
                    self.detect_nodes_parameters.lat_name,
                ]
            )
        if self.detect_nodes_parameters.lon_name is not None:
            dn_argslist.extend(
                [
                    "--lonname",
                    self.detect_nodes_parameters.lon_name,
                ]
            )
        if self.detect_nodes_parameters.min_lat is not None:
            dn_argslist.extend(
                [
                    "--minlat",
                    str(self.detect_nodes_parameters.min_lat),
                ]
            )
        if self.detect_nodes_parameters.max_lat is not None:
            dn_argslist.extend(
                [
                    "--maxlat",
                    str(self.detect_nodes_parameters.max_lat),
                ]
            )
        if self.detect_nodes_parameters.min_lon is not None:
            dn_argslist.extend(
                [
                    "--minlon",
                    str(self.detect_nodes_parameters.min_lon),
                ]
            )
        if self.detect_nodes_parameters.max_lon is not None:
            dn_argslist.extend(
                [
                    "--maxlon",
                    str(self.detect_nodes_parameters.max_lon),
                ]
            )
        if self.detect_nodes_parameters.regional:
            dn_argslist.extend(["--regional"])
        if self.detect_nodes_parameters.output_commands is not None:
            dn_argslist.extend(
                [
                    "--outputcmd",
                    lod_to_te(self.detect_nodes_parameters.output_commands),
                ]
            )

        return dn_argslist

    def detect_nodes(self):
        """
        Call the DetectNodes utility of Tempest Extremes.

        This will make a system call out to the DetectNodes method from Tempest Extremes
        (provided it has been installed as an external dependency). DetectNodes will be
        run according to the parameters in the :attr:`detect_nodes_parameters` attribute
        that were set when the :class:`TETracker` instance was created.

        The output file is a plain text file containing each of the TC candidates at
        each time from the input files. If :attr:`~DetectNodesParameters.out_header` is
        ``True`` the first two lines will be a header describing the structure of the
        data. After this each time is listed in the format:

        .. code-block:: text

           <year> <month> <day> <count> <hour>
                  <i> <j> <lon> <lat> <var1> <var2> ...
                  ...
                  <i> <j> <lon> <lat> <var1> <var2> ...

        - ``count`` is the number of nodes at that time.
        - ``i``, ``j`` are the grid indices of the node.
        - ``var1``, ``var2``, etc., are scalar variables as defined by
          :attr:`~DetectNodesParameters.output_commands` (typically, psl, orog).

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the DetectNodes executeable from TempestExtremes cannot be found.
        RuntimeError
            If Tempest Extremes DetectNodes returns a non-zero exit code.

        References
        ----------
        `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`_
        and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_

        Examples
        --------
        To set the parameters, instantiate a :class:`TETracker` instance and run
        DetectNodes:

        >>> my_params = DetectNodesParameters(...)
        >>> my_tracker = TETracker(detect_nodes_parameters=my_params)
        >>> result = TETracker.detect_nodes()
        """
        dn_call_list = self._make_detect_nodes_call()
        return self._run_te_process("DetectNodes", dn_call_list)

    def _make_stitch_nodes_call(self):
        """
        Construct a StitchNodes call based on options set in parameters.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a StitchNodes command
            based on the parameters set in self.stitch_nodes_parameters
        """
        sn_argslist = ["StitchNodes"]
        sn_params = self.stitch_nodes_parameters
        if sn_params.output_file is not None:
            sn_argslist.extend(["--out", sn_params.output_file])
        if sn_params.in_file is not None:
            sn_argslist.extend(["--in", sn_params.in_file])
        if sn_params.in_list is not None:
            sn_argslist.extend(["--in_list", sn_params.in_list])
        if sn_params.in_fmt is not None:
            sn_argslist.extend(["--in_fmt", sn_params.in_fmt])
        if sn_params.allow_repeated_times:
            sn_argslist.extend(["--allow_repeated_times"])
        sn_argslist.extend(["--caltype", str(sn_params.caltype)])
        if sn_params.time_begin is not None:
            sn_argslist.extend(["--time_begin", str(sn_params.time_begin)])
        if sn_params.time_end is not None:
            sn_argslist.extend(["--time_end", str(sn_params.time_end)])
        sn_argslist.extend(["--range", str(sn_params.max_sep)])
        sn_argslist.extend(["--maxgap", str(sn_params.max_gap)])
        sn_argslist.extend(["--mintime", str(sn_params.min_time)])
        sn_argslist.extend(["--min_endpoint_dist", str(sn_params.min_endpoint_dist)])
        sn_argslist.extend(["--min_path_dist", str(sn_params.min_path_dist)])
        if sn_params.threshold_filters is not None:
            sn_argslist.extend(["--threshold", lod_to_te(sn_params.threshold_filters)])
        if sn_params.prioritize is not None:
            sn_argslist.extend(["--prioritize", str(sn_params.prioritize)])
        if sn_params.add_velocity:
            sn_argslist.extend(["--add_velocity"])
        sn_argslist.extend(["--out_file_format", sn_params.out_file_format])
        if sn_params.out_seconds:
            sn_argslist.extend(["--out_seconds"])

        return sn_argslist

    def stitch_nodes(self):
        """Call the StitchNodes utility in Tempest Extremes.

        This will make a system call out to the StitchNodes method from Tempest Extremes
        (provided it has been installed as an external dependency).  StitchNodes will be
        run according to the parameters in the :attr:`stitch_nodes_parameters` attribute
        that were set when the :class:`TETracker` instance was created.

        The format of the output file containing the tracks depends on the
        :attr:`~StitchNodesParameters.out_file_format` parameter. The default ``"gfdl"``
        output is a plain-text "nodefile" format which contains a number of tracks, each
        of which in the form.

        .. code-block:: text

           start <N> <year> <month> <day> <hour>
                 <i> <j> <var1> <var2> ... <year> <month> <day> <hour>
                 ...
                 <i> <j> <var1> <var2> ... <year> <month> <day> <hour>

        - ``N`` is the number of nodes in the track (and number of lines below header).
        - ``i``, ``j`` are grid indices.
        - ``var1``, ``var2``, etc., are scalar variables as defined by
          :attr:`~StitchNodesParameters.in_fmt` (typically, lon, lat, psl, orog).
        - ``hour`` may instead be seconds if :attr:`~StitchNodesParameters.out_seconds`
          is ``True``.

        Returns
        -------
        dict
            dict of subprocess output corresponding to stdout, stderr, and returncode.

        Raises
        ------
        FileNotFoundError
            If the StitchNodes executeable from TempestExtremes cannot be found.
        RuntimeError
            If Tempest Extremes StitchNodes returns a non-zero exit code.

        References
        ----------
        `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`_
        and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_

        Examples
        --------
        To set the parameters, instantiate a :class:`TETracker` instance and run
        StitchNodes:

        >>> my_params = StitchNodesParameters(...)
        >>> my_tracker = TETracker(stitch_nodes_parameters=my_params)
        >>> result = TETracker.stitch_nodes()
        """
        sn_call_list = self._make_stitch_nodes_call()
        return self._run_te_process("StitchNodes", sn_call_list)

    def tracks(self) -> list[Track]:
        """
        Parse the tracks generated by Tempest Extremes into a list of :class:`Track`.

        The file to be read and its properties are based on the values in the
        :attr:`stitch_nodes_parameters` attribute.

        Returns
        -------
        list[Track]
            A list of :class:`Track` objects.
        """
        if self.stitch_nodes_parameters.out_file_format == "gfdl":
            tracks = self._parse_tracks_gfdl(self.stitch_nodes_parameters.output_file)
        elif self.stitch_nodes_parameters.out_file_format == "csv":
            tracks = self._parse_tracks_csv(
                self.stitch_nodes_parameters.output_file, has_header=True
            )
        elif self.stitch_nodes_parameters.out_file_format == "csvnohead":
            tracks = self._parse_tracks_csv(
                self.stitch_nodes_parameters.output_file, has_header=False
            )
        return tracks

    @staticmethod
    def _parse_gfdl_line_to_point(
        line: list[str], variable_names: list[str] | None = None
    ):
        """
        Parse a line from StitchNodes gfdl output into a :class:`Track` data point.

        Parameters
        ----------
        line : list[str]
            A list of strings representing the line split into parts.
        variable_names : list[str] | None
            List of variable names for the data columns. Defaults to None.

        Returns
        -------
        tuple
            A tuple containing (grid_i, grid_j, variables, year, month, day, hour).
        """
        grid_i = int(line[0])
        grid_j = int(line[1])
        if variable_names:
            variables = {
                name: float(value)
                for name, value in zip(variable_names, line[2:-4], strict=False)
            }
        else:
            variables = {
                f"var_{i}": float(value) for i, value in enumerate(line[2:-4], start=1)
            }
        year, month, day, hour = map(int, line[-4:])
        return grid_i, grid_j, variables, year, month, day, hour

    def _parse_tracks_gfdl(self, file_path):
        """
        Parse tracks from a gfdl file.

        Parameters
        ----------
        file_path : str
            Path to the input file.

        Returns
        -------
        list[Track]
            A list of :class:`Track` objects.
        """
        tracks = {}
        current_track_id = 0  # Initialize track ID

        with open(file_path, "r") as file:
            for line in file:
                items = line.split()
                if items[0] == "start":
                    # Start of new track, extract metadata and add Track to dict
                    current_track_id += 1
                    num_points = int(items[1])
                    year, month, day, hour = map(int, items[2:6])

                    tracks[current_track_id] = Track(
                        current_track_id,
                        num_points,
                        year,
                        month,
                        day,
                        hour,
                        calendar=self.stitch_nodes_parameters.caltype,
                    )

                # Continue processing ongoing track
                else:
                    tracks[current_track_id].add_point(
                        self._parse_gfdl_line_to_point(items)
                    )

        return list(tracks.values())

    def _parse_tracks_csv(self, file_path, has_header=False):
        """
        Generalized function to parse tracks from a csv file with or without header.

        Parameters
        ----------
        file_path : str
            Path to the input file.
        has_header : bool, optional
            Whether the file has a header. Defaults to False.

        Returns
        -------
        list[Track]
            A list of :class:`Track` objects.
        """
        tracks = {}

        with open(file_path, "r") as file:
            reader = (
                csv.DictReader(file, skipinitialspace=True)
                if has_header
                else csv.reader(file)
            )
            for row in reader:
                if has_header:
                    # Read from dict extracting variable names from keys/header
                    track_id = int(row["track_id"])
                    year, month, day, hour = map(
                        int, (row["year"], row["month"], row["day"], row["hour"])
                    )
                    grid_i = int(row["i"])
                    grid_j = int(row["j"])
                    variables = {
                        key: float(value)
                        for key, value in row.items()
                        if key
                        not in {"track_id", "year", "month", "day", "hour", "i", "j"}
                    }
                else:
                    # Read from csv assuming: id, y, m, d, h, i, j, var1, ..., varn
                    track_id = int(row[0])
                    year, month, day, hour = map(int, row[1:5])
                    grid_i = int(row[5])
                    grid_j = int(row[6])
                    variables = {
                        f"var_{i}": float(value)
                        for i, value in enumerate(row[7:], start=1)
                    }

                if track_id not in tracks:
                    tracks[track_id] = Track(
                        track_id=track_id,
                        num_points=0,
                        year=year,
                        month=month,
                        day=day,
                        hour=hour,
                        calendar=self.stitch_nodes_parameters.caltype,
                    )

                tracks[track_id].add_point(
                    (grid_i, grid_j, variables, year, month, day, hour)
                )
                tracks[track_id].num_points += 1

        return list(tracks.values())
