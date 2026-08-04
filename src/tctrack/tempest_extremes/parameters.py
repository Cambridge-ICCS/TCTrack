"""Module containing parameter classes for the Tempest Extremes tracker."""

from dataclasses import dataclass
from typing import TypedDict

from tctrack.core import TCTrackerParameters


class TEContour(TypedDict):
    """
    Data required for checking a contour of a single variable during detection.

    Points will be eliminated in a detection search if they fail this criterion.
    The closed contour is determined by breadth first search: if any paths exist from
    the candidate point (or nearby minima/maxima if minmaxdist is specified) that
    reach the specified distance before achieving the specified delta then we say no
    closed contour is present.
    Each contour takes the form of a ``dict`` with keys ``"var"``, ``"delta"``,
    ``"dist"``, and ``"minmaxdist"``.

    See Also
    --------
    TETracker.detect : The DetectNodes call from the TETracker object
    TEDetectParameters : The detection parameter class

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`__
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
    Lesser-circle radius (degrees) from the pivot within which the
    criteria must be satisfied.
    """

    minmaxdist: float
    """
    Lesser-circle radius away from the candidate to search for the minima/maxima.
    If delta is positive (negative), the pivot is a local minimum (maximum).
    """


class TEOutputCommand(TypedDict):
    """
    Data required to specify an additional column in the detection output.

    Each output command takes the form of a ``dict`` with keys ``"var"``, ``"op"``, and
    ``"dist"``.

    See Also
    --------
    TETracker.detect : The DetectNodes call from the TETracker object
    TEDetectParameters : The detection parameter class

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`__
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
    Lesser-circle radius (degrees) from the candidate within which the
    operator is applied.
    """


class TEThreshold(TypedDict):
    """Data required for a threshold filter for a track trajectory when stitching.

    Any storm track trajectories that do not satisfy the threshold value for a given
    number of points will be filtered out.
    Each condition is of the form "var,op,value,count" and multiple conditions are
    separated by ";".

    See Also
    --------
    TETracker.stitch : The StitchNodes call from the TETracker object.
    TEStitchParameters : The stitching parameter class.

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`__
    and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_

    Examples
    --------
    To add a filter requiring latitude ``"lat"`` to be less than 40 degrees for
    10 or more points in each track trajectory:

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


@dataclass(repr=False)
class TEDetectParameters(TCTrackerParameters):
    """
    Dataclass containing values used by the detection operation of TE.

    See Also
    --------
    TEContour : The class used to define contour criteria
    TEOutputCommand : The class used to define additional outputs

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`__
    and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_
    """

    out_header: bool = False
    """Include header at the top of the output file?"""

    output_dir: str = ""  # Not using None to keep type-checking simple
    """
    File path of the directory to output intermediate files. If left empty, it will use
    a temporary directory that will last only for the lifetime of the :class:`TETracker`
    instance.
    """

    output_file: str = "nodes.txt"
    """Name of output nodefile to write to in the :attr:`output_dir` directory."""

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
    DetectNodes merges candidate points with a distance (in degrees) shorter than the
    specified value. Among two candidates within the merge distance, only the candidate
    with the lowest value of the search_by_min field or highest value of the
    search_by_max field are retained.
    """

    time_filter: str | None = None
    """
    Filter for the input data frequency. Options are: `"3hr"`, `"6hr"`, `"daily"`.
    Alternatively, can be a regex for the datetime using format `"YYYY-MM-DD HH:MM:SS"`.
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


@dataclass(repr=False)
class TEStitchParameters(TCTrackerParameters):
    """Dataclass containing values used by the stitch operation of TE.

    References
    ----------
    `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`__
    and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_
    """

    output_dir: str = ""  # Not using None to keep type-checking simple
    """
    File path of the directory to output intermediate files. If left empty, it will use
    the :attr:`TEDetectParameters.output_dir` value if one is provided. Otherwise, it
    will use a temporary directory that will last only for the lifetime of the
    :class:`TETracker` instance.
    """

    output_file: str = "trajectories.txt"
    """
    The output filename to write the trajectories to in the :attr:`output_dir`
    directory.
    """

    in_file: str | None = None
    """
    Filepath of the DetectNodes output file. If this and ``in_list`` are ``None``, it
    will be determined from :attr:`TEDetectParameters.output_dir` and
    :attr:`TEDetectParameters.output_file`. Called "in" in TempestExtremes.
    """

    in_list: str | None = None
    """
    File containing a list of input files to be processed together. This is unadvised to
    use at present as it is likely to be changed.
    """

    in_fmt: list[str] | None = None
    """
    List of the variables in the order they appear in the input file.
    If ``None``, it will be ``["lon", "lat", ...]``, ending in variables defined in
    :attr:`TEDetectParameters.output_commands`.
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
    """Starting date/time for stitching trajectories. Earlier times will be ignored."""

    time_end: str | None = None
    """Ending date/time for stitching trajectories. Later times will be ignored."""

    max_sep: float = 5
    """
    The maximum distance allowed between candidates (degrees). Called "range" in
    TempestExtremes.
    """

    max_gap: int | str = 0
    """
    The number of missing points allowed between candidates, as an integer. Or as a
    string for the maximum time (inclusive) between points, e.g. ``"24hr"``.
    """

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
    Format of the output file. ``"gfdl"``, ``"csv"``, or ``"csvnohead"``.
    See :meth:`TETracker.stitch` for details.
    """

    out_seconds: bool = False
    """
    For GFDL output file types, determines whether to report the sub-daily time in
    seconds (``True``) or hours (``False``).
    """

    def __post_init__(self):
        """Validate parameters."""
        if self.out_file_format not in ("gfdl", "csv", "csvnohead"):
            msg = (
                f"Invalid out_file_format ({self.out_file_format}). "
                "Allowed values are 'gfdl', 'csv', or 'csvnohead'"
            )
            raise ValueError(msg)
        if self.caltype not in ("standard", "noleap", "360_day"):
            msg = (
                f"Invalid caltype ({self.caltype}). "
                "Allowed values are 'standard', 'noleap', or '360_day'"
            )
            raise ValueError(msg)
