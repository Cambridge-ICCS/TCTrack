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


class TEDetectThreshold(TypedDict):
    """Data for a threshold filter during the detection of candidate points.

    Any candidates that do not satisfy the threshold value within a certain distance are
    filtered out.

    When calling TempestExtremes each filter is mapped to the form "var,op,value,dist"
    and multiple conditions are separated by ";".

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
    To add a filter requiring surface windspeed ``"sfcWind"`` to be greater than 16 m/s
    somewhere within 1 degree of the centre:

    >>> TEDetectThreshold(var="sfcWind", op=">", value=16, dist=1)
    {"var": "sfcWind", "op": ">", "value": 16, "dist": 1},
    """

    var: str
    """Name of the variable being tested."""

    op: str
    """Operator used for the comparison (options include >,>=,<,<=,=,!=)."""

    value: float
    """Value on the right-hand-side of the comparison."""

    dist: float | str
    """Great-circle distance within which the threshold must be satisfied."""


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

    If using latitude or longitude these should be specified as ``"lat"`` or ``"lon"``
    regardless of their names in the input files.

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

    no_closed_contours: list[TEContour] | None = None
    """Criteria for candidates to be eliminated if they *do* have a closed contour
    as a list of separate :class:`TEContour` criteria.
    """

    thresholds: list[TEDetectThreshold] | None = None
    """Threshold criteria for candidates to be eliminated. I.e. if there is not a point
    within a given distance that satisfies the condition. Given as a list of separate
    :class:`TEDetectThreshold` criteria.
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

    lat_name: str = "latitude"
    """String for the latitude dimension in the NetCDF files."""

    lon_name: str = "longitude"
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


def _nc_names_defaults_uz(nc_names: dict[str, str]) -> None:
    """Set the default netcdf variable names for parameter_set_uz."""
    nc_names.setdefault("longitude", "longitude")
    nc_names.setdefault("latitude", "latitude")
    nc_names.setdefault("msl", "msl")
    nc_names.setdefault("orog", "orog")
    nc_names.setdefault("si10", "si10")
    gh_name = nc_names.get("gh", "gh")
    nc_names.setdefault("ghdiff", f"_DIFF({gh_name}(300hPa),{gh_name}(500hPa))")


def parameter_set_uz(
    nc_names: dict[str, str],
) -> tuple[TEDetectParameters, TEStitchParameters]:
    """Return parameters for the UZ algorithm for detecting tropical cyclones.

    Presented in [Ullrich2017]_.

    This uses sea-level pressure for detection

    Parameters
    ----------
    nc_names : dict[str, str]
        An optional dictionary of NetCDF variable names.

        The dictionary keys match the default variable names, which are shown in the
        table below. The difference in geopotential height uses the 300 and 500 hPa
        levels by default, but can be overridden, for example:
        ``{"ghdiff": "_DIFF(gh(250hPa),gh(500hPa))"}``.

        .. list-table:: NetCDF variable name defaults
           :header-rows: 1

           * - Variable
             - Key / |br| Default
             - Notes
           * - Longitude
             - ``longitude``
             -
           * - Latitude
             - ``latitude``
             -
           * - Sea-level pressure
             - ``msl``
             -
           * - Surface elevation
             - ``orog``
             -
           * - Surface wind speed
             - ``si10``
             -
           * - Geopotential height
             - ``gh``
             - Used to calculate ``ghdiff`` at 300 hPa and 500 hPa.
           * - Geopotential height difference
             - ``ghdiff``
             - Default: ``_DIFF(gh(300hPa),gh(500hPa))``.
    """
    nc_names = nc_names.copy()
    _nc_names_defaults_uz(nc_names)

    return (
        TEDetectParameters(
            search_by_min=nc_names["msl"],
            time_filter="6hr",
            merge_dist=6.0,
            closed_contours=[
                TEContour(var=nc_names["msl"], delta=200.0, dist=5.5, minmaxdist=0.0),
                TEContour(var=nc_names["ghdiff"], delta=-6.0, dist=6.5, minmaxdist=1.0),
            ],
            lon_name=nc_names["longitude"],
            lat_name=nc_names["latitude"],
            out_header=True,
            output_commands=[
                TEOutputCommand(var=nc_names["msl"], operator="min", dist=0.0),
                TEOutputCommand(var=nc_names["orog"], operator="max", dist=0.0),
                TEOutputCommand(var=nc_names["si10"], operator="max", dist=2.0),
            ],
        ),
        TEStitchParameters(
            caltype="360_day",
            max_sep=8.0,
            min_time="54h",
            max_gap="24h",
            min_endpoint_dist=8.0,
            threshold_filters=[
                TEThreshold(var=nc_names["latitude"], op="<=", value=50, count=10),
                TEThreshold(var=nc_names["latitude"], op=">=", value=-50, count=10),
                TEThreshold(var=nc_names["orog"], op="<=", value=150, count=10),
                TEThreshold(var=nc_names["si10"], op=">=", value=10, count=10),
            ],
        ),
    )


def _nc_names_defaults_owz(nc_names: dict[str, str]) -> None:
    """Set the default netcdf variable names for parameter_set_owz."""
    nc_names.setdefault("longitude", "longitude")
    nc_names.setdefault("latitude", "latitude")
    nc_names.setdefault("msl", "msl")
    nc_names.setdefault("u10", "u10")
    nc_names.setdefault("v10", "v10")

    for name, levels in (
        ("owz", (850, 500)),
        ("r", (950, 700)),
        ("q", (950,)),
        ("vo", (850,)),
        ("u", (925, 850, 200)),
        ("v", (925, 850, 200)),
    ):
        variable_name = nc_names.get(name, name)
        for level in levels:
            nc_names.setdefault(f"{name}{level}", f"{variable_name}({level}hPa)")

    vws_default = (
        "_VECMAG("
        f"_DIFF({nc_names['u850']},{nc_names['u200']}),"
        f"_DIFF({nc_names['v850']},{nc_names['v200']}))"
    )
    nc_names.setdefault("vws", vws_default)
    nc_names.setdefault("si10", f"_VECMAG({nc_names['u10']},{nc_names['v10']})")

    if ws := nc_names.get("ws"):
        ws_default = ws + "(925hPa)"
    else:
        ws_default = f"_VECMAG({nc_names['u925']},{nc_names['v925']})"
    nc_names.setdefault("ws925", ws_default)


def parameter_set_owz(
    nc_names: dict[str, str],
) -> tuple[TEDetectParameters, TEStitchParameters]:
    r"""Return parameters for the OWZ algorithm for detecting tropical cyclones.

    Presented in [Tory2013]_. Adapted to TempestExtremes in [Bourdin2022]_.

    This performs detection using the Okubo-Weiss-Zeta (OWZ) parameter which should be
    precomputed:

    .. math::

       OWZ = \eta \ \mathrm{sign}(f) \times
             \max\left(\frac{\zeta^2 - (E^2 + F^2)}{\zeta^2}, 0\right)

    where :math:`\eta` is the absolute vorticity, :math:`\zeta` is the sum of the
    relative vorticity, :math:`f` is the coriolis parameter, and :math:`E` and :math:`F`
    are the stretching and the shearing deformation, given by

    .. math::

       E = \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}, \quad
       F = \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}.


    Returns
    -------
    tuple[TEDetectParameters, TEStitchParameters]
        Detection and stitching parameters.

    Parameters
    ----------
    nc_names : dict[str, str]
        An optional dictionary of NetCDF variable names.

        The dictionary keys match the default variable names, which are shown in the
        table below. Some variables require specific vertical levels. These can also be
        overriden, for example if there is data for OWZ at 900 hPa but not 850 hPa:
        ``{"owz850": "owz(900hPa)"}``.

        .. list-table:: NetCDF variable name defaults
           :header-rows: 1

           * - Variable
             - Key / |br| Default
             - Levels
             - Notes
           * - Longitude
             - ``longitude``
             -
             -
           * - Latitude
             - ``latitude``
             -
             -
           * - OWZ parameter
             - ``owz``
             - ``owz850`` 850 hPa |br| ``owz500`` 500 hPa
             -
           * - Relative humidity
             - ``r``
             - ``r950`` 950 hPa |br| ``r700`` 700 hPa
             -
           * - Specific humidity
             - ``q``
             - ``q950`` 950 hPa
             -
           * - Relative vorticity
             - ``vo``
             - ``vo850`` 850 hPa
             -
           * - Wind components
             - ``u`` / ``v``
             - ``u925`` / ``v925`` 925 hPa |br|
               ``u850`` / ``v850`` 850 hPa |br|
               ``u200`` / ``v200`` 200 hPa |br|
             - Not neeeded if ``ws`` and ``vws`` are set.
           * - Wind speed
             - ``ws``
             - ``ws925`` 925 hPa
             - If not set, calculated from ``u925`` / ``v925``.
           * - Vertical wind shear
             - ``vws``
             -
             - If not set, calculated from ``u850`` / ``v850`` and ``u200`` / ``v200``.
           * - Surface wind |br| components (10m)
             - ``u10`` / ``v10``
             -
             - Not neeeded if ``si10`` is set.
           * - Surface wind speed |br| (10m)
             - ``si10``
             -
             - If not set, calculated from ``u10`` / ``v10``.
           * - Sea-level pressure
             - ``msl``
             -
             -
    """
    nc_names = nc_names.copy()
    _nc_names_defaults_owz(nc_names)

    return (
        TEDetectParameters(
            search_by_max=nc_names["owz850"],
            time_filter="6hr",
            merge_dist=5.0,
            thresholds=[
                TEDetectThreshold(var=nc_names["owz850"], op=">=", value=5e-5, dist=2),
                TEDetectThreshold(var=nc_names["owz500"], op=">=", value=4e-5, dist=2),
                TEDetectThreshold(var=nc_names["r950"], op=">=", value=70, dist=2),
                TEDetectThreshold(var=nc_names["r700"], op=">=", value=50, dist=2),
                TEDetectThreshold(var=nc_names["vws"], op="<=", value=25, dist=2),
                TEDetectThreshold(var=nc_names["q950"], op=">=", value=0.01, dist=2),
            ],
            lon_name=nc_names["longitude"],
            lat_name=nc_names["latitude"],
            out_header=True,
            output_commands=[
                TEOutputCommand(var=nc_names["owz500"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["r950"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["r700"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["vws"], operator="min", dist=1),
                TEOutputCommand(var=nc_names["q950"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["msl"], operator="min", dist=3),
                TEOutputCommand(var=nc_names["si10"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["ws925"], operator="max", dist=2),
                TEOutputCommand(var=nc_names["vo850"], operator="max", dist=2),
            ],
        ),
        TEStitchParameters(
            caltype="360_day",
            max_sep=5.0,
            min_time="48h",
            max_gap="24h",
            threshold_filters=[
                TEThreshold(var=nc_names["owz850"], op=">=", value=6e-5, count=9),
                TEThreshold(var=nc_names["owz500"], op=">=", value=5e-5, count=9),
                TEThreshold(var=nc_names["r950"], op=">=", value=85, count=9),
                TEThreshold(var=nc_names["r700"], op=">=", value=70, count=9),
                TEThreshold(var=nc_names["vws"], op="<=", value=12.5, count=9),
                TEThreshold(var=nc_names["q950"], op=">=", value=0.014, count=9),
                TEThreshold(var=nc_names["si10"], op=">=", value=16, count=1),
            ],
        ),
    )


def parameter_set(
    name: str,
    nc_names: dict[str, str] | None = None,
) -> tuple[TEDetectParameters, TEStitchParameters]:
    """Return a default TempestExtremes parameter set.

    Returns
    -------
    tuple[TEDetectParameters, TEStitchParameters]
        Detection and stitching parameters.

    Parameters
    ----------
    name : str
        The parameter set to return. Valid values are:

        - ``"default"`` or ``"UZ"`` for :func:`parameter_set_uz`.
        - ``"OWZ"`` for :func:`parameter_set_owz`.
    nc_names : dict[str, str]
        An optional dictionary of netcdf variable names to use instead of the defaults.
        The defaults typically correspond to ECMWF/ERA5 names or CMIP names where these
        don't exist. See the individual parameter_set functions for more details.

    Raises
    ------
    ValueError
        If ``name`` is not a supported parameter set.
    """
    if nc_names is None:
        nc_names = {}

    parameter_sets = [
        (["default", "UZ"], parameter_set_uz),
        (["OWZ"], parameter_set_owz),
    ]
    for aliases, parameter_set in parameter_sets:
        if name in aliases:
            return parameter_set(nc_names)

    available = ", ".join("/".join(aliases) for aliases, _ in parameter_sets)
    msg = f"Unknown parameter set: {name}. Available sets: {available}."
    raise ValueError(msg)
