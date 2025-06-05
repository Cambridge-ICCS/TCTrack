"""Module providing a class that binds the Tempest Extremes code.

References
----------
- Tempest Extremes code on GitHub: https://github.com/ClimateGlobalChange/tempestextremes
- Tempest Extremes User Guide: https://climate.ucdavis.edu/tempestextremes.php
- GMD paper on Tempest Extremes v2.1: https://doi.org/10.5194/gmd-14-5023-2021
"""

import subprocess
from dataclasses import asdict, dataclass
from typing import TypedDict


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

    Attributes
    ----------
    var : str
        Name of the variable being tested. "col" in TempestExtremes.
    op : str
        Operator being used for the comparison (options include >,>=,<,<=,=,!=,|>=,|<=).
    value : float
        Value on the right-hand-side of the comparison.
    cout : str
        Either the minimum number of points where the threshold must be satisfied or the
        instruction "all", "first", or "last". "all" for all points along the path,
        "first" for just the first point, and "last" for only the last point.

    See Also
    --------
    TETracker.stitch_nodes : The StitchNodes call from the TETracker object.
    StitchNodesParameters : The StitchNodes parameter class.

    References
    ----------
    TempestExtremes Documentation: https://climate.ucdavis.edu/tempestextremes.php#StitchNodes
    Source: https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp
    """

    var: str
    op: str
    value: float
    count: str


@dataclass
class DetectNodesParameters:
    """
    Dataclass containing values used by the DetectNodes operation of TE.

    Parameters
    ----------
    in_data : list[str] | None
        List of strings of NetCDF input files.
        Defaults to ``None``.
    out_header : bool
        Whether to include header at the top of the output file.
        Defaults to ``False``.
    output_file : str | None
        Output nodefile to write to from the detection procedure.
        Defaults to ``None``.
    search_by_min : str | None
        Input variable in NetCDF files for selecting candidate points (defined as local
        minima).
        Defaults to ``None`` in TCTrack (and then ``"PSL"`` in Tempest Extremes).
    search_by_max : str | None
        Input variable in NetCDF files for selecting candidate points (defined as local
        maxima).
        Defaults to ``None``.
    closed_contours : list[TEContour] | None
        Criteria for candidates to be eliminated if they do not have a closed contour.
        Criteria are provided as a list of separate ``TEContour`` criteria.
        Defaults to ``None``.
    merge_dist : float
        DetectNodes merges candidate points with a distance (in degrees
        great-circle-distance) shorter than the specified value. Among two candidates
        within the merge distance, only the candidate with the lowest value of the
        search_by_min field or highest value of the search_by_max field are retained.
        Defaults to ``0.0``.
    lat_name : str
        string for the longitude dimension in the NetCDF files.
        Defaults to ``"lat"``.
    lat_name : str
        string for the longitude dimension in the NetCDF files.
        Defaults to ``"lat"``.
    min_lat : float
        Minimum latitude for candidate points.
        Defaults to ``0.0``.
    min_lat : float
        Maximum latitude for candidate points.
        Defaults to ``0.0``.
        If max_lat and min_lat are equal then these arguments are ignored.
    min_lon : float
        Minimum longitude for candidate points.
        Defaults to ``0.0``.
    min_lon : float
        Maximum longitude for candidate points.
        Defaults to ``0.0``.
        If max_lon and min_lon are equal then these arguments are ignored.
    regional : bool
        should lat-lon grid be periodic in the longitudinal direction.
        Defaults to ``False``
    output_commands : list[TEOutputCommand] | None
        Criteria for any additional columns to be added to the output.
        Criteria are provided as a list of separate ``TEOutputCommand`` criteria.
        Defaults to ``None``.

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
    """List of strings of NetCDF input files. Defaults to ``None``."""

    out_header: bool = False
    """Include header at the top of the output file? Defaults to ``False``."""

    output_file: str | None = None
    """Output nodefile to write to. Defaults to ``None``."""

    search_by_min: str | None = None
    """
    Input variable in NetCDF files for selecting candidate points (defined as local
    minima). Defaults to ``None`` in TCTrack (and then ``"PSL"`` in Tempest Extremes).
    """

    search_by_max: str | None = None
    """
    Input variable in NetCDF files for selecting candidate points (defined as local
    maxima). Defaults to ``None``.
    """

    closed_contours: list[TEContour] | None = None
    """
    Criteria for candidates to be eliminated if they do not have a closed contour
    as a list of separate ``TEContour`` criteria.
    Defaults to ``None``.
    """

    merge_dist: float = 0.0
    """
    DetectNodes merges candidate points with a distance (in degrees
    great-circle-distance) shorter than the specified value. Among two candidates
    within the merge distance, only the candidate with the lowest value of the
    search_by_min field or highest value of the search_by_max field are retained.
    Defaults to ``0.0``.
    """

    lat_name: str = "lat"
    """String for the latitude dimension in the NetCDF files, defaults to ``"lat"``."""

    lon_name: str = "lon"
    """String for the longitude dimension in the NetCDF files, defaults to ``"lat"``."""

    min_lat: float = 0.0
    """Minimum latitude for candidate points. Defaults to ``0.0``."""

    max_lat: float = 0.0
    """
    Maximum latitude for candidate points. Defaults to ``0.0``.
    If max_lat and min_lat are equal then these arguments are ignored.
    """

    min_lon: float = 0.0
    """Minimum longitude for candidate points. Defaults to ``0.0``."""

    max_lon: float = 0.0
    """
    Maximum longitude for candidate points. Defaults to ``0.0``.
    If ``max_lon`` and ``min_lon`` are equal then these arguments are ignored.
    """

    regional: bool = False
    """Should lat-lon grid be periodic in longitude. Defaults to ``False``."""

    output_commands: list[TEOutputCommand] | None = None
    """
    Criteria for any additional columns to be added to the output.
    Criteria are provided as a list of separate ``TEOutputCommand`` criteria.
    Defaults to ``None``.
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

    Attributes
    ----------
    output_file : str | None
        The output filename to save the tracks. "out" in TempestExtremes.
    input_file : str | None
        Filename of the DetectNodes output file. If None, it will be taken from the
        DetectNodes parameters. "in" in TempestExtremes.
    in_fmt : str | None
        Comma-separated list of the variables in the input file. If None, it will be
        taken from the DetectNodes parameters.
    max_sep : float | None
        The maximum distance allowed between candidates (degrees). "range" in
        TempestExtremes.
    min_time : int | str
        The minimum required length of a path. Either as an integer for the number of
        candidates, or a string for total duration, e.g. "24h". Default: 1
    max_gap : int
        The number of missing points allowed between candidates. Default: 0
    min_endpoint_dist : float
        The minimum required distance between the first and last candidates (degrees).
        Default: 0
    min_path_dist : float
        The minimum required acumulated distance along the path (degrees). Default: 0
    threshold_filters: list[TEThreshold] | None
        Filters for paths based on the number of nodes that satisfy a threshold.
        "thresholdcmd" in TempestExtremes.
    timestride : int
        The frequency of the input times to consider. Default: 1
    out_file_format : str
        Format of the output file. "gfdl", "csv", or "csvnoheader". Default: "gfdl"
    """

    output_file: str | None = None
    input_file: str | None = None
    in_fmt: str | None = None
    max_sep: float | None = None
    min_time: int | str = 1
    max_gap: int = 0
    min_endpoint_dist: float = 0
    min_path_dist: float = 0
    threshold_filters: list[TEThreshold] | None = None
    timestride: int = 1
    out_file_format: str = "gfdl"

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
            self.stitch_nodes_parameters = StitchNodesParameters()

        # Set StitchNodes input arguments according to DetectNodes parameters,
        # if not provided
        sn_input_none = self.stitch_nodes_parameters.input_file is None
        dn_output = self.detect_nodes_parameters.output_file
        if sn_input_none and dn_output is not None:
            self.stitch_nodes_parameters.input_file = dn_output
        sn_infmt_none = self.stitch_nodes_parameters.in_fmt is None
        dn_outcmd = self.detect_nodes_parameters.output_commands
        if sn_infmt_none and dn_outcmd is not None:
            variables = [output["var"] for output in dn_outcmd]
            self.stitch_nodes_parameters.in_fmt = ",".join(["lon", "lat", *variables])

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

        This will make a system call out to the DetectNodes method from
        Tempest Extremes (provided it has been installed as an external dependency).
        DetectNodes will be run according to the parameters in the
        ``detect_nodes_parameters`` attribute that were set when the ``TETracker``
        instance was created.

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
        To set the parameters, instantiate a ``TETracker`` instance and run
        DetectNodes:

        >>> my_params = DetectNodesParameters(...)
        >>> my_tracker = TETracker(my_params)
        >>> result = TETracker.detect_nodes()
        """
        dn_call_list = self._make_detect_nodes_call()

        try:
            result = subprocess.run(  # noqa: S603 - no shell
                dn_call_list,
                check=True,
                capture_output=True,
                text=True,
            )
            print("DetectNodes completed successfully.")
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
                "DetectNodes failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc
        except subprocess.CalledProcessError as exc:
            msg = (
                f"DetectNodes failed with a non-zero exit code: ({exc.returncode}):\n"
                f"{exc.stderr}"
            )
            raise RuntimeError(msg) from exc

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
        if self.stitch_nodes_parameters.output_file is not None:
            sn_argslist.extend(
                [
                    "--out",
                    self.stitch_nodes_parameters.output_file,
                ]
            )
        if self.stitch_nodes_parameters.input_file is not None:
            sn_argslist.extend(
                [
                    "--in",
                    self.stitch_nodes_parameters.input_file,
                ]
            )
        if self.stitch_nodes_parameters.in_fmt is not None:
            sn_argslist.extend(
                [
                    "--in_fmt",
                    self.stitch_nodes_parameters.in_fmt,
                ]
            )
        if self.stitch_nodes_parameters.max_sep is not None:
            sn_argslist.extend(
                [
                    "--range",
                    str(self.stitch_nodes_parameters.max_sep),
                ]
            )
        if self.stitch_nodes_parameters.min_time is not None:
            sn_argslist.extend(
                [
                    "--mintime",
                    str(self.stitch_nodes_parameters.min_time),
                ]
            )
        if self.stitch_nodes_parameters.max_gap is not None:
            sn_argslist.extend(
                [
                    "--maxgap",
                    str(self.stitch_nodes_parameters.max_gap),
                ]
            )
        if self.stitch_nodes_parameters.min_endpoint_dist is not None:
            sn_argslist.extend(
                [
                    "--min_endpoint_dist",
                    str(self.stitch_nodes_parameters.min_endpoint_dist),
                ]
            )
        if self.stitch_nodes_parameters.min_path_dist is not None:
            sn_argslist.extend(
                [
                    "--min_path_dist",
                    str(self.stitch_nodes_parameters.min_path_dist),
                ]
            )
        if self.stitch_nodes_parameters.threshold_filters is not None:
            sn_argslist.extend(
                [
                    "--threshold",
                    lod_to_te(self.stitch_nodes_parameters.threshold_filters),
                ]
            )
        if self.stitch_nodes_parameters.timestride is not None:
            sn_argslist.extend(
                [
                    "--timestride",
                    str(self.stitch_nodes_parameters.timestride),
                ]
            )
        if self.stitch_nodes_parameters.out_file_format is not None:
            sn_argslist.extend(
                [
                    "--out_file_format",
                    self.stitch_nodes_parameters.out_file_format,
                ]
            )

        return sn_argslist

    def stitch_nodes(self):
        """Call the StitchNodes utility in Tempest Extremes."""
        sn_call_list = self._make_stitch_nodes_call()

        try:
            subprocess.run(sn_call_list, check=True)  # noqa: S603 - no shell
        except FileNotFoundError as exc:
            msg = (
                "StitchNodes failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc
