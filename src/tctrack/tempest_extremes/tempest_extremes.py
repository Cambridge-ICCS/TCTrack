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
    """
    return ";".join(",".join(str(value) for value in d.values()) for d in inputs)


class TEContour(TypedDict):
    """
    Data required for contouring a single variable in DetectNodes.

    Points will be eliminated in a DetectNodes search if they fail this criterion

    Attributes
    ----------
    var : str
        Name of the variable to contour in NetCDF files.
    delta : float
        Amount by which the field must change from the pivot value.
        If positive (negative) the field must increase (decrease) by this value along
        the contour.
    dist : float
         Great-circle distance (degrees) from the pivot within which the
         criteria must be satisfied.
    minmaxdist : float
        Great-circle distance away from the candidate to search for the minima/maxima.
        If delta is positive (negative), the pivot is a local minimum (maximum).

    See Also
    --------
    TETracker.DetectNodes : The Detect Nodes call from the TETracker object
    DetectNodesParameters : The DetectNodes parameter class

    References
    ----------
    TempestExtremes Documentation: https://climate.ucdavis.edu/tempestextremes.php#DetectNodes
    Source: https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp
    """

    var: str
    delta: float
    dist: float
    minmaxdist: float


class TEOutputCommand(TypedDict):
    """
    Data required to specify an additional column in the DetectNodes output.

    Each output command takes the form "var,op,dist".

    Attributes
    ----------
    var : str
        Name of the variable to write in output files.
    op : string
        Operator that is applied over all points within the specified distance of the
        candidate (options include max, min, avg, maxdist, mindist).
    dist : float
         Great-circle distance (degrees) from the candidate within which the
         operator is applied.

    See Also
    --------
    TETracker.DetectNodes : The Detect Nodes call from the TETracker object
    DetectNodesParameters : The DetectNodes parameter class

    References
    ----------
    TempestExtremes Documentation: https://climate.ucdavis.edu/tempestextremes.php#DetectNodes
    Source: https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp
    """

    var: str
    op: str
    dist: float


@dataclass
class DetectNodesParameters:
    """
    Dataclass containing values used by the DetectNodes operation of TE.

    Attributes
    ----------
    in_data : list[str] | None
        List of strings of NetCDF input files. Defaults to None.
    out_header : bool
        Whether to include header at the top of the output file. Defaults to False.
    output_file : str | None
        The output nodefile to write to from the detection procedure. Defaults to None.
    search_by_min : str | None
        Input variable in NetCDF files for selecting candidate points (defined as local
        minima). Defaults to None in TCTrack (and then "PSL" in Tempest Extremes).
    search_by_max : str | None
        Input variable in NetCDF files for selecting candidate points (defined as local
        maxima). Defaults to None.
    closed_contours : list[TEContour] | None
        Criteria for candidates to be eliminated if they do not have a closed contour.
        Defaults to None.
        Criteria are provided as a list of separate TEContour criteria.
    merge_dist : float
        DetectNodes merges candidate points with a distance (in degrees
        great-circle-distance) shorter than the specified value. Among two candidates
        within the merge distance, only the candidate with the lowest value of the
        search_by_min field or highest value of the search_by_max field are retained.
        Defaults to 0.0.
    lat_name : str
        string for the longitude dimension in the NetCDF files, defaults to "lat"
    lat_name : str
        string for the longitude dimension in the NetCDF files, defaults to "lat"
    min_lat : float
        Minimum latitude for candidate points. Defaults to `0.0`
    min_lat : float
        Maximum latitude for candidate points. Defaults to `0.0`
        If max_lat and min_lat are equal then these arguments are ignored.
    min_lon : float
        Minimum longitude for candidate points. Defaults to `0.0`
    min_lon : float
        Maximum longitude for candidate points. Defaults to `0.0`
        If max_lon and min_lon are equal then these arguments are ignored.
    regional : bool
        should lat-lon grid be periodic in the longitudinal direction.
        Defaults to False
    output_commands : list[TEOutputCommand] | None
        Criteria for any additional columns to be added to the output.
        Defaults to None.
        Criteria are provided as a list of separate TEOutputCommand criteria.

    References
    ----------
    TempestExtremes Documentation: https://climate.ucdavis.edu/tempestextremes.php#DetectNodes
    Source: https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp
    """

    in_data: list[str] | None = None
    out_header: bool = False
    output_file: str | None = None
    search_by_min: str | None = None
    search_by_max: str | None = None
    closed_contours: list[TEContour] | None = None
    merge_dist: float = 0.0
    lat_name: str = "lat"
    lon_name: str = "lon"
    min_lat: float = 0.0
    max_lat: float = 0.0
    min_lon: float = 0.0
    max_lon: float = 0.0
    regional: bool = False
    output_commands: list[TEOutputCommand] | None = None

    def __str__(self) -> str:
        """Improve the representation to users."""
        attributes = "\n\t".join(
            f"{key} \t = {value}" for key, value in asdict(self).items()
        )
        return f"DetectNodesParameters(\n\t{attributes}\n)"


class TETracker:
    """Class containing bindings to the Tempest Extremes code.

    Attributes
    ----------
    detect_nodes_parameters : DetectNodesParameters
        class containing the parameters for the DetectNodes algorithm

    """

    def __init__(
        self,
        detect_nodes_parameters: DetectNodesParameters | None = None,
    ):
        """
        Construct the TempestExtremes class.

        Parameters
        ----------
        detect_nodes_parameters : DetectNodesParameters
            class containing the parameters for the DetectNodes algorithm
            Defaults to the default values in DetectNodesParameters Class
        """
        if detect_nodes_parameters is not None:
            self.detect_nodes_parameters: DetectNodesParameters = (
                detect_nodes_parameters
            )
        else:
            self.detect_nodes_parameters = DetectNodesParameters()

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
        """Call the DetectNodes utility in Tempest Extremes."""
        dn_call_list = self._make_detect_nodes_call()

        try:
            subprocess.run(dn_call_list, check=True)  # noqa: S603 - no shell
        except FileNotFoundError as exc:
            msg = (
                "DetectNodes failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc
