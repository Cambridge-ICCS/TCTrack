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
    """Data required for contouring a single variable in DetectNodes."""

    var: str
    delta: float
    dist: float
    minmaxdist: float


class TEOutputComand(TypedDict):
    """Data required for output command of single variable in DetectNodes."""

    var: str
    operation: float  # `op` in Tempest Extremes
    dist: float


class TEThreshold(TypedDict):
    """Data required for a threshold filter for a track in StitchNodes."""

    var: str  # `col` in Tempest Extremes
    operation: str  # `op` in Tempest Extremes
    value: float
    count: str  # Number or 'all', 'first', or 'last'


@dataclass
class DetectNodesParameters:
    """Dataclass containing values used by the DetectNodes operation of TE."""

    in_data: list[str] | None = None
    out_header: bool = False
    output_file: str | None = None
    search_by_min: str = "PSL"
    merge_dist: float = 0.0  # TODO 6.0 used in paper
    closed_contours: list[TEContour] | None = None
    output_commands: list[TEOutputComand] | None = None

    def __str__(self) -> str:
        """Improve the representation to users."""
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
    detect_nodes_parameters : DetectNodesParameters | None
        class containing the parameters for the DetectNodes algorithm
    stitch_nodes_parameters : StitchNodesParameters | None
        class containing the parameters for the StitchNodes algorithm

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
            class containing the parameters for the DetectNodes algorithm
        stitch_nodes_parameters : StitchNodesParameters | None
            class containing the parameters for the StitchNodes algorithm
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

    def detect_nodes(self):
        """Call the DetectNodes utility in Tempest Extremes."""
        # Construct DetectNodes call based on options set in the TempestExtremes class
        dn_argslist = ["DetectNodes"]
        if self.detect_nodes_parameters.in_data is not None:
            dn_argslist.extend(
                [
                    "--in_data",
                    ";".join(self.detect_nodes_parameters.in_data),
                ]
            )
            dn_argslist.extend(
                [
                    "--searchbymin",
                    self.detect_nodes_parameters.search_by_min,
                ]
            )
            dn_argslist.extend(
                [
                    "--mergedist",
                    str(self.detect_nodes_parameters.merge_dist),
                ]
            )
            if self.detect_nodes_parameters.closed_contours is not None:
                dn_argslist.extend(
                    [
                        "--closedcontourcmd",
                        lod_to_te(self.detect_nodes_parameters.closed_contours),
                    ]
                )
            if self.detect_nodes_parameters.out_header:
                dn_argslist.extend(["--out_header"])
            if self.detect_nodes_parameters.output_commands is not None:
                dn_argslist.extend(
                    [
                        "--outputcmd",
                        lod_to_te(self.detect_nodes_parameters.output_commands),
                    ]
                )
            if self.detect_nodes_parameters.output_file is not None:
                dn_argslist.extend(["--out", self.detect_nodes_parameters.output_file])

        print(dn_argslist)
        print(" ".join(dn_argslist))

        try:
            subprocess.run(dn_argslist, check=True)  # noqa: S603 - no shell
        except FileNotFoundError as exc:
            msg = (
                "DetectNodes failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc

    def stitch_nodes(self):
        """Call the StitchNodes utility in Tempest Extremes."""
        # Construct StitchNodes call based on options set in the TempestExtremes class
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

        try:
            subprocess.run(sn_argslist, check=True)  # noqa: S603 - no shell
        except FileNotFoundError as exc:
            msg = (
                "StitchNodes failed because the executable could not be found.\n"
                "Did you provide the full executeable path or add it to $PATH?\n"
            )
            raise FileNotFoundError(msg) from exc
