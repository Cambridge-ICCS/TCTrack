"""Module providing a class that binds the Tempest Extremes code.

References
----------
- Tempest Extremes code on GitHub: https://github.com/ClimateGlobalChange/tempestextremes
- Tempest Extremes User Guide: https://climate.ucdavis.edu/tempestextremes.php
- GMD paper on Tempest Extremes v2.1: https://doi.org/10.5194/gmd-14-5023-2021
"""

import csv
import tempfile
from pathlib import Path
from typing import Iterable

import cf

from tctrack.core import TCTracker, TCTrackerMetadata, TCTrackerParameters, Trajectory
from tctrack.tempest_extremes import TEDetectParameters, TEStitchParameters


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


class TETracker(TCTracker):
    """Class containing bindings to the Tempest Extremes code.

    Attributes
    ----------
    detect_parameters : TEDetectParameters
        Class containing the parameters for the detection step
    stitch_parameters : TEStitchParameters | None
        Class containing the parameters for the stitching step
    """

    # Private attributes
    _tempdir: tempfile.TemporaryDirectory

    def __init__(
        self,
        detect_parameters: TEDetectParameters | None = None,
        stitch_parameters: TEStitchParameters | None = None,
    ):
        """
        Construct the TempestExtremes class.

        Parameters
        ----------
        detect_parameters : TEDetectParameters
            Class containing the parameters for the detection step
        stitch_parameters : TEStitchParameters | None
            Class containing the parameters for the stitching step
            Defaults to the default values in TEStitchParameters Class
        """
        if detect_parameters is not None:
            self.detect_parameters: TEDetectParameters = detect_parameters
        else:
            self.detect_parameters = TEDetectParameters()

        if stitch_parameters is not None:
            self.stitch_parameters: TEStitchParameters = stitch_parameters
        else:
            self.stitch_parameters = TEStitchParameters()

        dn_params = self.detect_parameters
        sn_params = self.stitch_parameters

        # Define the default output directories
        if dn_params.output_dir == "":
            # The temporary directory is deleted when the object goes out of scope.
            # It is stored in an attribute to persist for the lifetime of the tracker.
            self._tempdir = tempfile.TemporaryDirectory()
            dn_params.output_dir = self._tempdir.name
        if sn_params.output_dir == "":
            sn_params.output_dir = dn_params.output_dir

        # Set StitchNodes input arguments according to DetectNodes parameters,
        # if not provided
        if sn_params.in_file is None and sn_params.in_list is None:
            sn_params.in_file = str(Path(dn_params.output_dir) / dn_params.output_file)
        if sn_params.in_fmt is None and dn_params.output_commands is not None:
            variables = [output["var"] for output in dn_params.output_commands]
            sn_params.in_fmt = ["lon", "lat", *variables]

    @property
    def _parameters(self) -> list[TCTrackerParameters]:
        """A list of the parameter objects that is accessible from the base class."""
        return [self.detect_parameters, self.stitch_parameters]

    def _make_detect_nodes_call(self):  # noqa: PLR0912 - all branches same logic
        """
        Construct a DetectNodes call based on options set in parameters.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a DetectNodes command
            based on the parameters set in self.detect_parameters
        """
        dn_argslist = ["DetectNodes"]

        dn_argslist.extend(["--in_data", ";".join(self._input_files)])
        out_file = str(
            Path(self.detect_parameters.output_dir) / self.detect_parameters.output_file
        )
        dn_argslist.extend(["--out", out_file])
        if self.detect_parameters.out_header:
            dn_argslist.extend(["--out_header"])
        if self.detect_parameters.search_by_min is not None:
            dn_argslist.extend(
                [
                    "--searchbymin",
                    self.detect_parameters.search_by_min,
                ]
            )
        if self.detect_parameters.search_by_max is not None:
            dn_argslist.extend(
                [
                    "--searchbymax",
                    self.detect_parameters.search_by_max,
                ]
            )
        if self.detect_parameters.closed_contours is not None:
            dn_argslist.extend(
                [
                    "--closedcontourcmd",
                    lod_to_te(self.detect_parameters.closed_contours),
                ]
            )
        dn_argslist.extend(
            [
                "--mergedist",
                str(self.detect_parameters.merge_dist),
            ]
        )
        if self.detect_parameters.time_filter is not None:
            dn_argslist.extend(
                [
                    "--timefilter",
                    self.detect_parameters.time_filter,
                ]
            )
        if self.detect_parameters.lat_name is not None:
            dn_argslist.extend(
                [
                    "--latname",
                    self.detect_parameters.lat_name,
                ]
            )
        if self.detect_parameters.lon_name is not None:
            dn_argslist.extend(
                [
                    "--lonname",
                    self.detect_parameters.lon_name,
                ]
            )
        if self.detect_parameters.min_lat is not None:
            dn_argslist.extend(
                [
                    "--minlat",
                    str(self.detect_parameters.min_lat),
                ]
            )
        if self.detect_parameters.max_lat is not None:
            dn_argslist.extend(
                [
                    "--maxlat",
                    str(self.detect_parameters.max_lat),
                ]
            )
        if self.detect_parameters.min_lon is not None:
            dn_argslist.extend(
                [
                    "--minlon",
                    str(self.detect_parameters.min_lon),
                ]
            )
        if self.detect_parameters.max_lon is not None:
            dn_argslist.extend(
                [
                    "--maxlon",
                    str(self.detect_parameters.max_lon),
                ]
            )
        if self.detect_parameters.regional:
            dn_argslist.extend(["--regional"])
        if self.detect_parameters.output_commands is not None:
            dn_argslist.extend(
                [
                    "--outputcmd",
                    lod_to_te(self.detect_parameters.output_commands),
                ]
            )

        return dn_argslist

    def detect(self):
        """
        Call the DetectNodes utility of Tempest Extremes.

        This will make a system call out to the DetectNodes method from Tempest Extremes
        (provided it has been installed as an external dependency). DetectNodes will be
        run according to the parameters in the :attr:`detect_parameters` attribute
        that were set when the :class:`TETracker` instance was created.

        The output file is a plain text file containing each of the TC candidates at
        each time from the input files. If :attr:`~TEDetectParameters.output_dir` is
        ``None`` this will be a temporary file lasting the lifetime of the
        :class:`TETracker` instance. If :attr:`~TEDetectParameters.out_header` is
        ``True`` the first two lines of the file will be a header describing the
        structure of the data. After this each time is listed in the format:

        .. code-block:: text

           <year> <month> <day> <count> <hour>
                  <i> <j> <lon> <lat> <var1> <var2> ...
                  ...
                  <i> <j> <lon> <lat> <var1> <var2> ...

        - ``count`` is the number of nodes at that time.
        - ``i``, ``j`` are the grid indices of the node.
        - ``var1``, ``var2``, etc., are scalar variables as defined by
          :attr:`~TEDetectParameters.output_commands` (typically, psl, orog).

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
        `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#DetectNodes>`__
        and the `DetectNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/DetectNodes.cpp>`_

        Examples
        --------
        To set the parameters, instantiate a :class:`TETracker` instance and run
        DetectNodes:

        >>> my_params = TEDetectParameters(...)
        >>> my_tracker = TETracker(detect_parameters=my_params)
        >>> my_tracker.set_input_files("input.nc")
        >>> result = my_tracker.detect()
        """
        Path(self.detect_parameters.output_dir).mkdir(parents=True, exist_ok=True)
        dn_call_list = self._make_detect_nodes_call()
        return self.run_tracker_subprocess("DetectNodes", dn_call_list)

    def _make_stitch_nodes_call(self):
        """
        Construct a StitchNodes call based on options set in parameters.

        Returns
        -------
        list[str]
            list of strings that can be combined to form a StitchNodes command
            based on the parameters set in self.stitch_parameters
        """
        sn_argslist = ["StitchNodes"]
        sn_params = self.stitch_parameters
        out_file = str(Path(sn_params.output_dir) / sn_params.output_file)
        sn_argslist.extend(["--out", out_file])
        if sn_params.in_file is not None:
            sn_argslist.extend(["--in", sn_params.in_file])
        if sn_params.in_list is not None:
            sn_argslist.extend(["--in_list", sn_params.in_list])
        if sn_params.in_fmt is not None:
            sn_argslist.extend(["--in_fmt", ",".join(sn_params.in_fmt)])
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

    def stitch(self):
        """Call the StitchNodes utility in Tempest Extremes.

        This will make a system call out to the StitchNodes method from Tempest Extremes
        (provided it has been installed as an external dependency).  StitchNodes will be
        run according to the parameters in the :attr:`stitch_parameters` attribute
        that were set when the :class:`TETracker` instance was created.

        The output is a file containing the data for each node of each trajectory. If
        :attr:`~TEStitchParameters.output_dir` is ``None`` this will be a temporary
        file lasting the lifetime of the :class:`TETracker` instance. The format of the
        file depends on the :attr:`~TEStitchParameters.out_file_format` parameter.
        The default ``"gfdl"`` output is a plain-text "nodefile" format which contains a
        number of track trajectories, each of which in the form.

        .. code-block:: text

           start <N> <year> <month> <day> <hour>
                 <i> <j> <var1> <var2> ... <year> <month> <day> <hour>
                 ...
                 <i> <j> <var1> <var2> ... <year> <month> <day> <hour>

        - ``N`` is number of nodes in the trajectory (and number of lines below header).
        - ``i``, ``j`` are grid indices.
        - ``var1``, ``var2``, etc., are scalar variables as defined by
          :attr:`~TEStitchParameters.in_fmt` (typically, lon, lat, psl, orog).
        - ``hour`` may instead be seconds if :attr:`~TEStitchParameters.out_seconds`
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
        `TempestExtremes Documentation <https://climate.ucdavis.edu/tempestextremes.php#StitchNodes>`__
        and the `StitchNodes Source <https://github.com/ClimateGlobalChange/tempestextremes/blob/master/src/nodes/StitchNodes.cpp>`_

        Examples
        --------
        To set the parameters, instantiate a :class:`TETracker` instance and perform
        stitching:

        >>> my_params = TEStitchParameters(...)
        >>> my_tracker = TETracker(stitch_parameters=my_params)
        >>> result = my_tracker.stitch()
        """
        Path(self.stitch_parameters.output_dir).mkdir(parents=True, exist_ok=True)
        sn_call_list = self._make_stitch_nodes_call()
        return self.run_tracker_subprocess("StitchNodes", sn_call_list)

    def read_trajectories(self) -> list[Trajectory]:
        """
        Parse outputs from Tempest Extremes to list of :class:`tctrack.core.Trajectory`.

        The file to be read and its properties are based on the values in the
        :attr:`stitch_parameters` attribute.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        out_file = str(
            Path(self.stitch_parameters.output_dir) / self.stitch_parameters.output_file
        )
        if self.stitch_parameters.out_file_format == "gfdl":
            trajectories = self._parse_trajectories_gfdl(out_file)
        elif self.stitch_parameters.out_file_format == "csv":
            trajectories = self._parse_trajectories_csv(out_file, has_header=True)
        elif self.stitch_parameters.out_file_format == "csvnohead":
            trajectories = self._parse_trajectories_csv(out_file, has_header=False)
        return trajectories

    @staticmethod
    def _parse_gfdl_line_to_point(
        line: list[str], variable_names: list[str] | None = None
    ) -> tuple[list[int], dict[str, int | float]]:
        """
        Parse line from StitchNodes gfdl output into a trajectory data point.

        Data point format is that expected by a :class:`tctrack.core.Trajectory`.

        Parameters
        ----------
        line : list[str]
            A list of strings representing the line split into parts.
        variable_names : list[str] | None
            List of variable names for the data columns. Defaults to None.

        Returns
        -------
        tuple
            A tuple containing the time as an integer list of [year, day, month, hour]
            and a dict of variables.
        """
        return_vars: dict[str, int | float] = {}
        return_vars.update({"grid_i": int(line[0]), "grid_j": int(line[1])})
        if variable_names:
            return_vars.update(
                {
                    name: float(value)
                    for name, value in zip(variable_names, line[2:-4], strict=False)
                }
            )
        else:
            return_vars.update(
                {
                    f"var_{i}": float(value)
                    for i, value in enumerate(line[2:-4], start=1)
                }
            )
        time = list(map(int, line[-4:]))
        return time, return_vars

    def _parse_trajectories_gfdl(self, file_path):
        """
        Parse track trajectories from a gfdl file.

        Parameters
        ----------
        file_path : str
            Path to the input file.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        trajectories = {}
        current_trajectory_id = 0  # Initialize trajectory ID

        # Get variable names from in_fmt
        var_names = self.stitch_parameters.in_fmt or []

        with open(file_path, "r") as file:
            for line in file:
                items = line.split()
                if items[0] == "start":
                    # Start of new trajectory.
                    # Extract metadata and add Trajectory to dict
                    current_trajectory_id += 1
                    time = list(map(int, items[2:6]))

                    trajectories[current_trajectory_id] = Trajectory(
                        current_trajectory_id,
                        time,
                        calendar=self.stitch_parameters.caltype,
                    )

                # Continue processing ongoing trajectory
                else:
                    trajectories[current_trajectory_id].add_point(
                        *self._parse_gfdl_line_to_point(items, var_names)
                    )

        return list(trajectories.values())

    def _parse_trajectories_csv(self, file_path, has_header=False):
        """
        Generalized function to parse trajectories from csv file with/without header.

        Parameters
        ----------
        file_path : str
            Path to the input file.
        has_header : bool, optional
            Whether the file has a header. Defaults to False.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        trajectories = {}

        with open(file_path, "r") as file:
            reader = (
                csv.DictReader(file, skipinitialspace=True)
                if has_header
                else csv.reader(file)
            )
            for row in reader:
                if has_header:
                    # Read from dict extracting variable names from keys/header
                    trajectory_id = int(row["track_id"])
                    time = [int(row[k]) for k in ("year", "month", "day", "hour")]
                    variables_dict = {"grid_i": int(row["i"]), "grid_j": int(row["j"])}
                    variables_dict.update(
                        {
                            key: float(value)
                            for key, value in row.items()
                            if key
                            not in {
                                "track_id",
                                "year",
                                "month",
                                "day",
                                "hour",
                                "i",
                                "j",
                            }
                        }
                    )
                else:
                    # Read from csv assuming: id, y, m, d, h, i, j, var1, ..., varn
                    trajectory_id = int(row[0])
                    time = list(map(int, row[1:5]))
                    variables_dict = {"grid_i": int(row[5]), "grid_j": int(row[6])}
                    # Get variable names from in_fmt
                    var_names = self.stitch_parameters.in_fmt or [
                        f"var_{i + 1}" for i in range(len(row[7:]))
                    ]
                    variables_dict.update(
                        {
                            var_name: float(row[7 + i])
                            for i, var_name in enumerate(var_names)
                        }
                    )

                if trajectory_id not in trajectories:
                    trajectories[trajectory_id] = Trajectory(
                        trajectory_id=trajectory_id,
                        time=time,
                        calendar=self.stitch_parameters.caltype,
                    )

                trajectories[trajectory_id].add_point(time, variables_dict)

        return list(trajectories.values())

    def _set_metadata(self) -> None:
        """Set the time and variable metadata attributes by reading from input files.

        Reads metadata for each variable listed in
        :attr:`detect_parameters.output_commands` from the input NetCDF files that are
        set by :meth:`set_input_files` (matching the NetCDF variable
        name). These will be stored in the :attr:`variable_metadata` attribute as a
        dictionary of :class:`TCTrackerMetadata` objects. This will be called from the
        :meth:`set_metadata` method.

        Raises
        ------
        ValueError
            If a variable is not found in the input files.

        Examples
        --------
        To read in the metadata for ``psl`` from ``inputs.nc``:

        >>> detect_params = TEDetectParameters(
        >>>     output_commands=[TEOutputCommand(var="psl", operator="min", dist=0)],
        >>> )
        >>> tracker = TETracker(detect_params, stitch_params)
        >>> tracker.set_input_files("inputs.nc")
        >>> tracker.set_metadata()
        >>> tracker.variable_metadata
        {
            "psl": TCTrackerMetadata(
                properties={
                    "standard_name": "air_pressure_at_sea_level",
                    "long_name": "Sea Level Pressure",
                    "units": "Pa",
                },
                constructs=[<CF CellMethod: area: point>],
            ),
        }
        """
        # set time metadata
        variable_name = (
            self.detect_parameters.search_by_min
            or self.detect_parameters.search_by_max
            or "PSL"
        )
        fields = cf.read(
            self._input_files,
            select=f"ncvar%{variable_name}",  # type: ignore[operator]
            netcdf_backend="netCDF4",
        )
        if not fields:
            msg = f"Variable '{variable_name}' not found in input files."
            raise ValueError(msg)
        _, time_coord = fields[0].construct_item("time")
        time_arr = time_coord.datetime_array
        self._time_metadata = {
            "calendar": time_coord.get_property("calendar"),
            "units": time_coord.get_property("units"),
            "start_time": time_arr[0],
            "end_time": time_arr[-1],
        }

        # set variable metadata
        self._variable_metadata = {}

        # Set the variable metadata for the grid indices generated by Tempest Extremes
        self._variable_metadata["grid_i"] = TCTrackerMetadata(
            {"long_name": "longitudinal grid index"}
        )
        self._variable_metadata["grid_j"] = TCTrackerMetadata(
            {"long_name": "latitudinal grid index"}
        )

        # Set the variable metadata for the output variables
        var_outputs = self.detect_parameters.output_commands
        if var_outputs is None:
            return

        for var_output in var_outputs:
            var_name = var_output["var"]

            # Get the variable field from the netcdf file
            fields = cf.read(
                self._input_files,
                select=f"ncvar%{var_name}",  # type: ignore[operator]
                netcdf_backend="netCDF4",
            )
            if not fields:
                msg = f"Variable '{var_name}' not found in input files."
                raise ValueError(msg)
            field = fields[0]

            # Read and store the relevant metadata
            self._variable_metadata[var_name] = TCTrackerMetadata(
                {
                    "standard_name": field.get_property("standard_name", var_name),
                    "long_name": field.get_property("long_name", var_name),
                    "units": field.get_property("units", "unknown"),
                }
            )

            # Add information about how the value is determined using `output_commands`
            methods = {
                "max": "maximum",
                "min": "minimum",
                "avg": "mean",
            }
            method = methods.get(var_output["operator"])
            if method is not None:
                dist = var_output["dist"]
                if dist == 0:
                    cell_method = cf.CellMethod("area", "point")
                else:
                    qualifier = {"comment": f"lesser circle of radius {dist} degrees"}
                    cell_method = cf.CellMethod("area", method, qualifiers=qualifier)
                self._variable_metadata[var_name].constructs = [cell_method]

    def run_tracker(self, input_files: str | Iterable[str], output_file: str):
        """Run TempestExtremes tracker to obtain tropical cyclone track trajectories.

        This first runs :meth:`detect` to get TC candidates at each time. Then
        these are combined into trajectories using :meth:`stitch`.
        The output is then saved as a CF-compliant NetCDF trajectory file.

        Arguments
        ---------
        input_files : str | Iterable[str]
            A (list of) file path(s) containing NetCDF input data to use in the tracker.
        output_file : str
            Filename to which the tropical cyclone trajectories are saved.

        Raises
        ------
        FileNotFoundError
            - If the TempestExtremes executables cannot be found.
            - If the stitch output file does not exist.
        RuntimeError
            If the TempestExtremes commands return a non-zero exit code.

        Examples
        --------
        To create the tracker instance, then use run_tracker to perform the detection,
        stitching, and generate output.

        >>> detect_params = TEDetectParameters(...)
        >>> stitch_params = TEStitchParameters(...)
        >>> my_tracker = TETracker(detect_params, stitch_params)
        >>> my_tracker.run_tracker("input.nc", "output.nc")
        """
        self.set_input_files(input_files)
        self.detect()
        self.stitch()
        self.to_netcdf(output_file)
