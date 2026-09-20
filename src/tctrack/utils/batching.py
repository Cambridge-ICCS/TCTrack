"""Module providing functions for batching the tracking algorithms."""

import glob
import shutil
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, Literal, TypeAlias, TypedDict, cast

import cf

from tctrack.core import TCTracker
from tctrack.preprocessing import read_files, select_time_range


class InputArgs(TypedDict, total=False):
    """Additional arguments for an input file."""

    batch_file: str | Path | None
    store: str | Sequence[str | None]
    time_varying: bool


class _PreprocessRegistryArgs(TypedDict, total=False):
    use: str | Sequence[str]
    store: str | Sequence[str | None]


PreprocessResult: TypeAlias = cf.Field | Sequence[cf.Field] | None
PreprocessFn: TypeAlias = Callable[..., PreprocessResult]
PreprocessStep: TypeAlias = (
    tuple[PreprocessFn, dict[str, Any]]
    | tuple[PreprocessFn, dict[str, Any], _PreprocessRegistryArgs]
)


def _combine_trajectories(output_file_list: list[Path], output_file: Path) -> None:
    """Combine batched trajectory outputs into a single NetCDF file."""
    if not output_file_list:
        msg = "At least one output file is required to combine trajectories."
        raise ValueError(msg)

    fields_by_batch = [cf.read(str(file)) for file in output_file_list]  # type: ignore[operator]

    # Re-index the trajectory
    trajectory_offset = 0
    for batch_fields in fields_by_batch:
        n_trajectories = batch_fields[0].dimension_coordinate("trajectory").size
        trajectory_idx = range(trajectory_offset, trajectory_offset + n_trajectories)
        trajectory_offset = trajectory_idx[-1] + 1

        for field in batch_fields:
            field.dimension_coordinate("trajectory").set_data(trajectory_idx)

    fields = [field for fields in fields_by_batch for field in fields]

    # Pad the observation lengths to match
    max_observations = max(
        field.dimension_coordinate("observation").size for field in fields
    )
    for field in fields:
        field.pad_missing("observation", to_size=max_observations, inplace=True)
        field.dimension_coordinate("observation").set_data(range(max_observations))

    # Combine the fields
    combined_fields = cf.aggregate(fields, axes=["trajectory"], ncvar_identities=True)

    # Write the combined output file
    cf.write(combined_fields, str(output_file))  # type: ignore[operator]


def _get_fields(
    input_names: str | Sequence[str],
    fields: dict[str, cf.Field],
) -> list[cf.Field]:
    """Return the input fields requested by preprocessing."""
    if isinstance(input_names, str):
        input_names = [input_names]

    missing_names = [name for name in input_names if name not in fields]
    if missing_names:
        msg = (
            "The following fields are not available in the preprocessing registry: "
            + ", ".join(missing_names)
        )
        raise KeyError(msg)

    return [fields[name] for name in input_names]


def _store_fields(
    result: cf.Field | Sequence[cf.Field],
    store_names: str | Sequence[str | None],
    fields: dict[str, cf.Field],
    fn_name: str,
) -> None:
    """Store chosen preprocessing output fields in a registry."""
    # Ensure outputs / field names are arrays
    if isinstance(result, cf.Field):
        result = [result]
    if isinstance(store_names, str):
        store_names = [store_names]

    # Store all fields if the lengths are the same
    if len(store_names) == len(result):
        for name, field in zip(store_names, result, strict=True):
            if name is None:
                continue
            fields[name] = field

    # If not all output fields are to be stored then match by netcdf variable name
    elif len(store_names) < len(result):
        store_names = cast(
            list[str], [name for name in store_names if name is not None]
        )
        result_names = [field.nc_get_variable() for field in result]
        missing_names = [name for name in store_names if name not in result_names]
        if missing_names:
            msg = (
                f"Fields with the following names are not provided by {fn_name}: "
                + ", ".join(missing_names)
            )
            raise ValueError(msg)
        for name in store_names:
            fields[name] = result[result_names.index(name)]

    else:
        msg = f"Number of fields to store exceeds the number provided by {fn_name}."
        raise ValueError(msg)


def _run_preprocessing(
    step: PreprocessStep, batch: tuple[int, Path], fields: dict[str, cf.Field]
) -> PreprocessResult:
    """Run a batching step with supported keyword arguments."""
    fn = step[0]
    kwargs = step[1].copy()
    if len(step) == 2:  # noqa: PLR2004 - magic number
        input_names = None
        output_names = None
    elif len(step) == 3:  # noqa: PLR2004 - magic number
        input_names = step[2].get("use")
        output_names = step[2].get("store")
    else:
        msg = "Invalid preprocessing step format. It must be a tuple of length 2 or 3."
        raise ValueError(msg)

    # Replace any %BATCH% and %ITER% tags in string arguments
    i_iter, batch_dir = batch
    for k, v in kwargs.items():
        if isinstance(v, str):
            if v == "%ITER%":
                kwargs[k] = i_iter  # If just "%ITER%" use an integer instead
                continue
            kwargs[k] = v.replace("%BATCH%", str(batch_dir))
            kwargs[k] = kwargs[k].replace("%ITER%", str(i_iter))

    # Run the preprocessing function
    if input_names is None:
        result = fn(**kwargs)
    else:
        input_fields = _get_fields(input_names, fields)
        result = fn(input_fields, **kwargs)

    # (Optionally) Store the output field(s)
    if result is not None and output_names is not None:
        _store_fields(result, output_names, fields, fn.__name__)

    return result


def _parse_files(input_files: str | Iterable[str], directory: Path) -> list[str]:
    """Parse the input files to prepend the directory and expand wildcards."""
    input_files = [input_files] if isinstance(input_files, str) else input_files
    input_file_paths: list[str] = []

    for input_file in input_files:
        input_file_path = input_file

        # Prepend the batch directory for relative file paths
        if not Path(input_file).is_absolute():
            input_file_path = str(directory / input_file)

        # If it wildcards, expand these
        if glob.has_magic(input_file_path):
            matches = sorted(glob.glob(input_file_path))
            if matches:
                input_file_paths.extend(matches)
            else:
                msg = f"No files matched input pattern '{input_file_path}'."
                raise FileNotFoundError(msg)

        # Otherwise just append the path
        else:
            input_file_paths.append(input_file_path)

    return input_file_paths


def _batch_time_ranges(
    time_range: tuple[str, str],
    interval: Literal["month", "year"] | cf.TimeDuration,
    calendar: str,
) -> list[tuple[str, str]]:
    """Split a time range into time ranges for each batch."""
    time_start = cf.dt(time_range[0], calendar=calendar)
    time_end = cf.dt(time_range[1], calendar=calendar)

    if time_start >= time_end:
        msg = "The start of 'time_range' must be before its end."
        raise ValueError(msg)

    # Get the interval as a cf.TimeDuration
    if interval == "month":
        interval = cf.TimeDuration(1, "calendar_months")
    elif interval == "year":
        interval = cf.TimeDuration(1, "calendar_years")
    elif not isinstance(interval, cf.TimeDuration) or interval <= 0:
        msg = "Invalid batch interval. Use 'month', 'year' or a cf.TimeDuration object."
        raise ValueError(msg)

    # Get the list of time ranges
    batch_ranges: list[tuple[str, str]] = []
    batch_start = time_start
    while batch_start < time_end:
        batch_start, batch_end = interval.interval(batch_start)  # type: ignore
        batch_end = min(batch_end, time_end)
        batch_ranges.append((str(batch_start), str(batch_end)))
        batch_start = batch_end

    return batch_ranges


def _get_calendar(input_files: Iterable[str | tuple[str, InputArgs]]) -> str:
    """Return the calendar of the first time-varying input coordinate."""
    for input_ in input_files:
        filename, args = (input_, {}) if isinstance(input_, str) else input_
        if not args.get("time_varying", True):
            continue

        input_paths = _parse_files(filename, Path())

        # Get the calendar from the time dimension (if it exists)
        for field in cf.read(input_paths):  # type: ignore[operator]
            time_coord = field.dimension_coordinate("T", default=None)
            if time_coord is not None and time_coord.has_property("calendar"):
                return time_coord.get_property("calendar")

    # If no calendar in the input files use standard
    return "standard"


def _prepare_inputs(
    input_files: Iterable[str | tuple[str, InputArgs]],
    batch_dir: Path,
    time_range: tuple[str, str],
) -> dict[str, cf.Field]:
    """Write input files to the batch directory / store fields in the registry."""
    fields_registry: dict[str, cf.Field] = {}

    for input_ in input_files:
        # Note: filename may use wildcards to refer to multiple files
        filename, args = (input_, InputArgs()) if isinstance(input_, str) else input_
        input_file_paths = _parse_files(filename, Path())

        time_varying = args.get("time_varying", True)
        store_fields = args.get("store")

        # Load the fields
        if time_varying:
            input_fields = select_time_range(input_file_paths, time_range)
        else:
            input_fields = read_files(input_file_paths)

        # Store fields in the registry for preprocessing
        if store_fields:
            _store_fields(input_fields, store_fields, fields_registry, filename)

        # Save the fields to a file in the batch directory
        batch_file = args.get("batch_file", Path(filename).name)
        if batch_file is not None:
            cf.write(input_fields, str(batch_dir / batch_file))  # type: ignore[operator]

    return fields_registry


class BatchingConfig(TypedDict, total=False):
    """Additional arguments for :func:`batching` provided via the config argument."""

    output_dir: str | Path
    """The location to save the outputs. Default: ``"tctrack_outputs"``."""
    combine_outputs: bool
    """Whether to combine the outputs from each batch into a single ``tracks.nc`` file.
    Default: ``True``."""
    delete_batch_dirs: bool
    """Whether to delete the ``batch_[i]/`` directories. Default: ``True``."""


def batching(
    tracker: TCTracker,
    input_files: str | Iterable[str | tuple[str, InputArgs]],
    interval: Literal["month", "year"] | cf.TimeDuration,
    *,
    time_range: tuple[str, str],
    preprocessing: Sequence[PreprocessStep] | None = None,
    retrieve_data: Callable[[int, Path], None] | None = None,
    tracker_inputs: Iterable[str] = ["*"],
    config: BatchingConfig | None = None,
) -> None:
    """Perform tracking in batches with optional steps for retrieval and preprocessing.

    The outputs tracks are placed in the output directory with names ``tracks_[i].nc``.
    If ``combine_outputs`` is ``True`` then these will be combined into a single
    ``tracks.nc`` file.

    Parameters
    ----------
    tracker : TCTracker
        The tracker object used to perform the tropical cyclone tracking.
    input_files : str | Iterable[str | tuple[str, InputArgs]]
        Input netcdf files to use for each batch. These will be copied into the batch
        directory. The `*` wildcard can be used to match multiple files, in which case
        they will be combined into one file. A tuple can also be passed for each input
        file where the second value is a dictionary providing additional arguments.
        These can be:

        - ``batch_file``: The filename for the file in the batch directory. Use ``None``
          to not do so. By default it will use the same file name.
        - ``store``: Keys for storing the fields in-memory for use in the preprocessing.
          If there are multiple fields this must either match the full number of fields
          or match the netcdf variable names.
        - ``time_varying``: Whether the fields should be subset in time. Default:
          ``True``.
    interval : {"month", "year"} | cf.TimeDuration
        The calendar interval of each batch. The final batch may be shorter to end at
        the end of ``time_range``.
    time_range : tuple[str, str]
        The start and end datetimes for all batches. Must in YYYY-MM-DD format. The end
        time is open (not inclusive).
    preprocessing : Sequence[PreprocessStep] | None
        (optional) The list of preprocessing steps. These are each specified by a tuple.

        - The first entry in the tuple should be one of the preprocessing functions
          defined in :mod:`tctrack.preprocessing`. Alternatively, it can be a function
          that returns either a cf.Field, a list of cf.Field, or nothing. If it takes a
          field / fields as input this should be the first argument.
        - The second entry is a dictionary containing the arguments to pass to the
          function. String arguments can refer to the batch directory with ``%BATCH%``.
        - The optional third entry is a dictionary that can take ``store`` and/or
          ``use`` keys which allows fields to be stored and passed from memory to
          avoid unnecessary file IO. ``store`` behaves the same as in
          :attr:`input_files`.
    retrieve_data : Callable[[int, Path], None] | None
        (optional) A user-defined function that is called each iteration to retrieve
        data and put it in the batch directory. E.g. to download the data if it will not
        all fit on the filesystem in one go. The first argument of the function is the
        batch index, the second argument is the batch directory.
    tracker_inputs : Iterable[str]
        (optional) A list of filenames from the batch directory to pass to the tracker.
        By default it uses all the files (using ``["*"]``).
    config : BatchingConfig | None
        (optional) A dictionary of additional arguments. Valid keys:

        - output_dir: The location to save the outputs. Default: ``"tctrack_outputs"``.
        - combine_outputs: Whether to combine the outputs from each batch into a single
          ``tracks.nc`` file. Default: ``True``.
        - delete_batch_dirs: Whether to delete the ``batch_[i]/`` directories. Default:
          ``True``.

    Examples
    --------
    Track monthly data from 1950. The input file is loaded, selected to each monthly
    interval, and stored in memory. It is saved to the batch directory after
    preprocessing, which halves the latitude resolution and renames the netCDF variable.

    >>> from tctrack.utils import batching
    >>> from tctrack import tempest_extremes as te
    >>> from tctrack.preprocessing import subsample_field, set_nc_variable_name
    >>> preprocessing = [
    ...     (
    ...         subsample_field,
    ...         {"X": slice(0, None, 2)},
    ...         {"use": "psl", "store": "psl"},
    ...     ),
    ...     (
    ...         set_nc_variable_name,
    ...         {"field_name": "p", "output_file": "%BATCH%/psl_processed.nc"},
    ...         {"use": "psl"},
    ...     ),
    ... ]
    >>> tracker = te.TETracker()
    >>> batching(
    ...     tracker,
    ...     [("psl_*.nc", {"store": "psl", "batch_file": None})],
    ...     interval="month",
    ...     time_range=("1950-01-01", "1951-01-01"),
    ...     preprocessing=preprocessing,
    ... )
    """
    # Set the default config arguments
    if config is None:
        config = {}
    output_dir = Path(config.get("output_dir", "tctrack_outputs"))
    output_dir.mkdir(parents=True, exist_ok=True)
    combine_outputs = config.get("combine_outputs", True)
    delete_batch_dirs = config.get("delete_batch_dirs", True)

    # Get the batch time ranges
    input_files = [input_files] if isinstance(input_files, str) else list(input_files)
    calendar = _get_calendar(input_files)
    batch_time_ranges = _batch_time_ranges(time_range, interval, calendar)

    # List of output files to use for combining
    output_files: list[Path] = []

    # Perform the batching loop
    for i_iter, batch_time_range in enumerate(batch_time_ranges):
        # Create the directory for the input files
        batch_dir = output_dir / f"batch_{i_iter}"
        batch_dir.mkdir(parents=True, exist_ok=True)

        # Put input files in the batch directory & fields in the preprocessing registry
        fields = _prepare_inputs(input_files, batch_dir, batch_time_range)

        # Put additional files in the batch directory
        if retrieve_data is not None:
            retrieve_data(i_iter, batch_dir)

        # Preprocess the data
        for preprocessing_step in preprocessing or ():
            _run_preprocessing(
                preprocessing_step,
                batch=(i_iter, batch_dir),
                fields=fields,
            )

        # Run the tracker and keep track of the output files
        input_file_paths = _parse_files(tracker_inputs, batch_dir)
        output_file = output_dir / f"tracks_{i_iter}.nc"
        tracker.run_tracker(input_file_paths, str(output_file))
        output_files.append(output_file)

        # Optionally delete the batch directory
        if delete_batch_dirs:
            shutil.rmtree(batch_dir)

    if combine_outputs:
        _combine_trajectories(output_files, output_dir / "tracks.nc")
