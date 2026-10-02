"""Module providing functions for batching the tracking algorithms."""

import functools
import glob
import shutil
import warnings
from collections.abc import Callable, Iterable, Sequence
from datetime import timedelta
from pathlib import Path
from typing import Any, Literal, TypeAlias, TypedDict, cast

import cf
import cftime

from tctrack.core import TCTracker
from tctrack.preprocessing import read_files, select_time_range


class InputArgs(TypedDict, total=False):
    """Additional arguments for an input file."""

    batch_file: str | Path | None
    store: str | Sequence[str | None]
    time_varying: bool


PreprocessResult: TypeAlias = cf.Field | Sequence[cf.Field] | None
PreprocessFn: TypeAlias = (
    Callable[..., PreprocessResult] | functools.partial[PreprocessResult]
)
PreprocessStep: TypeAlias = (
    PreprocessFn | tuple[PreprocessFn, str | Sequence[str | None]]
)
TimeRange: TypeAlias = tuple[cftime.datetime, cftime.datetime]


def _combine_trajectories(
    output_file_list: list[Path],
    output_file: Path,
    batch_time_ranges: list[TimeRange],
) -> None:
    """Combine batched trajectory outputs into a single NetCDF file."""
    if not output_file_list:
        msg = "At least one output file is required to combine trajectories."
        raise ValueError(msg)

    # Get the fields for each batch
    # Selecting only the trajectories that start in the batch time period
    fields_by_batch = []
    for file, batch_time_range in zip(output_file_list, batch_time_ranges, strict=True):
        # Skip batches with no tracks, for which no output file was written
        if not file.is_file():
            continue

        batch_start_time, batch_end_time = batch_time_range
        batch_fields = cf.read(str(file))  # type: ignore[operator]

        trajectory_times = batch_fields[0].coordinate("time").datetime_array[:, 0]
        in_batch = (batch_start_time <= trajectory_times) & (
            trajectory_times < batch_end_time
        )

        if any(in_batch):
            fields_by_batch.append(
                [field.subspace(trajectory=in_batch) for field in batch_fields]
            )

    if not fields_by_batch:
        msg = "No tracks were found in any batch so no output file will be written."
        warnings.warn(msg, category=UserWarning, stacklevel=2)
        return

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


def _get_fields(arg: Any, fields: dict[str, cf.Field]) -> Any:
    """Return the argument with any registry field names replaced by fields.

    Any arguments that do not match the registry keys will be left unchanged.
    """
    if isinstance(arg, str):
        return fields.get(arg, arg)
    if isinstance(arg, (list, tuple)):
        return [
            fields.get(item, item) if isinstance(item, str) else item for item in arg
        ]
    return arg


def _substitute_tags(arg: Any, batch: tuple[int, Path]) -> Any:
    """Substitute %BATCH% and %ITER% tags in a string argument.

    Non-string arguments are returned unchanged. A bare "%ITER%" is replaced
    with an integer rather than a string.
    """
    if not isinstance(arg, str):
        return arg

    i_iter, batch_dir = batch
    if arg == "%ITER%":
        return i_iter  # If just "%ITER%" use an integer instead
    return arg.replace("%BATCH%", str(batch_dir)).replace("%ITER%", str(i_iter))


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

    if len(store_names) > len(result):
        msg = f"Number of fields to store exceeds the number provided by {fn_name}."
        raise ValueError(msg)

    # Prefer storing fields by netcdf variable name
    store_keys = cast(list[str], [name for name in store_names if name is not None])
    result_names = [field.nc_get_variable() for field in result]
    missing_names = [name for name in store_keys if name not in result_names]
    if not missing_names:
        for key in store_keys:
            fields[key] = result[result_names.index(key)]

    # Otherwise, if the number of keys matches the number of fields, store
    # positionally. Warn if any keys match netcdf variable names but not all do.
    elif len(store_names) == len(result):
        if len(missing_names) < len(store_keys):
            msg = (
                f"Some fields to store from {fn_name} match netcdf variable names but "
                f"others do not ({', '.join(missing_names)}), so the fields are "
                "being stored positionally. Check this is intentional."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)
        for name, field in zip(store_names, result, strict=True):
            if name is None:
                continue
            fields[name] = field

    else:
        msg = (
            f"Fields with the following names are not provided by {fn_name}: "
            + ", ".join(missing_names)
        )
        raise ValueError(msg)


def _run_preprocessing(
    step: PreprocessStep, batch: tuple[int, Path], fields: dict[str, cf.Field]
):
    """Run a preprocessing step.

    For partial functions, it first checks if the first argument should be replaced with
    field(s) already in memory. Also replaces any %BATCH% or %ITER% tags in string
    arguments.
    """
    if isinstance(step, tuple):
        if len(step) != 2:  # noqa: PLR2004 - magic number
            msg = (
                "Invalid preprocessing step format. It must be a callable, or a "
                "tuple of (callable, store)."
            )
            raise ValueError(msg)
        fn, store_names = step
    else:
        fn, store_names = step, None

    if not callable(fn):
        msg = "Invalid preprocessing step. It must be callable."
        raise ValueError(msg)

    if isinstance(fn, functools.partial):
        args = list(fn.args)
        # Replace the first argument with any fields from the registry
        if len(args) > 0:
            args[0] = _get_fields(args[0], fields)
        # Substitute the batch tags in any string arguments
        args = [_substitute_tags(arg, batch) for arg in args]
        kwargs = {key: _substitute_tags(arg, batch) for key, arg in fn.keywords.items()}
        # Run the function
        result = fn.func(*args, **kwargs)

    # For a bare function simply run it
    else:
        result = fn()

    # (Optionally) Store the output field(s)
    if result is not None and store_names is not None:
        if isinstance(fn, functools.partial):
            fn_name = getattr(fn.func, "__name__", repr(fn.func))
        else:
            fn_name = fn.__name__
        _store_fields(result, store_names, fields, fn_name)


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
    time_range: TimeRange, interval: Literal["month", "year"] | cf.TimeDuration
) -> list[TimeRange]:
    """Split a time range into time ranges for each batch."""
    if time_range[0] >= time_range[1]:
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
    batch_ranges: list[TimeRange] = []
    batch_start = time_range[0]
    while batch_start < time_range[1]:
        batch_start, batch_end = interval.interval(batch_start)  # type: ignore
        batch_end = min(batch_end, time_range[1])
        batch_ranges.append((batch_start, batch_end))
        batch_start = batch_end

    return batch_ranges


def _expand_batch_time_range(
    time_range: TimeRange, full_time_range: TimeRange, config: "BatchingConfig"
) -> TimeRange:
    """Expand the batch time range, without extending beyond the overall range."""
    buffer_period = config.get("buffer_period")
    start_buffer_period = config.get("start_buffer_period", timedelta(days=1))

    if buffer_period is None:
        return time_range

    batch_start, batch_end = time_range
    overall_start, overall_end = full_time_range

    return (
        max(overall_start, batch_start - start_buffer_period),
        min(overall_end, batch_end + buffer_period),
    )


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
    time_range: TimeRange,
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
        batch_file_default = None if store_fields else Path(filename).name
        batch_file = args.get("batch_file", batch_file_default)
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
    buffer_period: timedelta | None
    """Extra time included at the end of each batch. Default: ``None``."""
    start_buffer_period: timedelta
    """Extra time included at the start of each batch. Defaults to one day when
    ``buffer_period`` is set, otherwise not used."""


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
          to not do so. By default it will use the same file name, or ``None`` if
          ``store`` is set.
        - ``store``: Keys for storing the fields in-memory for use in the preprocessing.
          If there are multiple fields these will be matched by netcdf variable name if
          possible, otherwise they must match the full number of fields, in which case
          they are matched positionally.
        - ``time_varying``: Whether the fields should be subset in time. Default:
          ``True``.
    interval : {"month", "year"} | cf.TimeDuration
        The calendar interval of each batch. The final batch may be shorter to end at
        the end of ``time_range``.
    time_range : tuple[str, str]
        The start and end datetimes for all batches. Must in YYYY-MM-DD format. The end
        time is open (not inclusive).
    preprocessing : Sequence[PreprocessStep] | None
        (optional) The list of preprocessing functions. These can either be simple
        functions or provided using :class:`functools.partial`. The functions in
        :mod:`tctrack.preprocessing` can be used for this.

        Each step can also be provided as a 2-tuple where the second value is
        the key / list of keys for storing the output fields in memory to prevent
        unnecessary file IO. This behaves the same as ``store`` in :attr:`input_files`.
        To use this the function must return a cf.Field, a list of cf.Fields.

        If using `partial`, the first function argument can specify the key(s) of the
        stored fields and it will be replaced by the cf.Field object(s). Any ``%BATCH%``
        or ``%ITER%`` tags in string arguments will also be replaced by the batch
        directory and iteration number respectively.
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
        - buffer_period: ``timedelta`` for extra time at the end of each batch so tracks
          reaching a batch boundary are not cut off. Default: ``None``.
        - start_buffer_period: ``timedelta`` for extra time included at the start of
          each batch. This is to ensure that only tracks starting in the batch period
          are retained. Defaults to one day, but is only used when ``buffer_period`` is
          set.

    Examples
    --------
    Track monthly data from 1950. The input file is loaded, selected to each monthly
    interval, and stored in memory. It is saved to the batch directory after
    preprocessing, which halves the latitude resolution and renames the netCDF variable.

    >>> from tctrack.utils import batching
    >>> from tctrack import tempest_extremes as te
    >>> from tctrack.preprocessing import subsample_field, set_netcdf_info
    >>> from functools import partial
    >>> input_files = [("psl_*.nc", {"store": "psl"})]
    >>> preprocessing = [
    ...     (partial(subsample_field, "psl", {"X": slice(0, None, 2)}), "psl"),
    ...     partial(
    ...         set_netcdf_info,
    ...         "psl",
    ...         nc_name="p",
    ...         output_file="%BATCH%/psl_processed.nc",
    ...     ),
    ... ]
    >>> tracker = te.TETracker()
    >>> batching(
    ...     tracker,
    ...     input_files,
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
    full_time_range: TimeRange = (
        cf.dt(time_range[0], calendar=calendar),
        cf.dt(time_range[1], calendar=calendar),
    )
    batch_time_ranges = _batch_time_ranges(full_time_range, interval)

    # List of output files to use for combining
    output_files: list[Path] = []

    # Perform the batching loop
    for i_iter, batch_time_range in enumerate(batch_time_ranges):
        # Create the directory for the input files
        batch_dir = output_dir / f"batch_{i_iter}"
        batch_dir.mkdir(parents=True, exist_ok=True)

        # Put input files in the batch directory & fields in the preprocessing registry
        buffer_time_range = _expand_batch_time_range(
            batch_time_range, full_time_range, config
        )
        fields = _prepare_inputs(input_files, batch_dir, buffer_time_range)

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
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="There are no trajectories in this period.*"
            )
            tracker.run_tracker(input_file_paths, str(output_file))
        output_files.append(output_file)

        # Optionally delete the batch directory
        if delete_batch_dirs:
            shutil.rmtree(batch_dir)

    if combine_outputs:
        _combine_trajectories(output_files, output_dir / "tracks.nc", batch_time_ranges)
