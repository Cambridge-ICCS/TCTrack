"""Module providing functions for batching the tracking algorithms."""

import glob
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, TypeAlias, TypedDict, cast

import cf

from tctrack.core import TCTracker

PreprocessResult: TypeAlias = cf.Field | Sequence[cf.Field] | None


class _PreprocessRegistryArgs(TypedDict, total=False):
    use: str | Sequence[str]
    store: str | Sequence[str | None]


PreprocessFn: TypeAlias = Callable[..., PreprocessResult]
PreprocessStep: TypeAlias = (
    tuple[PreprocessFn, dict[str, Any]]
    | tuple[PreprocessFn, dict[str, Any], _PreprocessRegistryArgs]
)


def _combine_trajectories(output_file_list: list[Path], output_file: Path) -> None:
    """Combine batched trajectory outputs into a single NetCDF file."""


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
                f"Fields with the following names are not returned by {fn_name}: "
                + ", ".join(missing_names)
            )
            raise ValueError(msg)
        for name in store_names:
            fields[name] = result[result_names.index(name)]

    else:
        msg = f"Number of fields to store exceeds the number returned by {fn_name}."
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
        if v == "%ITER%":
            kwargs[k] = i_iter  # If just "%ITER%" use an integer instead
        elif isinstance(v, str):
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


def _parse_input_files(input_files: Iterable[str], batch_dir: Path) -> list[str]:
    """Parse the input files to prepend the batch directory and expand wildcards."""
    input_file_paths: list[str] = []

    for input_file in input_files:
        input_file_path = input_file

        # Prepend the batch directory for relative file paths
        if not Path(input_file).is_absolute():
            input_file_path = str(batch_dir / input_file)

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
    n_iter: int,
    input_files: str | Iterable[str],
    *,
    preprocessing: Sequence[PreprocessStep] | None = None,
    retrieve_data: Callable[[int, Path], None] | None = None,
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
    n_iter : int
        The number of batching iterations to perform.
    input_files : str | Iterable[str]
        The input file names to pass to the tracker. These will be taken relative to the
        batch directory unless absolute file paths are provided. The `*` wildcard can be
        used to match multiple files.
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
          avoid unnecessary file IO. If ``store`` is used and the function returns more
          than one field it must either match the full number of fields or match the
          netcdf variable names.
    retrieve_data : Callable[[int, Path], None] | None
        (optional) A user-defined function that is called each iteration to retrieve the
        appropriate data and put it in. The first argument is the batch index, the
        second argument is the batch directory.
    config : BatchingConfig | None
        (optional) A dictionary of additional arguments. Valid keys:

        - output_dir: The location to save the outputs. Default: ``"tctrack_outputs"``.
        - combine_outputs: Whether to combine the outputs from each batch into a single
          ``tracks.nc`` file. Default: ``True``.
        - delete_batch_dirs: Whether to delete the ``batch_[i]/`` directories. Default:
          ``True``.

    Examples
    --------
    Do 10 iterations. The input files are copied to the batch directory. The data is
    then preprocessed to halve the latitude resolution and rename the netcdf variable.

    >>> from tctrack.utils import batching
    >>> from tctrack import tempest_extremes as te
    >>> from tctrack.preprocessing import subsample_field, set_nc_variable_name
    >>> def copy_files(i_iter, batch_dir):
    ...     year = 1950 + i_iter
    ...     file = Path(f"psl_{year}.nc")
    ...     file.copy_into(batch_dir)
    >>> preprocessing = [
    ...     (
    ...         subsample_field,
    ...         {"input": "%BATCH%/psl_*.nc", "X": slice(0, None, 2)},
    ...         {"store": "psl"},
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
    ...     10,
    ...     "psl_processed.nc",
    ...     retrieve_data=copy_files,
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

    if isinstance(input_files, str):
        input_files = [input_files]

    # List of output files to use for combining
    output_files: list[Path] = []

    # Perform the batching loop
    for i_iter in range(n_iter):
        # Create the directory for the input files
        batch_dir = output_dir / f"batch_{i_iter}"
        batch_dir.mkdir(parents=True, exist_ok=True)

        # Put the relevant data in the batch directory
        if retrieve_data is not None:
            retrieve_data(i_iter, batch_dir)

        # Preprocess the data
        fields: dict[str, cf.Field] = {}  # Fields to keep in memory across steps
        for preprocessing_step in preprocessing or ():
            _run_preprocessing(
                preprocessing_step,
                batch=(i_iter, batch_dir),
                fields=fields,
            )

        # Run the tracker and keep track of the output files
        input_file_paths = _parse_input_files(input_files, batch_dir)
        output_file = output_dir / f"tracks_{i_iter}.nc"
        tracker.run_tracker(input_file_paths, str(output_file))
        output_files.append(output_file)

        # Optionally delete the batch directory
        if delete_batch_dirs:
            batch_dir.rmdir()

    if combine_outputs:
        _combine_trajectories(output_files, output_dir / "tracks.nc")
