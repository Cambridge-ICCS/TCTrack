"""Module providing functions for batching the tracking algorithms."""

import glob
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, TypeAlias, TypedDict

import cf

from tctrack.core import TCTracker

PreprocessResult: TypeAlias = cf.Field | Sequence[cf.Field] | None


PreprocessFn: TypeAlias = Callable[..., PreprocessResult]
PreprocessStep: TypeAlias = tuple[PreprocessFn, dict[str, Any]]


def _combine_trajectories(output_file_list: list[Path], output_file: Path) -> None:
    """Combine batched trajectory outputs into a single NetCDF file."""


def _run_preprocessing(step: PreprocessStep, batch_dir: Path) -> PreprocessResult:
    """Run a batching step with supported keyword arguments."""
    fn = step[0]
    kwargs = step[1].copy()

    # Replace any %BATCH% tags in string arguments
    for k, v in kwargs.items():
        if isinstance(v, str):
            kwargs[k] = v.replace("%BATCH%", str(batch_dir))

    # Run the preprocessing function
    result = fn(**kwargs)

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
        for preprocessing_step in preprocessing or ():
            _run_preprocessing(preprocessing_step, batch_dir=batch_dir)

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
