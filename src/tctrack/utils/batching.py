"""Module providing functions for batching the tracking algorithms."""

from collections.abc import Callable
from pathlib import Path

from tctrack.core import TCTracker


def _combine_trajectories(output_file_list: list[str], output_file: str):
    pass


def batching(
    tracker: TCTracker,
    n_iter: int,
    preprocessing: Callable[[int, TCTracker], TCTracker],
    retrieve_data: Callable[[int], None] | None = None,
    output_dir: str | Path = "tctrack_outputs",
    combine_outputs: bool = True,
):
    """Utility to batch the tracking and combine the outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i_iter in range(n_iter):
        if retrieve_data is not None:
            retrieve_data(i_iter)
        tracker = preprocessing(i_iter, tracker)

        output_file = output_dir / f"tracks_{i_iter}.nc"
        tracker.run_tracker(str(output_file))

    if combine_outputs:
        _combine_trajectories(
            [output_dir / f"tracks_{i}.nc" for i in range(n_iter)],
            output_dir / "tracks.nc",
        )
