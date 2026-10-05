"""Integration test for the MLTracker pipeline.

Runs the small ERA5 sample in ``data/machine_learning`` through preprocess() ->
detect() -> stitch(), then writes the results with both output writers and reads
them back. Shapes and invariants are checked at each step; how *well* the model
detects storms is deliberately not asserted on, since it is trained on full
721x1440 global fields and behaves unreliably on small crops.

The sample is an 80x80 crop over 10 consecutive 6-hourly timesteps. The normalisation
statistics bundled with the package are used. The model is downloaded from the
HuggingFace Hub, and the tests are skipped if it cannot be obtained, for example
without network access or without the ``HF_TOKEN`` environment variable if
the repository requires a token.

Run from the repo root with:
    pytest tests/integration
"""

# The test inspects the tracker's internal state on purpose - that is what it
# is verifying - so private-member access is expected throughout.
# ruff: noqa: SLF001

from pathlib import Path

import cf
import numpy as np
import pytest
from huggingface_hub.errors import HfHubHTTPError

from tctrack.machine_learning import MLParameters, MLTracker

SAMPLE_FILE = (
    Path(__file__).parents[2]
    / "data"
    / "machine_learning"
    / "era5_dikeledi_2025-01-10.nc"
)

N_TIME = 10
N_CHANNELS = 17
N_CLASSES = 5


@pytest.fixture
def tracker() -> MLTracker:
    """Build a tracker on the sample file, skipping the test if the model is missing.

    The tracker downloads the model from the HuggingFace Hub when it is constructed.
    """
    try:
        return MLTracker(MLParameters(input_file=str(SAMPLE_FILE)))
    except (OSError, HfHubHTTPError) as error:
        pytest.skip(f"Requires the model from the HuggingFace Hub: {error}")


def test_pipeline(tracker: MLTracker) -> None:
    """Run the pipeline at default settings and check each step's output."""
    tensor = tracker.preprocess()
    n_lat, n_lon = len(tracker._lats), len(tracker._lons)
    assert tuple(tensor.shape) == (N_CHANNELS, N_TIME, n_lat, n_lon)
    assert not np.isnan(tensor.numpy()).any(), "input tensor contains NaNs"
    assert len(tracker._times) == N_TIME, "wrong number of timesteps read"

    tracker.detect()
    assert len(tracker._scores) == N_TIME * n_lat * n_lon
    # detect() builds every score in one loop, so checking one covers the shape.
    assert set(tracker._scores[0]) == {"time", "lat", "lon", "probs"}
    assert len(tracker._scores[0]["probs"]) == N_CLASSES
    sums = np.array([sum(score["probs"]) for score in tracker._scores])
    assert np.abs(sums - 1.0).max() < 1e-4, "softmax outputs must sum to 1"

    trajectories = tracker.stitch()
    assert trajectories, "stitch() found no trajectories in the sample"
    assert tracker.read_trajectories() == trajectories, (
        "read_trajectories() must return what stitch() produced, else to_netcdf() "
        "would write nothing"
    )
    for trajectory in trajectories:
        assert trajectory.observations >= tracker.stitch_parameters.min_length
        for key in ("time", "lat", "lon"):
            assert key in trajectory.data, f"trajectory missing '{key}'"


def test_output_writers(tracker: MLTracker, tmp_path: Path) -> None:
    """Check run_tracker() and both writers produce readable CF-NetCDF files.

    Drives the pipeline through run_tracker(), the public entry point, rather
    than calling the steps by hand, so that is covered too.
    """
    # run_tracker() runs detect() and stitch() and writes the trajectories.
    tracks_file = tmp_path / "tracks.nc"
    tracker.run_tracker(str(tracks_file))

    trajectories = tracker.read_trajectories()
    assert tracker._candidates, (
        "no detections in the sample, so the writers below cannot be exercised"
    )
    assert trajectories, "run_tracker() produced no trajectories"

    # -- detections: CF point layout ---------------------------------------
    detections_file = tmp_path / "detections.nc"
    tracker.detections_to_netcdf(str(detections_file))
    fields = cf.read(str(detections_file))
    written = {field.nc_get_variable("?"): field for field in fields}

    assert len(written) == len(fields), "netCDF variable names collided"
    for required in ("class_index", "score", "sea_surface_temperature"):
        assert required in written, f"{required} missing from the detections file"

    sst = written["sea_surface_temperature"]
    assert sst.get_property("standard_name") == "sea_surface_temperature", (
        "co-located variables should carry a CF standard name"
    )
    assert written["class_index"].get_property("flag_meanings"), (
        "the categorical class should carry CF flag attributes"
    )
    assert sst.get_property("featureType") == "point", (
        "detections are points, not trajectories"
    )
    for coordinate in ("time", "latitude", "longitude"):
        assert sst.coordinate(coordinate) is not None, f"missing {coordinate}"

    in_memory = [c["sea_surface_temperature"] for c in tracker._candidates]
    assert np.allclose(sst.array.tolist(), in_memory), "values changed on write/read"

    # -- trajectories: CF trajectory layout, written by run_tracker() ------
    track_fields = cf.read(str(tracks_file))
    assert track_fields, "trajectory file has no variables"
    assert track_fields[0].get_property("featureType") == "trajectory"
