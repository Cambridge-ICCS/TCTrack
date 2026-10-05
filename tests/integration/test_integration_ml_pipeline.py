"""Integration test for the MLTracker pipeline.

Runs the small ERA5 sample in ``data/machine_learning`` through preprocess() ->
detect() -> stitch(), then writes the results with both output writers and reads
them back. Shapes and invariants are checked at each step; how *well* the model
detects storms is deliberately not asserted on, since it is trained on full
721x1440 global fields and behaves unreliably on small crops.

The sample is an 80x80 crop over 10 consecutive 6-hourly timesteps, so it is large
enough to link detections across timesteps. The normalisation statistics bundled
with the package are used. The model is downloaded from the HuggingFace Hub.

Run from the repo root with:
    conda run -n tctrack-env python tests/integration/test_integration_ml_pipeline.py
"""

# The test inspects the tracker's internal state on purpose - that is what it
# is verifying - so private-member access is expected throughout.
# ruff: noqa: SLF001

import tempfile
from pathlib import Path

import cf
import numpy as np

from tctrack.machine_learning.cyclone_track_ml import MLParameters, MLTracker

SAMPLE_FILE = (
    Path(__file__).parents[2]
    / "data"
    / "machine_learning"
    / "era5_dikeledi_2025-01-10.nc"
)

N_TIME = 10
N_CHANNELS = 17
N_CLASSES = 5

RULE = "-" * 66


def section(title: str) -> None:
    """Print a titled section header."""
    print(f"\n{RULE}\n{title}\n{RULE}")


def make_tracker() -> MLTracker:
    """Build a tracker on the sample file with the default parameters."""
    return MLTracker(MLParameters(input_file=str(SAMPLE_FILE)))


def test_pipeline() -> None:
    """Run the pipeline at default settings and check each step's output."""
    section("sample ERA5 data through preprocess -> detect -> stitch")

    tracker = make_tracker()

    tensor = tracker.preprocess()
    n_lat, n_lon = len(tracker._lats), len(tracker._lons)
    assert tuple(tensor.shape) == (N_CHANNELS, N_TIME, n_lat, n_lon), (
        f"unexpected tensor shape {tuple(tensor.shape)}"
    )
    assert not np.isnan(tensor.numpy()).any(), "input tensor contains NaNs"
    assert len(tracker._times) == N_TIME, "wrong number of timesteps read"
    print(f"  preprocess(): {tuple(tensor.shape)}, no NaNs                    OK")

    tracker.detect()
    expected = N_TIME * n_lat * n_lon
    assert len(tracker._scores) == expected, (
        f"expected {expected} scores, got {len(tracker._scores)}"
    )
    # detect() builds every score in one loop, so checking one covers the shape.
    assert set(tracker._scores[0]) == {"time", "lat", "lon", "probs"}
    assert len(tracker._scores[0]["probs"]) == N_CLASSES
    sums = np.array([sum(score["probs"]) for score in tracker._scores])
    assert np.abs(sums - 1.0).max() < 1e-4, "softmax outputs must sum to 1"
    detected = sum(
        1 for score in tracker._scores if int(np.argmax(score["probs"])) != 0
    )
    print(f"  detect():     {len(tracker._scores)} scores, softmax valid       OK")
    print(f"                {detected} non-background pixels (not asserted on)")

    trajectories = tracker.stitch()
    assert tracker.read_trajectories() == trajectories, (
        "read_trajectories() must return what stitch() produced, else to_netcdf() "
        "would write nothing"
    )
    for trajectory in trajectories:
        assert trajectory.observations >= tracker.stitch_parameters.min_length
        for key in ("time", "lat", "lon"):
            assert key in trajectory.data, f"trajectory missing '{key}'"
    print(f"  stitch():     {len(trajectories)} trajectories, all well-formed   OK")
    if not trajectories:
        print("                (none found - fine, model quality is not under test)")


def test_output_writers(work_dir: str) -> None:
    """Check run_tracker() and both writers produce readable CF-NetCDF files.

    Drives the pipeline through run_tracker(), the public entry point, rather
    than calling the steps by hand, so that is covered too.
    """
    section("run_tracker() and the output writers")

    tracker = make_tracker()

    # run_tracker() runs detect() and stitch() and writes the trajectories.
    tracks_file = f"{work_dir}/tracks.nc"
    tracker.run_tracker(tracks_file)

    trajectories = tracker.read_trajectories()
    assert tracker._candidates, (
        "no detections in the sample, so the writers below cannot be exercised"
    )
    assert trajectories, "run_tracker() produced no trajectories"
    print(
        f"  run_tracker(): {len(tracker._candidates)} detections, "
        f"{len(trajectories)} trajectories                OK"
    )

    # -- detections: CF point layout ---------------------------------------
    detections_file = f"{work_dir}/detections.nc"
    tracker.detections_to_netcdf(detections_file)
    fields = cf.read(detections_file)
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
    print(f"  detections_to_netcdf(): {len(fields)} variables, point layout    OK")

    # -- trajectories: CF trajectory layout, written by run_tracker() ------
    track_fields = cf.read(tracks_file)
    assert track_fields, "trajectory file has no variables"
    assert track_fields[0].get_property("featureType") == "trajectory"
    print(
        f"  to_netcdf():            {len(track_fields)} variables, "
        "trajectory layout OK"
    )


def main() -> None:
    """Run the integration test on the sample data."""
    banner = "=" * 66
    print(f"{banner}\nMLTracker pipeline integration test")
    print(f"(checks the code works; does NOT judge model accuracy)\n{banner}")

    # Fresh temporary directory per run, removed on exit.
    with tempfile.TemporaryDirectory(prefix="tctrack_ml_test_") as work_dir:
        test_pipeline()
        test_output_writers(work_dir)

    print(f"\n{banner}\nALL CHECKS PASSED\n{banner}")


if __name__ == "__main__":
    main()
