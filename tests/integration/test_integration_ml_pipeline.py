"""Integration test for the MLTracker pipeline.

Runs a small crop of real ERA5 data through preprocess() -> detect() ->
stitch(), then writes the results with both output writers and reads them
back. Shapes and invariants are checked at each step; how *well* the model
detects storms is deliberately not asserted on, since it is trained on full
721x1440 global fields and behaves unreliably on crops this small.

Everything is built fresh on every run - the ERA5 crop is re-cut from the
source files and the statistics re-fetched from the reference repo, into a
temporary directory that is discarded afterwards. Nothing is cached, so the
test can never pass against stale inputs.

Kept small deliberately: a 32x32 crop over 2 timesteps. Two timesteps is the
minimum that exercises cross-timestep track linking, and 32x32 is the smallest
safe size given the U-Net downsamples by 16.

Run from the repo root with:
    conda run -n tctrack-env python tests/integration/test_integration_ml_pipeline.py
"""

# The test inspects the tracker's internal state on purpose - that is what it
# is verifying - so private-member access is expected throughout.
# ruff: noqa: SLF001

import subprocess
import tempfile

import cf
import numpy as np

from tctrack.machine_learning.cyclone_track_ml import MLParameters, MLTracker

RDS_DIR = "/home/sg2147/rds/rds-inspire-tc-TqEGHMWTn8A/sg2147"
PRESSURE_FILE = f"{RDS_DIR}/era5_pressure_2025_1.nc"
SURFACE_FILE = f"{RDS_DIR}/era5_surface_2025_1.nc"

REFERENCE_REPO = "/home/sg2147/tc-track/cyclone-track-ml"
REFERENCE_STATS_REF = "origin/main:data/normalisation_parameters.nc"

MODEL_PATH = (
    "/home/sg2147/.cache/huggingface/hub/models--surbhigoel456--cyclone-TC-ML/"
    "snapshots/529eccbf13ee27056e2f01ec6bd06163133e8923/cyclone-detect-ml-scripted.pt"
)

# Small and cheap: 32x32 grid points over 2 consecutive 6-hourly timesteps.
CENTRE_LAT, CENTRE_LON = -12.6, 51.0
HALF_BOX = 16
TIME_START, N_TIME = 40, 2

N_CHANNELS = 17
N_CLASSES = 5

# The default threshold of 0.5 finds nothing in this window - with five classes
# a winning probability can be as low as ~0.21 - so the output writers are
# exercised at a lower threshold, purely so that there is something to write.
LOW_THRESHOLD = 0.2

RULE = "-" * 66


def section(title: str) -> None:
    """Print a titled section header."""
    print(f"\n{RULE}\n{title}\n{RULE}")


def make_tracker(input_file: str, stats_file: str, **overrides) -> MLTracker:
    """Build a tracker against the prepared inputs, overriding any parameter."""
    return MLTracker(
        MLParameters(
            input_file=input_file,
            model_path=MODEL_PATH,
            normalisation_stats_path=stats_file,
            **overrides,
        )
    )


def build_era5_input(input_file: str) -> None:
    """Build the tracker's input file from the real ERA5 source files.

    Cuts a small window out of the pressure and surface files and merges them
    into one file, since MLTracker expects a single input containing every
    variable it needs, while the ERA5 downloads keep them separate.
    """
    print("  reading real ERA5 files...")
    params = MLParameters()

    # One read of each source file, reused for both the coordinate lookup and
    # the field extraction - these files are several GB each.
    pressure_fields = cf.read(PRESSURE_FILE)
    surface_fields = cf.read(SURFACE_FILE)

    reference = pressure_fields.select_field(params.pressure_variables[0])
    lats = reference.coordinate("latitude").array
    lons = reference.coordinate("longitude").array
    lat_idx = int(np.argmin(np.abs(lats - CENTRE_LAT)))
    lon_idx = int(np.argmin(np.abs(lons - CENTRE_LON)))

    lat_slice = slice(lat_idx - HALF_BOX, lat_idx + HALF_BOX)
    lon_slice = slice(lon_idx - HALF_BOX, lon_idx + HALF_BOX)
    time_slice = slice(TIME_START, TIME_START + N_TIME)

    fields = [
        pressure_fields.select_field(variable)[time_slice, :, lat_slice, lon_slice]
        for variable in params.pressure_variables
    ]
    fields += [
        surface_fields.select_field(variable)[time_slice, lat_slice, lon_slice]
        for variable in (params.sst_variable, params.t2m_variable)
    ]

    cf.write(fields, input_file)
    print(f"  cut {N_TIME}x{2 * HALF_BOX}x{2 * HALF_BOX} window from real ERA5 data")


def fetch_norm_stats(stats_file: str) -> None:
    """Pull the real normalisation statistics from the reference repo."""
    git = ["git", "-C", REFERENCE_REPO]
    subprocess.run([*git, "fetch", "origin"], check=True)  # noqa: S603
    with open(stats_file, "wb") as handle:
        subprocess.run(  # noqa: S603
            [*git, "show", REFERENCE_STATS_REF], stdout=handle, check=True
        )
    print(f"  fetched statistics from {REFERENCE_STATS_REF}")


def test_pipeline(input_file: str, stats_file: str) -> None:
    """Run the pipeline at default settings and check each step's output."""
    section("real ERA5 data through preprocess -> detect -> stitch")

    tracker = make_tracker(input_file, stats_file)

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
        assert trajectory.observations >= tracker.parameters.stitch_min_length
        for key in ("time", "lat", "lon"):
            assert key in trajectory.data, f"trajectory missing '{key}'"
    print(f"  stitch():     {len(trajectories)} trajectories, all well-formed   OK")
    if not trajectories:
        print("                (none found - fine, model quality is not under test)")


def test_output_writers(work_dir: str, input_file: str, stats_file: str) -> None:
    """Check run_tracker() and both writers produce readable CF-NetCDF files.

    Drives the pipeline through run_tracker(), the public entry point, rather
    than calling the steps by hand, so that is covered too.

    Uses a lowered threshold so that detections actually exist: at the default
    the model finds nothing in this small window, and the writers would take
    their empty-input path without the file-writing code ever running.
    """
    section("run_tracker() and the output writers")

    tracker = make_tracker(
        input_file,
        stats_file,
        threshold=LOW_THRESHOLD,
        stitch_min_length=1,  # keep single-timestep tracks, so tracks exist
    )

    # run_tracker() runs detect() and stitch() and writes the trajectories.
    tracks_file = f"{work_dir}/tracks.nc"
    tracker.run_tracker(tracks_file)

    trajectories = tracker.read_trajectories()
    assert tracker._candidates, (
        f"no detections even at threshold={LOW_THRESHOLD}, so the writers "
        "below cannot be exercised"
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
    """Run the integration test against freshly prepared real data."""
    banner = "=" * 66
    print(f"{banner}\nMLTracker pipeline integration test")
    print(f"(checks the code works; does NOT judge model accuracy)\n{banner}")

    # Fresh temporary directory per run, removed on exit.
    with tempfile.TemporaryDirectory(prefix="tctrack_ml_test_") as work_dir:
        input_file = f"{work_dir}/era5_window.nc"
        stats_file = f"{work_dir}/normalisation_parameters.nc"
        build_era5_input(input_file)
        fetch_norm_stats(stats_file)

        test_pipeline(input_file, stats_file)
        test_output_writers(work_dir, input_file, stats_file)

    print(f"\n{banner}\nALL CHECKS PASSED\n{banner}")


if __name__ == "__main__":
    main()
