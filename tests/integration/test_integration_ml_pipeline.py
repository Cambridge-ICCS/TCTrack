"""Integration test for the MLTracker pipeline: preprocess -> detect -> stitch.

A small crop of real ERA5 data is run through preprocess() -> detect() ->
stitch(), checking shapes and invariants at each step. Finding zero storms is
an acceptable outcome.

Everything is built fresh on every run - the ERA5 crop is re-cut from the
source files and the statistics re-fetched from the reference repo, into a
temporary directory that is discarded afterwards. Kept small deliberately: a 32x32 crop over 2 timesteps. Two timesteps is the
minimum that exercises cross-timestep track linking, and 32x32 is the smallest
safe size given the U-Net downsamples by 16.

Run from the repo root with:
    conda run -n tctrack-env python tests/integration/test_integration_ml_pipeline.py
"""

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
# Two timesteps is the minimum that exercises cross-timestep track linking.
CENTRE_LAT, CENTRE_LON = -12.6, 51.0
HALF_BOX = 16
TIME_START, N_TIME = 40, 2

N_CHANNELS = 17
N_CLASSES = 5

# The default threshold of 0.5 finds nothing in this window - with five classes
# a winning probability can be as low as ~0.21 - so the output writers are
# exercised at a lower threshold, purely so that there is something to write.
LOW_THRESHOLD = 0.2


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
    subprocess.run(["git", "fetch", "origin"], cwd=REFERENCE_REPO, check=True)
    with open(stats_file, "wb") as handle:
        subprocess.run(
            ["git", "show", REFERENCE_STATS_REF],
            cwd=REFERENCE_REPO,
            stdout=handle,
            check=True,
        )
    print(f"  fetched statistics from {REFERENCE_STATS_REF}")


def test_real_data(input_file: str, stats_file: str) -> None:
    """Run the real pipeline and check every step is structurally sound."""
    print("\n" + "-" * 66)
    print("real ERA5 data through preprocess -> detect -> stitch")
    print("-" * 66)

    params = MLParameters(
        input_file=input_file,
        model_path=MODEL_PATH,
        normalisation_stats_path=stats_file,
    )
    tracker = MLTracker(params)

    # -- preprocess --------------------------------------------------------
    tensor = tracker.preprocess()
    n_lat, n_lon = len(tracker._lats), len(tracker._lons)
    assert tuple(tensor.shape) == (N_CHANNELS, N_TIME, n_lat, n_lon), (
        f"unexpected tensor shape {tuple(tensor.shape)}"
    )
    assert not np.isnan(tensor.numpy()).any(), "input tensor contains NaNs"
    assert len(tracker._times) == N_TIME, "wrong number of timesteps read"
    print(f"  preprocess(): {tuple(tensor.shape)}, no NaNs                    OK")

    # -- detect ------------------------------------------------------------
    tracker.detect()
    expected = N_TIME * n_lat * n_lon
    assert len(tracker._scores) == expected, (
        f"expected {expected} scores, got {len(tracker._scores)}"
    )
    for score in tracker._scores[:100]:
        assert set(score) == {"time", "lat", "lon", "probs"}, "unexpected score keys"
        assert len(score["probs"]) == N_CLASSES, "expected 5 class probabilities"
    sums = np.array([sum(s["probs"]) for s in tracker._scores])
    assert np.abs(sums - 1.0).max() < 1e-4, "softmax outputs must sum to 1"
    detected = sum(1 for s in tracker._scores if int(np.argmax(s["probs"])) != 0)
    print(f"  detect():     {len(tracker._scores)} scores, softmax valid       OK")
    print(f"                {detected} non-background pixels (not asserted on)")

    # -- stitch ------------------------------------------------------------
    trajectories = tracker.stitch()
    assert isinstance(trajectories, list), "stitch() must return a list"
    assert tracker.read_trajectories() == trajectories, (
        "read_trajectories() must return what stitch() produced, else to_netcdf() "
        "would write nothing"
    )
    for trajectory in trajectories:
        assert trajectory.observations >= params.stitch_min_length
        for key in ("time", "lat", "lon"):
            assert key in trajectory.data, f"trajectory missing '{key}'"
    print(f"  stitch():     {len(trajectories)} trajectories, all well-formed   OK")
    if not trajectories:
        print("                (none found - fine, model quality is not under test)")


def test_outputs(work_dir: str, input_file: str, stats_file: str) -> None:
    """Check both output writers produce readable CF-NetCDF files.

    Uses a lowered threshold so that detections actually exist: at the default
    the model finds nothing in this small window, and the writers would take
    their empty-input path without the file-writing code ever running.
    """
    print("\n" + "-" * 66)
    print("output writers: detections_to_netcdf() and to_netcdf()")
    print("-" * 66)

    tracker = MLTracker(
        MLParameters(
            input_file=input_file,
            model_path=MODEL_PATH,
            normalisation_stats_path=stats_file,
            threshold=LOW_THRESHOLD,
            stitch_min_length=1,  # keep single-timestep tracks, so tracks exist
        )
    )
    tracker.detect()
    trajectories = tracker.stitch()
    assert tracker._candidates, (
        f"no detections even at threshold={LOW_THRESHOLD}; the writers below "
        "cannot be exercised"
    )
    assert trajectories, "no trajectories to write"
    print(f"  {len(tracker._candidates)} detections, {len(trajectories)} trajectories to write")

    # -- detections: CF point layout ---------------------------------------
    detections_file = f"{work_dir}/detections.nc"
    tracker.detections_to_netcdf(detections_file)
    fields = cf.read(detections_file)

    written = {field.nc_get_variable("?"): field for field in fields}
    assert len(written) == len(fields), "netCDF variable names collided"
    for required in ("class_index", "score", "sea_surface_temperature"):
        assert required in written, f"{required} missing from the detections file"
    assert written["sea_surface_temperature"].get_property("standard_name") == (
        "sea_surface_temperature"
    ), "co-located variables should carry a CF standard name"
    assert written["class_index"].get_property("flag_meanings"), (
        "the categorical class should carry CF flag attributes"
    )

    example = next(iter(written.values()))
    assert example.get_property("featureType") == "point", (
        "detections are points, not trajectories"
    )
    for coordinate in ("time", "latitude", "longitude"):
        assert example.coordinate(coordinate) is not None, f"missing {coordinate}"

    # values should survive the round trip unchanged
    from_file = written["sea_surface_temperature"].array.tolist()
    in_memory = [c["sea_surface_temperature"] for c in tracker._candidates]
    assert np.allclose(from_file, in_memory), "values changed on write/read"
    print(f"  detections_to_netcdf(): {len(fields)} variables, point layout    OK")

    # -- trajectories: CF trajectory layout --------------------------------
    tracks_file = f"{work_dir}/tracks.nc"
    tracker.to_netcdf(tracks_file)
    track_fields = cf.read(tracks_file)
    assert track_fields, "trajectory file has no variables"
    assert track_fields[0].get_property("featureType") == "trajectory"
    print(f"  to_netcdf():            {len(track_fields)} variables, trajectory layout OK")


def main() -> None:
    """Run the integration test against freshly prepared real data."""
    print("=" * 66)
    print("MLTracker pipeline integration test")
    print("(checks the code works; does NOT judge model accuracy)")
    print("=" * 66)

    # Fresh temporary directory per run, removed on exit.
    with tempfile.TemporaryDirectory(prefix="tctrack_ml_test_") as work_dir:
        input_file = f"{work_dir}/era5_window.nc"
        stats_file = f"{work_dir}/normalisation_parameters.nc"
        build_era5_input(input_file)
        fetch_norm_stats(stats_file)

        test_real_data(input_file, stats_file)
        test_outputs(work_dir, input_file, stats_file)

    print("\n" + "=" * 66)
    print("ALL CHECKS PASSED")
    print("=" * 66)


if __name__ == "__main__":
    main()
