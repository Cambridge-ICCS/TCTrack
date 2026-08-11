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


def test_real_data(work_dir: str) -> None:
    """Run the real pipeline and check every step is structurally sound."""
    print("\n" + "-" * 66)
    print("real ERA5 data through preprocess -> detect -> stitch")
    print("-" * 66)

    input_file = f"{work_dir}/era5_window.nc"
    stats_file = f"{work_dir}/normalisation_parameters.nc"
    build_era5_input(input_file)
    fetch_norm_stats(stats_file)

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


def main() -> None:
    """Run the integration test against freshly prepared real data."""
    print("=" * 66)
    print("MLTracker pipeline integration test")
    print("(checks the code works; does NOT judge model accuracy)")
    print("=" * 66)

    # Fresh temporary directory per run, removed on exit - nothing is cached,
    # so the test can never pass against inputs left over from a previous run.
    with tempfile.TemporaryDirectory(prefix="tctrack_ml_test_") as work_dir:
        test_real_data(work_dir)

    print("\n" + "=" * 66)
    print("ALL CHECKS PASSED")
    print("=" * 66)


if __name__ == "__main__":
    main()
