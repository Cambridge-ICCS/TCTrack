"""Run MLTracker.run_tracker() on real ERA5 data and save the results.

Builds a merged input file for window (2025-01-10 to
2025-01-12 (6-hourly), the same region used in the tracking GIF, where the model is known
to produce detections), then runs the tracker end to end and writes both
output files:

Merged inputs are written to the directory (era5_run/).

Detection and Stitched files written to the current working directory:
 - cyclone_trajectories.nc: the linked storm tracks (CF trajectory layout)
 - cyclone_detections.nc: every individual detection (CF point layout)

Run the script from the repo root - as a slurm job,
    `sbatch run_tracker_on_real_data.sh`
"""

import os

import cf
import numpy as np

from tctrack.machine_learning.cyclone_track_ml import MLParameters, MLTracker

RDS_DIR = "/home/sg2147/rds/rds-inspire-tc-TqEGHMWTn8A/sg2147"
PRESSURE_FILE = f"{RDS_DIR}/era5_pressure_2025_1.nc"
SURFACE_FILE = f"{RDS_DIR}/era5_surface_2025_1.nc"
STATS_FILE = f"{RDS_DIR}/normalisation_parameters.nc"
MODEL_PATH = (
    "/home/sg2147/.cache/huggingface/hub/models--surbhigoel456--cyclone-TC-ML/"
    "snapshots/529eccbf13ee27056e2f01ec6bd06163133e8923/cyclone-detect-ml-scripted.pt"
)

RUN_DIR = "era5_run"
INPUT_FILE = f"{RUN_DIR}/input.nc"
TRAJECTORIES_FILE = "cyclone_trajectories.nc"
DETECTIONS_FILE = "cyclone_detections.nc"

CENTRE_LAT, CENTRE_LON = -12.6, 51.0
HALF_BOX = 40
TIME_START, N_TIME = 36, 10

THRESHOLD = 0.25


def build_input() -> None:
    """Crop and merge the run's ERA5 window into a single input file."""
    print("reading ERA5 files...", flush=True)
    params = MLParameters()
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
        pressure_fields.select_field(v)[time_slice, :, lat_slice, lon_slice]
        for v in params.pressure_variables
    ]
    fields += [
        surface_fields.select_field(v)[time_slice, lat_slice, lon_slice]
        for v in (params.sst_variable, params.t2m_variable)
    ]
    cf.write(fields, INPUT_FILE)
    print(f"  wrote {INPUT_FILE}: {N_TIME} timesteps, {2 * HALF_BOX}x{2 * HALF_BOX} grid")


def main() -> None:
    os.makedirs(RUN_DIR, exist_ok=True)
    build_input()

    print("\nrunning tracker...", flush=True)
    tracker = MLTracker(
        MLParameters(
            input_file=INPUT_FILE,
            model_path=MODEL_PATH,
            normalisation_stats_path=STATS_FILE,
            threshold=THRESHOLD,
        )
    )
    tracker.run_tracker(TRAJECTORIES_FILE)
    tracker.detections_to_netcdf(DETECTIONS_FILE)

    print(
        f"\n{len(tracker._candidates)} detections across "  # noqa: SLF001
        f"{N_TIME} timesteps"
    )
    print(f"{len(tracker._trajectories)} trajectories written to {TRAJECTORIES_FILE}")  # noqa: SLF001
    for trajectory in tracker._trajectories:  # noqa: SLF001
        lats = [round(v, 2) for v in trajectory.data["lat"]]
        lons = [round(v, 2) for v in trajectory.data["lon"]]
        print(f"  id={trajectory.trajectory_id}  obs={trajectory.observations}  "
            f"lat={lats}  lon={lons}")
    print(f"detections written to {DETECTIONS_FILE}")


if __name__ == "__main__":
    main()
