"""Script to run TSTORMS as part of the TCTrack tutorial."""

import os

from tctrack import tstorms

input_files = [
    "data_processed/u_ref_day_ASO50.nc",
    "data_processed/v_ref_day_ASO50.nc",
    "data_processed/vort850_day_ASO50.nc",
    "data_processed/tm_day_ASO50.nc",
    "data_processed/slp_day_ASO50.nc",
]

# Get the default TSTORMS parameters
# This returns three parameter classes:
# - TSTORMSBaseParameters
# - TSTORMSDetectParameters
# - TSTORMSStitchParameters
tstorms_params, detect_params, stitch_params = tstorms.parameter_set(
    "default",
    tstorms_dir=f"{os.getcwd()}/TSTORMS",
    output_dir="tstorms_outputs",
)

# Define and run the tracker
tstorms_tracker = tstorms.TSTORMSTracker(tstorms_params, detect_params, stitch_params)
tstorms_tracker.run_tracker(input_files, "tracks_tstorms.nc")
