"""Script to run TSTORMS as part of the TCTrack tutorial."""

import os

from tctrack import tstorms

tstorms_params = tstorms.TSTORMSBaseParameters(
    tstorms_dir=f"{os.getcwd()}/TSTORMS",
    output_dir="tstorms_outputs",
)

detect_params = tstorms.TSTORMSDetectParameters(
    vort_crit=3.5e-5,
    tm_crit=0.0,
    thick_crit=50.0,
    dist_crit=4.0,
    lat_bound_n=70.0,
    lat_bound_s=-70.0,
    do_spline=False,
    do_thickness=False,
    use_sfc_wind=True,
)

stitch_params = tstorms.TSTORMSStitchParameters(
    r_crit=900.0,
    wind_crit=15.0,
    vort_crit=3.5e-5,
    tm_crit=0.0,
    n_day_crit=1,
    do_filter=True,
    lat_bound_n=70.0,
    lat_bound_s=-70.0,
)

input_files = [
    "data_processed/u_ref_day_ASO50.nc",
    "data_processed/v_ref_day_ASO50.nc",
    "data_processed/vort850_day_ASO50.nc",
    "data_processed/tm_day_ASO50.nc",
    "data_processed/slp_day_ASO50.nc",
]

tstorms_tracker = tstorms.TSTORMSTracker(tstorms_params, detect_params, stitch_params)
tstorms_tracker.run_tracker(input_files, "tracks_tstorms.nc")

# tstorms_tracker.stitch(verbosity=2)
# tstorms_tracker.to_netcdf("tracks_tstorms.nc")
