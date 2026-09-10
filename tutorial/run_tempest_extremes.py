"""Script to run Tempest Extremes as part of the TCTrack tutorial."""

import tctrack.tempest_extremes as te

data_dir = "data_processed/"

input_files = [
    f"{data_dir}/zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
]

# Get the default TempestExtremes parameters, changing the default netcdf variable names
# to match the inputs. This returns two parameter classes:
# - TEDetectParameters
# - TEStitchParameters
nc_names = {
    "msl": "psl",
    "si10": "sfcWind",
    "ghdiff": "_DIFF(zg7h(250hPa),zg7h(500hPa))",
    "latitude": "lat",
    "longitude": "lon",
}
detect_params, stitch_params = te.parameter_set("default", nc_names)

# Define and run the tracker
te_tracker = te.TETracker(detect_params, stitch_params)
te_tracker.run_tracker(input_files, "tracks_tempest_extremes.nc")
