"""Script to run Tempest Extremes as part of the TCTrack tutorial."""

import tctrack.tempest_extremes as te

data_dir = "data_processed/"

input_files = [
    f"{data_dir}/z.nc",
    f"{data_dir}/msl.nc",
    f"{data_dir}/si10.nc",
    f"{data_dir}/orog.nc",
]

# ======== Tempest Extremes Parameters ========
closed_contours = [
    te.TEContour(var="msl", delta=200.0, dist=5.5, minmaxdist=0.0),
    te.TEContour(
        var="_DIFF(z(300hPa),z(500hPa))", delta=-58.8, dist=6.5, minmaxdist=1.0
    ),
]

output_commands = [
    te.TEOutputCommand(var="msl", operator="min", dist=0.0),
    te.TEOutputCommand(var="orog", operator="max", dist=0.0),
    te.TEOutputCommand(var="si10", operator="max", dist=2.0),
]

threshold_filters = [
    te.TEThreshold(var="lat", op="<=", value=50, count=10),
    te.TEThreshold(var="lat", op=">=", value=-50, count=10),
    te.TEThreshold(var="orog", op="<=", value=150, count=10),
    te.TEThreshold(var="si10", op=">=", value=10, count=10),
]

dn_params = te.TEDetectParameters(
    in_data=input_files,
    search_by_min="msl",
    time_filter="6hr",
    merge_dist=6.0,
    closed_contours=closed_contours,
    lat_name="latitude",
    lon_name="longitude",
    out_header=True,
    output_commands=output_commands,
    output_dir="te_outputs",
)

sn_params = te.TEStitchParameters(
    caltype="standard",
    max_sep=8.0,
    min_time="54h",
    max_gap="24h",
    min_endpoint_dist=8.0,
    threshold_filters=threshold_filters,
)

# ======== Run Tempest extremes ========
te_tracker = te.TETracker(dn_params, sn_params)

te_tracker.run_tracker("tracks_tempest_extremes.nc")
