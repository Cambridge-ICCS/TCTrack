"""Script to run Tempest Extremes as part of the TCTrack tutorial."""

import tctrack.tempest_extremes as te

data_dir = "data_processed/"

input_files = [
    f"{data_dir}/zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
    f"{data_dir}/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
]

closed_contours = [
    te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
    te.TEContour(
        var="_DIFF(zg7h(250hPa),zg7h(500hPa))", delta=-6.0, dist=6.5, minmaxdist=1.0
    ),
]

output_commands = [
    te.TEOutputCommand(var="psl", operator="min", dist=0.0),
    te.TEOutputCommand(var="orog", operator="max", dist=0.0),
    te.TEOutputCommand(var="sfcWind", operator="max", dist=2.0),
]

threshold_filters = [
    te.TEThreshold(var="lat", op="<=", value=50, count=10),
    te.TEThreshold(var="lat", op=">=", value=-50, count=10),
    te.TEThreshold(var="orog", op="<=", value=150, count=10),
    te.TEThreshold(var="sfcWind", op=">=", value=10, count=10),
]

dn_params = te.DetectNodesParameters(
    in_data=input_files,
    search_by_min="psl",
    time_filter="6hr",
    merge_dist=6.0,
    closed_contours=closed_contours,
    out_header=True,
    output_commands=output_commands,
    output_dir="te_outputs",
)

sn_params = te.StitchNodesParameters(
    caltype="360_day",
    max_sep=8.0,
    min_time="54h",
    max_gap="24h",
    min_endpoint_dist=8.0,
    threshold_filters=threshold_filters,
)

te_tracker = te.TETracker(dn_params, sn_params)

te_tracker.run_tracker("tracks_tempest_extremes.nc")
