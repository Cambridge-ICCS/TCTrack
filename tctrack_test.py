"""Short script to test TCTrack functionality as we develop."""

from tctrack import tempest_extremes as te

input_files = [
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/psl_E3hr_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001010300-195006302100.nc",
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/zgdiff_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001-195003.nc",
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
]

closed_contours = [
    {"var": "psl", "delta": 200.0, "dist": 5.5, "minmaxdist": 0},
    {"var": "zgdiff", "delta": -6.0, "dist": 6.5, "minmaxdist": 1.0},
]

output_commands = [
    {"var": "psl", "operation": "min", "dist": 0.0},
    {"var": "orog", "operation": "max", "dist": 0.0},
]

dn_params = te.DetectNodesParameters(
    in_data=input_files,
    search_by_min="psl",
    merge_dist=6.0,
    closed_contours=closed_contours,
    out_header=True,
    output_commands=output_commands,
    output_file="/rds/project/rds-TqEGHMWTn8A/test_data/test_out/nodes_jack.dat",
)

threshold_filters = [
    {"var": "lat", "operation": "<=", "value": 40, "count": 10},
    {"var": "lat", "operation": ">=", "value": -40, "count": 10},
    {"var": "orog", "operation": "<=", "value": 1500, "count": 10},
    {"var": "orog", "operation": "<=", "value": 10, "count": 4},
]

sn_params = te.StitchNodesParameters(
    input_file="/rds/project/rds-TqEGHMWTn8A/test_data/test_out/nodes_jack.dat",
    in_fmt="lon,lat,psl,orog",
    output_file="/rds/project/rds-TqEGHMWTn8A/test_data/test_out/tracks_tctrack.txt",
    max_sep=8.0,
    min_time=10,
    max_gap=3,
    min_endpoint_dist=8.0,
    threshold_filters=threshold_filters,
)

print(dn_params)
print(sn_params)

tracker = te.TETracker(dn_params, sn_params)

result = tracker.detect_nodes()
tracker.stitch_nodes()
