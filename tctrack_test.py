"""Short script to test TCTrack functionality as we develop."""

from tctrack import tempest_extremes as te

input_files = [
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/psl_E3hr_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001010300-195006302100.nc",
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/zgdiff_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001-195003.nc",
    "/rds/project/rds-TqEGHMWTn8A/test_data/hist-1950/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
]

closed_contours = [
    te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
    te.TEContour(var="zgdiff", delta=-6.0, dist=6.5, minmaxdist=1.0),
]

output_commands = [
    te.TEOutputCommand(var="psl", operator="min", dist=0.0),
    te.TEOutputCommand(var="orog", operator="max", dist=0.0),
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

print(dn_params)

tracker = te.TETracker(dn_params)

result = tracker.detect_nodes()

print(result["returncode"])
