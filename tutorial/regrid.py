"""Script to pre-process and regrid CMIP6 data for use in TCTrack analysis."""

import os
import shutil

import cf

from tctrack import preprocessing

# Set up file structure
data_dir = "data/"
data_out = "data_processed/"
os.makedirs(data_out, exist_ok=True)


# Define time window for data - ASO 1950
time_bounds = ("1950-08-01", "1950-11-01")
time_window = cf.wi(cf.dt(time_bounds[0]), cf.dt(time_bounds[1]), open_upper=True)


# ======== Tempest Extremes ========
# Extract ASO from annual data files
print("Extracting subspace from psl...", end="", flush=True)
field_psl = preprocessing.combine_time(
    f"{data_dir}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc",
    time_bounds,
    f"{data_out}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)[0]
print("done.")

print("Extracting subspace from sfcWind...", end="", flush=True)
preprocessing.combine_time(
    f"{data_dir}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc",
    time_bounds,
    f"{data_out}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)
print("done.")

# Combine the monthly zg files into one, subspacing in time
print("Combining zg files into one...", end="", flush=True)
# zg is 3hr data from 00:00 but we want daily at 12:00, so set time of lower bound
preprocessing.combine_time(
    f"{data_dir}/zg7h_*.nc",
    (time_bounds[0] + " 12:00", time_bounds[1]),
    f"{data_out}/zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)
print("done.")

# Copy orography across directly:
print("Copying orography file...", end="", flush=True)
shutil.copy(
    f"{data_dir}/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc", f"{data_out}"
)
print("done.")


# ======== TSTORMS ========
print("Renaming slp...", end="", flush=True)
preprocessing.set_netcdf_variable_name(
    field_psl,
    "slp",
    f"{data_out}/slp_day_ASO50.nc",
)
del field_psl
print("done.")

print("Extracting subspace from uas and renaming...", end="", flush=True)
field_uas = preprocessing.combine_time(
    f"{data_dir}/uas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc",
    time_bounds,
)[0]
field_uas = preprocessing.set_netcdf_variable_name(
    field_uas,
    "u_ref",
    f"{data_out}/u_ref_day_ASO50.nc",
)
print("done.")

print("Extracting subspace from vas and renaming...", end="", flush=True)
field_vas = preprocessing.combine_time(
    f"{data_dir}/vas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc",
    time_bounds,
)[0]
field_vas = preprocessing.regrid_to_field(field_vas, field_uas, method="linear")
preprocessing.set_netcdf_variable_name(
    field_vas,
    "v_ref",
    f"{data_out}/v_ref_day_ASO50.nc",
)
del field_vas
print("done.")


print("Extracting subspace from ua and va to calculate vorticity and renaming...")
print("\tExtracting subspace from ua for u850...", end="", flush=True)
field_u850 = preprocessing.subsample_field(
    f"{data_dir}/ua_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc",
    {"T": time_window, "Z": [1]},
    squeeze=True,
)
field_u850 = preprocessing.regrid_to_field(field_u850, field_uas, method="linear")
field_u850 = preprocessing.set_netcdf_variable_name(
    field_u850,
    "u850",
    f"{data_out}/u850_day_ASO50.nc",
)
print("done.")

print("\tExtracting subspace from va for v850...", end="", flush=True)
field_v850 = preprocessing.subsample_field(
    f"{data_dir}/va_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc",
    {"T": time_window, "Z": [1]},
    squeeze=True,
)
field_v850 = preprocessing.regrid_to_field(field_v850, field_uas, method="linear")
del field_uas
field_v850 = preprocessing.set_netcdf_variable_name(
    field_v850,
    "v850",
    f"{data_out}/v850_day_ASO50.nc",
)
print("done.")

print("\tCalculating vorticity for vort850...", end="", flush=True)
field_vort850 = preprocessing.calculate_vorticity(field_u850, field_v850)
del field_u850
del field_v850
field_vort850 = preprocessing.replace_fill_value(field_vort850, 0.0)
preprocessing.set_netcdf_variable_name(
    field_vort850,
    "vort850",
    f"{data_out}/vort850_day_ASO50.nc",
)
del field_vort850
print("done.")

print("done.")

print("Extracting subspace and taking mean of ta and renaming...", end="", flush=True)
field_ta = preprocessing.subsample_field(
    f"{data_dir}/ta_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc",
    {"T": time_window, "Z": slice(3, -3)},
)
field_ta = preprocessing.collapse_field(field_ta, "mean", axes="Z")
preprocessing.set_netcdf_variable_name(
    field_ta,
    "tm",
    f"{data_out}/tm_day_ASO50.nc",
)
del field_ta
print("done.")
