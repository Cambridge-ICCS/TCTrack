"""Script to pre-process and regrid CMIP6 data for use in TCTrack analysis."""

import os
import shutil

import cf

# Set up file structure
data_dir = "data/"
data_out = "data_processed/"
os.makedirs(data_out, exist_ok=True)


# Define time window for data - ASO 1950
ASO_1950 = cf.wi(cf.dt("1950-08-01"), cf.dt("1950-11-01"), open_upper=True)


# ======== Tempest Extremes ========
# Extract ASO from annual data files
print("Extracting subspace from psl...", end="", flush=True)
field_psl = cf.read(
    f"{data_dir}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc"
)[0]
field_psl = field_psl.subspace(T=ASO_1950)
print("writing data...", end="", flush=True)
cf.write(
    field_psl,
    f"{data_out}/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)
print("done.")

print("Extracting subspace from sfcWind...", end="", flush=True)
field_sfcWind = cf.read(
    f"{data_dir}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc"
)[0]
field_sfcWind = field_sfcWind.subspace(T=ASO_1950)
print("writing data...", end="", flush=True)
cf.write(
    field_sfcWind,
    f"{data_out}/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)
del field_sfcWind
print("done.")

# Combine the monthly zg files into one, subspacing in time to the above resolutions:
print("Combining zg files into one...", end="", flush=True)
input_files = [
    f"{data_dir}/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195008010000-195008302100.nc",
    f"{data_dir}/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195009010000-195009302100.nc",
    f"{data_dir}/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195010010000-195010302100.nc",
]

field_zg_in = cf.read(input_files)[0]
# zg is 3hr data from 00:00 but we want daily at 12:00, so subspace with a slice
field_zg = field_zg_in.subspace(T=slice(4, None, 8))
del field_zg_in
print("writing data...", end="", flush=True)
cf.write(
    field_zg,
    f"{data_out}/zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
)
del field_zg
print("done.")

# Copy orography across directly:
print("Copying orography file...", end="", flush=True)
shutil.copy(
    f"{data_dir}/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc", f"{data_out}"
)
print("done.")


# ======== TSTORMS ========
print("Renaming slp...", end="", flush=True)
field_psl.nc_set_variable("slp")
print("writing data...", end="", flush=True)
cf.write(
    field_psl,
    f"{data_out}/slp_day_ASO50.nc",
)
del field_psl
print("done.")

print("Extracting subspace from uas and renaming...", end="", flush=True)
field_uas_in = cf.read(
    f"{data_dir}/uas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc"
)[0]
field_uas = field_uas_in.subspace(T=ASO_1950)
del field_uas_in
field_uas.nc_set_variable("u_ref")
print("writing data...", end="", flush=True)
cf.write(
    field_uas,
    f"{data_out}/u_ref_day_ASO50.nc",
)
print("done.")

print("Extracting subspace from vas and renaming...", end="", flush=True)
field_vas_in = cf.read(
    f"{data_dir}/vas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc"
)[0]
field_vas = field_vas_in.subspace(T=ASO_1950)
field_vas.regrids(field_uas, method="linear", inplace=True)
field_vas.nc_clear_dataset_chunksizes()
del field_vas_in
field_vas.nc_set_variable("v_ref")
print("writing data...", end="", flush=True)
cf.write(
    field_vas,
    f"{data_out}/v_ref_day_ASO50.nc",
)
del field_vas
print("done.")


print("Extracting subspace from ua and va to calculate vorticity and renaming...")
print("\tExtracting subspace from ua for u850...", end="", flush=True)
field_ua = cf.read(
    f"{data_dir}/ua_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc"
)[0]
field_u850 = field_ua.subspace(T=ASO_1950, Z=[1])
del field_ua
field_u850.squeeze(inplace=True)
field_u850_rg = field_u850.regrids(field_uas, method="linear")
del field_u850
field_u850_rg.nc_clear_dataset_chunksizes()
field_u850_rg.nc_set_variable("u850")
print("writing data...", end="", flush=True)
cf.write(
    field_u850_rg,
    f"{data_out}/u850_day_ASO50.nc",
)
print("done.")

print("\tExtracting subspace from va for v850...", end="", flush=True)
field_va = cf.read(
    f"{data_dir}/va_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc"
)[0]
field_v850 = field_va.subspace(T=ASO_1950, Z=[1])
del field_va
field_v850.squeeze(inplace=True)
field_v850_rg = field_v850.regrids(field_uas, method="linear")
del field_uas
del field_v850
field_v850_rg.nc_clear_dataset_chunksizes()
field_v850_rg.nc_set_variable("v850")
print("writing data...", end="", flush=True)
cf.write(
    field_v850_rg,
    f"{data_out}/v850_day_ASO50.nc",
)
print("done.")

print("\tCalculating vorticity for vort850...", end="", flush=True)
field_vort850 = cf.curl_xy(field_u850_rg, field_v850_rg, radius="earth")
del field_u850_rg
del field_v850_rg
field_vort850.filled(fill_value=0.0, inplace=True)
field_vort850.nc_set_variable("vort850")
field_vort850.set_property("standard_name", "atmosphere_upward_absolute_vorticity")
field_vort850.set_property("units", "s-1")
print("writing data...", end="", flush=True)
cf.write(
    field_vort850,
    f"{data_out}/vort850_day_ASO50.nc",
)
del field_vort850
print("done.")

print("done.")

# # Combine the monthly vortmean files into one, subspacing in time to the above resolutions:
# print("Combining vortmean files into one...", end="", flush=True)
# input_files = [
#     f"{data_dir}/vortmean_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195008010000-195008302100.nc",
#     f"{data_dir}/vortmean_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195009010000-195009302100.nc",
#     f"{data_dir}/vortmean_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195010010000-195010302100.nc",
# ]
#
# field_vort = cf.read(input_files)[0]
# # vortmean is 3hr data from 00:00 but we want daily at 12:00, so subspace with a slice
# field_vort = field_vort.subspace(T=slice(4, None, 8))
# print("writing data...", end="", flush=True)
# cf.write(
#     field_zg,
#     f"{data_out}/zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
# )
# del field_zg
# print("done.")

print("Extracting subspace and taking mean of ta and renaming...", end="", flush=True)
field_ta_full = cf.read(
    f"{data_dir}/ta_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc"
)[0]
# Extract 500 and 250 pressure levels and take mean
field_ta = field_ta_full.subspace(T=ASO_1950, Z=slice(3, -3))
del field_ta_full
field_ta.collapse("mean", axes="Z", inplace=True)
field_ta.squeeze(inplace=True)
field_ta.nc_set_variable("tm")
print("writing data...", end="", flush=True)
cf.write(
    field_ta,
    f"{data_out}/tm_day_ASO50.nc",
)
del field_ta
print("done.")
