"""Script to pre-process ERA5 data for use in TCTrack analysis.

Fields are manually deleted when they are no longer needed to reduce memory usage.
"""

import os

import cf

from tctrack import preprocessing

# Set up file structure
data_dir = "data"
data_out = "data_processed"
os.makedirs(data_out, exist_ok=True)

GRAVITY = 9.80665  # Standard gravitational acceleration [m s-2]

# ======== Tempest Extremes ========
# Copy the geopotential and sea-level pressure across unchanged
preprocessing.read_files(f"{data_dir}/era5_z.nc", output_file=f"{data_out}/z.nc")

field_u10, field_v10, field_msl = preprocessing.separate_variables(
    f"{data_dir}/era5_sfc.nc",
    output_files={"msl": f"{data_out}/msl.nc"},
    return_order=["u10", "v10", "msl"],
)

print("Calculating 10m wind speed from components...", end="", flush=True)
preprocessing.calculate_wind_speed(
    field_u10, field_v10, nc_name="si10", output_file=f"{data_out}/si10.nc"
)
print("done.")

print("Converting surface geopotential to orography...", end="", flush=True)
# Remove time dimension and convert from geopotential to metres.
field_orog = preprocessing.squeeze_field(f"{data_dir}/era5_sfc_z.nc")
field_orog = preprocessing.multiply_field(field_orog, 1 / GRAVITY)
field_orog = preprocessing.set_netcdf_info(
    field_orog,
    nc_name="orog",
    properties={
        "standard_name": "surface_altitude",
        "long_name": "Surface Altitude",
        "units": "m",
    },
    output_file=f"{data_out}/orog.nc",
)
del field_orog
print("done.")


# ======== TSTORMS ========
# TSTORMS requires specific netcdf names, an unlimited time dimension, a change to the
# time units as it cannot handle negative times, and ascending latitude
def preprocess_tstorms_input(
    field: str | cf.Field, nc_name: str, output_file: str
) -> None:
    """Remove extra dims, convert time units, and write a TSTORMS input."""
    print(f"Preprocessing {nc_name}...", end="", flush=True)
    field = preprocessing.squeeze_field(field)
    field = preprocessing.set_time_units(field, "days since 1950-01-01")
    field = preprocessing.flip_axis(field, "Y")
    preprocessing.set_netcdf_info(
        field,
        nc_name=nc_name,
        output_file=output_file,
        coord_nc_names={"time": "time", "latitude": "lat", "longitude": "lon"},
        axis_unlimited="T",
    )
    print("done.")


preprocess_tstorms_input(field_msl, "slp", f"{data_out}/slp.nc")
del field_msl

preprocess_tstorms_input(field_u10, "u_ref", f"{data_out}/u_ref.nc")
del field_u10

preprocess_tstorms_input(field_v10, "v_ref", f"{data_out}/v_ref.nc")
del field_v10

preprocess_tstorms_input(f"{data_dir}/era5_vo.nc", "vort850", f"{data_out}/vort850.nc")

field_t = preprocessing.collapse_field(f"{data_dir}/era5_t.nc", "mean", axes="Z")
preprocess_tstorms_input(field_t, "tm", f"{data_out}/tm.nc")
del field_t
