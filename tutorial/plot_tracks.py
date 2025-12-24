"""
Script to plot tracks from TCTrack output.

Change the value of `TCTRACK_DATA` to use different input files.
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

TCTRACK_DATA = "tracks_tempest_extremes.nc"

# Open the NetCDF file
with netCDF4.Dataset(TCTRACK_DATA) as ncfile:
    # Read variables
    lat_var = ncfile.variables["lat"]
    lon_var = ncfile.variables["lon"]
    time_var = ncfile.variables["time"]
    intensity_var = ncfile.variables["wind_speed"]
    traj_var = ncfile.variables["trajectory"]

    lats = lat_var[:]
    lons = lon_var[:]
    intensity = intensity_var[:]
    traj_labels = traj_var[:]
    times = time_var[:]

    # Convert times to datetime objects
    missing_time = getattr(time_var, "missing_value", np.nan)
    times = np.ma.masked_where(times == missing_time, times)
    time_units = time_var.units
    time_calendar = time_var.calendar
    times_dt = netCDF4.num2date(times, units=time_units, calendar=time_calendar)

    # Get intensity metadata for labels
    intensity_name = intensity_var.long_name
    intensity_units = getattr(intensity_var, "units", "")
    min_intensity = np.nanmin(intensity)
    max_intensity = np.nanmax(intensity)

    plt.figure(figsize=(10, 3))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True

    # Plot each trajectory
    for i in traj_labels:
        times_i = times_dt[i, :].compressed()
        label = (
            f"{times_i[0].strftime('%Y-%m-%d %H:%M')} to "
            f"{times_i[-1].strftime('%Y-%m-%d %H:%M')}"
        )
        pl = ax.plot(
            lons[i], lats[i], "--", transform=ccrs.PlateCarree(), label=f"{label}"
        )
        sc = ax.scatter(
            lons[i],
            lats[i],
            c=intensity[i],
            cmap="viridis",
            s=40,
            vmin=min_intensity,
            vmax=max_intensity,
            transform=ccrs.PlateCarree(),
        )

    plt.colorbar(sc, label=f"{intensity_name} ({intensity_units})")
    plt.title(f"All Trajectories Colored by {intensity_name}")
    plt.legend()
    plt.tight_layout()
    plt.savefig("my_tracks.png")
