"""Integration tests for the TempestExtremes tracker pipeline."""

import cftime
import numpy as np
import xarray as xr


def synthetic_data_file(tmp_path):
    """NetCDF file with a Gaussian psl minimum planted at (15°N, 45°E) over 5 days.

    The depression has amplitude -1500 Pa and sigma=3°, so at 5.5° from the
    centre psl rises by ~1220 Pa — well above the 200 Pa closed-contour threshold
    used in the detection parameters.
    """
    nlat, nlon, ntimes = 90, 180, 5
    lat = np.linspace(-89, 89, nlat)
    lon = np.linspace(0, 358, nlon)
    times = [
        cftime.datetime(1950, 8, i + 1, 0, calendar="360_day") for i in range(ntimes)
    ]

    lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")
    dist2 = (lat2d - 15.0) ** 2 + (lon2d - 45.0) ** 2
    psl_anomaly = -1500.0 * np.exp(-dist2 / 18.0)  # sigma = 3 degrees
    psl = np.broadcast_to(101325.0 + psl_anomaly, (ntimes, nlat, nlon)).copy()

    ds = xr.Dataset(
        {
            "psl": (["time", "lat", "lon"], psl.astype(np.float32), {"units": "Pa"}),
            "sfcWind": (
                ["time", "lat", "lon"],
                np.full((ntimes, nlat, nlon), 5.0, dtype=np.float32),
                {"units": "m s-1"},
            ),
        },
        coords={"time": times, "lat": lat, "lon": lon},
    )
    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})
    path = tmp_path / "synthetic.nc"
    ds.to_netcdf(
        path,
        encoding={"time": {"calendar": "360_day", "units": "days since 1950-01-01"}},
    )
    return str(path)
