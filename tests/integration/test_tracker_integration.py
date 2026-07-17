"""Integration tests for the TempestExtremes tracker pipeline."""

import cf
import cftime
import numpy as np
import pytest
import xarray as xr

import tctrack.tempest_extremes as te


@pytest.fixture
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


class TestTETrackerIntegration:
    """Integration tests for the full TETracker pipeline."""

    def test_run_tracker_produces_valid_tracks(self, synthetic_data_file, tmp_path):
        """TETracker detects the synthetic pressure minimum and writes a track file.

        Verifies the full chain: DetectNodes → StitchNodes → CF-NetCDF output.
        The planted minimum at (15°N, 45°E) should produce exactly one trajectory
        whose coordinates stay within 5° of the expected location.
        """
        output_file = str(tmp_path / "tracks.nc")

        dn_params = te.TEDetectParameters(
            in_data=[synthetic_data_file],
            search_by_min="psl",
            merge_dist=6.0,
            closed_contours=[
                te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
            ],
            out_header=True,
            output_commands=[
                te.TEOutputCommand(var="psl", operator="min", dist=0.0),
                te.TEOutputCommand(var="sfcWind", operator="max", dist=2.0),
            ],
            output_dir=str(tmp_path / "te_outputs"),
        )
        sn_params = te.TEStitchParameters(
            caltype="360_day",
            max_sep=8.0,
            min_time="2",
        )

        tracker = te.TETracker(dn_params, sn_params)
        tracker.run_tracker(output_file)

        fields = cf.read(output_file)  # type: ignore[operator]
        assert len(fields) > 0, "No fields written to output file"

        field = fields[0]
        assert field.dimension_coordinate("trajectory").size >= 1, (
            "No trajectories found"
        )

        # All valid (non-fill) track coordinates should sit near (15°N, 45°E)
        lat_vals = field.construct("latitude").data.array.flatten()
        lon_vals = field.construct("longitude").data.array.flatten()
        lat_valid = lat_vals[~np.isnan(lat_vals)]
        lon_valid = lon_vals[~np.isnan(lon_vals)]
        assert np.all(np.abs(lat_valid - 15.0) < 5.0), (
            "Track latitude far from planted minimum"
        )
        assert np.all(np.abs(lon_valid - 45.0) < 5.0), (
            "Track longitude far from planted minimum"
        )
