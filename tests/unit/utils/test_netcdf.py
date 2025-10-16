"""Unit tests for netcdf.py in the TCTrack utils Python package."""

import pytest
from netCDF4 import Dataset

from tctrack.utils import netcdf


def create_netcdf_file(path, dims):
    """Create a dummy NetCDF file with given dimensions."""
    with Dataset(path, "w") as ds:
        for name, size in dims.items():
            ds.createDimension(name, size)


@pytest.mark.parametrize(
    "dims, lon_lat_kwargs, expected, raises",
    [
        pytest.param({"lon": 10, "lat": 20}, {}, (10, 20), None, id="standard names"),
        pytest.param({"lon": 0, "lat": 0}, {}, (0, 0), None, id="empty dimensions"),
        pytest.param(
            {"longitude": 7, "latitude": 9},
            {"lon": "longitude", "lat": "latitude"},
            (7, 9),
            None,
            id="custom_names",
        ),
        pytest.param({"lon": 10}, {}, (10, 20), KeyError, id="missing dimension"),
    ],
)
def test_lon_lat_sizes_valid(tmp_path, dims, lon_lat_kwargs, expected, raises):
    """Test lon_lat_sizes returns the expected values / raises the expected error."""
    ncfile = tmp_path / "test.nc"
    create_netcdf_file(ncfile, dims)
    if raises:
        with pytest.raises(raises):
            netcdf.lon_lat_sizes(str(ncfile), **lon_lat_kwargs)
    else:
        assert netcdf.lon_lat_sizes(str(ncfile), **lon_lat_kwargs) == expected
