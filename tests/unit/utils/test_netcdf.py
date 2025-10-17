"""Unit tests for netcdf utility functions in the TCTrack utils Python package."""

import pytest
from netCDF4 import Dataset

from tctrack.utils import lat_lon_sizes


class TestNetCDF:
    """Tests for the netcdf utility functions."""

    def create_netcdf_file(self, path, dims):
        """Create a dummy NetCDF file with given dimensions."""
        with Dataset(path, "w") as ds:
            for name, size in dims.items():
                ds.createDimension(name, size)

    @pytest.mark.parametrize(
        "dims, lat_lon_names, expected, raises",
        [
            pytest.param({"lon": 10, "lat": 20}, {}, (20, 10), None, id="defaults"),
            pytest.param({"lon": 0, "lat": 0}, {}, (0, 0), None, id="zero dimensions"),
            pytest.param(
                {"longitude": 7, "latitude": 9},
                {"lon": "longitude", "lat": "latitude"},
                (9, 7),
                None,
                id="custom names",
            ),
            pytest.param({"lon": 10}, {}, (20, 10), KeyError, id="missing dimension"),
        ],
    )
    def test_lat_lon_sizes_valid(self, tmp_path, dims, lat_lon_names, expected, raises):
        """Test lat_lon_sizes returns the expected values / the expected error."""
        ncfile = tmp_path / "test.nc"
        self.create_netcdf_file(ncfile, dims)
        if raises:
            with pytest.raises(raises):
                lat_lon_sizes(str(ncfile), **lat_lon_names)
        else:
            assert lat_lon_sizes(str(ncfile), **lat_lon_names) == expected
