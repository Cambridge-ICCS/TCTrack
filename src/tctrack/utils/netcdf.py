"""Module providing utility functions for working with netcdf files."""

from netCDF4 import Dataset


def lon_lat_sizes(filename: str, lon: str = "lon", lat: str = "lat") -> tuple[int, int]:
    """Return the sizes of the longitude and latitude for a netcdf file."""
    dataset = Dataset(filename, "r")
    lon_size = len(dataset.dimensions[lon])
    lat_size = len(dataset.dimensions[lat])
    dataset.close()
    return lon_size, lat_size
