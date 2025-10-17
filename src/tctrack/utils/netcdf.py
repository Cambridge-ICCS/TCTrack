"""Module providing utility functions for working with netcdf files."""

from netCDF4 import Dataset


def lat_lon_sizes(filename: str, lat: str = "lat", lon: str = "lon") -> tuple[int, int]:
    """Get the sizes of the longitude and latitude for a netcdf file."""
    dataset = Dataset(filename, "r")
    lon_size = len(dataset.dimensions[lon])
    lat_size = len(dataset.dimensions[lat])
    dataset.close()
    return lat_size, lon_size
