"""Module providing utility functions for working with netcdf files."""

from netCDF4 import Dataset


def lat_lon_sizes(filename: str, lat: str = "lat", lon: str = "lon") -> tuple[int, int]:
    """Get the sizes of the longitude and latitude for a netcdf file.

    Parameters
    ----------
    filename : str
        The path to the netcdf file to read from.
    lat : str, default="lat"
        The netcdf variable name of the latitude.
    lon : str, default="lon"
        The netcdf variable name of the longitude.

    Returns
    -------
    int
        The size of the latitude dimension.
    int
        The size of the longitude dimension.
    """
    dataset = Dataset(filename, "r")
    lon_size = len(dataset.dimensions[lon])
    lat_size = len(dataset.dimensions[lat])
    dataset.close()
    return lat_size, lon_size
