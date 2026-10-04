"""Module providing classification of trajectory points against a world map.

Points are matched to the cells of a bundled quarter-degree world map grid
generated from the Natural Earth 1:50m physical vector data (see
`generate_map.py`). This makes classification fast, with an accuracy on the
order of the grid resolution, and free of any runtime map data downloads.
"""

from functools import cache
from pathlib import Path

import numpy as np

DATA_PATH = Path(__file__).parent / "world_map.npz"

#: Maximum absolute latitude in degrees.
MAX_ABS_LATITUDE = 90.0

#: Number of dimensions of the latitude and longitude arrays of a dataset.
EXPECTED_NDIM = 2


@cache
def _world_map() -> dict:
    """
    Load the bundled world map grid data, cached after the first call.

    Returns
    -------
    Dictionary with the grid arrays and their resolution in degrees:

    - ``land``: boolean grid, True where a cell is over land
    - ``basin``: grid of ocean name indices
    - ``distance``: grid of distances in km to the nearest cell of the
      opposite surface (coast distance for land cells, land distance for sea
      cells)
    - ``names``: ocean names, indexed by the ``basin`` grid
    - ``resolution``: grid resolution in degrees
    """
    with np.load(DATA_PATH) as data:
        world_map = {
            "land": data["land"],
            "basin": data["basin"],
            "distance": data["distance"].astype(np.float64),
            "names": data["names"],
        }

    world_map["resolution"] = 360.0 / world_map["land"].shape[1]

    return world_map


def _to_nan_array(values: np.typing.ArrayLike) -> np.ndarray:
    """
    Return values as a float array with masked entries set to NaN.

    Parameters
    ----------
    values
        Array of numbers, possibly a masked array from netCDF4.

    Returns
    -------
    Float array with the same shape, masked entries replaced by NaN.
    """
    array: np.ma.MaskedArray = np.ma.asanyarray(values)

    if not np.ma.is_masked(array):
        return np.asarray(array, dtype=np.float64)

    return np.asarray(array.filled(np.nan), dtype=np.float64)


def _grid_indices(
    lat: np.ndarray, lon: np.ndarray, world_map: dict
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the (row, column) of the grid cell containing each point.

    Parameters
    ----------
    lat, lon
        Latitudes and longitudes of the points in degrees. All points must be
        finite with latitudes within [-90, 90].
    world_map
        World map dictionary from `_world_map`.

    Returns
    -------
    Row and column index arrays with one entry per point.
    """
    n_rows, n_cols = world_map["land"].shape

    # Wrap longitudes to [-180, 180) so either convention can be used
    wrapped = np.mod(lon + 180.0, 360.0) - 180.0
    cols = np.floor((wrapped + 180.0) / world_map["resolution"]).astype(np.int64)
    rows = np.floor((lat + 90.0) / world_map["resolution"]).astype(np.int64)

    # Points at the edge of the grid fall in the last row or column
    np.clip(rows, 0, n_rows - 1, out=rows)
    np.clip(cols, 0, n_cols - 1, out=cols)

    return rows, cols


def classify_points(
    latitude: np.typing.ArrayLike,
    longitude: np.typing.ArrayLike,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Classify points as land or sea, with coast distance and ocean basin.

    Each point is matched to the cell of a quarter-degree world map grid
    containing it, so results are accurate to the order of the grid
    resolution (about 25 km).

    Parameters
    ----------
    latitude
        Latitude values in degrees north. NaN entries are skipped.
    longitude
        Longitude values in degrees east, in either the [-180, 180) or
        [0, 360) convention. NaN entries are skipped.

    Returns
    -------
    landfall
        Boolean array, True where a point is over land.
    distance_to_coast_km
        Float array of distances in km to the nearest coastline: distance
        inland from the coast for land points, and distance to the nearest
        land for sea points. NaN entries are NaN.
    ocean_basin
        Array of ocean names, of the containing cell for sea points and of
        the nearest sea cell for land points. NaN entries are None.

    Raises
    ------
    ValueError
        If latitude and longitude have different shapes, or if any finite
        latitude is outside the range [-90, 90].
    """
    lat = _to_nan_array(latitude)
    lon = _to_nan_array(longitude)

    if lat.shape != lon.shape:
        msg = f"latitude shape {lat.shape} does not match longitude shape {lon.shape}"
        raise ValueError(msg)

    valid = np.isfinite(lat) & np.isfinite(lon)

    if (valid & (np.abs(lat) > MAX_ABS_LATITUDE)).any():
        msg = "latitude values outside the range [-90, 90] degrees"
        raise ValueError(msg)

    landfall = np.zeros(lat.shape, dtype=bool)
    distance = np.full(lat.shape, np.nan)
    basin = np.full(lat.shape, None, dtype=object)

    if not valid.any():
        return landfall, distance, basin

    world_map = _world_map()
    rows, cols = _grid_indices(lat[valid], lon[valid], world_map)

    landfall[valid] = world_map["land"][rows, cols]
    distance[valid] = world_map["distance"][rows, cols]
    basin[valid] = world_map["names"][world_map["basin"][rows, cols]]

    return landfall, distance, basin


def classify_tracks(netcdf_data: dict) -> dict[int, dict]:
    """
    Classify the points of every trajectory in a read_netcdf dataset.

    Each point is tested against a bundled quarter-degree world map grid to
    find whether it is over land, its distance to the nearest coastline and
    the ocean it is in.

    Parameters
    ----------
    netcdf_data
        Dictionary from `tctrack.build_db.build_db.read_netcdf`.

    Returns
    -------
    Dictionary with one entry per trajectory, keyed by trajectory index and
    therefore aligned with the rows of the input arrays. Each entry holds
    arrays over the observations of that trajectory:

    - ``landfall``: boolean array, True where a point is over land
    - ``distance_to_coast_km``: float array of distances in km to the nearest
      coastline, NaN where the input point is NaN
    - ``ocean_basin``: array of ocean names, None where the input point is
      NaN

    The arrays span the full observation dimension of the input, with NaN
    padded observations left unclassified.

    Raises
    ------
    ValueError
        If latitude or longitude data is missing, is not two-dimensional or
        the two have different shapes.
    """
    for key in ("latitude", "longitude"):
        if key not in netcdf_data:
            msg = f"missing {key!r} key in netcdf data"
            raise ValueError(msg)

    lat = _to_nan_array(netcdf_data["latitude"])
    lon = _to_nan_array(netcdf_data["longitude"])

    if lat.ndim != EXPECTED_NDIM or lon.ndim != EXPECTED_NDIM:
        msg = "latitude and longitude must be two-dimensional arrays"
        raise ValueError(msg)

    landfall, distance, basin = classify_points(lat, lon)

    return {
        traj_idx: {
            "landfall": landfall[traj_idx],
            "distance_to_coast_km": distance[traj_idx],
            "ocean_basin": basin[traj_idx],
        }
        for traj_idx in range(lat.shape[0])
    }
