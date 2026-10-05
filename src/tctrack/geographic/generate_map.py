"""Script to regenerate the world map data file used by the geographic module.

Reads land and ocean polygons from the Natural Earth 1:50m physical vector
data (https://www.naturalearthdata.com/) and rasterises them onto a regular
quarter-degree grid. For every grid cell the following values are precomputed:

- land: whether the cell is over land.
- basin: the ocean the cell belongs to, one of the seven oceans of
  OCEAN_NAMES, with land cells assigned the ocean of their nearest sea cell.
- distance: the great-circle distance in km to the nearest cell of the
  opposite surface (distance to the coast for land cells, distance to land
  for sea cells).

The result is written to a compressed NumPy data file (world_map.npz)
and bundled with the package.

This is a development script: it requires matplotlib and scipy, which are not
runtime dependencies of TCTrack. They are available via the `mapdata` optional
dependency: `pip install .[mapdata]`.

Execute: python -m tctrack.geographic.generate_map
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath
from scipy.ndimage import binary_dilation, distance_transform_edt

LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_50m_land.geojson"
)
MARINE_POLYS_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_50m_geography_marine_polys.geojson"
)

OUTPUT_PATH = Path(__file__).parent / "world_map.npz"

#: Resolution of the generated grid in degrees.
RESOLUTION = 0.25

#: Mean Earth radius in km used for great-circle distances.
EARTH_RADIUS_KM = 6371.0

#: Land fraction expected in the generated map, as a sanity range.
EXPECTED_LAND_FRACTION = (0.20, 0.40)

#: Maximum number of basin names storable in the uint8 basin grid.
MAX_BASINS = 256

#: Spot checks used to sanity check the generated map: (lat, lon, on land).
SPOT_CHECKS = [
    (39.7, -105.0, True),  # Denver
    (51.5, -0.1, True),  # London
    (30.0, -40.0, False),  # mid-Atlantic
    (0.0, -140.0, False),  # mid-Pacific
]

#: The seven oceans the map classifies cells into (from the named sea
#: polygons dataset). All seas take the ocean they connect to.
OCEAN_NAMES = (
    "North Atlantic Ocean",
    "South Atlantic Ocean",
    "North Pacific Ocean",
    "South Pacific Ocean",
    "Indian Ocean",
    "Arctic Ocean",
    "Southern Ocean",
)


def read_geojson(source: str) -> dict:
    """
    Read a GeoJSON dataset from a local file path or a URL.

    Parameters
    ----------
    source
        Path of a GeoJSON file, or an http(s) URL to download it from.

    Returns
    -------
    GeoJSON dataset as a dictionary.
    """
    if not source.startswith(("http://", "https://")):
        with open(source, encoding="utf-8") as file:
            return json.load(file)

    with urllib.request.urlopen(source) as response:  # noqa: S310 - scheme checked above
        return json.loads(response.read().decode("utf-8"))


def outer_rings(feature: dict) -> list[np.ndarray]:
    """
    Return the outer rings of the (multi-)polygon geometry of a feature.

    Interior rings are skipped so that lakes and rivers count as land.

    Parameters
    ----------
    feature
        GeoJSON feature with a Polygon or MultiPolygon geometry.

    Returns
    -------
    List of outer rings as (N, 2) arrays of (longitude, latitude) vertices.
    """
    geometry = feature["geometry"]
    if geometry is None:
        return []

    if geometry["type"] == "Polygon":
        polygons = [geometry["coordinates"]]
    else:
        polygons = geometry["coordinates"]

    return [np.asarray(polygon[0]) for polygon in polygons]


def rasterise(
    rings: list[np.ndarray],
    lon_centres: np.ndarray,
    lat_centres: np.ndarray,
) -> np.ndarray:
    """
    Return a boolean grid of cells whose centre lies inside any of the rings.

    Parameters
    ----------
    rings
        Polygon outer rings as (N, 2) arrays of (longitude, latitude).
    lon_centres
        Longitudes of the grid cell centres.
    lat_centres
        Latitudes of the grid cell centres.

    Returns
    -------
    Boolean array of shape (n_lat, n_lon).
    """
    mask = np.zeros((lat_centres.size, lon_centres.size), dtype=bool)

    for ring in rings:
        lon_min, lat_min = ring.min(axis=0)
        lon_max, lat_max = ring.max(axis=0)

        # Only cell centres within the bounding box of the ring can be inside
        cols = np.nonzero((lon_centres >= lon_min) & (lon_centres <= lon_max))[0]
        rows = np.nonzero((lat_centres >= lat_min) & (lat_centres <= lat_max))[0]
        if cols.size == 0 or rows.size == 0:
            continue

        lon_mesh, lat_mesh = np.meshgrid(lon_centres[cols], lat_centres[rows])
        points = np.column_stack([lon_mesh.ravel(), lat_mesh.ravel()])
        inside = MplPath(ring).contains_points(points).reshape(rows.size, cols.size)
        mask[np.ix_(rows, cols)] |= inside

    return mask


def haversine_km(
    lat1: np.ndarray,
    lon1: np.ndarray,
    lat2: np.ndarray,
    lon2: np.ndarray,
) -> np.ndarray:
    """
    Return the great-circle distance in km between pairs of points.

    Parameters
    ----------
    lat1, lon1
        Latitudes and longitudes of the first points in degrees.
    lat2, lon2
        Latitudes and longitudes of the second points in degrees.

    Returns
    -------
    Distances in km, with one entry per point pair.
    """
    lat1, lon1, lat2, lon2 = (np.radians(array) for array in (lat1, lon1, lat2, lon2))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    h = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(h))


def build(land_source: str, marine_polys_source: str, output_path: Path) -> None:  # noqa: PLR0912, PLR0915
    """
    Rasterise the map data and write the world map data file.

    Only the ocean polygons of the named seas dataset are used, so that cells
    are classified into the seven oceans of OCEAN_NAMES. Sea cells covered by
    no ocean polygon take the ocean nearest through connected sea cells, or
    for waters cut off by straits narrower than a grid cell, the nearest
    ocean by straight-line distance.

    Parameters
    ----------
    land_source
        Path or URL of the land polygons GeoJSON dataset.
    marine_polys_source
        Path or URL of the named sea polygons GeoJSON dataset, of which only
        the ocean polygons are used.
    output_path
        Path of the compressed NumPy data file to write.
    """
    lat_centres = -90.0 + RESOLUTION * (np.arange(180 // RESOLUTION) + 0.5)
    lon_centres = -180.0 + RESOLUTION * (np.arange(360 // RESOLUTION) + 0.5)
    lat_grid, lon_grid = np.meshgrid(lat_centres, lon_centres, indexing="ij")

    # Land cells: centres inside any land polygon (outer rings only)
    land_rings = [
        ring
        for feature in read_geojson(land_source)["features"]
        for ring in outer_rings(feature)
    ]
    land = rasterise(land_rings, lon_centres, lat_centres)

    # Basin labels: rasterise each of the seven ocean polygons separately
    names: list[str] = []
    entries: list[tuple[int, np.ndarray]] = []
    for feature in read_geojson(marine_polys_source)["features"]:
        name = feature["properties"]["name"]
        if name.isupper():  # e.g. "INDIAN OCEAN" in the source data
            name = name.title()

        if name not in OCEAN_NAMES:
            continue

        mask = rasterise(outer_rings(feature), lon_centres, lat_centres)
        if not mask.any():
            continue

        if name not in names:
            names.append(name)
        entries.append((names.index(name), mask))

    # Paint labels from the largest polygon to the smallest, so that smaller
    # polygons take priority where they overlap
    entries.sort(key=lambda entry: -int(entry[1].sum()))
    basin = np.full(land.shape, -1, dtype=np.int16)
    for basin_id, mask in entries:
        basin[mask] = basin_id

    sea = ~land

    # Fill sea cells covered by no ocean polygon from the ocean nearest
    # through connected sea cells, so that land does not cut a sea off from
    # the ocean it flows into: each ocean grows by one cell per round,
    # claiming whichever cells it reaches first
    unlabelled = sea & (basin < 0)
    labelled_sea = sea & (basin >= 0)
    disconnected = 0
    if unlabelled.any():
        if not labelled_sea.any():
            msg = "no sea cells received a basin label"
            raise RuntimeError(msg)

        # Neighbours are 8-connected, so straits crossing the grid diagonally
        # still connect their waters
        structure = np.ones((3, 3), dtype=bool)
        unreached = unlabelled.copy()
        while unreached.any():
            progressed = False

            for basin_id in range(len(names)):
                newly = (
                    binary_dilation(sea & (basin == basin_id), structure=structure)
                    & unreached
                )
                if newly.any():
                    basin[newly] = basin_id
                    unreached &= ~newly
                    progressed = True

            if not progressed:
                break

        # Waters cut off from every ocean by straits narrower than a grid
        # cell (e.g. the Bosphorus) take the basin of the nearest reached
        # sea cell by straight-line distance
        if unreached.any():
            disconnected = int(unreached.sum())
            _, (rows, cols) = distance_transform_edt(
                ~(sea & (basin >= 0)), return_indices=True
            )
            basin[unreached] = basin[rows[unreached], cols[unreached]]

    # Index of the nearest sea cell for every cell, and of the nearest land
    # cell for every cell, from the grid distances to the opposite surface
    _, (sea_rows, sea_cols) = distance_transform_edt(land, return_indices=True)
    _, (land_rows, land_cols) = distance_transform_edt(sea, return_indices=True)

    # Land cells take the basin of their nearest sea cell
    basin[land] = basin[sea_rows[land], sea_cols[land]]

    # Distance to the nearest cell of the opposite surface for every cell
    distance = np.empty(land.shape, dtype=np.float64)
    distance[land] = haversine_km(
        lat_grid[land],
        lon_grid[land],
        lat_grid[sea_rows[land], sea_cols[land]],
        lon_grid[sea_rows[land], sea_cols[land]],
    )
    distance[sea] = haversine_km(
        lat_grid[sea],
        lon_grid[sea],
        lat_grid[land_rows[sea], land_cols[sea]],
        lon_grid[land_rows[sea], land_cols[sea]],
    )

    _sanity_check(land, basin, names)

    np.savez_compressed(
        output_path,
        land=land,
        basin=basin.astype(np.uint8),
        distance=np.rint(distance).astype(np.uint16),
        names=np.array(names),
    )

    print(
        f"World map written to {output_path}: "
        f"{land.mean():.1%} land, {len(names)} oceans, "
        f"{int(unlabelled.sum())} sea cells filled from connected oceans, "
        f"{disconnected} only by straight-line distance"
    )


def _sanity_check(land: np.ndarray, basin: np.ndarray, names: list[str]) -> None:
    """Raise RuntimeError if the generated map fails basic checks."""
    if not EXPECTED_LAND_FRACTION[0] < land.mean() < EXPECTED_LAND_FRACTION[1]:
        msg = f"land fraction {land.mean():.2f} outside {EXPECTED_LAND_FRACTION}"
        raise RuntimeError(msg)

    if (basin < 0).any():
        msg = "some cells were not assigned an ocean basin"
        raise RuntimeError(msg)

    if len(names) > MAX_BASINS:
        msg = f"{len(names)} basins exceed the {MAX_BASINS} storable in the grid"
        raise RuntimeError(msg)

    for lat, lon, expected_land in SPOT_CHECKS:
        row = int((lat + 90.0) / RESOLUTION)
        col = int((lon + 180.0) / RESOLUTION)
        if land[row, col] != expected_land:
            msg = f"spot check ({lat}, {lon}) expected land={expected_land}"
            raise RuntimeError(msg)

    for expected in OCEAN_NAMES:
        if expected not in names:
            msg = f"expected ocean {expected!r} not found in {names}"
            raise RuntimeError(msg)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Regenerate the world map data file for the geographic module.",
    )
    parser.add_argument(
        "--land",
        default=LAND_URL,
        help="Land polygons GeoJSON file path or URL",
    )
    parser.add_argument(
        "--marine-polys",
        default=MARINE_POLYS_URL,
        help="Named sea polygons GeoJSON file path or URL",
    )
    parser.add_argument(
        "--output",
        default=str(OUTPUT_PATH),
        help="Path of the data file to write",
    )
    args = parser.parse_args(argv)

    build(args.land, args.marine_polys, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
