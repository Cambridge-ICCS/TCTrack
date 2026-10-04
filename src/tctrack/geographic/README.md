# Geographic

Classification of TCTrack trajectory points against a bundled world map.
	
Each lat/lon point of a dataset from `read_netcdf` is matched to a cell of a quarter-degree world map grid to find whether it is over land, its distance to the nearest coastline in km, and the ocean it is in. The map data is generated from the Natural Earth 1:50m physical vectors, so classification is fast, has an accuracy on the order of the grid resolution (about 25 km), and requires no runtime downloads.

## Usage

Presuming TCTrack has been installed via `pip install`:

```python
from tctrack.build_db.build_db import read_netcdf
from tctrack.geographic import geographic

netcdf_data = read_netcdf("tracks.nc")
result = geographic.classify_tracks(netcdf_data)
```

`result` is a dictionary with one entry per trajectory, keyed by trajectory index and therefore aligned with the rows of the input arrays. Each entry holds arrays over the observations of that trajectory:

| Key                      | Description                                     |
|--------------------------|-------------------------------------------------|
| `landfall`               | Boolean array, True where a point is over land. |
| `distance_to_coast_km`   | Float array of distances in km to the nearest coastline, NaN where the input point is NaN. Distance inland from the coast for land points, distance to the nearest land for sea points. |
| `ocean_basin`            | Array of ocean names, None where the input point is NaN. Name of the containing cell for sea points, of the nearest sea cell for land points. |

Single arrays of points can be classified directly:

```python
landfall, distance, basin = geographic.classify_points([51.5], [-0.1])
```

Longitudes may be in either the [-180, 180) or [0, 360) convention. Lakes and rivers count as land. Points are classified into the seven oceans of Natural Earth; marginal seas such as the Gulf of Mexico take the ocean they connect to.

The arrays span the full observation dimension of the input, with NaN padded observations left unclassified.

## World map data

The grid data is bundled as [`world_map.npz`](world_map.npz). It can be regenerated with [generate_map.py](generate_map.py), a development script whose extra requirements are available via the `mapdata` optional dependency:

```bash
pip install .[mapdata]
python -m tctrack.geographic.generate_map
```
