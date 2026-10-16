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
| `gate_id`                | Object array of gate codes for gated landfalls, None otherwise (including when gates are disabled). |

Single arrays of points can be classified directly:

```python
landfall, distance, basin = geographic.classify_points([51.5], [-0.1])
```

Longitudes may be in either the [-180, 180) or [0, 360) convention. Lakes and rivers count as land. Points are classified into the seven oceans of Natural Earth; marginal seas such as the Gulf of Mexico take the ocean they connect to.

The arrays span the full observation dimension of the input, with NaN padded observations left unclassified.

## Optional landfall gates

Landfall gates add named rectangular zones to the existing classification. Supply an explicit CSV file with the installed `build-db` command (the Python component is named `build_db`). From the repository root, after installing TCTrack:

```shell
build-db --output east_coast.db --gates src/tctrack/geographic/east_coast_gates.csv tests/unit/build_db/6hr_track_mm_2013_natl.nc
```

Use a fresh output database for the gate-enabled schema. The [example CSV](src/tctrack/geographic/east_coast_gates.csv) is a user-supplied configuration, not an automatically loaded package resource or built-in preset. No additional dependency or map download is needed for gates.

For programmatic classification, load definitions once and pass them to `classify_tracks`; `classify_points` retains its existing three-value return:

```python
from tctrack.geographic import classify_tracks, load_gates

netcdf_data = read_netcdf("tracks.nc")
gates = load_gates("src/tctrack/geographic/east_coast_gates.csv")
result = classify_tracks(netcdf_data, gates=gates)
```

### Matching and summary semantics

- An observation receives a gate code only if its existing `landfall` result is true and its latitude/longitude is inside that gate. Sea points, land outside every gate, and invalid, masked, or padded observations have no gate.
- Gates do not filter observations, create landfalls, or change landfall, ocean, coast-distance, or stored coordinate values. The existing detector calls any point over land a landfall, including subsequent inland observations; it does not require a sea-to-land transition. Gates do not add coastline crossings or interpolation.
- Gates do not check direction of travel or whether a storm approached from the Atlantic. A land observation inside a gate can match regardless of the storm's previous position or approach direction.
- Rectangles use WGS84 decimal degrees, north-positive latitude and east-positive longitude. Both axes use inclusive minima and exclusive maxima: `min_lat <= latitude < max_lat` and `min_lon <= longitude < max_lon`. Observation longitudes are wrapped to `[-180, 180)` for matching only, so equivalent negative and 0–360 longitudes match the same gate.
- Each observation can match at most one gate. Overlapping rectangle interiors are rejected, not resolved by CSV row order. Shared edges are allowed and the half-open bounds make their assignment unambiguous.
- A trajectory's `gate_id` is the gate of its first observation with a non-NULL gate ID, in persisted sequence order. An earlier landfall outside all gates does not prevent a later gated landfall from supplying the summary. Later observations retain their own gate IDs even if they differ from the summary. No gated landfall means NULL. Unlike the gate summary, the trajectory ocean is taken from its first observation (genesis); trajectory landfall remains true if any selected observation is over land.
- Without `--gates`, no definitions are loaded from a file and newly imported observations and trajectories receive NULL gate IDs, even if the database already contains gate definitions. NULL means unassigned, not necessarily at sea; it is not an empty string or a replacement landfall code.

### CSV contract

The CSV requires these six columns and at least one data row:

```csv
code,name,min_lat,max_lat,min_lon,max_lon
FLA,Atlantic Florida,24.4,30.7,-81.8,-79.5
```

Codes and names have surrounding whitespace trimmed and must be nonempty. Codes must be unique, case-sensitive explicit identifiers; they are not constructed from names. The example uses three-character, project-defined codes without a country prefix; these are not official identifiers. Names need not be unique.

Bounds must be finite numbers satisfying:

- `-90 <= min_lat < max_lat <= 90`
- `-180 <= min_lon < max_lon <= 180`

Wrapped or antimeridian-crossing intervals are not supported. Missing or extra row fields, malformed numbers (including NaN and infinity), invalid bounds, and overlapping interiors are errors. Overlap means positive intersection on both axes, and errors identify both conflicting codes. Missing/unreadable files fail clearly; malformed content is reported with actionable filename/row information. Definitions are validated once before opening or writing the output database, including for programmatic `build` calls, so invalid configuration does not create or mutate it.

The complete example defines six approximate Atlantic coastal sectors:

| Code | Sector |
|------|--------|
| `FLA` | Atlantic Florida |
| `GSC` | Georgia-South Carolina coast |
| `NCA` | North Carolina coast |
| `MAT` | Mid-Atlantic coast |
| `NYL` | New York Bight-Long Island |
| `NEN` | New England coast |

These are illustrative storm-tracking zones, **not state rectangles, official NHC regions, or a validated coastal mask**. Boundaries are approximate and use shared latitude edges to avoid overlapping interiors. The New York Bight-Long Island sector includes nearby coastal areas, not just New York. Florida targets the Atlantic coast, not the whole state or Keys/Gulf coast; a coarse rectangle cannot exclude every Gulf-side location. Rectangles include water and inland areas; the existing landfall result remains authoritative. Scientific use should include a coastline/track review.

### Database and import constraints

CSV `code` maps to the text primary key `gates.id`; the lookup table stores the name and four numeric bounds. Both `observations.gate_id` and `trajectories.gate_id` are nullable text foreign keys to `gates(id)`, with `obs_gate_idx` and `traj_gate_idx` indexes. There is no numeric surrogate ID or association table.

`points` exposes the observation gate ID and `tracks` the trajectory summary. The corresponding layered views are `points_layer_collection`, `points_layer_file`, `tracks_layer_collection`, and `tracks_layer_file`, using `layer_collection` or `layer_file` respectively. All expose `gate_id`.

Gate definitions are stored before referencing observations and trajectories. When reopening a gate-enabled database, identical definitions are reused, comparing codes and values rather than CSV row order. An empty gate table may be populated, but a configuration differing from an existing nonempty gate set fails before importing or changing definitions; old gate IDs are never silently reinterpreted.

Already imported file paths remain skipped. Supplying gates later does not reclassify skipped files or backfill earlier NULL assignments. Only fresh-database DDL is supported for this feature: there are no migrations, automatic schema upgrades, or backfills for old-schema databases.

## World map data

The grid data is bundled as [`world_map.npz`](world_map.npz). It can be regenerated with [generate_map.py](generate_map.py), a development script whose extra requirements are available via the `mapdata` optional dependency:

```bash
pip install .[mapdata]
python -m tctrack.geographic.generate_map
```
