# Build Database

A CLI tool to import TCTrack NetCDF files into SQLite databases.

The resulting database forms the foundation for a dashboard. However, it is also valuable as a organisational tool for TCTrack output files. The data within can be queried and manipulated easily using local SQLite viewers such as [DB Browser](https://sqlitebrowser.org/).

Imported files can be grouped together into collections, e.g. geographic areas or a particular research hypothesis. Collections are defined using the `--collection` option. Collections are optional.


## Usage

Presuming TCTrack has been installed via `pip install`:

```bash
build-db [--help] --output OUTPUT [--collection COLLECTION] [--hadsst PATH] [--roni PATH] files...
```

To use the tool without installation, bypassing the need for all other TCTrack dependencies:

```bash
cd /src
python -m tctrack.build_db.build_db --output OUTPUT files...
```

### Arguments

| Flag                 | Required | Description                                                                                              |
|----------------------|----------|----------------------------------------------------------------------------------------------------------|
| `--output`, `-o`     | yes      | SQLite database path and name. Created if it doesn't exist, appended to if it does.                      |
| `--collection`, `-c` | no       | Collection name. Creates or reuses an existing collection. Will use a default collection if unspecified. |
| `--hadsst PATH`      | no       | HadSST median monthly NetCDF file used to add sea surface temperature anomaly values to observations. |
| `--roni PATH`        | no       | Monthly NetCDF file used to add matching-month Relative Oceanic Niño Index values to observations. |
| `files`              | yes      | One or more NetCDF files from TCTrack.                                                                   |

### Examples

Import a single file into a new database:
```bash
build-db --output tracks.db tracks.nc
```

Import into a named collection:
```bash
build-db --output tracks.db --collection hadgem3 tracks.nc
```

Append files to an existing database and collection:
```bash
build-db --output tracks.db --collection hadgem3 more_tracks.nc
```

Import tracks with sea surface temperature anomalies from the required HadSST median monthly dataset:
```bash
build-db --output tracks.db --hadsst HadSST.4.2.0.0_median.nc tracks.nc
```

Import tracks with matching-month values from a monthly RONI dataset:
```bash
build-db --output tracks.db --roni roni.nc tracks.nc
```

## Database Schema

```
collections              Named groups of track files.
└─ files                 Individual NetCDF files with metadata.
   └─ trajectories       Cyclone tracks.
      ├─ observations    Individual observation rows.
      └─ statistics      Summary statistics for each track.

oceans                   Ocean basins.
```

See [`schema.sql`](schema.sql) for the full schema definition.


## Trajectory Statistics

Each track has summary statistics calculated from its observations.

| Columns | Definition |
| --- | --- |
| `date_start`, `date_end`, `date_median` | Earliest, latest and median observation timestamps. An even-count median is the midpoint of the two central timestamps. |
| `duration_hours` | Elapsed hours between the start and end. |
| `observation_count` | Number of stored observations. |
| `wind_speed_min`, `wind_speed_max`, `wind_speed_median`, `wind_speed_mean` | Wind summaries in m/s. |
| `air_pressure_at_sea_level_min` | Minimum sea-level pressure in the stored metric units. |
| `ace` | Accumulated Cyclone Energy, calculated as described below. |
| `sea_surface_temp_anomaly_mean`, `roni_mean` | Means of the stored environmental values. |

ACE uses the researchers' approved calculation over all wind observations:

```text
ace = sum(1e-4 * (wind_speed_ms * 1.944)^2)
```

Every stored observation with a valid wind speed contributes; there is no wind threshold, storm-status filtering, cadence selection, interpolation or time weighting.


## Notes

- The database is SQLite with foreign keys enabled.
- Trajectories in NetCDF files may be padded with NaN values. Only valid observations are imported.
- Duplicate files are skipped (by filename) on attempted re-import.
- Every observation is classified using the geographic module: whether it is over land, its distance to the nearest coastline and its ocean basin.
- Trajectories record the ocean basin of their first (genesis) observation and whether any observation made landfall.
- `--hadsst` requires the HadSST median monthly anomalies NetCDF, e.g. [`HadSST.4.2.0.0_median.nc`](https://www.metoffice.gov.uk/hadobs/hadsst4/data/download.html) via [Met Office Hadley Centre](https://www.metoffice.gov.uk/hadobs/hadsst4/). Values are selected by calendar year and month, so `360_day` model dates such as February 30 match their corresponding monthly HadSST field.
- `--roni` requires a NetCDF file with one-dimensional `time` and `value` variables, e.g. [`roni.nc`](https://psl.noaa.gov/data/timeseries/month/data/roni.nc) via [NOAA Physical Sciences Laboratory](https://psl.noaa.gov/data/timeseries/month/DS/RONI/). Each observation receives the value for its calendar year and month; months without a value remain null.
