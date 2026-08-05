# Build Database

A CLI tool to import TCTrack NetCDF files into SQLite databases.

The resulting database forms the foundation for a dashboard. However, it is also valuable as a organisational tool for TCTrack output files. The data within can be queried and manipulated easily using local SQLite viewers such as [DB Browser](https://sqlitebrowser.org/).

Imported files can be grouped together into collections, e.g. geographic areas or a particular research hypothesis. Collections are defined using the `--collection` option. Collections are optional.

Longitude coordinates are wrapped from 0-360 to -180-180 for compatibility with the GeoJSON standard (RFC 7946).

Storm tracks and properties are converted into GeoJSON format in the `trajectories` table. These representations can be viewed directly in GIS tools.


## Usage

Presuming TCTrack has been installed via `pip install`:

```bash
build-db [--help] --output OUTPUT [--collection COLLECTION] files...
```

To use the tool without installation, bypassing the need for all other TCTrack dependencies:

```bash
cd /src/tctrack/build_db
python -m build_db --output OUTPUT files...
```

### Arguments

| Flag           | Required | Description |
|----------------|----------|----------------------------------------------------------------------------------------------------------|
| `--output`     | yes      | SQLite database path and name. Created if it doesn't exist, appended to if it does.                      |
| `--collection` | no       | Collection name. Creates or reuses an existing collection. Will use a default collection if unspecified. |
| `files`        | yes      | One or more NetCDF files from TCTrack.                                                                   |

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

## Database Schema

```
collections              Named groups of track files.
└─ files                 Individual NetCDF files with metadata.
   └─ trajectories       Cyclone tracks stored as GeoJSON FeatureCollections.
      └─ observations    Individual observation rows.
```

See [`schema.sql`](schema.sql) for the full schema definition.


## Notes

- The database is SQLite with foreign keys enabled.
- Trajectories in NetCDF files may be padded with `NaN` values. Only valid observations are imported.
- Duplicate files are skipped (by filename) on attempted re-import.
