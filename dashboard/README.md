# TCTrack Dashboard

A web-based dashboard for viewing and comparing cyclone trajectories. Tracks and observations are viewed on a map along with other environmental data for context.

The dashboard runs on the [Datasette platform](https://datasette.io/).

## Installation

The custom Datasette plugin, `datasette-maplibre`, renders a [MapLibre GL JS](https://maplibre.org/projects/gl-js/) map when one or more columns prefixed *geojson* are detected, e.g. trajectories.geojson_track.

Datasette and the `datasette-maplibre` plugin are installed together via the [dashboard] extra. Run from the project root:

	pip install -e .[dashboard]

For plugin development only, install the plugin directly (the Datasette dependency will also be installed):

	pip install -e dashboard/plugins/datasette_maplibre


## Starting the Dashboard

A database of track files must be built for loading into the dashboard. See the [build_db instructions](src\tctrack\build_db\README.md) for instructions.

Dashboard configuration files for metadata and style are located in the `dashboard` directory. If run from the project root:

	datasette serve <your_database.db> --metadata dashboard/metadata.yaml --static static:dashboard/static --setting max_returned_rows 6000

- `metadata.yaml` applies the dashboard theme, metadata, table configuration and server settings.
- `--static static:dashboard/static` serves the dashboard/static directory at location `/static/`, so that files linked in the metadata can be loaded.
- `--setting max_returned_rows` sets the maximum number of rows that can be loaded in one request. This should be set to accommodate the total number of trajectories in your dataset.

View a running dashboard at http://localhost:8001
