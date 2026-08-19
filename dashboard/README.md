# TCTrack Dashboard

A web-based dashboard for viewing and comparing cyclone trajectories. Tracks and observations are viewed on a map along with other environmental data for context.


## Installation

A database of track files must be built for loading into the dashboard. See the [build_db instructions](src\tctrack\build_db\README.md) for instructions.

The dashboard runs on the [Datasette platform](https://datasette.io/). Datasette is Python-based and [easily installed](https://docs.datasette.io/en/stable/installation.html).

Datasette plugins are also required for mapping capabilities. Plugin, [datasette-cluster-map](https://datasette.io/plugins/datasette-cluster-map) is shown for tables with columns named `latitude` and `longitude`. The [datasette-geojson-map](https://datasette.io/plugins/datasette-geojson-map) is shown when a column named `geometry` is present.
The custom plugin, `datasette-maplibre`, renders a [MapLibre GL JS](https://maplibre.org/) map when one or more `geojson`-prefixed columns exist (e.g. `trajectories.geojson_track`).

Plugin installation:

	datasette install datasette-geojson-map datasette-cluster-map
	pip install -e dashboard/plugins/datasette_maplibre


## Starting the Dashboard

Dashboard configuration files for metadata and style are located in the `dashboard` directory.

	datasette serve <your_database.db> --metadata dashboard/metadata.yaml --static static:dashboard/static --setting max_returned_rows 6000

- `metadata.yaml` applies the dashboard theme, metadata, table configuration and server settings.
- `--static static:dashboard/static` serves the dashboard/static directory at location `/static/`, so that files linked in the metadata can be loaded.
- `--setting max_returned_rows` sets the maximum number of rows that can be loaded in one request. This should be set to accommodate the total number of trajectories in your dataset.

View a running dashboard at http://localhost:8001
