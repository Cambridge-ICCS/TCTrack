# TCTrack Dashboard

A web-based dashboard for viewing and comparing cyclone trajectories. Tracks and observations are viewed on a map along with other environmental data for context.

The dashboard runs on the [Datasette platform](https://datasette.io/).

## Installation

The custom Datasette plugin, `datasette-maplibre`, renders a [MapLibre GL JS](https://maplibre.org/projects/gl-js/) map when columns named `latitude` and `longitude` are detected.

Datasette and the `datasette-maplibre` plugin are installed together via the [dashboard] extra. Run from the project root:

	pip install -e .[dashboard]

For plugin development only, install the plugin directly (the Datasette dependency will also be installed):

	pip install -e dashboard/plugins/datasette_maplibre


## Starting the Dashboard

A database of track files must be built for loading into the dashboard. See the [build_db instructions](src\tctrack\build_db\README.md) for instructions.

Dashboard configuration files for metadata and style are located in the `dashboard` directory. If run from the project root:

	datasette serve <your_database.db> --metadata dashboard/metadata.yaml --static static:dashboard/static --setting max_returned_rows 6000 --setting sql_time_limit_ms 3000

- `metadata.yaml` applies the dashboard theme, metadata, table configuration and server settings.
- `--static static:dashboard/static` serves the dashboard/static directory at location `/static/`, so that files linked in the metadata can be loaded.
- `--setting max_returned_rows` [optional, default 1000] sets the maximum number of rows that can be loaded in one request. This should be set to accommodate the total number of observations in your dataset *if all data is needed to be seen on the map at any time*. If not, this can be set to a lower number to reduce load times.
- `--setting sql_time_limit_ms 3000` [optional, default 1000] sets the timeout for SQL query execution. Increase this if complex queries on large datasets return HTTP 400 errors.

View a running dashboard at http://localhost:8001


## General Usage

Datasette is a powerful tool for exploring data. Preset views of the data are available as part of the schema. Custom SQL queries can also be run.

Each view shows the data in table form along with filters for refining the dataset. These are especially useful for reducing data shown on the map. To help with filtering, columns with common values are shown as facets that can be selected quickly.

Resulting data from any view can be exported to JSON or CSV files. There is also a JSON API to allow access to external software.


## Map Display

The MapLibre map is shown automatically when a Datasette query, table or view has columns named `latitude` and `longitude`. All rows and columns are dynamically converted into GeoJSON for use as the map data source. Just as with the table view, data shown on the map is determined by the underlying dataset filters. There are also controls on the map to hide and show layers and switch between flat map and globe views. Map projection, zoom level and position are persisted in local storage.

Data is automatically split across map layers when a column name beginning `layer_` is found. Distinct values from this layer column determine the name of each new layer. See `map_file_layers_view` in the database and the dashboard for an example of this in action.


### Map Control

Scrolling via scroll-wheel or trackpad defaults to scrolling the whole page so that filters and the table view can be easily navigated to. Instead, the map can be zoomed in a few ways:

	Ctrl + scroll-wheel/trackpad
	+/- keys
	Double-click/tap
	Shift + Double-click/tap


### Configuration

The map can be configured in the `datasette-maplibre` section of the `metadata.yaml` file. The following settings are available:

- `basemap` URL of the MapLibre-compatible basemap to project onto the map. See the [basemap gallery](https://madewithmaplibre.com/basemaps/gallery) for alternatives.
- `group_by` A collection of settings for defining groups in the source data. This is primarily used to create map lines from the observation points of each track. Data is added to the group in the order given by the query so it should be in sequence, e.g. order by trajectory_id, sequence.
	- `column` The column name to group by. Use `trajectory_id` to define a group for each track.
	- `properties` An array of column names that will be used as properties for each group if they exist in the dataset. These are linked to the group and displayed when selected. The columns are presumed to contain unique values across each group.
- `layer_palette` Any number of hex colours to use for layers (in sequence). The default accessible set is taken from [Qualitative Colour Schemes](https://sronpersonalpages.nl/~pault/).
- `max_layers` The maximum number of layers allowed (the layer pair of lines and points is counted as one). Data in layers beyond the maximum is not shown; a warning is sent to the console log.


### Ideas for the Future

If the quantity or complexity of track data increases, it might be better served via the [MapLibre Tile Specification](https://maplibre.org/maplibre-tile-spec/) (MLT) instead of through standard Datasette delivery. This would increase performance and allow for much larger datasets. However, for this first version, leaning on the flexibility of Datasette gave the most immediate opportunities. Using MLT as the data source would require new map filtering capabilities.
