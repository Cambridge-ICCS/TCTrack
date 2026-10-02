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

	datasette serve <your_database.db> --metadata dashboard/metadata.yaml --config .\dashboard\datasette.yaml --static static:dashboard/static

- `metadata.yaml` applies the database and table metadata.
- `datasette.yaml` applies the dashboard theme, map configuration and server settings.
- `static:dashboard/static` serves the `dashboard/static` directory at location `/static/`, so that files linked in the metadata can be loaded.

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

The map can be configured in the `datasette-maplibre` section of the `metadata.yaml` file. *The `datasette` server process must be restarted for changes to take effect*.

The following settings are available:

- `basemap` URL of the MapLibre-compatible basemap to project onto the map. See the [basemap gallery](https://madewithmaplibre.com/basemaps/gallery) for alternatives. The MapTiler dataviz map ([light](https://www.maptiler.com/maps/#style=dataviz-v4)/[dark](https://www.maptiler.com/maps/#style=dataviz-v4-dark)) has proven effective (API key required).
- `group_by` A collection of settings for defining groups in the source data. This is primarily used to create map lines from the observation points of each track. Data is added to the group in the order given by the query so it should be in sequence, e.g. `order by trajectory_id, sequence`.
	- `column` The column name to group by. Use `trajectory_id` to define a group for each track.
	- `properties` An array of column names that will be used as properties for each group if they exist in the dataset. These are linked to the group and displayed when selected. The columns are presumed to contain unique values across each group.
- `point` Point rendering options:
	- `radius` Standard point radius.
	- `radius_property` Link radius size to a property value:
		- `column`: Column name to use as the linked property.
		- `min`: Value/radius pair for smallest point, e.g. [0, 2] when property = 0, radius = 2.
		- `max`: Value/radius pair for largest point, e.g. [30, 12] when property >= 30, radius = 12.
		- `exponential`: Exponential base for radius scaling.
	- `opacity_property`: Link opacity (0.0–1.0) to a property value.
		- `column`: Column name to use as the linked property.
		- `min`: Value/opacity pair for smallest value, e.g. [0, 0.1] when property = 0, opacity = 0.1.
		- `max`: Value/opacity pair for largest value, e.g. [25, 1] when property >= 25, opacity = 1.0.
- `line` Line rendering options:
	- `thickness`: Line thickness.
	- `opacity_property`: Link opacity (0.0–1.0) to a property value.
		- `column`: Column name to use as the linked property.
		- `min`: Value/opacity pair for smallest value, e.g. [0, 0.1] when property = 0, opacity = 0.1.
		- `max`: Value/opacity pair for largest value, e.g. [25, 1] when property >= 25, opacity = 1.0.
- `heatmap` Configuration to create a density map layer:
	- `weight_property`: Link a property value to heatmap weights (0.0–1.0):
		- `column`: Column name to use as the linked property.
		- `min`: Value/weight pair for lowest weight, e.g. [0, 0] when property = 0, weight = 0.
		- `max`: Value/weight pair for strongest weight, e.g. [30, 1] when property >= 30, weight = 1.
		- `exponential`: Exponential base for weight scaling.
	- `palette` Colour ramp for heatmap density (`#rrggbb` hex format, no alpha).
		- `low`: Starting colour.
		- `mid`: Colour at midpoint.
		- `high`: Colour at maximum weight.
	- `max_zoom`: Maximum zoom level. The heatmap will be faded out by this level of zoom.
	- `opacity`: Heatmap opacity.
- `palette` Colours for rendering features (`#rrggbb[aa]` hex format). The default accessible set is taken from [Qualitative Colour Schemes](https://sronpersonalpages.nl/~pault/):
	- `single` Feature colour when there is only one layer.
	- `layers` Any number of colours to use for features on each layer. Colours are applied to layers in sequence. When all colours are used, subsequent layers will use the final colour in the list.
- `max_layers` The maximum number of layers allowed (the layer pair of lines and points is counted as one). Data in layers beyond the maximum is not shown; a warning is sent to the console log.


### Ideas for the Future

If the quantity or complexity of track data increases, it might be better served via the [MapLibre Tile Specification](https://maplibre.org/maplibre-tile-spec/) (MLT) instead of through standard Datasette delivery. This would increase performance and allow for much larger datasets. However, for this first version, leaning on the flexibility of Datasette gave the most immediate opportunities. Using MLT as the data source would require new map filtering capabilities.
