import json

from datasette import hookimpl

PLUGIN = "datasette-maplibre"

MAPLIBRE_CSS = "https://unpkg.com/maplibre-gl@^6.4.1/dist/maplibre-gl.css"
DEFAULT_BASEMAP = "https://demotiles.maplibre.org/style.json"


def find_geojson_columns(columns):
    """Return columns with names beginning "geojson"."""
    return [c for c in (columns or []) if c.lower().startswith("geojson")]


@hookimpl
def extra_css_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre CSS URLs when plugin triggered."""
    if view_name not in ("table", "query"):
        return []
    if not find_geojson_columns(columns):
        return []

    return [
        MAPLIBRE_CSS,
        datasette.urls.static_plugins(PLUGIN, "map.css"),
    ]


@hookimpl
def extra_js_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre JS module URL when plugin triggered."""
    if view_name not in ("table", "query"):
        return []
    if not find_geojson_columns(columns):
        return []

    return [
        {"url": datasette.urls.static_plugins(PLUGIN, "map.js"), "module": True},
    ]


@hookimpl
def extra_body_script(
    template, database, table, columns, view_name, request, datasette
):
    """Pass configuration to the JavaScript layer."""
    if view_name not in ("table", "query"):
        return ""

    geojson_columns = find_geojson_columns(columns)
    if not geojson_columns:
        return ""

    config = datasette.plugin_config(PLUGIN, database=database, table=table) or {}
    basemap_style = config.get("basemap", DEFAULT_BASEMAP)

    # Pass this configuration to map.js via the JavaScript window object
    return (
        f"window.DATASETTE_MAPLIBRE_GEOJSON_COLUMNS = {json.dumps(geojson_columns)};\n"
        f"window.DATASETTE_MAPLIBRE_STYLE = {json.dumps(basemap_style)};\n"
    )
