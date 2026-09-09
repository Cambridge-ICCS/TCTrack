"""Datasette hooks for the MapLibre plugin."""

import json

from datasette import hookimpl

PLUGIN = "datasette-maplibre"

MAPLIBRE_CSS = "https://unpkg.com/maplibre-gl@^6.8.0/dist/maplibre-gl.css"
LAYER_CONTROL_CSS = "https://unpkg.com/maplibre-gl-layer-control@^0.17.4/dist/maplibre-gl-layer-control.css"
DEFAULT_BASEMAP = "https://demotiles.maplibre.org/style.json"


def find_geojson_columns(columns):
    """Return columns with names beginning "geojson"."""
    return [c for c in (columns or []) if c.lower().startswith("geojson")]


# ruff: noqa: PLR0913 ARG001

@hookimpl
def extra_css_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre CSS URLs when plugin triggered."""
    if view_name not in ("database", "table"):
        return []
    if not find_geojson_columns(columns):
        return []

    return [
        MAPLIBRE_CSS,
        LAYER_CONTROL_CSS,
        datasette.urls.static_plugins(PLUGIN, "map.css"),
    ]


@hookimpl
def extra_js_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre JS module URL when plugin triggered."""
    if view_name not in ("database", "table"):
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
    if view_name not in ("database", "table"):
        return ""

    geojson_columns = find_geojson_columns(columns)
    if not geojson_columns:
        return ""

    # Get optional layer column - name format: layer_<geojson property field>
    layer_column = next(
        (c for c in (columns or []) if c.lower().startswith("layer_")), None
    )

    config = datasette.plugin_config(PLUGIN, database=database, table=table) or {}
    basemap_style = config.get("basemap", DEFAULT_BASEMAP)
    layer_palette = config.get("layer_palette",
        ("#ee7733", "#0077bb", "#33bbee", "#ee3377", "#cc3311", "#009988", "#bbbbbb"))

    # Pass configuration to map.js via the JavaScript window object
    return (
        f"window.DATASETTE_MAPLIBRE_STYLE = {json.dumps(basemap_style)};\n"
        f"window.DATASETTE_MAPLIBRE_GEOJSON_COLUMNS = {json.dumps(geojson_columns)};\n"
        f"window.DATASETTE_MAPLIBRE_LAYER_COLUMN = {json.dumps(layer_column)};\n"
        f"window.DATASETTE_MAPLIBRE_LAYER_PALETTE = {json.dumps(layer_palette)};\n"
    )
