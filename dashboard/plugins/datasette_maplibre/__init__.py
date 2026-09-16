"""Datasette hooks for the MapLibre plugin."""

import json

from datasette import hookimpl

PLUGIN = "datasette-maplibre"

MAPLIBRE_CSS = "https://unpkg.com/maplibre-gl@6.8.0/dist/maplibre-gl.css"
LAYER_CONTROL_CSS = "https://unpkg.com/maplibre-gl-layer-control@0.17.4/dist/maplibre-gl-layer-control.css"
DEFAULT_BASEMAP = "https://demotiles.maplibre.org/style.json"


def has_lat_lon_columns(columns):
    """Return true when columns named latitude and longitude exist."""
    names = {c.lower() for c in (columns or [])}
    return "latitude" in names and "longitude" in names


# ruff: noqa: PLR0913 PLR0917 ARG001

@hookimpl
def extra_css_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre CSS URLs when plugin triggered."""
    if view_name not in ("database", "table") or not has_lat_lon_columns(columns):
        return []

    return [
        MAPLIBRE_CSS,
        LAYER_CONTROL_CSS,
        datasette.urls.static_plugins(PLUGIN, "map.css"),
    ]


@hookimpl
def extra_js_urls(template, database, table, columns, view_name, request, datasette):
    """Return MapLibre JS module URL when plugin triggered."""
    if view_name not in ("database", "table") or not has_lat_lon_columns(columns):
        return []

    return [
        {"url": datasette.urls.static_plugins(PLUGIN, "map.js"), "module": True},
    ]


@hookimpl
def extra_body_script(
    template, database, table, columns, view_name, request, datasette
):
    """Pass configuration to the JavaScript layer."""
    if view_name not in ("database", "table") or not has_lat_lon_columns(columns):
        return ""

    # Get optional layer column
    layer_column = next(
        (c for c in (columns or []) if c.lower().startswith("layer_")), None
    )

    config = datasette.plugin_config(PLUGIN, database=database, table=table) or {}
    basemap_style = config.get("basemap", DEFAULT_BASEMAP)
    group_by = config.get("group_by")
    layer_palette = config.get("layer_palette",
        ("#ee7733", "#0077bb", "#33bbee", "#ee3377", "#cc3311", "#009988", "#bbbbbb"))

    # Pass configuration to map.js via the JavaScript window object
    return (
        f"window.DATASETTE_MAPLIBRE_STYLE = {json.dumps(basemap_style)};\n"
        f"window.DATASETTE_MAPLIBRE_GROUP_BY = {json.dumps(group_by)};\n"
        f"window.DATASETTE_MAPLIBRE_LAYER_COLUMN = {json.dumps(layer_column)};\n"
        f"window.DATASETTE_MAPLIBRE_LAYER_PALETTE = {json.dumps(layer_palette)};\n"
    )
