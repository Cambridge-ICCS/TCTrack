"""Test MapLibre hook activation for observation and trajectory geometry."""

import runpy
from pathlib import Path
from unittest.mock import Mock

import pytest

pytest.importorskip("datasette")

PLUGIN_PATH = (
    Path(__file__).resolve().parents[3]
    / "dashboard"
    / "plugins"
    / "datasette_maplibre"
    / "__init__.py"
)
PLUGIN = runpy.run_path(str(PLUGIN_PATH))


@pytest.mark.parametrize(
    ("columns", "active"),
    [
        (None, False),
        ([], False),
        (["latitude"], False),
        (["longitude"], False),
        (["trajectory_id"], False),
        (["LATITUDE", "Longitude"], True),
        (["LAT_LON"], True),
        (["lat_lon", "latitude", "longitude", "layer_collection"], True),
    ],
)
def test_geometry_activates_map_hooks(columns, active):
    """Enable all map hooks for either geometry format, preserving configuration."""
    datasette = Mock()
    datasette.static.side_effect = lambda name, plugin: f"/{plugin}/{name}"
    datasette.plugin_config.return_value = {
        "group_by": {"column": "trajectory_id"},
        "line": {"thickness": 3},
    }
    args = {
        "template": None,
        "database": "tracks",
        "table": "track_stats",
        "columns": columns,
        "view_name": "table",
        "request": None,
        "datasette": datasette,
    }

    assert PLUGIN["has_lat_lon_columns"](columns) is active
    assert bool(PLUGIN["extra_css_urls"](**args)) is active
    assert bool(PLUGIN["extra_js_urls"](**args)) is active
    script = PLUGIN["extra_body_script"](**args)
    assert bool(script) is active
    if active:
        assert (
            'window.DATASETTE_MAPLIBRE_GROUP_BY = {"column": "trajectory_id"};'
            in script
        )
        assert 'window.DATASETTE_MAPLIBRE_LINE = {"thickness": 3};' in script

    args["view_name"] = "index"
    assert PLUGIN["extra_css_urls"](**args) == []
    assert PLUGIN["extra_js_urls"](**args) == []
    assert PLUGIN["extra_body_script"](**args) == ""
