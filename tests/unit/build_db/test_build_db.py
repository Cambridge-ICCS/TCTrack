"""Tests for the build_db CLI module."""

import json
import math
import re
import sqlite3
from itertools import groupby
from pathlib import Path

import netCDF4
import pytest

from tctrack.build_db import build_db

NETCDF_FILE = str(Path(__file__).parent / "test_tracks.nc")
NETCDF_SE_FILE = str(Path(__file__).parent / "test_tracks_se.nc")


@pytest.fixture
def tmp_db(tmp_path):
    """Return a unique temporary database path."""
    return str(tmp_path / "test.db")


def open_db(db_path):
    """Open SQLite database."""
    db = sqlite3.connect(db_path)
    db.execute("pragma foreign_keys = on")
    return db


class TestFileImport:
    """Test import of a single NetCDF file."""

    def test_creates_database(self, tmp_db):
        """Test that build_db.main creates a database file."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        assert Path(tmp_db).exists()

    def test_file_record(self, tmp_db):
        """Test record inserted into the files table."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute(
            "select filepath, filename, trajectories, observations from files"
        ).fetchone()

        assert row is not None
        assert Path(row[0]).is_absolute()
        assert Path(row[0]) == Path(NETCDF_FILE)
        assert row[1] == "test_tracks.nc"
        assert row[2] == 6
        assert row[3] == 19

    def test_duplicate_file_skipped(self, tmp_db):
        """Test that specifiying the same file again does nothing."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)
        count = db.execute("select count(*) from files").fetchone()[0]
        assert count == 1

    def test_file_metadata(self, tmp_db):
        """Test file metadata (tracker, version, time units)."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute(
            "select tracker_name, tctrack_version, time_units, time_calendar from files"
        ).fetchone()

        assert row[0] == "TETracker"
        assert row[1] == "0.1.dev253+ga96ab077a"
        assert row[2] == "days since 1950-01-01"
        assert row[3] == "360_day"

    def test_tracker_parameters(self, tmp_db):
        """Test that tracker_parameters contains valid JSON."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute("select tracker_parameters from files").fetchone()
        assert row[0] is not None
        params = json.loads(row[0])
        assert isinstance(params, dict)


class TestTrajectoryImport:
    """Test import of trajectories."""

    def test_trajectories_inserted(self, tmp_db):
        """Test that all trajectories are inserted into the database."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        count = db.execute("select count(*) from trajectories").fetchone()[0]
        assert count == 6

    def test_trajectories_start_end_flag_none(self, tmp_db):
        """Test start_end flag is unset in all trajectory rows."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        # Start/end flags are unset for all trajectories in test NetCDF file
        rows = db.execute("select start_end from trajectories").fetchall()
        for row in rows:
            assert row[0] is None

    def test_trajectories_start_end_flag_some(self, tmp_db):
        """Test start_end flag is set in some trajectory rows."""
        build_db.main(["--output", tmp_db, NETCDF_SE_FILE])
        db = open_db(tmp_db)

        rows = db.execute("select start_end from trajectories").fetchall()

        # First few rows have a starting flag set, the last few have none
        assert rows[0][0] == "S"
        assert rows[1][0] == "S"
        assert rows[-2][0] is None
        assert rows[-1][0] is None


class TestGeoJSONConversion:
    """Test track conversion to GeoJSON."""

    def test_geojson_structure(self, tmp_db):
        """Test that geojson contains a LineString path and Point features."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute("select geojson from trajectories limit 1").fetchone()
        geojson = json.loads(row[0])

        assert geojson["type"] == "FeatureCollection"

        # First feature is the LineString path
        features = geojson["features"]
        assert features[0]["geometry"]["type"] == "LineString"
        assert features[0]["properties"]["source_file"] == "test_tracks.nc"

        # Remaining features are Points
        for f in features[1:]:
            assert f["geometry"]["type"] == "Point"
            assert "date" in f["properties"]
            assert "air_pressure_at_sea_level" in f["properties"]
            assert "surface_altitude" in f["properties"]
            assert "wind_speed" in f["properties"]

    def test_geojson_coordinate_count(self, tmp_db):
        """Test that geojson has the correct number of coordinates and features."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute("select geojson from trajectories limit 1").fetchone()
        geojson = json.loads(row[0])

        # First trajectory has 1 LineString and 11 valid observations
        features = geojson["features"]
        assert len(features) == 12
        assert len(features[0]["geometry"]["coordinates"]) == 11

    def test_geojson_start_end_flag(self, tmp_db):
        """Test that start_end is present in geojson."""
        build_db.main(["--output", tmp_db, NETCDF_SE_FILE])
        db = open_db(tmp_db)

        row = db.execute("select geojson from trajectories limit 1").fetchone()
        geojson = json.loads(row[0])

        linestring_props = geojson["features"][0]["properties"]
        assert "start_end" in linestring_props
        assert linestring_props["start_end"] == "S"


class TestObservationImport:
    """Test import of observations."""

    def test_observation_count_per_trajectory(self, tmp_db):
        """Test that the correct number of observations exist per trajectory."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        rows = db.execute(
            "select trajectory_id, count(*) "
            "from observations group by trajectory_id order by trajectory_id"
        ).fetchall()

        # Valid counts per trajectory
        assert rows[0][1] == 11
        assert rows[1][1] == 10
        assert rows[2][1] == 10
        assert rows[3][1] == 19
        assert rows[4][1] == 13
        assert rows[5][1] == 13

    def test_observation_values(self, tmp_db):
        """Test that observation values are real numbers."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        row = db.execute(
            "select latitude, longitude, "
            "air_pressure_at_sea_level, surface_altitude, wind_speed "
            "from observations where trajectory_id = 1 limit 1"
        ).fetchone()

        assert row[0] == pytest.approx(13.476562)
        assert row[1] == pytest.approx(183.339844)

        assert row[2] == pytest.approx(100472.1)
        assert row[3] == pytest.approx(0.0)
        assert row[4] == pytest.approx(10.42255)

    def test_observation_date(self, tmp_db):
        """Test that dates are in ISO format and monotonically increasing."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        # Observations ordered by rowid to retain insertion order
        rows = db.execute(
            "select trajectory_id, date from observations order by trajectory_id, rowid"
        ).fetchall()

        for _traj_id, date_str in rows:
            assert re.match(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$", date_str)

        # Dates within each trajectory should increase monotonically
        for _traj_id, group in groupby(rows, key=lambda r: r[0]):
            dates = [r[1] for r in group]
            assert dates == sorted(dates)

    def test_coords_match_source_netcdf(self, tmp_db):
        """Test that database coordinates match the source NetCDF file."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)

        ds = netCDF4.Dataset(NETCDF_FILE, "r")
        lats = ds.variables["latitude"][:]
        lons = ds.variables["longitude"][:]

        # Collect valid (non-NaN) coords per trajectory, same as build_db
        expected: list = []
        for traj_idx in range(lats.shape[0]):
            for obs_idx in range(lats.shape[1]):
                lat, lon = (
                    float(lats[traj_idx, obs_idx]),
                    float(lons[traj_idx, obs_idx]),
                )
                if not math.isnan(lat):
                    expected.append((traj_idx + 1, lat, lon))
        ds.close()

        rows = db.execute(
            "select trajectory_id, latitude, longitude from observations "
            "order by trajectory_id, date"
        ).fetchall()

        assert len(rows) == len(expected)
        for (db_traj, db_lat, db_lon), (exp_traj, exp_lat, exp_lon) in zip(
            rows, expected, strict=False
        ):
            assert db_traj == exp_traj
            assert db_lat == pytest.approx(exp_lat)
            assert db_lon == pytest.approx(exp_lon)


class TestCollections:
    """Test --collection option."""

    def test_default_collection(self, tmp_db):
        """Test that default collection is created when --collection is unspecified."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        db = open_db(tmp_db)
        row = db.execute(
            "select name from collections where name = 'default'"
        ).fetchone()
        assert row is not None

    def test_named_collection(self, tmp_db):
        """Test that a named collection is created when --collection is specified."""
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_FILE])
        db = open_db(tmp_db)
        row = db.execute("select name from collections where name = 'test'").fetchone()
        assert row is not None

    def test_files_associated_with_collection(self, tmp_db):
        """Test that files are correctly associated with the specified collection."""
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_FILE])
        db = open_db(tmp_db)
        row = db.execute(
            "select c.name, f.filename from files f "
            "join collections c on f.collection_id = c.id"
        ).fetchone()
        assert row is not None
        assert row[0] == "test"

    def test_append_to_existing_collection(self, tmp_db):
        """Test that appending to an existing database use the correct collection."""
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_FILE])
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_SE_FILE])
        db = open_db(tmp_db)

        rows = db.execute(
            "select c.name, count(f.id) as files "
            "from collections c left join files f on f.collection_id = c.id "
            "group by c.name"
        ).fetchall()

        assert len(rows) == 1
        assert rows[0][0] == "test"
        assert rows[0][1] == 2

    def test_multiple_collections(self, tmp_db):
        """Test use of two collections."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])  # use default collection
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_SE_FILE])
        db = open_db(tmp_db)

        rows = db.execute(
            "select c.name, count(f.id) as files "
            "from collections c left join files f on f.collection_id = c.id "
            "group by c.name"
        ).fetchall()

        assert len(rows) == 2
        assert rows[0][0] == "default"
        assert rows[0][1] == 1
        assert rows[1][0] == "test"
        assert rows[1][1] == 1


class TestWrapLongitude:
    """Test the --wrap flag for longitude coordinate wrapping."""

    def test_wrap_longitude(self, tmp_db):
        """Test that --wrap stores longitudes in the -180 to 180 range."""
        build_db.main(["--output", tmp_db, "--wrap", NETCDF_FILE])
        db = open_db(tmp_db)
        rows = db.execute("select longitude from observations").fetchall()
        for (lon,) in rows:
            assert -180 <= lon <= 180

    def test_wrap_in_geojson(self, tmp_db):
        """Test that --wrap affects longitude values in geojson output."""
        build_db.main(["--output", tmp_db, "--wrap", NETCDF_FILE])
        db = open_db(tmp_db)
        rows = db.execute("select geojson from trajectories").fetchall()
        for (geojson_str,) in rows:
            geojson = json.loads(geojson_str)
            for feature in geojson["features"]:
                geom = feature["geometry"]
                if geom["type"] == "LineString":
                    for lon, _lat in geom["coordinates"]:
                        assert -180 <= lon <= 180
                elif geom["type"] == "Point":
                    lon, _lat = geom["coordinates"]
                    assert -180 <= lon <= 180


class TestCliExitCode:
    """Test CLI exit behaviour."""

    def test_returns_zero(self, tmp_db):
        """Test that successful CLI invocation returns exit code 0."""
        rc = build_db.main(["--output", tmp_db, NETCDF_FILE])
        assert rc == 0

    def test_missing_output_errors(self):
        """Test that CLI exits with an error when required arguments are missing."""
        with pytest.raises(SystemExit):
            build_db.main([])
