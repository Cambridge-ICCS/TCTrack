"""Tests for the build_db CLI module."""

import json
import math
import re
import sqlite3
from itertools import groupby
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tctrack.build_db import build_db
from tctrack.geographic import geographic

NETCDF_FILE = str(Path(__file__).parent / "test_tracks.nc")
NETCDF_SE_FILE = str(Path(__file__).parent / "test_tracks_se.nc")
NON_STANDARD_FILE = str(Path(__file__).parent / "test_non_standard.nc")
NETCDF_NATL_FILE = str(Path(__file__).parent / "6hr_track_mm_2013_natl.nc")
NETCDF_360_DAY_FILE = str(
    Path(__file__).parent / "HadGEM3-GC31-HH_r1i1p1f1_tracks_2000.nc"
)
HADSST_2000_FILE = str(Path(__file__).parent / "HadSST.4.2.0.0_median_2000.nc")
RONI_FILE = str(Path(__file__).parent / "roni.nc")


@pytest.fixture
def tmp_db(tmp_path):
    """Return a unique temporary database path."""
    return str(tmp_path / "test.db")


def open_db(db_path):
    """Open SQLite database."""
    db = sqlite3.connect(db_path)
    db.execute("pragma foreign_keys = on")
    return db


GATE_HEADER = "code,name,min_lat,max_lat,min_lon,max_lon\n"
GATE_ROWS = [
    "S,South,0,30,-180,180\n",
    "N,North,30,90,-180,180\n",
]


@pytest.fixture
def gates_csv(tmp_path):
    """Provide two shared-edge gates covering northern hemisphere landfalls."""
    path = tmp_path / "gates.csv"
    path.write_text(GATE_HEADER + "".join(GATE_ROWS), encoding="utf-8")
    return path


class TestGates:
    """Test gate validation, persistence, summaries and CLI integration."""

    def test_invalid_configuration_does_not_touch_database(self, tmp_db, tmp_path):
        """Validate gates before creating or modifying any database."""
        missing = tmp_path / "missing.csv"
        invalid = tmp_path / "invalid.csv"
        invalid.write_text(GATE_HEADER + "A,Invalid,0,nan,0,10\n", encoding="utf-8")

        for path, error in [(missing, OSError), (invalid, ValueError)]:
            with pytest.raises(error):
                build_db.build(tmp_db, NETCDF_FILE, gates_filepath=path)
            assert not Path(tmp_db).exists()

        build_db.build(tmp_db, NETCDF_FILE)
        before = Path(tmp_db).read_bytes()
        with pytest.raises(ValueError, match="finite"):
            build_db.build(tmp_db, NETCDF_SE_FILE, gates_filepath=invalid)
        assert Path(tmp_db).read_bytes() == before

    def test_import_assigns_gates_and_trajectory_summaries(self, tmp_db, gates_csv):
        """Store gate definitions, observation assignments and first-gate summaries."""
        build_db.build(tmp_db, NETCDF_NATL_FILE, gates_filepath=gates_csv)

        with open_db(tmp_db) as db:
            assert db.execute("select * from gates order by id").fetchall() == [
                ("N", "North", 30, 90, -180, 180),
                ("S", "South", 0, 30, -180, 180),
            ]
            assert db.execute("pragma foreign_key_check").fetchall() == []
            assert set(db.execute("select distinct gate_id from points")) == {
                (None,),
                ("S",),
                ("N",),
            }
            assert db.execute(
                "select count(*) from points where landfall = 0 and gate_id is not null"
            ).fetchone() == (0,)

            for trajectory_id, gate_id in db.execute(
                "select trajectory_id, gate_id from tracks"
            ):
                observations = db.execute(
                    "select gate_id from observations "
                    "where trajectory_id = ? order by sequence",
                    (trajectory_id,),
                ).fetchall()
                first_gate = next((g for (g,) in observations if g is not None), None)
                assert gate_id == first_gate

    def test_reopening_reuse_conflict_skipping_and_disabled(self, tmp_db, gates_csv):
        """Require identical definitions, never backfill or implicitly enable gates."""
        build_db.build(tmp_db, NETCDF_FILE)
        build_db.build(tmp_db, NETCDF_FILE, gates_filepath=gates_csv)
        build_db.build(tmp_db, NETCDF_NATL_FILE, gates_filepath=gates_csv)
        gates_csv.write_text(
            GATE_HEADER + "".join(reversed(GATE_ROWS)), encoding="utf-8"
        )
        build_db.build(tmp_db, NETCDF_NATL_FILE, gates_filepath=gates_csv)
        build_db.build(tmp_db, NETCDF_SE_FILE)
        with open_db(tmp_db) as db:
            assert db.execute("select count(*) from files").fetchone() == (3,)
            assert db.execute(
                "select count(*) from observations ob join trajectories tr "
                "on tr.id = ob.trajectory_id join files f on f.id = tr.file_id "
                "where f.filename in ('test_tracks.nc', 'test_tracks_se.nc') "
                "and (ob.gate_id is not null or tr.gate_id is not null)"
            ).fetchone() == (0,)
        before = Path(tmp_db).read_bytes()
        gates_csv.write_text(GATE_HEADER + GATE_ROWS[0], encoding="utf-8")
        with pytest.raises(ValueError, match="differs"):
            build_db.build(
                tmp_db, NETCDF_SE_FILE, "new collection", gates_filepath=gates_csv
            )
        assert Path(tmp_db).read_bytes() == before

    def test_cli_gates_and_errors(self, tmp_db, tmp_path, gates_csv, capsys):
        """Accept valid optional gates and report configuration errors clearly."""
        assert (
            build_db.main(
                ["--output", tmp_db, "--gates", str(gates_csv), NETCDF_NATL_FILE]
            )
            == 0
        )
        with open_db(tmp_db) as db:
            assert db.execute("select count(*) from gates").fetchone() == (2,)
        invalid_db = tmp_path / "invalid.db"
        gates_csv.write_text("invalid header", encoding="utf-8")
        with pytest.raises(SystemExit) as exc:
            build_db.main(
                ["--output", str(invalid_db), "--gates", str(gates_csv), NETCDF_FILE]
            )
        assert exc.value.code == 2
        assert str(gates_csv) in capsys.readouterr().err
        assert not invalid_db.exists()


class TestStatisticsIntegration:
    """Test that database building produces trajectory statistics."""

    def test_statistics_match_observation_counts(self, tmp_db):
        """Build one statistics row per observed trajectory with the correct count."""
        build_db.build(tmp_db, NETCDF_FILE)

        with open_db(tmp_db) as db:
            expected = db.execute(
                "select trajectory_id, count(*) from observations "
                "group by trajectory_id order by trajectory_id"
            ).fetchall()
            actual = db.execute(
                "select trajectory_id, observation_count from statistics "
                "order by trajectory_id"
            ).fetchall()

        assert expected
        assert actual == expected

        with open_db(tmp_db) as db:
            view_rows = db.execute(
                "select trajectory_id, lat_lon, observation_count "
                "from tracks order by trajectory_id"
            ).fetchall()

        assert len(view_rows) == len(expected)
        for (trajectory_id, lat_lon, count), expected_row in zip(
            view_rows, expected, strict=True
        ):
            assert (trajectory_id, count) == expected_row
            assert len(json.loads(lat_lon)) == count


class TestFileImport:
    """Test import of a single NetCDF file."""

    def test_creates_database(self, tmp_db):
        """Test that build_db.main creates a database file."""
        build_db.build(tmp_db, NETCDF_FILE)
        assert Path(tmp_db).exists()

    def test_file_record(self, tmp_db):
        """Test record inserted into the files table."""
        build_db.build(tmp_db, NETCDF_FILE)
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
        build_db.build(tmp_db, NETCDF_FILE)
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)
        count = db.execute("select count(*) from files").fetchone()[0]
        assert count == 1

    def test_invalid_hadsst_fails_before_database_creation(self, tmp_db):
        """Test an unusable HadSST path fails before importing any files."""
        with pytest.raises(ValueError, match="HadSST dataset does not exist"):
            build_db.build(tmp_db, NETCDF_FILE, hadsst_filepath="missing-hadsst.nc")

        assert not Path(tmp_db).exists()

    def test_invalid_roni_fails_before_database_creation(self, tmp_db):
        """Test an unusable RONI path fails before importing any files."""
        with pytest.raises(ValueError, match="RONI dataset does not exist"):
            build_db.build(tmp_db, NETCDF_FILE, roni_filepath="missing-roni.nc")

        assert not Path(tmp_db).exists()

    def test_file_metadata(self, tmp_db):
        """Test file metadata (tracker, version, time units)."""
        build_db.build(tmp_db, NETCDF_FILE)
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
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        row = db.execute("select tracker_parameters from files").fetchone()
        assert row[0] is not None
        params = json.loads(row[0])
        assert isinstance(params, dict)


class TestTrajectoryImport:
    """Test import of trajectories."""

    @pytest.mark.parametrize("longitude", [np.nan, np.inf, -np.inf])
    def test_invalid_coordinates_rejected(self, tmp_db, longitude):
        """Reject non-finite geometry before inserting a trajectory."""
        build_db.build(tmp_db, NETCDF_FILE)
        netcdf_data = build_db.read_netcdf(NETCDF_FILE)
        traj = build_db.extract_trajectory(netcdf_data, 0)
        traj["longitude"][traj["indices"][0]] = longitude

        with open_db(tmp_db) as db:
            count = db.execute("select count(*) from trajectories").fetchone()[0]
            with pytest.raises(ValueError, match="Out of range float values"):
                build_db.insert_trajectory(db, 1, traj, build_db.get_ocean_ids(db))

            assert (
                db.execute("select count(*) from trajectories").fetchone()[0] == count
            )

    @pytest.mark.parametrize("source", [NETCDF_FILE, NON_STANDARD_FILE])
    def test_lat_lon_matches_observations(self, tmp_db, source):
        """Store ordered JSON coordinate pairs for every trajectory."""
        build_db.build(tmp_db, source)

        with open_db(tmp_db) as db:
            trajectories = db.execute(
                "select id, lat_lon from trajectories order by id"
            ).fetchall()
            assert trajectories
            for trajectory_id, lat_lon in trajectories:
                observations = db.execute(
                    "select latitude, longitude from observations "
                    "where trajectory_id = ? order by sequence",
                    (trajectory_id,),
                ).fetchall()
                pairs = json.loads(lat_lon)
                assert pairs == [list(row) for row in observations]
                assert all(math.isfinite(value) for pair in pairs for value in pair)

    def test_trajectories_inserted(self, tmp_db):
        """Test that all trajectories are inserted into the database."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        count = db.execute("select count(*) from trajectories").fetchone()[0]
        assert count == 6

    def test_trajectories_start_end_flag_none(self, tmp_db):
        """Test start_end flag is unset in all trajectory rows."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        # Start/end flags are unset for all trajectories in test NetCDF file
        rows = db.execute("select start_end from trajectories").fetchall()
        for row in rows:
            assert row[0] is None

    def test_trajectories_start_end_flag_some(self, tmp_db):
        """Test start_end flag is set in some trajectory rows."""
        build_db.build(tmp_db, NETCDF_SE_FILE)
        db = open_db(tmp_db)

        rows = db.execute("select start_end from trajectories").fetchall()

        # First few rows have a starting flag set, the last few have none
        assert rows[0][0] == "S"
        assert rows[1][0] == "S"
        assert rows[-2][0] is None
        assert rows[-1][0] is None


class TestObservationImport:
    """Test import of observations."""

    def test_observation_count_per_trajectory(self, tmp_db):
        """Test that the correct number of observations exist per trajectory."""
        build_db.build(tmp_db, NETCDF_FILE)
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
        build_db.build(tmp_db, NETCDF_FILE)
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

    def test_non_standard_data_types(self, tmp_db):
        """Test that NumPy scalar values are stored as SQLite numbers."""
        build_db.build(tmp_db, NON_STANDARD_FILE)
        db = open_db(tmp_db)

        row = db.execute(
            "select latitude, longitude, wind_speed from observations "
            "where wind_speed is not null limit 1"
        ).fetchone()

        assert row is not None
        assert all(isinstance(value, (int, float)) for value in row)

    def test_observation_date(self, tmp_db):
        """Test that dates are in ISO format and monotonically increasing."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        rows = db.execute(
            "select trajectory_id, date from observations "
            "order by trajectory_id, sequence"
        ).fetchall()

        for _traj_id, date_str in rows:
            assert re.match(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$", date_str)

        # Dates within each trajectory should increase monotonically
        for _traj_id, group in groupby(rows, key=lambda r: r[0]):
            dates = [r[1] for r in group]
            assert dates == sorted(dates)

    def test_coords_match_source_netcdf(self, tmp_db):
        """Test that database coordinates match the source NetCDF file."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        ds = netCDF4.Dataset(NETCDF_FILE, "r")
        lats = ds.variables["latitude"][:]
        lons = ds.variables["longitude"][:]

        # Collect valid (non-NaN) coords per trajectory, same as build_db.
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
            "order by trajectory_id, sequence"
        ).fetchall()

        assert len(rows) == len(expected)
        for (db_traj, db_lat, db_lon), (exp_traj, exp_lat, exp_lon) in zip(
            rows, expected, strict=False
        ):
            assert db_traj == exp_traj
            assert db_lat == pytest.approx(exp_lat)
            assert db_lon == pytest.approx(exp_lon)

    def test_observation_sequence(self, tmp_db):
        """Test that sequence is 0-based, contiguous, and unique per trajectory."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        rows = db.execute(
            "select trajectory_id, sequence from observations "
            "order by trajectory_id, sequence"
        ).fetchall()

        # Group by trajectory and verify each group is 0-based and contiguous
        for _traj_id, group in groupby(rows, key=lambda r: r[0]):
            sequences = [r[1] for r in group]
            assert sequences == list(range(len(sequences)))

    def test_sea_surface_temp_anomaly_without_hadsst_is_null(self, tmp_db):
        """Test the anomaly remains null when a HadSST dataset is not supplied."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        values = db.execute(
            "select distinct sea_surface_temp_anomaly from observations"
        ).fetchall()

        assert values == [(None,)]

    def test_roni_without_dataset_is_null(self, tmp_db):
        """Test RONI remains null when a RONI dataset is not supplied."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        values = db.execute("select distinct roni from observations").fetchall()

        assert values == [(None,)]

    def test_roni_from_dataset(self, tmp_db):
        """Test observations receive values from the supplied RONI dataset."""
        build_db.build(
            tmp_db,
            NETCDF_360_DAY_FILE,
            roni_filepath=RONI_FILE,
        )
        db = open_db(tmp_db)

        roni = db.execute(
            "select roni from observations where trajectory_id = 1 and sequence = 0"
        ).fetchone()[0]

        assert roni == pytest.approx(-1.59)

    def test_sea_surface_temp_anomaly_from_hadsst(self, tmp_db):
        """Test observations receive an anomaly from the supplied HadSST dataset."""
        build_db.build(
            tmp_db,
            NETCDF_360_DAY_FILE,
            hadsst_filepath=HADSST_2000_FILE,
        )
        db = open_db(tmp_db)

        ssta = db.execute(
            "select sea_surface_temp_anomaly from observations "
            "where trajectory_id = 1 and sequence = 0"
        ).fetchone()[0]

        assert ssta is not None

    def test_360_day_observations_match_hadsst_month(self, tmp_db):
        """Test 360_day dates, including February 30, use February HadSST."""
        build_db.build(
            tmp_db,
            NETCDF_360_DAY_FILE,
            hadsst_filepath=HADSST_2000_FILE,
        )
        db = open_db(tmp_db)

        rows = db.execute(
            "select sea_surface_temp_anomaly from observations "
            "where trajectory_id = 16 and sequence between 102 and 104 "
            "order by sequence"
        ).fetchall()

        assert len(rows) == 3
        assert all(value is not None for (value,) in rows)

    def test_missing_hadsst_grid_cell_remains_null(self, tmp_db):
        """Test missing HadSST source data is imported as null anomalies."""
        with netCDF4.Dataset(HADSST_2000_FILE) as dataset:
            # February 2000, latitude -17.5, longitude 142.5. The observations
            # below all resolve to this 5-degree HadSST cell.
            source_ssta = dataset.variables["tos"][1, 14, 64]

        assert np.ma.is_masked(source_ssta)

        build_db.build(
            tmp_db,
            NETCDF_360_DAY_FILE,
            hadsst_filepath=HADSST_2000_FILE,
        )
        db = open_db(tmp_db)

        values = db.execute(
            "select sea_surface_temp_anomaly from observations "
            "where trajectory_id = 16 and sequence between 27 and 34 "
            "order by sequence"
        ).fetchall()

        assert values == [(None,)] * 8


class TestClassification:
    """Test geographic classification of imported trajectories."""

    def test_open_ocean_tracks(self, tmp_db):
        """Test classification of tracks that stay over the open ocean."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        north_pacific = db.execute(
            "select id from oceans where name = 'North Pacific Ocean'"
        ).fetchone()[0]

        # Every observation is at sea, far from land, in the North Pacific
        rows = db.execute(
            "select landfall, distance_to_coast_km, ocean_id from observations"
        ).fetchall()
        for landfall, distance_to_coast_km, ocean_id in rows:
            assert landfall == 0
            assert distance_to_coast_km > 400
            assert ocean_id == north_pacific

        # No trajectory makes landfall, all form in the North Pacific
        rows = db.execute("select landfall, ocean_id from trajectories").fetchall()
        for landfall, ocean_id in rows:
            assert landfall == 0
            assert ocean_id == north_pacific

    def test_landfalling_tracks(self, tmp_db):
        """Test classification of tracks that cross land."""
        build_db.build(tmp_db, NETCDF_NATL_FILE)
        db = open_db(tmp_db)

        oceans = {
            name: ocean_id
            for ocean_id, name in db.execute("select id, name from oceans")
        }
        atlantic = oceans["North Atlantic Ocean"]
        pacific = oceans["North Pacific Ocean"]

        trajectories = db.execute(
            "select landfall, ocean_id from trajectories"
        ).fetchall()
        observations = db.execute(
            "select landfall, ocean_id from observations"
        ).fetchall()

        # Some but not all tracks make landfall (Central America, Yucatan,
        # the US East Coast and Nova Scotia in 2013), and only some
        # observations are over land
        assert {row[0] for row in trajectories} == {0, 1}
        assert {row[0] for row in observations} == {0, 1}

        # Tracks form in the North Atlantic and, on the Pacific side of
        # Central America, the North Pacific
        assert {row[1] for row in trajectories} == {atlantic, pacific}
        assert {row[1] for row in observations} == {atlantic, pacific}


class TestOceans:
    """Test the seeded oceans type table."""

    def test_oceans_seeded_from_world_map(self, tmp_db):
        """Test that the oceans table holds the world map basin names."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)

        rows = {row[0] for row in db.execute("select name from oceans")}

        with np.load(geographic.DATA_PATH) as data:
            names = {str(name) for name in data["names"]}

        assert rows == names


class TestCollections:
    """Test file collections."""

    def test_default_collection(self, tmp_db):
        """Test that default collection is created when unspecified."""
        build_db.build(tmp_db, NETCDF_FILE)
        db = open_db(tmp_db)
        row = db.execute(
            "select name from collections where name = 'default'"
        ).fetchone()
        assert row is not None

    def test_named_collection(self, tmp_db):
        """Test that a named collection is created when specified."""
        build_db.build(tmp_db, NETCDF_FILE, collection_name="test")
        db = open_db(tmp_db)
        row = db.execute("select name from collections where name = 'test'").fetchone()
        assert row is not None

    def test_files_associated_with_collection(self, tmp_db):
        """Test that files are correctly associated with the specified collection."""
        build_db.build(tmp_db, NETCDF_FILE, collection_name="test")
        db = open_db(tmp_db)
        row = db.execute(
            "select c.name, f.filename from files f "
            "join collections c on f.collection_id = c.id"
        ).fetchone()
        assert row is not None
        assert row[0] == "test"

    def test_append_to_existing_collection(self, tmp_db):
        """Test that appending to an existing database use the correct collection."""
        build_db.build(tmp_db, NETCDF_FILE, collection_name="test")
        build_db.build(tmp_db, NETCDF_SE_FILE, collection_name="test")
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
        build_db.build(tmp_db, NETCDF_FILE)  # use default collection
        build_db.build(tmp_db, NETCDF_SE_FILE, collection_name="test")
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


class TestCli:
    """Test CLI behaviour."""

    def test_returns_zero(self, tmp_db):
        """Test that successful CLI invocation returns exit code 0."""
        rc = build_db.main(["--output", tmp_db, NETCDF_FILE])
        assert rc == 0

    def test_missing_output_errors(self):
        """Test that CLI exits with an error when required arguments are missing."""
        with pytest.raises(SystemExit):
            build_db.main([])

    def test_creates_database(self, tmp_db):
        """Test that build_db.main creates a database file."""
        build_db.main(["--output", tmp_db, NETCDF_FILE])
        assert Path(tmp_db).exists()

    def test_hadsst_option(self, tmp_db):
        """Test that --hadsst passes the lookup dataset to the importer."""
        build_db.main(
            [
                "--output",
                tmp_db,
                "--hadsst",
                HADSST_2000_FILE,
                NETCDF_360_DAY_FILE,
            ]
        )
        db = open_db(tmp_db)

        ssta = db.execute(
            "select sea_surface_temp_anomaly from observations "
            "where trajectory_id = 1 and sequence = 0"
        ).fetchone()[0]

        assert ssta is not None

    def test_roni_option(self, tmp_db):
        """Test that --roni passes the lookup dataset to the importer."""
        build_db.main(
            [
                "--output",
                tmp_db,
                "--roni",
                RONI_FILE,
                NETCDF_360_DAY_FILE,
            ]
        )
        db = open_db(tmp_db)

        roni = db.execute(
            "select roni from observations where trajectory_id = 1 and sequence = 0"
        ).fetchone()[0]

        assert roni == pytest.approx(-1.59)

    def test_named_collection(self, tmp_db):
        """Test that a named collection is created when --collection is specified."""
        build_db.main(["--output", tmp_db, "--collection", "test", NETCDF_FILE])
        db = open_db(tmp_db)
        row = db.execute("select name from collections where name = 'test'").fetchone()
        assert row is not None
