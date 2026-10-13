"""Tests for trajectory statistics calculated from stored observations."""

import sqlite3
from pathlib import Path

import pytest

from tctrack.build_db import build_db, statistics

NETCDF_FILE = str(Path(__file__).parent / "test_tracks.nc")
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
    """Open the database with foreign keys enabled."""
    db = sqlite3.connect(db_path)
    db.execute("pragma foreign_keys = on")
    return db


@pytest.fixture
def observation():
    """Read one existing NetCDF observation without building a database."""
    data = build_db.read_netcdf(NETCDF_FILE)
    track = build_db.extract_trajectory(data, 0)
    index = track["indices"][0]
    return {
        "date": track["times"][index].isoformat(" "),
        "wind_speed": float(track["wind_speed"][index]),
        "air_pressure_at_sea_level": float(track["air_pressure_at_sea_level"][index]),
        "sea_surface_temp_anomaly": None,
        "roni": None,
    }


class TestCalculations:
    """Test summary values and important input edge cases."""

    def test_invalid_input(self, observation):
        """Reject empty input, missing fields, invalid calendars and invalid dates."""
        for records in (None, [], [None]):
            with pytest.raises((TypeError, ValueError)):
                statistics.summarise_observations(records, "360_day")

        for field in observation:
            record = observation.copy()
            del record[field]
            with pytest.raises(ValueError, match=field):
                statistics.summarise_observations([record], "360_day")

        for calendar in (None, "", "not-a-calendar"):
            with pytest.raises(ValueError):
                statistics.summarise_observations([observation], calendar)

        for date in (None, "not-a-date", "2000-02-31 00:00:00"):
            with pytest.raises(ValueError, match="date"):
                statistics.summarise_observations(
                    [{**observation, "date": date}], "360_day"
                )

        with pytest.raises(ValueError, match="date"):
            statistics.summarise_observations(
                [{**observation, "date": "2000-02-30 00:00:00"}], "standard"
            )

    def test_missing_values(self, observation):
        """Missing and nonfinite measurements produce unavailable summaries."""
        records = [
            {
                **observation,
                "wind_speed": value,
                "air_pressure_at_sea_level": value,
                "sea_surface_temp_anomaly": value,
                "roni": value,
            }
            for value in (None, float("nan"), float("inf"), float("-inf"))
        ]
        summary = statistics.summarise_observations(records, "360_day")

        assert summary["observation_count"] == 4
        assert summary["wind_speed_min"] is None
        assert summary["wind_speed_max"] is None
        assert summary["wind_speed_median"] is None
        assert summary["wind_speed_mean"] is None
        assert summary["ace"] is None
        assert summary["air_pressure_at_sea_level_min"] is None
        assert summary["sea_surface_temp_anomaly_mean"] is None
        assert summary["roni_mean"] is None

    def test_summary_values(self, observation):
        """Summarise each field independently and count every valid wind sample."""
        record = {**observation, "date": "2000-02-30 01:00:00"}
        records = [
            {
                **record,
                "wind_speed": 0,
                "air_pressure_at_sea_level": 100000,
                "sea_surface_temp_anomaly": -2,
                "roni": -1,
            },
            {
                **record,
                "wind_speed": 5,
                "air_pressure_at_sea_level": None,
                "sea_surface_temp_anomaly": 0,
                "roni": -1,
            },
            {
                **record,
                "wind_speed": 10,
                "air_pressure_at_sea_level": 98000,
                "sea_surface_temp_anomaly": 4,
                "roni": -1,
            },
            {
                **record,
                "wind_speed": 20,
                "air_pressure_at_sea_level": float("inf"),
                "sea_surface_temp_anomaly": None,
                "roni": 3,
            },
            {
                **record,
                "wind_speed": float("nan"),
                "air_pressure_at_sea_level": 99000,
                "sea_surface_temp_anomaly": None,
                "roni": None,
            },
        ]
        summary = statistics.summarise_observations(records, "360_day")

        assert summary == pytest.approx(
            {
                "date_start": "2000-02-30 01:00:00",
                "date_end": "2000-02-30 01:00:00",
                "date_median": "2000-02-30 01:00:00",
                "duration_hours": 0.0,
                "observation_count": 5,
                "wind_speed_min": 0.0,
                "wind_speed_max": 20.0,
                "wind_speed_median": 7.5,
                "wind_speed_mean": 8.75,
                "air_pressure_at_sea_level_min": 98000.0,
                "ace": 0.19840464,
                "sea_surface_temp_anomaly_mean": 2 / 3,
                "roni_mean": 0.0,
            }
        )

    def test_zero_wind(self, observation):
        """A valid zero wind gives zero ACE, not a missing value."""
        summary = statistics.summarise_observations(
            [{**observation, "wind_speed": 0}], "360_day"
        )
        assert summary["wind_speed_mean"] == 0.0
        assert summary["ace"] == 0.0

    @pytest.mark.parametrize(
        ("calendar", "dates", "hours"),
        [
            (
                "360_day",
                ("2000-02-29 00:00:00", "2000-03-01 00:00:00", "2000-02-30 00:00:00"),
                48,
            ),
            (
                "standard",
                ("2000-02-28 00:00:00", "2000-03-01 00:00:00", "2000-02-29 00:00:00"),
                48,
            ),
            (
                "noleap",
                ("2000-02-28 00:00:00", "2000-03-01 00:00:00", "2000-02-28 12:00:00"),
                24,
            ),
        ],
    )
    def test_dates(self, observation, calendar, dates, hours):
        """Unsorted dates give a calendar-aware midpoint and elapsed duration."""
        start, end, median = dates
        records = [{**observation, "date": date} for date in (end, start)]
        summary = statistics.summarise_observations(records, calendar)
        assert summary["date_start"] == start
        assert summary["date_end"] == end
        assert summary["date_median"] == median
        assert summary["duration_hours"] == hours

        single = statistics.summarise_observations([records[0]], calendar)
        assert single["date_start"] == end
        assert single["date_end"] == end
        assert single["date_median"] == end
        assert single["duration_hours"] == 0


class TestDatabaseStatistics:
    """Test saving, updating and rolling back statistics."""

    @pytest.mark.parametrize("include_environment", [False, True])
    def test_saved_values(self, tmp_db, include_environment):
        """Save a known fixture track with rounded values and optional environment."""
        build_db.build(
            tmp_db,
            NETCDF_360_DAY_FILE,
            hadsst_filepath=HADSST_2000_FILE if include_environment else None,
            roni_filepath=RONI_FILE if include_environment else None,
        )
        with open_db(tmp_db) as db:
            db.row_factory = sqlite3.Row
            row = db.execute(
                "select * from statistics where trajectory_id = 1"
            ).fetchone()

        assert dict(row) == {
            "trajectory_id": 1,
            "date_start": "2000-01-01 00:00:00",
            "date_end": "2000-01-05 15:00:00",
            "date_median": "2000-01-03 07:30:00",
            "duration_hours": 111.0,
            "observation_count": 38,
            "wind_speed_min": 11.77,
            "wind_speed_max": 24.87,
            "wind_speed_median": 20.93,
            "wind_speed_mean": 19.94,
            "air_pressure_at_sea_level_min": 98669.44,
            "ace": 5.8736,
            "sea_surface_temp_anomaly_mean": 0.6 if include_environment else None,
            "roni_mean": -1.59 if include_environment else None,
        }

    def test_updates(self, tmp_db):
        """Regeneration updates existing statistics without adding duplicate rows."""
        build_db.build(tmp_db, NETCDF_FILE)
        with open_db(tmp_db) as db:
            db.execute("update observations set wind_speed = 5 where trajectory_id = 1")
            db.commit()
            statistics.generate_statistics(db)

            row = db.execute(
                "select wind_speed_mean, ace from statistics where trajectory_id = 1"
            ).fetchone()
            assert row == (5.0, 0.1039)
            assert db.execute("select count(*) from statistics").fetchone()[0] == 6

    def test_empty_tracks(self, tmp_db):
        """A trajectory without observations gets no statistics row."""
        build_db.build(tmp_db, NETCDF_FILE)
        with open_db(tmp_db) as db:
            cursor = db.execute("insert into trajectories (file_id) values (1)")
            trajectory_id = cursor.lastrowid
            db.commit()
            statistics.generate_statistics(db)

            assert (
                db.execute(
                    "select * from statistics where trajectory_id = ?", (trajectory_id,)
                ).fetchone()
                is None
            )

    def test_deleting_tracks(self, tmp_db):
        """Deleting a trajectory also deletes its statistics."""
        build_db.build(tmp_db, NETCDF_FILE)
        with open_db(tmp_db) as db:
            db.execute("delete from trajectories where id = 1")
            assert (
                db.execute(
                    "select * from statistics where trajectory_id = 1"
                ).fetchone()
                is None
            )

    def test_rollback(self, tmp_db):
        """A failure on a later track leaves all previous statistics unchanged."""
        build_db.build(tmp_db, NETCDF_FILE)
        with open_db(tmp_db) as db:
            before = db.execute(
                "select * from statistics order by trajectory_id"
            ).fetchall()
            db.execute("update observations set wind_speed = 1 where trajectory_id = 1")
            db.execute(
                "update observations set date = 'invalid-date' "
                "where trajectory_id = 6 and sequence = 0"
            )
            db.commit()

            with pytest.raises(ValueError, match="trajectory 6"):
                statistics.generate_statistics(db)

            after = db.execute(
                "select * from statistics order by trajectory_id"
            ).fetchall()
            assert after == before
