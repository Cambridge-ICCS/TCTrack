"""Calculate and persist trajectory summaries from stored observations."""

import re
import sqlite3
from datetime import timedelta
from itertools import groupby

import cftime
import numpy as np


def _summarise_dates(date_strings: list[str], time_calendar: str) -> dict:
    """Summarise stored ISO dates using their original NetCDF calendar."""
    if not isinstance(time_calendar, str) or not time_calendar:
        msg = "time_calendar must be a nonempty calendar name"
        raise ValueError(msg)

    dates = []
    for value in date_strings:
        match = re.fullmatch(
            r"(-?\d+)-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})(?:\.(\d{1,6}))?",
            value if isinstance(value, str) else "",
        )
        if match is None:
            msg = f"invalid observation date {value!r}"
            raise ValueError(msg)

        parts = [int(part) for part in match.groups()[:6]]
        microsecond = int((match.group(7) or "0").ljust(6, "0"))
        try:
            date = cftime.datetime(
                *parts, microsecond=microsecond, calendar=time_calendar
            )
        except (TypeError, ValueError) as exc:
            msg = f"invalid observation date {value!r} for calendar {time_calendar!r}"
            raise ValueError(msg) from exc
        dates.append(date)

    dates.sort()
    seconds = np.asarray([(date - dates[0]).total_seconds() for date in dates])
    median = dates[0] + timedelta(seconds=float(np.median(seconds)))

    return {
        "date_start": dates[0].isoformat(" "),
        "date_end": dates[-1].isoformat(" "),
        "date_median": median.isoformat(" "),
        "duration_hours": float((seconds[-1] - seconds[0]) / 3600),
    }


def summarise_observations(observations: list, time_calendar: str) -> dict:
    """
    Calculate trajectory statistics from a nonempty group of stored observations.

    Parameters
    ----------
    observations
        Mapping records with date, wind_speed, air_pressure_at_sea_level,
        sea_surface_temp_anomaly and roni fields. NULL and nonfinite measurements
        are excluded independently for each field.
    time_calendar
        Original NetCDF calendar from the trajectory's file record.

    Returns
    -------
    Dictionary of scalar statistics, with None for unavailable measurements.
    Environmental means are observation-weighted. ACE includes every valid wind
    using the researchers' approved factor of 1.944, without cadence filtering.

    Raises
    ------
    ValueError
        If observations are empty, required fields are missing or dates/calendar
        are invalid.
    """
    if not observations:
        msg = "observations must be a nonempty group of records"
        raise ValueError(msg)

    records = [dict(record) for record in observations]
    fields = (
        "wind_speed",
        "air_pressure_at_sea_level",
        "sea_surface_temp_anomaly",
        "roni",
    )
    for record in records:
        for field in ("date", *fields):
            if field not in record:
                msg = f"missing {field!r} in observation record"
                raise ValueError(msg)

    summary = _summarise_dates([record["date"] for record in records], time_calendar)
    summary["observation_count"] = len(records)
    values = {}
    for field in fields:
        array = np.asarray([record[field] for record in records], dtype=float)
        values[field] = array[np.isfinite(array)]

    wind = values["wind_speed"]
    pressure = values["air_pressure_at_sea_level"]
    ssta = values["sea_surface_temp_anomaly"]
    roni = values["roni"]
    summary.update(
        wind_speed_min=float(np.min(wind)) if wind.size else None,
        wind_speed_max=float(np.max(wind)) if wind.size else None,
        wind_speed_median=float(np.median(wind)) if wind.size else None,
        wind_speed_mean=float(np.mean(wind)) if wind.size else None,
        air_pressure_at_sea_level_min=float(np.min(pressure))
        if pressure.size
        else None,
        ace=float(np.sum(1e-4 * (wind * 1.944) ** 2)) if wind.size else None,
        sea_surface_temp_anomaly_mean=float(np.mean(ssta)) if ssta.size else None,
        roni_mean=float(np.mean(roni)) if roni.size else None,
    )
    return summary


def generate_statistics(db: sqlite3.Connection) -> None:
    """
    Recalculate and upsert statistics from stored observations in one transaction.

    Floating-point summaries are rounded before insertion: four decimal places
    for ACE and two for all other values. Dates, counts and None are unchanged.

    Parameters
    ----------
    db
        Database connection with completed imports committed. Observation groups
        are streamed by trajectory ID; trajectories without observations are ignored.

    Raises
    ------
    ValueError
        If a trajectory's observations cannot be summarised. The entire statistics
        pass is rolled back, leaving completed imports and previous statistics intact.
    """
    read_cursor = db.cursor()
    read_cursor.row_factory = sqlite3.Row
    try:
        with db:
            read_cursor.execute(
                """select ob.trajectory_id, ob.date, ob.wind_speed,
                          ob.air_pressure_at_sea_level, ob.sea_surface_temp_anomaly,
                          ob.roni, f.time_calendar
                   from observations ob
                   join trajectories tr on tr.id = ob.trajectory_id
                   join files f on f.id = tr.file_id
                   order by ob.trajectory_id, ob.sequence"""
            )
            for trajectory_id, group in groupby(
                read_cursor, key=lambda row: row["trajectory_id"]
            ):
                records = list(group)
                try:
                    summary = summarise_observations(
                        records, records[0]["time_calendar"]
                    )
                except (TypeError, ValueError) as exc:
                    msg = f"cannot summarise trajectory {trajectory_id}: {exc}"
                    raise ValueError(msg) from exc

                for column, value in summary.items():
                    if isinstance(value, float):
                        summary[column] = round(value, 4 if column == "ace" else 2)

                db.execute(
                    """insert into statistics
                       (trajectory_id, date_start, date_end, date_median,
                        duration_hours, observation_count, wind_speed_min,
                        wind_speed_max, wind_speed_median, wind_speed_mean,
                        air_pressure_at_sea_level_min, ace,
                        sea_surface_temp_anomaly_mean, roni_mean)
                       values (:trajectory_id, :date_start, :date_end, :date_median,
                               :duration_hours, :observation_count, :wind_speed_min,
                               :wind_speed_max, :wind_speed_median, :wind_speed_mean,
                               :air_pressure_at_sea_level_min, :ace,
                               :sea_surface_temp_anomaly_mean, :roni_mean)
                       on conflict(trajectory_id) do update set
                           date_start = excluded.date_start,
                           date_end = excluded.date_end,
                           date_median = excluded.date_median,
                           duration_hours = excluded.duration_hours,
                           observation_count = excluded.observation_count,
                           wind_speed_min = excluded.wind_speed_min,
                           wind_speed_max = excluded.wind_speed_max,
                           wind_speed_median = excluded.wind_speed_median,
                           wind_speed_mean = excluded.wind_speed_mean,
                           air_pressure_at_sea_level_min =
                               excluded.air_pressure_at_sea_level_min,
                           ace = excluded.ace,
                           sea_surface_temp_anomaly_mean =
                               excluded.sea_surface_temp_anomaly_mean,
                           roni_mean = excluded.roni_mean""",
                    {"trajectory_id": trajectory_id, **summary},
                )
    finally:
        read_cursor.close()
