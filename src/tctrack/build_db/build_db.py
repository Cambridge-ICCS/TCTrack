"""CLI tool to import NetCDF trajectory files into SQLite databases."""

import argparse
import json
import logging
import os
import sqlite3
import sys
from pathlib import Path

import netCDF4
import numpy as np

logger = logging.getLogger(__name__)

SCHEMA_PATH = Path(__file__).parent / "schema.sql"

def init_db(db_path: str) -> sqlite3.Connection:
    """
    Create a new database from the schema, or open an existing one.

    Parameters
    ----------
    db_path
        Path and name of database to open.

    Returns
    -------
    An open connection with foreign keys enabled.
    """
    db = sqlite3.connect(db_path)
    db.execute("pragma foreign_keys = on")

    existing_tables = {
        row[0]
        for row in db.execute(
            "select name from sqlite_master where type='table'"
        ).fetchall()
    }
    if not {"collections", "files", "trajectories", "observations"} <= existing_tables:
        logger.info("Creating new database at %s", db_path)
        schema = SCHEMA_PATH.read_text(encoding="utf-8")
        db.executescript(schema)
        db.commit()

    return db


def get_collection(db: sqlite3.Connection, name: str) -> int:
    """
    Get or create a collection.

    Parameters
    ----------
    db
        Database connection.
    name
        Collection name.

    Returns
    -------
    Row id of collection in the database.
    """
    row = db.execute("select id from collections where name = ?", (name,)).fetchone()
    if row is not None:
        return row[0]

    cur = db.execute("insert into collections (name) values (?)", (name,))
    if cur.lastrowid is None:
        msg = "Insert into collections table failed"
        raise RuntimeError(msg)

    return cur.lastrowid


def file_exists(db: sqlite3.Connection, filepath: str) -> bool:
    """
    Check if a file record already exists.

    Parameters
    ----------
    db
        Database connection.
    filepath
        Full path and name of tracker file.

    Returns
    -------
    True if file exists in the database.
    """
    return db.execute(
        "select 1 from files where filepath = ?", (filepath,)
    ).fetchone() is not None


def read_netcdf(netcdf_filepath: str, wrap_longitude: bool = False) -> dict:
    """
    Read relevant data from a NetCDF trajectory file.

    Parameters
    ----------
    netcdf_filepath
        Full path to a NetCDF file.
    wrap_longitude
        True to wrap longitude coordinates from 0-360 to the -180 to 180 range.

    Returns
    -------
    Dictionary with dimensions, global attributes and variable arrays.
    """
    ds = netCDF4.Dataset(netcdf_filepath, "r")

    n_trajectories = len(ds.dimensions["trajectory"])
    n_observations = len(ds.dimensions["observation"])

    data = {
        "filepath": netcdf_filepath,
        "n_trajectories": n_trajectories,
        "n_observations": n_observations,
        "attrs": {k: ds.getncattr(k) for k in ds.ncattrs()},
        "latitude": ds.variables["latitude"][:],
        "longitude": ds.variables["longitude"][:],
        "grid_i": ds.variables["grid_i"][:],
        "grid_j": ds.variables["grid_j"][:],
        "air_pressure_at_sea_level": ds.variables["air_pressure_at_sea_level"][:],
        "surface_altitude": ds.variables["surface_altitude"][:],
        "wind_speed": ds.variables["wind_speed"][:],
    }

    time_var = ds.variables["time"]
    data["time_units"] = time_var.getncattr("units")
    data["time_calendar"] = time_var.getncattr("calendar")
    data["time"] = ds.variables["time"][:]

    # Read start/end flags if present (per-trajectory booleans)
    if "start_flag" in ds.variables:
        data["start_flag"] = ds.variables["start_flag"][:]
    if "end_flag" in ds.variables:
        data["end_flag"] = ds.variables["end_flag"][:]

    if wrap_longitude:
        data["longitude"] = ((data["longitude"] + 180) % 360) - 180

    ds.close()
    return data


def extract_trajectory(netcdf_data: dict, traj_idx: int) -> dict:
    """
    Extract and precompute all fields for a single trajectory.

    Observations padded with NaN are skipped.

    Parameters
    ----------
    netcdf_data
        Dictionary from `read_netcdf`.
    traj_idx
        Index of trajectory to extract.

    Returns
    -------
    Dictionary of attributes for the specified trajectory.
    """
    # Get indices of valid (non-NaN) observations for the trajectory
    indices = list(np.where(~np.isnan(netcdf_data["latitude"][traj_idx]))[0])

    times = netCDF4.num2date(
        netcdf_data["time"][traj_idx],
        units=netcdf_data["time_units"],
        calendar=netcdf_data["time_calendar"],
    )

    # Compute the trajectory start/end indicator (S, E, SE, None)
    start = "start_flag" in netcdf_data and bool(netcdf_data["start_flag"][traj_idx])
    end = "end_flag" in netcdf_data and bool(netcdf_data["end_flag"][traj_idx])
    start_end = ("S" if start else "") + ("E" if end else "") or None

    return {
        "filepath": netcdf_data["filepath"],
        "start_end": start_end,
        "indices": indices,
        "times": times,
        "latitude": netcdf_data["latitude"][traj_idx],
        "longitude": netcdf_data["longitude"][traj_idx],
        "grid_i": netcdf_data["grid_i"][traj_idx],
        "grid_j": netcdf_data["grid_j"][traj_idx],
        "air_pressure_at_sea_level": netcdf_data["air_pressure_at_sea_level"][traj_idx],
        "surface_altitude": netcdf_data["surface_altitude"][traj_idx],
        "wind_speed": netcdf_data["wind_speed"][traj_idx],
    }


def build_geojson(traj: dict) -> str:
    """
    Build a GeoJSON FeatureCollection for a single trajectory.

    Contains a LineString feature for the full path and Point features
    for each observation along with attributes.

    Parameters
    ----------
    traj
        Dictionary from `extract_trajectory`.

    Returns
    -------
    Trajectory vector line and points in a GeoJSON string.
    """
    features = []

    # 1. Vector track (LineString)
    coordinates = []
    for i in traj["indices"]:
        coordinates.append([float(traj["longitude"][i]), float(traj["latitude"][i])])

    features.append({
        "type": "Feature",
        "geometry": {
            "type": "LineString",
            "coordinates": coordinates,
        },
        "properties": {
            "feature_type": "track",
            "start_end": traj["start_end"],
            "source_file": os.path.basename(traj["filepath"]),
        },
    })

    # 2. Point features for each observation (Point)
    for idx, i in enumerate(traj["indices"]):
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [
                    float(traj["longitude"][i]), float(traj["latitude"][i])],
            },
            "properties": {
                "feature_type": "observation",
                "index": idx,
                "date": traj["times"][i].isoformat(" "),
                "air_pressure_at_sea_level":
                    float(traj["air_pressure_at_sea_level"][i]),
                "surface_altitude": float(traj["surface_altitude"][i]),
                "wind_speed": float(traj["wind_speed"][i]),
            },
        })

    return json.dumps({"type": "FeatureCollection", "features": features})


def insert_file(
    db: sqlite3.Connection, collection_id: int, netcdf_data: dict,
) -> int:
    """
    Insert a file record into the database.

    Parameters
    ----------
    db
        Database connection.
    collection_id
        ID of collection row in the database.
    netcdf_data
        Dictionary from `read_netcdf`.

    Returns
    -------
    Row id of inserted file in the database.
    """
    attrs = netcdf_data["attrs"]

    tracker_name = attrs.get("tctrack_tracker", "unknown")

    # Combine detect and stitch parameters
    raw_params: dict = {}
    if "detect_parameters" in attrs:
        raw_params["detect_parameters"] = attrs["detect_parameters"]
    if "stitch_parameters" in attrs:
        raw_params["stitch_parameters"] = attrs["stitch_parameters"]
    tracker_parameters = json.dumps(raw_params) if raw_params else "{}"

    tctrack_version = attrs.get("tctrack_version", None)

    cur = db.execute(
        """insert into files
           (collection_id, filepath, filename, tctrack_version, tracker_name,
            tracker_parameters, trajectories, observations,
            time_units, time_calendar)
           values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            collection_id,
            os.path.abspath(netcdf_data["filepath"]),
            os.path.basename(netcdf_data["filepath"]),
            tctrack_version,
            tracker_name,
            tracker_parameters,
            netcdf_data["n_trajectories"],
            netcdf_data["n_observations"],
            netcdf_data["time_units"],
            netcdf_data["time_calendar"],
        ),
    )
    if cur.lastrowid is None:
        msg = "Insert into files table failed"
        raise RuntimeError(msg)

    return cur.lastrowid


def insert_trajectory(db: sqlite3.Connection, file_id: int, traj: dict) -> int:
    """
    Insert a trajectory and its observations.

    Parameters
    ----------
    db
        Database connection.
    file_id
        ID of file row in the database.
    traj
        Dictionary from `extract_trajectory`.

    Returns
    -------
    Row id of trajectory in the database.
    """
    geojson_str = build_geojson(traj)

    cur = db.execute(
        "insert into trajectories (file_id, start_end, geojson) values (?, ?, ?)",
        (file_id, traj["start_end"], geojson_str),
    )
    if cur.lastrowid is None:
        msg = "Insert into trajectories table failed"
        raise RuntimeError(msg)

    traj_id = cur.lastrowid

    rows = [
        (
            traj_id,
            traj["times"][i].isoformat(" "),
            float(traj["latitude"][i]),
            float(traj["longitude"][i]),
            float(traj["grid_i"][i]),
            float(traj["grid_j"][i]),
            float(traj["air_pressure_at_sea_level"][i]),
            float(traj["surface_altitude"][i]),
            float(traj["wind_speed"][i]),
        )
        for i in traj["indices"]
    ]

    if rows:
        cur = db.executemany(
            "insert into observations"
            " (trajectory_id, date, latitude, longitude, grid_i, grid_j,"
            "  air_pressure_at_sea_level, surface_altitude, wind_speed)"
            " values (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        if cur.rowcount != len(rows):
            msg = "Inserts into observations table failed"
            raise RuntimeError(msg)

    return traj_id


def import_file(
    db: sqlite3.Connection,
    collection_id: int,
    netcdf_filepath: str,
    wrap_longitude: bool = False,
) -> int:
    """
    Import a single NetCDF file into the database.

    Files already imported are skipped.

    Parameters
    ----------
    db
        Database connection.
    collection_id
        Row id of collection to use in database.
    netcdf_filepath
        Full path to a NetCDF file to import.
    wrap_longitude
        True to wrap longitude coordinates from 0-360 to the -180 to 180 range.

    Returns
    -------
    Row id of inserted file in the database.
    """
    filepath = os.path.abspath(netcdf_filepath)

    if file_exists(db, filepath):
        logger.info("File %s already imported, skipping.", filepath)
        return

    filename = os.path.basename(netcdf_filepath)
    logger.info("Importing %s", filename)
    netcdf_data = read_netcdf(netcdf_filepath, wrap_longitude)

    with db:
        file_id = insert_file(db, collection_id, netcdf_data)

        for traj_idx in range(netcdf_data["n_trajectories"]):
            traj = extract_trajectory(netcdf_data, traj_idx)
            insert_trajectory(db, file_id, traj)

    logger.info(
        "Imported %d trajectories with %d observations each from %s",
        netcdf_data["n_trajectories"],
        netcdf_data["n_observations"],
        filename,
    )


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Import NetCDF trajectory files into SQLite databases.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="SQLite database filename to open or create",
    )
    parser.add_argument(
        "--collection",
        default="default",
        metavar="NAME",
        help="Collection name to use or create",
    )
    parser.add_argument(
        "--wrap",
        action="store_true",
        help="Wrap longitude coordinates from 0 to 360 to the -180 to 180 range",
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="NetCDF files to import",
    )

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s", force=True
    )

    with init_db(args.output) as db:

        collection_id = get_collection(db, args.collection)

        for nc_file in args.files:
            import_file(db, collection_id, nc_file, args.wrap)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
