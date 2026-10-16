"""Load landfall gate rectangles and match already classified observations."""

import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import TypeAlias

import numpy as np

_COLUMNS = ("code", "name", "min_lat", "max_lat", "min_lon", "max_lon")
_LAT_LIMIT = 90
_LON_LIMIT = 180


@dataclass(frozen=True)
class Gate:
    """A named WGS84 rectangle with inclusive minima and exclusive maxima."""

    code: str
    name: str
    min_lat: float
    max_lat: float
    min_lon: float
    max_lon: float

    def __post_init__(self):
        """Validate a gate's identifiers and bounds once at construction."""
        if not self.code.strip() or not self.name.strip():
            msg = "code and name must be nonempty"
            raise ValueError(msg)

        bounds = (self.min_lat, self.max_lat, self.min_lon, self.max_lon)
        if not all(isfinite(value) for value in bounds):
            msg = "bounds must be finite numbers"
            raise ValueError(msg)

        if not -_LAT_LIMIT <= self.min_lat < self.max_lat <= _LAT_LIMIT:
            msg = "latitude bounds must satisfy -90 <= min_lat < max_lat <= 90"
            raise ValueError(msg)

        if not -_LON_LIMIT <= self.min_lon < self.max_lon <= _LON_LIMIT:
            msg = "longitude bounds must satisfy -180 <= min_lon < max_lon <= 180"
            raise ValueError(msg)


Gates: TypeAlias = tuple[Gate, ...]


def load_gates(path: str | Path) -> Gates:
    """Read a CSV into an immutable tuple of validated, nonoverlapping gates.

    Parameters
    ----------
    path
        CSV file with code, name, min_lat, max_lat, min_lon and max_lon columns.

    Returns
    -------
    tuple of Gate
        Definitions in CSV order, with surrounding code/name whitespace removed.

    Raises
    ------
    OSError
        If the file cannot be opened or read; includes the filename.
    ValueError
        If columns, rows, identifiers, bounds or overlaps are invalid. Errors
        include the filename and, for data rows, the CSV line number.
    """
    gates: list[Gate] = []
    codes: set[str] = set()

    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, strict=True)
        try:
            if reader.fieldnames is None or sorted(reader.fieldnames) != sorted(
                _COLUMNS
            ):
                msg = f"{path}: row 1: expected columns {', '.join(_COLUMNS)}"
                raise ValueError(msg)

            for row in reader:
                context = f"{path}: row {reader.line_num}"
                if None in row or any(value is None for value in row.values()):
                    msg = f"{context}: missing or extra row fields"
                    raise ValueError(msg)

                try:
                    gate = Gate(
                        row["code"].strip(),
                        row["name"].strip(),
                        *(float(row[column]) for column in _COLUMNS[2:]),
                    )
                except ValueError as error:
                    msg = f"{context}: {error}"
                    raise ValueError(msg) from error

                if gate.code in codes:
                    msg = f"{context}: duplicate code {gate.code!r}"
                    raise ValueError(msg)

                for previous in gates:
                    if max(gate.min_lat, previous.min_lat) < min(
                        gate.max_lat, previous.max_lat
                    ) and max(gate.min_lon, previous.min_lon) < min(
                        gate.max_lon, previous.max_lon
                    ):
                        msg = (
                            f"{context}: overlapping gates "
                            f"{previous.code!r} and {gate.code!r}"
                        )
                        raise ValueError(msg)

                gates.append(gate)
                codes.add(gate.code)
        except csv.Error as error:
            msg = f"{path}: row {reader.line_num}: {error}"
            raise ValueError(msg) from error

    if not gates:
        msg = f"{path}: expected at least one gate data row"
        raise ValueError(msg)

    return tuple(gates)


def match_gates(
    latitude: np.typing.ArrayLike,
    longitude: np.typing.ArrayLike,
    landfall: np.typing.ArrayLike,
    gates: Gates,
) -> np.ndarray:
    """Match existing landfalls to validated gates without reclassifying points.

    Parameters
    ----------
    latitude, longitude
        Same-shaped coordinates in degrees; masked/nonfinite values are skipped.
        Longitudes are wrapped to [-180, 180) without modifying the input.
    landfall
        Same-shaped boolean results from the existing landfall detector.
        Masked entries are skipped.
    gates
        Validated definitions returned by `load_gates`; not revalidated here.

    Returns
    -------
    numpy.ndarray
        Object array of the input shape containing gate codes or None.

    Raises
    ------
    ValueError
        If coordinates and landfall have different shapes.
    """
    lat = np.ma.asarray(latitude, dtype=float).filled(np.nan)
    lon = np.ma.asarray(longitude, dtype=float).filled(np.nan)
    land = np.ma.asarray(landfall, dtype=bool).filled(False)

    if lat.shape != lon.shape or lat.shape != land.shape:
        msg = "latitude, longitude and landfall must have the same shape"
        raise ValueError(msg)

    result = np.full(lat.shape, None, dtype=object)
    valid = land & np.isfinite(lat) & np.isfinite(lon)
    wrapped = np.full(lon.shape, np.nan)
    wrapped[valid] = np.mod(lon[valid] + 180.0, 360.0) - 180.0

    for gate in gates:
        inside = (
            valid
            & (lat >= gate.min_lat)
            & (lat < gate.max_lat)
            & (wrapped >= gate.min_lon)
            & (wrapped < gate.max_lon)
        )
        result[inside] = gate.code

    return result
