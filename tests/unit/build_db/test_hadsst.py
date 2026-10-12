"""Tests for monthly, nearest-grid-cell HadSST lookup."""

from datetime import datetime
from pathlib import Path

import cftime
import netCDF4
import numpy as np
import pytest

from tctrack.build_db.hadsst import HadSST

HADSST_FILE = Path(__file__).parent / "HadSST.4.2.0.0_median_2000.nc"
INVALID_NC_FILE = Path(__file__).parent / "roni.nc"


def test_invalid_datasets(tmp_path):
    """Reject missing paths, directories, unreadable files and missing variables."""
    with pytest.raises(ValueError, match="does not exist or is not a file"):
        HadSST(tmp_path / "missing.nc")

    with pytest.raises(ValueError, match="does not exist or is not a file"):
        HadSST(tmp_path)

    invalid_file = tmp_path / "invalid.nc"
    invalid_file.write_text("not a NetCDF dataset")
    with pytest.raises(ValueError, match="Could not open HadSST dataset"):
        HadSST(invalid_file)

    with pytest.raises(ValueError, match="missing required variable 'latitude'"):
        HadSST(INVALID_NC_FILE)


def test_monthly_grid_lookup():
    """Match source cells, missing values, calendar months and wrapped longitudes."""
    hadsst = HadSST(HADSST_FILE)

    with netCDF4.Dataset(HADSST_FILE) as dataset:
        dates = netCDF4.num2date(
            dataset.variables["time"][:], dataset.variables["time"].units
        )
        latitudes = dataset.variables["latitude"][:]
        longitudes = dataset.variables["longitude"][:]
        values = dataset.variables["tos"][:]

    for time_index, date in enumerate(dates):
        for latitude_index in (0, 14, 35):
            for longitude_index in (0, 35, 64, 71):
                source = values[time_index, latitude_index, longitude_index]
                expected = None
                if not np.ma.is_masked(source) and np.isfinite(source):
                    expected = round(float(source), 3)

                latitude = float(latitudes[latitude_index])
                longitude = float(longitudes[longitude_index])
                for offset in (-360, 0, 360):
                    assert hadsst.lookup(date, latitude, longitude + offset) == expected

    february = datetime(2000, 2, 1)
    february_30 = cftime.Datetime360Day(2000, 2, 30)
    assert hadsst.lookup(february_30, -17.5, 142.5) is None
    assert hadsst.lookup(february_30, 2.5, -177.5) == hadsst.lookup(
        february, 2.5, -177.5
    )

    # Outside the dataset's range, use the nearest available calendar month.
    assert hadsst.lookup(datetime(1999, 12, 1), 2.5, -177.5) == hadsst.lookup(
        dates[0], 2.5, -177.5
    )
    assert hadsst.lookup(datetime(2001, 1, 1), 2.5, -177.5) == hadsst.lookup(
        dates[-1], 2.5, -177.5
    )

    # Points across the dateline can be nearest to the first longitude cell.
    assert hadsst.lookup(february, 3.0, 181.0) == hadsst.lookup(february, 2.5, -177.5)
