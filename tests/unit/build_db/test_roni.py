"""Tests for monthly Relative Oceanic Niño Index lookup."""

from datetime import datetime
from pathlib import Path

import pytest

from tctrack.build_db.roni import RONI

RONI_FILE = Path(__file__).parent / "roni.nc"
INVALID_NC_FILE = Path(__file__).parent / "HadSST.4.2.0.0_median_2000.nc"


class TestRONIValidation:
    """Test RONI source validation."""

    def test_missing_path_raises_value_error(self):
        """Test a missing RONI dataset is rejected."""
        with pytest.raises(ValueError, match="RONI dataset does not exist"):
            RONI("missing-roni.nc")

    def test_malformed_dataset_raises_value_error(self, tmp_path):
        """Test a file that is not NetCDF is rejected."""
        filepath = tmp_path / "malformed.nc"
        filepath.write_text("not a NetCDF file", encoding="utf-8")

        with pytest.raises(ValueError, match="Could not open RONI dataset"):
            RONI(filepath)

    def test_missing_value_variable_raises_value_error(self):
        """Test a NetCDF dataset without the RONI value variable is rejected."""
        with pytest.raises(ValueError, match="missing required variable 'value'"):
            RONI(INVALID_NC_FILE)


class TestRONILookup:
    """Test monthly lookup against the supplied RONI fixture."""

    def test_returns_values_for_known_fixture_dates(self):
        """Test exact monthly source dates return their matching RONI values."""
        roni = RONI(RONI_FILE)

        assert roni.lookup(datetime(2000, 1, 1)) == -1.59
        assert roni.lookup(datetime(2000, 2, 1)) == -1.31

    def test_uses_value_from_observation_month(self):
        """Test late-month observations do not take next month's value."""
        roni = RONI(RONI_FILE)

        assert roni.lookup(datetime(2000, 1, 16)) == -1.59
        assert roni.lookup(datetime(2000, 1, 31, 23)) == -1.59

    def test_month_outside_source_range_returns_none(self):
        """Test a month without source data does not use a distant value."""
        roni = RONI(RONI_FILE)

        assert roni.lookup(datetime(2100, 1, 1)) is None
        assert roni.lookup(datetime(1948, 12, 31)) is None

    def test_preserves_missing_source_value(self):
        """Test a masked RONI source value is returned as null."""
        roni = RONI(RONI_FILE)

        assert roni.lookup(datetime(1949, 1, 1)) is None
        assert roni.lookup(datetime(1949, 1, 31)) is None
