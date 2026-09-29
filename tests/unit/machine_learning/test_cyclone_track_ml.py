"""Unit tests for cyclone_track_ml.py of the TCTrack machine_learning package."""

# The tests inspect module-private helpers on purpose.
# ruff: noqa: SLF001

import pytest

from tctrack.core.ml_tracker import TCMLParameters
from tctrack.machine_learning import (
    MLParameters,
    MLStitchParameters,
)
from tctrack.machine_learning.cyclone_track_ml import (
    _angular_distance_deg,
    _point_variables,
)


class TestAngularDistanceDeg:
    """Tests for the _angular_distance_deg helper."""

    def test_identical_points(self):
        """Test that the distance between identical points is zero."""
        assert _angular_distance_deg(10.0, 50.0, 10.0, 50.0) == 0.0

    def test_symmetric(self):
        """Test that swapping the two points does not change the distance."""
        forward = _angular_distance_deg(10.0, 50.0, 15.0, 55.0)
        backward = _angular_distance_deg(15.0, 55.0, 10.0, 50.0)
        assert forward == pytest.approx(backward)

    def test_pure_latitude_difference(self):
        """Test that a latitude-only separation is returned unscaled."""
        assert _angular_distance_deg(10.0, 50.0, 13.0, 50.0) == pytest.approx(3.0)

    def test_longitude_difference_at_equator(self):
        """Test that a longitude separation at the equator is unscaled."""
        assert _angular_distance_deg(0.0, 50.0, 0.0, 53.0) == pytest.approx(3.0)

    def test_longitude_scaled_by_latitude(self):
        """Test that longitude differences shrink away from the equator."""
        # cos(60 degrees) = 0.5
        assert _angular_distance_deg(60.0, 50.0, 60.0, 54.0) == pytest.approx(2.0)

    @pytest.mark.parametrize(
        "lon1,lon2",
        [(179.0, -179.0), (-179.0, 179.0), (359.0, 1.0), (1.0, 359.0)],
    )
    def test_wraps_across_antimeridian(self, lon1, lon2):
        """Test that points either side of the date line are close together."""
        assert _angular_distance_deg(0.0, lon1, 0.0, lon2) == pytest.approx(2.0)

    def test_equivalent_longitudes(self):
        """Test that longitudes given on different conventions agree."""
        assert _angular_distance_deg(5.0, -170.0, 5.0, 190.0) == pytest.approx(0.0)

    def test_returns_float(self):
        """Test that the distance is a plain Python float."""
        assert isinstance(_angular_distance_deg(1.0, 2.0, 3.0, 4.0), float)


class TestPointVariables:
    """Tests for the _point_variables helper."""

    def test_removes_time(self):
        """Test that the time key is dropped from a candidate."""
        candidate = {"time": "2025-01-01", "lat": 1.0, "lon": 2.0}
        assert "time" not in _point_variables(candidate)

    def test_keeps_other_keys(self):
        """Test that every non-time key is kept, including extra channels."""
        candidate = {
            "time": "2025-01-01",
            "lat": 1.0,
            "lon": 2.0,
            "class_index": 3.0,
            "score": 0.9,
            "sea_surface_temperature": 300.0,
        }
        assert _point_variables(candidate) == {
            "lat": 1.0,
            "lon": 2.0,
            "class_index": 3.0,
            "score": 0.9,
            "sea_surface_temperature": 300.0,
        }

    def test_does_not_modify_candidate(self):
        """Test that the original candidate is left unchanged."""
        candidate = {"time": "2025-01-01", "lat": 1.0, "lon": 2.0}
        _point_variables(candidate)
        assert candidate == {"time": "2025-01-01", "lat": 1.0, "lon": 2.0}


class TestMLParameters:
    """Tests for the MLParameters dataclass."""

    def test_subclass_of_base_parameters(self):
        """Test that MLParameters extends the shared ML parameters."""
        assert issubclass(MLParameters, TCMLParameters)

    def test_defaults(self):
        """Test the default values of MLParameters."""
        params = MLParameters()
        assert params.input_file == ""
        assert params.hf_repo_id == "surbhigoel456/cyclone-TC-ML"
        assert params.normalisation_stats_path is None
        assert params.merge_distance_deg == 2.0
        assert params.pressure_levels == (1000, 750, 500)

    def test_overrides(self):
        """Test that values passed at construction override the defaults."""
        params = MLParameters(input_file="in.nc", merge_distance_deg=0.0)
        assert params.input_file == "in.nc"
        assert params.merge_distance_deg == 0.0


class TestMLStitchParameters:
    """Tests for the MLStitchParameters dataclass."""

    def test_defaults(self):
        """Test the default values of MLStitchParameters."""
        params = MLStitchParameters()
        assert params.max_distance_deg == 3.0
        assert params.max_gap == 1
        assert params.min_length == 2

    def test_overrides(self):
        """Test that values passed at construction override the defaults."""
        params = MLStitchParameters(max_distance_deg=5.0, max_gap=0, min_length=1)
        assert (params.max_distance_deg, params.max_gap, params.min_length) == (
            5.0,
            0,
            1,
        )
