"""Unit tests for cyclone_track_ml.py of the TCTrack machine_learning package.

These tests cover the logic of the machine-learning tracker that turns model output
into tracks.

No test needs a real model or a download. A fixture creates trackers with model loading
disabled, the model itself is replaced by a mock that returns chosen output, and the
inputs are small hand-built grids and candidate lists. The tests of reading the input
file use the small ERA5 sample in ``data/machine_learning``. The tests therefore check
the tracker's logic.

The data used for tests is sample subset of ERA5 data in the directory
/data/machine_learning, and is licensed under the Open Government Licence v3.0 (OGL).
"""

# The tests inspect module-private helpers on purpose.
# ruff: noqa: SLF001

import json
from dataclasses import asdict
from datetime import timedelta
from pathlib import Path
from unittest.mock import MagicMock

import cf
import h5py
import numpy as np
import pytest
import torch
from cftime import datetime
from tctrack.machine_learning.cyclone_track_ml import (
    _angular_distance_deg,
    _point_variables,
)

from tctrack.machine_learning import (
    MLParameters,
    MLStitchParameters,
    MLTracker,
)

SAMPLE_FILE = (
    Path(__file__).parents[3]
    / "data"
    / "machine_learning"
    / "era5_dikeledi_2025-01-10.nc"
)
FAKE_TOKEN = "fake-token"  # noqa: S105 - not a real credential
FAKE_ENV_TOKEN = "fake-env-token"  # noqa: S105 - not a real credential
N_CHANNELS = 17  # 5 pressure variables x 3 levels, land mask, sea surface temperature


@pytest.fixture
def make_tracker(monkeypatch):
    """Create an MLTracker without loading a model."""
    monkeypatch.setattr(MLTracker, "_load_model", lambda *_args, **_kwargs: None)

    def _make(parameters=None, stitch_parameters=None):
        return MLTracker(parameters or MLParameters(), stitch_parameters)

    return _make


def _time(hour: int) -> datetime:
    """Return a 6-hourly timestep on a fixed day."""
    return datetime(2025, 1, 1, calendar="gregorian") + timedelta(hours=6 * hour)


def _candidate(
    hour: int, lat: float, lon: float, score: float = 0.9, data: dict | None = None
) -> dict:
    """Return a candidate detection at the given timestep and location."""
    return {
        "time": _time(hour),
        "lat": lat,
        "lon": lon,
        "class_index": 2.0,
        "score": score,
        "data": data or {},
    }


def _logits(storms=None, n_lat=6, n_lon=6):
    """Model output where every pixel is background except those in ``storms``.

    ``storms`` maps ``(y, x)`` to ``(class_index, logit)``.
    """
    logits = torch.zeros(1, 5, n_lat, n_lon)
    logits[0, 0] = 10.0
    for (y, x), (class_index, value) in (storms or {}).items():
        logits[0, :, y, x] = 0.0
        logits[0, class_index, y, x] = value
    return logits


class TestAngularDistanceDeg:
    """Tests for the _angular_distance_deg helper."""

    def test_pure_latitude_difference(self):
        """Test that a latitude-only separation is returned unscaled."""
        assert _angular_distance_deg(10.0, 50.0, 13.0, 50.0) == pytest.approx(3.0)

    def test_longitude_scaled_by_latitude(self):
        """Test that longitude differences shrink away from the equator."""
        # cos(60 degrees) = 0.5
        assert _angular_distance_deg(60.0, 50.0, 60.0, 54.0) == pytest.approx(2.0)

    @pytest.mark.parametrize(
        "lon1,lon2",
        [(179.0, -179.0), (359.0, 1.0)],
    )
    def test_wraps_across_antimeridian(self, lon1, lon2):
        """Test that points either side of the date line are close together."""
        assert _angular_distance_deg(0.0, lon1, 0.0, lon2) == pytest.approx(2.0)


class TestPointVariables:
    """Tests for the _point_variables helper."""

    def test_keeps_other_keys(self):
        """Test that time and data are dropped and the data values are merged in."""
        candidate = {
            "time": "2025-01-01",
            "lat": 1.0,
            "lon": 2.0,
            "class_index": 3.0,
            "score": 0.9,
            "data": {"sea_surface_temperature": 300.0},
        }
        assert _point_variables(candidate) == {
            "lat": 1.0,
            "lon": 2.0,
            "class_index": 3.0,
            "score": 0.9,
            "sea_surface_temperature": 300.0,
        }


class TestMLParametersThreshold:
    """Tests for the validation of the confidence threshold."""

    @pytest.mark.parametrize("threshold", [-0.1, 1.1])
    def test_invalid_threshold_raises(self, threshold):
        """Thresholds outside 0 to 1 raise a ValueError."""
        with pytest.raises(ValueError, match="threshold must be in"):
            MLParameters(threshold=threshold)


class TestMLTrackerInit:
    """Tests for MLTracker setup."""

    @pytest.fixture
    def fake_hf(self, monkeypatch, tmp_path):
        """Replace the model download and loading with mocks."""
        model_file = tmp_path / "downloaded.pt"
        model_file.touch()
        download = MagicMock(return_value=str(model_file))
        load = MagicMock()
        monkeypatch.setattr("tctrack.core.ml_tracker.hf_hub_download", download)
        monkeypatch.setattr(torch.jit, "load", load)
        monkeypatch.delenv("HF_TOKEN", raising=False)
        return download, load, model_file

    @pytest.mark.usefixtures("fake_hf")
    def test_token_cleared_from_parameters(self):
        """Test that the token is kept privately and removed from parameters."""
        parameters = MLParameters(hf_token=FAKE_TOKEN)
        tracker = MLTracker(parameters)
        assert parameters.hf_token is None
        assert tracker._hf_token == FAKE_TOKEN

    def test_token_passed_to_download(self, fake_hf):
        """Test that the token given in the parameters is used to download."""
        download, _, _ = fake_hf
        MLTracker(MLParameters(hf_token=FAKE_TOKEN))
        assert download.call_args.kwargs["token"] == FAKE_TOKEN

    @pytest.mark.usefixtures("fake_hf")
    def test_token_not_in_serialised_parameters(self):
        """Test that the token cannot end up in the output metadata."""
        tracker = MLTracker(MLParameters(hf_token=FAKE_TOKEN))
        serialised = json.dumps([asdict(p) for p in tracker._parameters])
        assert FAKE_TOKEN not in serialised

    def test_env_token_fallback(self, fake_hf, monkeypatch):
        """Test that fake HF_TOKEN is used when no token is given."""
        download, _, _ = fake_hf
        monkeypatch.setenv("HF_TOKEN", FAKE_ENV_TOKEN)
        MLTracker(MLParameters())
        assert download.call_args.kwargs["token"] == FAKE_ENV_TOKEN

    def test_local_model_path_skips_download(self, fake_hf, tmp_path):
        """Test that a local model file is loaded without downloading."""
        download, load, _ = fake_hf
        model_path = tmp_path / "local.pt"
        model_path.touch()
        MLTracker(MLParameters(model_path=str(model_path), device="cuda"))
        download.assert_not_called()
        load.assert_called_once_with(str(model_path), map_location="cuda")

    @pytest.mark.usefixtures("fake_hf")
    def test_missing_model_path_raises(self, tmp_path):
        """Test that a model_path that does not exist raises an OSError."""
        parameters = MLParameters(model_path=str(tmp_path / "missing.pt"))
        with pytest.raises(OSError, match="Model file not found"):
            MLTracker(parameters)


class TestNormalisationStats:
    """Tests for MLTracker._load_normalisation_stats."""

    def test_bundled_stats_used_by_default(self, make_tracker):
        """The statistics shipped with the package cover every channel."""
        tracker = make_tracker()
        mean, value_range = tracker._load_normalisation_stats()
        assert mean.shape == value_range.shape == (N_CHANNELS,)
        assert np.isfinite(mean).all()
        assert (value_range > 0).all()

    def test_custom_stats_file_used(self, make_tracker, tmp_path):
        """The mean and range are read from normalisation_stats_path when set."""
        stats_file = tmp_path / "stats.nc"
        with h5py.File(stats_file, "w") as stats:
            stats["mean"] = np.arange(N_CHANNELS, dtype=float)
            stats["range"] = np.arange(1, N_CHANNELS + 1, dtype=float)
        tracker = make_tracker(MLParameters(normalisation_stats_path=str(stats_file)))
        mean, value_range = tracker._load_normalisation_stats()
        assert np.array_equal(mean, np.arange(N_CHANNELS))
        assert np.array_equal(value_range, np.arange(1, N_CHANNELS + 1))


class TestSetMetadata:
    """Tests for MLTracker._set_metadata, run on the ERA5 sample file."""

    @pytest.fixture
    def tracker(self, make_tracker):
        """Tracker on the sample file with its metadata set."""
        tracker = make_tracker(MLParameters(input_file=str(SAMPLE_FILE)))
        tracker._set_metadata()
        return tracker

    def test_time_metadata(self, tracker):
        """The calendar, units and first and last times come from the file."""
        metadata = tracker._time_metadata
        assert metadata["calendar"] == "proleptic_gregorian"
        assert metadata["units"] == "seconds since 1970-01-01"
        assert (metadata["start_time"].month, metadata["start_time"].day) == (1, 10)
        assert (metadata["end_time"].day, metadata["end_time"].hour) == (12, 6)

    def test_variable_metadata_matches_channels(self, tracker):
        """Every input channel and both model outputs have metadata."""
        expected = {*tracker._channel_names, "class_index", "score"}
        assert set(tracker._variable_metadata) == expected

    def test_file_without_time_raises(self, make_tracker, tmp_path):
        """An input file with no time coordinate raises a ValueError."""
        field = cf.Field(properties={"standard_name": "air_temperature", "units": "K"})
        lat_axis = field.set_construct(cf.DomainAxis(3))
        lon_axis = field.set_construct(cf.DomainAxis(4))
        for name, axis, size, units in (
            ("latitude", lat_axis, 3, "degrees_north"),
            ("longitude", lon_axis, 4, "degrees_east"),
        ):
            coordinate = cf.DimensionCoordinate(
                data=cf.Data(np.arange(size, dtype=float), units=units),
                properties={"standard_name": name},
            )
            field.set_construct(coordinate, axes=axis)
        field.set_data(cf.Data(np.zeros((3, 4)), units="K"), axes=(lat_axis, lon_axis))
        input_file = tmp_path / "no_time.nc"
        cf.write(field, str(input_file))

        tracker = make_tracker(MLParameters(input_file=str(input_file)))
        with pytest.raises(ValueError, match="time"):
            tracker._set_metadata()


class TestPreprocess:
    """Tests for MLTracker.preprocess, run on the ERA5 sample file."""

    @pytest.fixture
    def tracker(self, make_tracker, monkeypatch):
        """Tracker on the sample file with normalisation switched off."""
        tracker = make_tracker(MLParameters(input_file=str(SAMPLE_FILE)))
        monkeypatch.setattr(
            tracker,
            "_load_normalisation_stats",
            lambda: (np.zeros(N_CHANNELS), np.ones(N_CHANNELS)),
        )
        return tracker

    @pytest.fixture(scope="class")
    def fields(self):
        """Return the fields in the sample file, for comparison with the tensor."""
        return cf.read(str(SAMPLE_FILE))

    def test_tensor_shape(self, tracker):
        """The tensor is (channel, time, lat, lon) float32 with no NaNs."""
        data = tracker.preprocess()
        assert tuple(data.shape) == (N_CHANNELS, 10, 80, 80)
        assert data.dtype == torch.float32
        assert not torch.isnan(data).any()

    def test_grid_and_times_stored(self, tracker):
        """The latitudes, longitudes and times of the file are stored on the tracker."""
        tracker.preprocess()
        assert len(tracker._lats) == 80
        assert len(tracker._lons) == 80
        assert len(tracker._times) == 10
        assert (tracker._lats.min(), tracker._lats.max()) == (-22.25, -2.5)
        assert (tracker._lons.min(), tracker._lons.max()) == (41.0, 60.75)
        assert (tracker._times[0].month, tracker._times[0].day) == (1, 10)
        assert (tracker._times[-1].day, tracker._times[-1].hour) == (12, 6)

    def test_pressure_channels_ordered_variable_then_level(self, tracker, fields):
        """Channels run through each variable's levels in turn: 1000, 750, 500 hPa."""
        data = tracker.preprocess().numpy()
        temperature = fields.select_field("air_temperature")
        for channel, level in zip((3, 4, 5), (1000, 750, 500), strict=True):
            expected = temperature.subspace(Z=level).squeeze("Z").array
            assert np.allclose(data[channel], expected)

    def test_normalisation_applied(self, tracker, monkeypatch):
        """Each channel is normalised as (x - mean) / range."""
        raw = tracker.preprocess()
        mean = np.arange(N_CHANNELS, dtype=float)
        value_range = np.full(N_CHANNELS, 2.0)
        monkeypatch.setattr(
            tracker, "_load_normalisation_stats", lambda: (mean, value_range)
        )
        expected = (raw - torch.from_numpy(mean).float()[:, None, None, None]) / 2.0
        assert torch.allclose(tracker.preprocess(), expected)

    def test_land_mask_channel(self, tracker, fields):
        """The land mask is 1 where sea surface temperature is undefined, all times."""
        data = tracker.preprocess().numpy()
        land_mask = data[15]
        sst_mask = np.ma.getmaskarray(fields.select_field("ncvar%sst").array)
        assert set(np.unique(land_mask)) == {0.0, 1.0}
        assert np.array_equal(land_mask[0], sst_mask[0].astype(np.float32))
        assert all(np.array_equal(land_mask[0], frame) for frame in land_mask)

    def test_sea_surface_temperature_filled_over_land(self, tracker, fields):
        """Sea surface temperature is replaced by 2 m temperature over land."""
        data = tracker.preprocess().numpy()
        sst = fields.select_field("ncvar%sst").array
        t2m = fields.select_field("ncvar%t2m").array
        land = np.ma.getmaskarray(sst)
        assert land.any()
        assert np.allclose(data[16][land], np.asarray(t2m)[land])
        assert np.allclose(data[16][~land], np.ma.getdata(sst)[~land])


class TestDetect:
    """Tests for MLTracker.detect, using a fake model and a synthetic input grid."""

    n_lat, n_lon, n_time, n_channels = 6, 6, 2, 17

    @pytest.fixture
    def tracker(self, make_tracker, monkeypatch):
        """Tracker with a small grid and mocked preprocessing."""
        tracker = make_tracker()
        tracker._lats = np.arange(self.n_lat, dtype=float)
        tracker._lons = np.arange(10, 10 + self.n_lon, dtype=float)
        tracker._times = [_time(i) for i in range(self.n_time)]
        # Channel c holds the constant value c, so the physical value reported for
        # it is c * range + mean and can be predicted.
        data = (
            torch.arange(self.n_channels, dtype=torch.float32)
            .reshape(-1, 1, 1, 1)
            .expand(self.n_channels, self.n_time, self.n_lat, self.n_lon)
            .contiguous()
        )
        monkeypatch.setattr(tracker, "preprocess", lambda: data)
        monkeypatch.setattr(
            tracker,
            "_load_normalisation_stats",
            lambda: (np.full(self.n_channels, 10.0), np.full(self.n_channels, 2.0)),
        )
        return tracker

    def test_storm_pixel_becomes_candidate(self, tracker):
        """Test that a confident storm pixel is reported with its location and class."""
        tracker.model = MagicMock(side_effect=[_logits({(2, 3): (2, 10.0)}), _logits()])
        tracker.detect()
        assert len(tracker._candidates) == 1
        candidate = tracker._candidates[0]
        assert (candidate["lat"], candidate["lon"]) == (2.0, 13.0)
        assert candidate["class_index"] == 2.0
        assert candidate["score"] > 0.99
        assert candidate["time"] == tracker._times[0]
        assert set(candidate) == {"time", "lat", "lon", "class_index", "score", "data"}

    def test_low_confidence_discarded(self, tracker):
        """Test that a storm class below the threshold is not a candidate."""
        # A winning probability of about 0.4, below the default threshold of 0.5.
        tracker.model = MagicMock(return_value=_logits({(2, 3): (2, 1.0)}))
        tracker.detect()
        assert tracker._candidates == []

    def test_scores_hold_every_pixel(self, tracker):
        """Test that the class probabilities are kept for every pixel and timestep."""
        tracker.model = MagicMock(return_value=_logits())
        tracker.detect()
        assert len(tracker._scores) == self.n_time * self.n_lat * self.n_lon
        assert set(tracker._scores[0]) == {"time", "lat", "lon", "probs"}
        assert len(tracker._scores[0]["probs"]) == 5
        assert sum(tracker._scores[0]["probs"]) == pytest.approx(1.0)

    def test_colocated_variables_in_physical_units(self, tracker):
        """Test that input variables at the storm are reported un-normalised."""
        tracker.model = MagicMock(side_effect=[_logits({(2, 3): (2, 10.0)}), _logits()])
        tracker.detect()
        candidate = tracker._candidates[0]
        assert set(candidate["data"]) == set(tracker._channel_names)
        # Last channel (sea surface temperature) is 16 -> 16 * 2 + 10.
        assert candidate["data"]["sea_surface_temperature"] == pytest.approx(42.0)

    def test_grid_mismatch_raises(self, tracker):
        """Test that a model output on a different grid is rejected."""
        tracker.model = MagicMock(return_value=_logits(n_lat=4, n_lon=4))
        with pytest.raises(ValueError, match="does not match"):
            tracker.detect()


class TestClusterCandidates:
    """Tests for MLTracker._cluster_candidates."""

    @staticmethod
    def _grid_tracker(make_tracker):
        """Build a tracker with an example 6x6 grid of lats 0-5 and lons 10-15."""
        tracker = make_tracker()
        tracker._lats = np.arange(6, dtype=float)
        tracker._lons = np.arange(10, 16, dtype=float)
        return tracker

    @staticmethod
    def _arrays(pixels):
        """Build the input arrays from pixel values."""
        is_storm = np.zeros((6, 6), dtype=bool)
        class_idx = np.zeros((6, 6))
        class_prob = np.full((6, 6), 0.1)
        for (y, x), (cls, prob) in pixels.items():
            is_storm[y, x] = True
            class_idx[y, x] = cls
            class_prob[y, x] = prob
        return is_storm, class_idx, class_prob

    def test_separate_blobs(self, make_tracker):
        """Test that non-adjacent pixels give separate candidates."""
        tracker = self._grid_tracker(make_tracker)
        pixels = {(0, 0): (1, 0.7), (5, 5): (2, 0.9)}
        candidates = tracker._cluster_candidates(*self._arrays(pixels), _time(0))
        assert len(candidates) == 2
        by_lat = sorted(candidates, key=lambda c: c["lat"])
        assert (by_lat[0]["lat"], by_lat[0]["lon"]) == (0.0, 10.0)
        assert (by_lat[1]["lat"], by_lat[1]["lon"]) == (5.0, 15.0)
        # Every candidate carries the frame time and an empty data dictionary.
        assert all(c["time"] == _time(0) and c["data"] == {} for c in candidates)

    def test_class_boundary_gives_one_cluster(self, make_tracker):
        """Test that adjacent pixels of different classes form one cluster."""
        tracker = self._grid_tracker(make_tracker)
        pixels = {(2, 2): (1, 0.6), (2, 3): (3, 0.9)}
        candidates = tracker._cluster_candidates(*self._arrays(pixels), _time(0))
        assert len(candidates) == 1

    def test_class_and_score_from_peak_pixel(self, make_tracker):
        """Test that class and score are read off the most confident pixel."""
        tracker = self._grid_tracker(make_tracker)
        pixels = {(2, 2): (1, 0.6), (2, 3): (3, 0.9), (2, 4): (2, 0.7)}
        candidate = tracker._cluster_candidates(*self._arrays(pixels), _time(0))[0]
        assert candidate["class_index"] == 3.0
        assert candidate["score"] == pytest.approx(0.9)

    def test_centroid_weighted_by_confidence(self, make_tracker):
        """Test that the centroid is pulled towards the more confident pixel."""
        tracker = self._grid_tracker(make_tracker)
        pixels = {(2, 2): (1, 0.9), (2, 3): (1, 0.3)}
        candidate = tracker._cluster_candidates(*self._arrays(pixels), _time(0))[0]
        assert candidate["lat"] == pytest.approx(2.0)
        assert candidate["lon"] == pytest.approx((12 * 0.9 + 13 * 0.3) / 1.2)


class TestMergeNearbyCandidates:
    """Tests for MLTracker._merge_nearby_candidates."""

    @staticmethod
    def _cand(lat, lon, score):
        return {"lat": lat, "lon": lon, "class_index": 2.0, "score": score}

    def test_close_candidates_keep_strongest(self, make_tracker):
        """Test that of two close candidates only the higher score remains."""
        candidates = [self._cand(0.0, 50.0, 0.6), self._cand(1.0, 50.0, 0.9)]
        merged = make_tracker()._merge_nearby_candidates(candidates)
        assert merged == [self._cand(1.0, 50.0, 0.9)]

    def test_distant_candidates_all_kept_strongest_first(self, make_tracker):
        """Test that far-apart candidates are kept, ordered by score."""
        candidates = [self._cand(0.0, 50.0, 0.6), self._cand(10.0, 50.0, 0.9)]
        merged = make_tracker()._merge_nearby_candidates(candidates)
        assert [c["score"] for c in merged] == [0.9, 0.6]

    def test_zero_distance_disables_merging(self, make_tracker):
        """Test that merge_distance_deg=0 returns the candidates untouched."""
        tracker = make_tracker(MLParameters(merge_distance_deg=0.0))
        candidates = [self._cand(0.0, 50.0, 0.6), self._cand(0.0, 50.0, 0.9)]
        assert tracker._merge_nearby_candidates(candidates) == candidates

    def test_merge_across_antimeridian(self, make_tracker):
        """Test that candidates either side of the date line are merged."""
        candidates = [self._cand(0.0, 179.5, 0.6), self._cand(0.0, -179.5, 0.9)]
        merged = make_tracker()._merge_nearby_candidates(candidates)
        assert merged == [self._cand(0.0, -179.5, 0.9)]

    def test_only_compared_with_kept_candidates(self, make_tracker):
        """Test that a merged-away candidate does not suppress others."""
        candidates = [
            self._cand(0.0, 50.0, 0.9),
            self._cand(1.5, 50.0, 0.8),  # within 2 degrees of the first: dropped
            self._cand(3.0, 50.0, 0.7),  # within 2 of the dropped one only: kept
        ]
        merged = make_tracker()._merge_nearby_candidates(candidates)
        assert [c["score"] for c in merged] == [0.9, 0.7]


class TestStitch:
    """Tests for MLTracker.stitch and MLTracker._nearest_candidate."""

    @staticmethod
    def _stitch(tracker, n_times, candidates):
        """Run stitch() on the given candidates over ``n_times`` timesteps."""
        tracker._times = [_time(i) for i in range(n_times)]
        tracker._candidates = candidates
        return tracker.stitch()

    def test_single_storm_linked(self, make_tracker):
        """Test that a slowly moving storm becomes a single trajectory."""
        candidates = [
            _candidate(
                i,
                10.0 + 0.5 * i,
                50.0 + 0.5 * i,
                data={"sea_surface_temperature": 300.0 + i},
            )
            for i in range(3)
        ]
        trajectories = self._stitch(make_tracker(), 3, candidates)
        assert len(trajectories) == 1
        assert trajectories[0].observations == 3
        assert trajectories[0].data["lat"] == [10.0, 10.5, 11.0]
        assert trajectories[0].data["lon"] == [50.0, 50.5, 51.0]
        assert trajectories[0].data["time"] == [_time(i) for i in range(3)]
        # The values sampled from the input are carried into the trajectory points.
        assert trajectories[0].data["sea_surface_temperature"] == [300.0, 301.0, 302.0]

    def test_jump_beyond_max_distance_not_linked(self, make_tracker):
        """Test that candidates further than max_distance_deg start a new track."""
        stitch_parameters = MLStitchParameters(max_distance_deg=3.0, min_length=1)
        tracker = make_tracker(stitch_parameters=stitch_parameters)
        candidates = [_candidate(0, 10.0, 50.0), _candidate(1, 10.0, 54.0)]
        assert len(self._stitch(tracker, 2, candidates)) == 2

    def test_link_across_antimeridian(self, make_tracker):
        """Test that a storm crossing the date line stays a single track."""
        candidates = [_candidate(0, 0.0, 179.5), _candidate(1, 0.0, -179.5)]
        trajectories = self._stitch(make_tracker(), 2, candidates)
        assert len(trajectories) == 1
        assert trajectories[0].observations == 2

    def test_gap_within_max_gap_bridged(self, make_tracker):
        """Test that a track survives being unmatched for up to max_gap steps."""
        stitch_parameters = MLStitchParameters(max_gap=1)
        tracker = make_tracker(stitch_parameters=stitch_parameters)
        candidates = [_candidate(i, 10.0, 50.0) for i in (0, 1, 3)]
        trajectories = self._stitch(tracker, 4, candidates)
        assert len(trajectories) == 1
        assert trajectories[0].observations == 3

    def test_gap_beyond_max_gap_splits_track(self, make_tracker):
        """Test that a longer gap closes the track and a new one starts."""
        stitch_parameters = MLStitchParameters(max_gap=1, min_length=1)
        tracker = make_tracker(stitch_parameters=stitch_parameters)
        candidates = [_candidate(i, 10.0, 50.0) for i in (0, 1, 4)]
        trajectories = self._stitch(tracker, 5, candidates)
        assert sorted(t.observations for t in trajectories) == [1, 2]

    def test_short_tracks_removed(self, make_tracker):
        """Test that tracks below min_length are dropped."""
        candidates = [_candidate(0, 10.0, 50.0)]
        assert self._stitch(make_tracker(), 2, candidates) == []

    def test_candidate_only_claimed_once(self, make_tracker):
        """Test that one candidate cannot extend two tracks at the same time."""
        tracker = make_tracker(stitch_parameters=MLStitchParameters(min_length=1))
        candidates = [
            _candidate(0, 0.0, 50.0),
            _candidate(0, 2.0, 50.0),
            _candidate(1, 1.0, 50.0),
        ]
        trajectories = self._stitch(tracker, 2, candidates)
        assert sorted(t.observations for t in trajectories) == [1, 2]

    def test_nearest_candidate_picks_closest(self, make_tracker):
        """Test that the closest unmatched candidate within range is chosen."""
        tracker = make_tracker()
        track = {"last": _candidate(0, 10.0, 50.0)}
        candidates = [_candidate(1, 12.0, 50.0), _candidate(1, 10.5, 50.0)]
        assert tracker._nearest_candidate(track, candidates, {0, 1}) == 1

    def test_nearest_candidate_none_in_range(self, make_tracker):
        """Test that None is returned if every candidate is too far away."""
        tracker = make_tracker()
        track = {"last": _candidate(0, 10.0, 50.0)}
        candidates = [_candidate(1, 20.0, 50.0)]
        assert tracker._nearest_candidate(track, candidates, {0}) is None


class TestDetectionsToNetcdf:
    """Tests for MLTracker.detections_to_netcdf."""

    @pytest.fixture
    def tracker(self, make_tracker):
        """Tracker on the sample file holding two detections."""
        tracker = make_tracker(MLParameters(input_file=str(SAMPLE_FILE)))
        start = datetime(2025, 1, 10, calendar="proleptic_gregorian")
        tracker._candidates = [
            {
                "time": start,
                "lat": -13.0,
                "lon": 55.0,
                "class_index": 3.0,
                "score": 0.6,
                "data": {
                    "sea_surface_temperature": 301.0,
                    "air_temperature_500": 265.0,
                },
            },
            {
                "time": start + timedelta(hours=6),
                "lat": -12.5,
                "lon": 54.0,
                "class_index": 4.0,
                "score": 0.7,
                "data": {
                    "sea_surface_temperature": 302.0,
                    "air_temperature_500": 266.0,
                },
            },
        ]
        return tracker

    @staticmethod
    def _read(path):
        """Read a written file into {netCDF variable name: field}."""
        return {field.nc_get_variable("?"): field for field in cf.read(str(path))}

    def test_no_detections_warns_and_writes_nothing(self, make_tracker, tmp_path):
        """With no detections a warning is given and no file is created."""
        tracker = make_tracker(MLParameters(input_file=str(SAMPLE_FILE)))
        output_file = tmp_path / "detections.nc"
        with pytest.warns(UserWarning, match="no detections"):
            tracker.detections_to_netcdf(str(output_file))
        assert not output_file.exists()

    def test_values_and_coordinates_round_trip(self, tracker, tmp_path):
        """The values and the latitude and longitude read back unchanged."""
        tracker.detections_to_netcdf(str(tmp_path / "detections.nc"))
        score = self._read(tmp_path / "detections.nc")["score"]
        assert np.allclose(score.array, [0.6, 0.7])
        assert np.allclose(score.coordinate("latitude").array, [-13.0, -12.5])
        assert np.allclose(score.coordinate("longitude").array, [55.0, 54.0])
        assert score.coordinate("time").array.size == 2

    def test_cf_attributes(self, tracker, tmp_path):
        """The fields are CF point data with the metadata of each variable."""
        tracker.detections_to_netcdf(str(tmp_path / "detections.nc"))
        written = self._read(tmp_path / "detections.nc")
        assert all(f.get_property("featureType") == "point" for f in written.values())
        assert written["class_index"].get_property("flag_meanings")
        sst = written["sea_surface_temperature"]
        assert sst.get_property("standard_name") == "sea_surface_temperature"
        assert np.allclose(sst.array, [301.0, 302.0])


class TestRunTracker:
    """Tests for MLTracker.run_tracker, with a mocked model on the ERA5 sample."""

    def test_storm_becomes_trajectory_in_file(self, make_tracker, tmp_path):
        """A storm in the first 3 timesteps becomes one trajectory, written to file."""
        tracker = make_tracker(MLParameters(input_file=str(SAMPLE_FILE)))
        storm = _logits({(40, 40): (4, 10.0)}, 80, 80)
        quiet = _logits(None, 80, 80)
        tracker.model = MagicMock(side_effect=[storm] * 3 + [quiet] * 7)
        output_file = tmp_path / "tracks.nc"

        tracker.run_tracker(str(output_file))

        trajectories = tracker.read_trajectories()
        assert len(trajectories) == 1
        assert trajectories[0].observations == 3
        assert trajectories[0].data["lat"] == [float(tracker._lats[40])] * 3
        written = cf.read(str(output_file))
        assert written
        assert written[0].get_property("featureType") == "trajectory"
