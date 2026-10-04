"""Tests for the geographic module."""

from pathlib import Path

import numpy as np
import pytest

from tctrack.build_db.build_db import read_netcdf
from tctrack.geographic import geographic

NETCDF_NATL_FILE = str(Path(__file__).parent / "6hr_track_mm_2013_natl.nc")
NETCDF_PACIFIC_FILE = str(Path(__file__).parent / "test_tracks.nc")


class TestClassifyPointsErrors:
    """Test error handling of classify_points."""

    def test_shape_mismatch(self):
        """Test that mismatched latitude and longitude shapes error."""
        with pytest.raises(ValueError, match="shape"):
            geographic.classify_points([10.0, 20.0], [-140.0])

    def test_latitude_out_of_range(self):
        """Test that a latitude outside [-90, 90] errors."""
        with pytest.raises(ValueError, match="latitude"):
            geographic.classify_points([95.0], [0.0])


class TestClassifyTracksErrors:
    """Test error handling of classify_tracks."""

    def test_missing_keys(self):
        """Test that a dataset without latitude or longitude data errors."""
        with pytest.raises(ValueError, match="latitude"):
            geographic.classify_tracks({"longitude": [[0.0]]})

        with pytest.raises(ValueError, match="longitude"):
            geographic.classify_tracks({"latitude": [[0.0]]})

    def test_not_two_dimensional(self):
        """Test that one-dimensional coordinates error."""
        data = {"latitude": [10.0, 20.0], "longitude": [-140.0, -150.0]}

        with pytest.raises(ValueError, match="two-dimensional"):
            geographic.classify_tracks(data)


class TestClassifyPoints:
    """Test classification of individual points."""

    def test_land_points(self):
        """Test points well inside land are land, inland of the coast."""
        landfall, distance, basin = geographic.classify_points(
            [51.5, 39.7, 21.3], [-0.1, -105.0, -157.9]
        )

        assert landfall.tolist() == [True, True, True]
        assert (distance > 0).all()

        # London, Denver and Honolulu: about 100, 1200 and 25 km from the sea
        assert 50 < distance[0] < 200
        assert 1000 < distance[1] < 1400
        assert distance[2] < 50

        # Honolulu sits in the middle of the North Pacific Ocean
        assert basin[2] == "North Pacific Ocean"

    def test_sea_points_and_basins(self):
        """Test mid-ocean and sea points get sensible basins."""
        landfall, distance, basin = geographic.classify_points(
            [30.0, -30.0, 10.0, -20.0, -20.0, 85.0, 25.0, 15.0, 35.0],
            [-40.0, -20.0, -140.0, -170.0, 80.0, 0.0, -90.0, -75.0, 18.0],
        )

        assert landfall.tolist() == [False] * 9
        assert (distance > 0).all()

        # Marginal seas take the ocean they connect to: the Gulf of Mexico,
        # Caribbean Sea and Mediterranean points are all North Atlantic
        assert basin.tolist() == [
            "North Atlantic Ocean",
            "South Atlantic Ocean",
            "North Pacific Ocean",
            "South Pacific Ocean",
            "Indian Ocean",
            "Arctic Ocean",
            "North Atlantic Ocean",
            "North Atlantic Ocean",
            "North Atlantic Ocean",
        ]

    def test_nan_points_skipped(self):
        """Test that NaN points are unclassified rather than erroring."""
        landfall, distance, basin = geographic.classify_points(
            [10.0, np.nan], [-140.0, np.nan]
        )

        assert landfall.tolist() == [False, False]
        assert distance[0] > 0
        assert np.isnan(distance[1])
        assert basin[0] == "North Pacific Ocean"
        assert basin[1] is None

    def test_longitude_conventions(self):
        """Test that both longitude conventions give the same result."""
        for lon in (-140.0, 220.0):
            landfall, _, basin = geographic.classify_points([10.0], [lon])

            assert not landfall[0]
            assert basin[0] == "North Pacific Ocean"


class TestClassifyTracks:
    """Test classification of whole NetCDF datasets."""

    def test_output_links_to_trajectory_index(self):
        """Test the output is keyed by trajectory index, one entry per point."""
        netcdf_data = read_netcdf(NETCDF_NATL_FILE)
        result = geographic.classify_tracks(netcdf_data)

        assert set(result) == set(range(netcdf_data["n_trajectories"]))

        lat = np.asarray(netcdf_data["latitude"], dtype=float)
        lon = np.asarray(netcdf_data["longitude"], dtype=float)

        for traj_idx, entry in result.items():
            assert set(entry) == {"landfall", "distance_to_coast_km", "ocean_basin"}

            expected_shape = (netcdf_data["n_observations"],)
            assert entry["landfall"].shape == expected_shape
            assert entry["distance_to_coast_km"].shape == expected_shape
            assert entry["ocean_basin"].shape == expected_shape

            # Entries align positionally with the input observation dimension
            valid = np.nonzero(~np.isnan(lat[traj_idx]))[0]
            landfall, distance, basin = geographic.classify_points(
                lat[traj_idx, valid], lon[traj_idx, valid]
            )

            assert landfall.tolist() == entry["landfall"][valid].tolist()
            assert distance.tolist() == entry["distance_to_coast_km"][valid].tolist()
            assert basin.tolist() == entry["ocean_basin"][valid].tolist()

    def test_padding_points_unclassified(self):
        """Test that NaN padded observations are left unclassified."""
        netcdf_data = read_netcdf(NETCDF_NATL_FILE)
        result = geographic.classify_tracks(netcdf_data)

        for traj_idx, entry in result.items():
            valid = ~np.isnan(
                np.asarray(netcdf_data["latitude"][traj_idx], dtype=float)
            )
            padding = ~valid

            assert entry["landfall"][padding].sum() == 0
            assert np.isnan(entry["distance_to_coast_km"][padding]).all()
            assert all(entry["ocean_basin"][i] is None for i in np.nonzero(padding)[0])

        # An entirely NaN trajectory yields full-length unclassified arrays
        # rather than erroring
        data = {"latitude": [[np.nan, np.nan]], "longitude": [[np.nan, np.nan]]}
        entry = geographic.classify_tracks(data)[0]

        assert entry["landfall"].shape == (2,)
        assert not entry["landfall"].any()
        assert np.isnan(entry["distance_to_coast_km"]).all()
        assert all(basin is None for basin in entry["ocean_basin"])

    def test_north_atlantic_tracks(self):
        """Test tracks with landfall in the North Atlantic region."""
        netcdf_data = read_netcdf(NETCDF_NATL_FILE)
        result = geographic.classify_tracks(netcdf_data)

        lat = np.asarray(netcdf_data["latitude"], dtype=float)
        valid = ~np.isnan(lat)

        landfall = np.array([entry["landfall"] for entry in result.values()])
        distance = np.stack(
            [entry["distance_to_coast_km"] for entry in result.values()]
        )
        basins = {
            basin
            for entry in result.values()
            for basin in entry["ocean_basin"]
            if basin is not None
        }

        # Some but not all points are over land (Central America, Yucatan,
        # the US East Coast and Nova Scotia in 2013)
        assert 0 < landfall.sum() < valid.sum()

        # Land points are all within 600 km of a coast, sea points are not
        # all far from land
        assert (distance[landfall] > 0).all()
        assert (distance[landfall] < 600).all()
        assert (distance[~landfall & valid] < 100).any()

        # Most points are North Atlantic, with some on the Pacific side of
        # Central America
        assert {"North Atlantic Ocean", "North Pacific Ocean"} <= basins

    def test_pacific_tracks(self):
        """Test tracks that stay over the open western North Pacific."""
        netcdf_data = read_netcdf(NETCDF_PACIFIC_FILE)
        result = geographic.classify_tracks(netcdf_data)

        landfall = np.array([entry["landfall"] for entry in result.values()])
        distance = np.stack(
            [entry["distance_to_coast_km"] for entry in result.values()]
        )
        basins = {
            basin
            for entry in result.values()
            for basin in entry["ocean_basin"]
            if basin is not None
        }

        # No point makes landfall, and every point is far from land
        assert landfall.sum() == 0
        assert (distance[~np.isnan(distance)] > 400).all()

        assert basins == {"North Pacific Ocean"}
