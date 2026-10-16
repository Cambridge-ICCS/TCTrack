"""Tests for the geographic module."""

from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

from tctrack.build_db.build_db import read_netcdf
from tctrack.geographic import Gate, geographic, load_gates, match_gates

NETCDF_NATL_FILE = str(Path(__file__).parent / "6hr_track_mm_2013_natl.nc")
NETCDF_PACIFIC_FILE = str(Path(__file__).parent / "test_tracks.nc")


GATE_HEADER = "code,name,min_lat,max_lat,min_lon,max_lon\n"
GATE_ROW = "A,First,0,10,-20,-10\n"


class TestGatesErrors:
    """Reject invalid definitions before any observation matching."""

    @pytest.mark.parametrize(
        "content",
        [
            "",
            GATE_HEADER,
            "code,name,min_lat,max_lat,min_lon\n",
            GATE_HEADER.replace("max_lon", "code"),
            GATE_HEADER.rstrip() + ",extra\n",
            GATE_HEADER + "A,First,0,10,-20\n",
            GATE_HEADER + GATE_ROW.rstrip() + ",extra\n",
            GATE_HEADER + GATE_ROW.replace("A,", " ,"),
            GATE_HEADER + GATE_ROW + GATE_ROW,
            GATE_HEADER + GATE_ROW.replace("First", " "),
            GATE_HEADER + GATE_ROW.replace("0,10", "bad,10"),
            GATE_HEADER + GATE_ROW.replace("0,10", "nan,10"),
            GATE_HEADER + GATE_ROW.replace("0,10", "0,inf"),
            GATE_HEADER + GATE_ROW.replace("0,10", "-inf,10"),
            GATE_HEADER + GATE_ROW.replace("0,10", "10,0"),
            GATE_HEADER + GATE_ROW.replace("0,10", "10,10"),
            GATE_HEADER + GATE_ROW.replace("0,10", "-91,10"),
            GATE_HEADER + GATE_ROW.replace("0,10", "0,91"),
            GATE_HEADER + GATE_ROW.replace("-20,-10", "170,-170"),
            GATE_HEADER + GATE_ROW.replace("-20,-10", "-10,-10"),
            GATE_HEADER + GATE_ROW.replace("-20,-10", "-181,-10"),
            GATE_HEADER + GATE_ROW.replace("-20,-10", "-20,181"),
            GATE_HEADER + GATE_ROW.replace("First", '"unclosed'),
        ],
    )
    def test_invalid_csv(self, tmp_path, content):
        """Group malformed headers, fields, identifiers and bounds."""
        path = tmp_path / "invalid.csv"
        path.write_text(content, encoding="utf-8")

        with pytest.raises(ValueError) as caught:
            load_gates(path)

        assert str(path) in str(caught.value)
        assert "row" in str(caught.value)

    def test_file_errors(self, tmp_path):
        """Missing and unreadable paths report their filename."""
        for path in (tmp_path / "missing.csv", tmp_path):
            with pytest.raises(OSError) as caught:
                load_gates(path)

            assert caught.value.filename == str(path)

    def test_overlap_and_duplicate_errors(self, tmp_path):
        """Reject positive-area overlaps and duplicates with useful context."""
        path = tmp_path / "invalid.csv"
        for row, message in (
            ("B,Second,5,15,-15,-5\n", "overlapping gates 'A' and 'B'"),
            ("A,Second,10,20,-20,-10\n", "duplicate code 'A'"),
        ):
            path.write_text(GATE_HEADER + GATE_ROW + row, encoding="utf-8")
            with pytest.raises(ValueError, match=message) as caught:
                load_gates(path)

            assert str(path) in str(caught.value)
            assert "row 3" in str(caught.value)

    def test_match_shape_errors(self):
        """Do not silently broadcast misaligned observation arrays."""
        for lat, lon, land in (([1, 2], [1], [True, True]), ([1], [1], [True, False])):
            with pytest.raises(ValueError, match="same shape"):
                match_gates(lat, lon, land, ())


class TestGates:
    """Test immutable configuration and half-open rectangle matching."""

    def test_load_valid_gates(self, tmp_path):
        """Trim identifiers, preserve case, and allow shared edges and names."""
        path = tmp_path / "gates.csv"
        path.write_text(
            GATE_HEADER
            + " A , First ,0,10,-20,-10\n"
            + "a,First,10,20,-20,-10\n"
            + "B,First,0,10,-10,0\n"
            + "WORLD,Limits,-90,90,0,180\n",
            encoding="utf-8",
        )
        gates = load_gates(str(path))

        assert isinstance(gates, tuple)
        assert gates[0] == Gate("A", "First", 0, 10, -20, -10)
        assert [gate.code for gate in gates] == ["A", "a", "B", "WORLD"]
        with pytest.raises(FrozenInstanceError):
            gates[0].name = "Changed"

    def test_match_edges_and_invalid_observations(self, tmp_path):
        """Match only landfalls, wrapping longitude and excluding masked padding."""
        path = tmp_path / "gates.csv"
        path.write_text(
            GATE_HEADER + GATE_ROW + "B,Second,10,20,-20,-10\n",
            encoding="utf-8",
        )
        gates = load_gates(path)
        lat = np.ma.array(
            [0, 10, 20, 5, 5, 5, 5, 5, np.nan, 5, 91, 5, 5, 5],
            mask=[False] * 11 + [True, False, False],
        )
        lon = np.ma.array(
            [-20, -20, -20, -10, -21, 345, -15, -15, -15, np.inf, -15, -15, -15, -15],
            mask=[False] * 12 + [True, False],
        )
        land = np.ma.array([True] * 14, mask=[False] * 13 + [True])
        land[6] = False
        lat_before, lon_before, land_before = lat.copy(), lon.copy(), land.copy()

        result = match_gates(lat, lon, land, gates)

        assert result.dtype == object
        assert result.shape == lat.shape
        assert result.tolist() == [
            "A",
            "B",
            None,
            None,
            None,
            "A",
            None,
            "A",
            None,
            None,
            None,
            None,
            None,
            None,
        ]
        for actual, expected in (
            (lat, lat_before),
            (lon, lon_before),
            (land, land_before),
        ):
            np.testing.assert_equal(actual, expected)
            np.testing.assert_array_equal(actual.mask, expected.mask)

        assert match_gates(lat, lon, land, ()).tolist() == [None] * len(lat)
        assert match_gates([], [], [], gates).shape == (0,)

        # Shared longitude boundaries belong only to the gate on their east.
        path.write_text(
            GATE_HEADER + GATE_ROW + "E,East,0,10,-10,0\n",
            encoding="utf-8",
        )
        result = match_gates(
            [[0, 0, 0], [9, 10, 5]],
            [[-20, -10, 0], [350, -10, -15]],
            [[True] * 3] * 2,
            load_gates(path),
        )
        assert result.tolist() == [["A", "E", None], ["E", None, "A"]]


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
            assert set(entry) == {
                "landfall",
                "distance_to_coast_km",
                "ocean_basin",
                "gate_id",
            }
            assert entry["gate_id"].shape == (netcdf_data["n_observations"],)
            assert entry["gate_id"].dtype == object
            assert all(code is None for code in entry["gate_id"])

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

    def test_gates_preserve_classification(self, tmp_path):
        """Gate existing fixture landfalls without changing any existing output."""
        path = tmp_path / "gates.csv"
        path.write_text(
            GATE_HEADER + "WEST,West,-90,90,-180,0\nEAST,East,-90,90,0,180\n",
            encoding="utf-8",
        )
        gates = load_gates(path)
        data = read_netcdf(NETCDF_NATL_FILE)
        original = geographic.classify_tracks(data)
        gated = geographic.classify_tracks(data, gates=gates)

        for index, entry in gated.items():
            for key in ("landfall", "distance_to_coast_km", "ocean_basin"):
                np.testing.assert_equal(entry[key], original[index][key])

            assert all(code is None for code in entry["gate_id"][~entry["landfall"]])
            assert all(code == "WEST" for code in entry["gate_id"][entry["landfall"]])

        # A broad rectangle contains London and open Atlantic water, but not Denver.
        path.write_text(GATE_HEADER + "UK,UK,20,60,-50,5\n", encoding="utf-8")
        data = {
            "latitude": np.ma.array(
                [[51.5, 30, 39.7, np.nan, 51.5, 51.5]],
                mask=[[False, False, False, False, True, False]],
            ),
            "longitude": np.ma.array(
                [[-0.1, -40, -105, np.nan, -0.1, -0.1]],
                mask=[[False, False, False, False, False, True]],
            ),
        }
        before = {key: value.copy() for key, value in data.items()}
        original = geographic.classify_tracks(data)[0]
        entry = geographic.classify_tracks(data, load_gates(path))[0]

        assert entry["landfall"].tolist() == [True, False, True, False, False, False]
        assert entry["gate_id"].tolist() == ["UK", None, None, None, None, None]
        for key in ("landfall", "distance_to_coast_km", "ocean_basin"):
            np.testing.assert_equal(entry[key], original[key])
        for key in data:
            np.testing.assert_equal(data[key], before[key])
            np.testing.assert_array_equal(data[key].mask, before[key].mask)

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
