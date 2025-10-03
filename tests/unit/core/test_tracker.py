"""Unit tests for tracker.py of the TCTrack Core Python package."""

from dataclasses import dataclass

import cf
import numpy as np
import pytest

from tctrack.core import TCTracker, TCTrackerParameters, Trajectory


class TestTCTrackerParameters:
    """Tests for the TCTrackerParameters base class."""

    @dataclass(repr=False)
    class ExampleParameters(TCTrackerParameters):
        """Concrete implementation of TCTrackerParameters for testing purposes."""

        param_a: int
        param_b: str

    def test_str_representation(self):
        """Test the string representation of TCTrackerParameters."""
        params = self.ExampleParameters(param_a=42, param_b="test")
        expected_output = (
            "ExampleParameters(\n\tparam_a \t = 42\n\tparam_b \t = test\n)"
        )
        assert str(params) == expected_output

    def test_repr_representation(self):
        """Test the repr representation of TCTrackerParameters."""
        params = self.ExampleParameters(param_a=42, param_b="test")
        expected_output = (
            "ExampleParameters(\n\tparam_a \t = 42\n\tparam_b \t = test\n)"
        )
        assert repr(params) == expected_output


class TestTCTracker:
    """Tests for the TCTracker abstract base class."""

    class ExampleTracker(TCTracker):
        """Concrete implementation of TCTracker for testing purposes."""

        def __init__(self, example_trajectories):
            self._example_trajectories = example_trajectories

        def read_variable_metadata(self) -> None:
            """Implement a dummy of the read_variable_metadata abstractmethod."""
            self._variable_metadata = {
                "test_var": {
                    "standard_name": "test_standard_name",
                    "long_name": "Test Long Name",
                    "units": "test_units",
                },
                "lat": {"long_name": "latitude", "units": "degrees_north"},
                "lon": {"long_name": "longitude", "units": "degrees_east"},
            }

        def trajectories(self) -> list[Trajectory]:
            """Implement a dummy of the trajectories abstractmethod."""
            return self._example_trajectories

        def run_tracker(self, output_file: str) -> None:  # noqa:ARG002 unused arg
            """Implement a dummy of the run_tracker abstractmethod."""
            return None

    def test_abstract_methods_instantiation(self):
        """Ensure TCTracker ABC cannot be instantiated directly."""
        with pytest.raises(
            TypeError, match="Can't instantiate abstract class TCTracker"
        ):
            TCTracker()

    def test_variable_metadata_uninitialized(self):
        """Test accessing `variable_metadata` before initialization raises an error."""
        tracker = self.ExampleTracker(example_trajectories=None)
        with pytest.raises(
            AttributeError, match="_variable_metadata has not been initialized"
        ):
            _ = tracker.variable_metadata

    def test_variable_metadata_initialized(self):
        """Test that `variable_metadata` is correctly initialized by the subclass."""
        tracker = self.ExampleTracker(example_trajectories=None)
        tracker.read_variable_metadata()
        expected_metadata = {
            "test_var": {
                "standard_name": "test_standard_name",
                "long_name": "Test Long Name",
                "units": "test_units",
            },
            "lat": {"long_name": "latitude", "units": "degrees_north"},
            "lon": {"long_name": "longitude", "units": "degrees_east"},
        }
        assert tracker.variable_metadata == expected_metadata

    def test_to_netcdf(self, tmp_path):
        """
        Test writing trajectories to NetCDF.

        We will take some predefined Trajectories (matching the variable_metadata of
        ExampleTracker), write to NetCDF, then read teh resulting file back in to
        confirm it matches what was expected.
        """
        # Simple example trajectories matching variable_metadata of ExampleTracker
        trajectory1 = Trajectory(
            trajectory_id=1,
            observations=2,
            year=2023,
            month=10,
            day=1,
            hour=0,
            calendar="standard",
        )
        trajectory1.add_multiple_points(
            years=[2023, 2023],
            months=[10, 10],
            days=[1, 1],
            hours=[0, 6],
            variables={
                "lat": [10.0, 15.0],
                "lon": [20.0, 25.0],
                "test_var": [5.0, 10.0],
            },
        )

        trajectory2 = Trajectory(
            trajectory_id=2,
            observations=3,
            year=2023,
            month=10,
            day=2,
            hour=0,
            calendar="standard",
        )
        trajectory2.add_multiple_points(
            years=[2023, 2023, 2023],
            months=[10, 10, 10],
            days=[2, 2, 2],
            hours=[0, 18, 21],
            variables={
                "lat": [30.0, 35.0, 35.0],
                "lon": [40.0, 45.0, 50.0],
                "test_var": [15.0, 20.0, 15.0],
            },
        )

        # Instantiate the dummy tracker with valid trajectories
        tracker = self.ExampleTracker([trajectory1, trajectory2])

        # Write to NetCDF
        output_file = tmp_path / "trajectories.nc"
        tracker.read_variable_metadata()
        tracker.to_netcdf(str(output_file))

        # Read back with cf.read
        fields = cf.read(output_file)

        # Validate the structure and content - one variable field
        assert len(fields) == 1
        field = fields[0]

        # Validate trajectory and observation dimensions
        trajectory_ids = field.dimension_coordinate("trajectory").size
        observation_count = field.dimension_coordinate("observation").size
        assert trajectory_ids == 2
        assert observation_count == 3

        # Validate key coordinates like time, latitude, and longitude exist
        time_coord = field.construct("time")
        lat_coord = field.construct("lat")
        lon_coord = field.construct("lon")
        assert time_coord.shape == (trajectory_ids, observation_count)
        assert lat_coord.shape == lon_coord.shape == (trajectory_ids, observation_count)

        # Check the metadata - not we need to treat constructs (coordinates) and
        # fields (variables) separately
        for variable, metadata in tracker.variable_metadata.items():
            if variable in {"lat", "lon", "time"}:  # Coordinate constructs
                construct = field.construct(variable)
                assert construct is not None, f"Construct for {variable} is missing"
                for key, value in metadata.items():
                    assert construct.get_property(key) == value, (
                        f"Metadata mismatch for {variable}: {key}"
                    )
            else:  # Field data (variables)
                for key, value in metadata.items():
                    assert field.get_property(key) == value, (
                        f"Metadata mismatch for field data: {variable} - {key}"
                    )

        # Validate data in arrays
        # Note that the first trajectory has a nan appended due to observation mismatch
        var_data = field.data.array
        assert np.allclose(var_data[0, :], [5.0, 10.0, float("nan")], equal_nan=True)
        assert np.allclose(var_data[1, :], [15.0, 20.0, 15.0])
