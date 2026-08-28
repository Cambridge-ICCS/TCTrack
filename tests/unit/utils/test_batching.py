"""Unit tests for batching utility functions."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import cf
import pytest

from tctrack.core import TCTracker
from tctrack.utils.batching import batching

### Dummy functions for testing the preprocessing


def dummy_step(comment="dummy step", log: list | None = None) -> None:
    """Perform a dummy preproprocessing step that just updates the log."""
    if log is not None:
        log.append(comment)


def make_field(name: str, log: list[str] | None = None) -> cf.Field:
    """Create a field with a netcdf variable name."""
    field = cf.example_field(0).copy()
    field.nc_set_variable(name)
    if log is not None:
        log.append(f"created {name}")
    return field


### Dummy tracker to use in the batching


class DummyTracker(TCTracker):
    """Concrete tracker used to test batching utilities."""

    def __init__(self):
        self.inputs: list[list[str]] = []
        self.outputs: list[str] = []

    @property
    def _parameters(self) -> list:
        """Return no parameters for the test tracker."""
        return []

    def _set_metadata(self) -> None:
        """Implement the abstract method for testing."""

    def read_trajectories(self) -> list:
        """Implement the abstract method for testing."""
        return []

    def run_tracker(self, input_files: str | Iterable[str], output_file: str) -> None:
        """Write a minimal trajectory NetCDF file for combination tests."""
        if isinstance(input_files, str):
            input_files = [input_files]
        self.inputs.append(list(input_files))
        self.outputs.append(output_file)


### Tests for the batching function


@pytest.fixture
def config(tmp_path: Path):
    """Batching config for the tests."""
    return {
        "output_dir": tmp_path,
        "combine_outputs": False,
        "delete_batch_dirs": False,
    }


class TestBatchingPreprocessing:
    """Tests for the preprocessing stage of the batching utility."""

    def test_preprocessing_steps(self, config) -> None:
        """Test the preprocessing functions are called by batching."""
        tracker = DummyTracker()
        log: list[str] = []

        batching(
            tracker,
            n_iter=1,
            input_files="input.nc",
            preprocessing=[
                (dummy_step, {"log": log}),
                (make_field, {"name": "p", "log": log}),
            ],
            config=config,
        )

        assert log == ["dummy step", "created p"]

    def test_tag_expansion(self, config) -> None:
        """Test the %ITER% and %BATCH% tags in any arguments get expanded."""
        tracker = DummyTracker()
        log: list[str] = []

        batching(
            tracker,
            n_iter=2,
            input_files="input.nc",
            preprocessing=[(dummy_step, {"comment": "%ITER%%BATCH%", "log": log})],
            config=config,
        )

        batch_dirs = [str(Path(config["output_dir"]) / f"batch_{i}") for i in range(2)]
        assert log == [f"{i}{batch_dirs[i]}" for i in range(2)]

    def test_tag_expansion_iter_only(self, config) -> None:
        """Test %ITER% tag on its own gets replaced with an integer."""
        tracker = DummyTracker()
        log: list = []

        batching(
            tracker,
            n_iter=2,
            input_files="input.nc",
            preprocessing=[(dummy_step, {"comment": "%ITER%", "log": log})],
            config=config,
        )

        assert log == [0, 1]
        assert isinstance(log[0], int)


class TestBatching:
    """Tests for the batching utility."""

    def test_batch_directories(self, config) -> None:
        """Test batching creates the expected batch directories."""
        tracker = DummyTracker()
        n_iter = 2

        batching(tracker, n_iter=n_iter, input_files="input.nc", config=config)

        # Check the directories have been created
        batch_dirs = [config["output_dir"] / f"batch_{i}" for i in range(n_iter)]
        for batch_dir in batch_dirs:
            assert batch_dir.exists()

    def test_batch_directories_deletion(self, config) -> None:
        """Test batching deletes the batch directories when delete_batch_dirs=True."""
        tracker = DummyTracker()
        n_iter = 2

        del config["delete_batch_dirs"]  # Default is True

        batching(tracker, n_iter=n_iter, input_files="input.nc", config=config)

        # Check the directories have been deleted
        batch_dirs = [config["output_dir"] / f"batch_{i}" for i in range(n_iter)]
        for batch_dir in batch_dirs:
            assert not batch_dir.exists()

    def test_retrieve_data(self, config) -> None:
        """Test batching correctly calls a provided retrieve_data function."""
        tracker = DummyTracker()
        n_iter = 2
        log: list[str] = []

        def retrieve_data(i_iter: int, _: Path) -> None:
            log.append(f"retrieve_data {i_iter}")

        batching(
            tracker,
            n_iter=n_iter,
            input_files="input.nc",
            retrieve_data=retrieve_data,
            config=config,
        )

        assert log == [f"retrieve_data {i}" for i in range(n_iter)]

    def test_input_output_files(self, config) -> None:
        """Test batching correctly sets the per-batch input and output files."""
        tracker = DummyTracker()
        n_iter = 2
        input_files = ["file1", "file2"]

        batching(tracker, n_iter=n_iter, input_files=input_files, config=config)

        # Check the input and output files were set correctly
        batch_dirs = [config["output_dir"] / f"batch_{i}" for i in range(n_iter)]
        assert tracker.inputs == [
            [str(batch_dir / file) for file in input_files] for batch_dir in batch_dirs
        ]
        assert tracker.outputs == [
            str(config["output_dir"] / f"tracks_{i}.nc") for i in range(n_iter)
        ]

    def test_input_file_wildcards(self, config) -> None:
        """Test batching expands input file wildcards in each batch directory."""
        tracker = DummyTracker()

        def retrieve_data(_: int, batch_dir: Path) -> None:
            (batch_dir / "input_2.nc").touch()
            (batch_dir / "input_1.nc").touch()

        batching(
            tracker,
            n_iter=2,
            input_files="input_*.nc",
            retrieve_data=retrieve_data,
            config=config,
        )

        assert tracker.inputs == [
            [
                str(config["output_dir"] / f"batch_{i}" / "input_1.nc"),
                str(config["output_dir"] / f"batch_{i}" / "input_2.nc"),
            ]
            for i in range(2)
        ]

    def test_input_file_wildcards_missing(self, config) -> None:
        """Test batching fails if there are no input files that match a wildcard."""
        tracker = DummyTracker()

        with pytest.raises(FileNotFoundError, match="No files matched"):
            batching(tracker, n_iter=1, input_files="input_*.nc", config=config)
