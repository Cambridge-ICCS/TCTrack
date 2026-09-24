"""Unit tests for batching utility functions."""

from __future__ import annotations

from contextlib import nullcontext
from datetime import timedelta
from pathlib import Path
from typing import Iterable
from unittest.mock import call

import cf
import pytest

from tctrack.core import TCTracker
from tctrack.utils.batching import (
    _batch_time_ranges,
    _expand_batch_time_range,
    _get_calendar,
    batching,
)

BATCH_INTERVAL = cf.TimeDuration(1, "days")
BATCH_TIME_RANGE = ("2000-01-01", "2000-01-02")
TWO_BATCH_TIME_RANGE = ("2000-01-01", "2000-01-03")

### Dummy functions for testing the preprocessing


def dummy_step(comment="dummy step", log: list | None = None) -> None:
    """Perform a dummy preproprocessing step that just updates the log."""
    if log is not None:
        log.append(comment)


def make_field(
    name: str, log: list[str] | None = None, calendar: str = "standard"
) -> cf.Field:
    """Create a field with a netcdf variable name."""
    field = cf.example_field(0).copy()
    field.nc_set_variable(name)
    field.dimension_coordinate("T").set_data([cf.dt("2000-01-01", calendar=calendar)])
    if log is not None:
        log.append(f"created {name}")
    return field


def make_fields(names: list[str], log: list[str] | None = None) -> list[cf.Field]:
    """Create multiple fields with specific netcdf variable names."""
    fields = []
    for name in names:
        field = cf.example_field(0).copy()
        field.nc_set_variable(name)
        fields.append(field)
    if log is not None:
        log.append("created " + " ".join(names))
    return fields


def load_field(fields: cf.Field | list[cf.Field], log: list[str] | None = None) -> None:
    """Emulate a preprocessing step that takes a field / fields."""
    if not isinstance(fields, list):
        fields = [fields]
    if log is not None:
        log.append("loaded " + " ".join([field.nc_get_variable() for field in fields]))


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
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[],
            preprocessing=[
                (dummy_step, {"log": log}),
                (make_field, {"name": "p", "log": log}),
            ],
            tracker_inputs=[],
            config=config,
        )

        assert log == ["dummy step", "created p"]

    def test_tag_expansion(self, config) -> None:
        """Test the %ITER% and %BATCH% tags in any arguments get expanded."""
        tracker = DummyTracker()
        log: list[str] = []

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            input_files=[],
            preprocessing=[(dummy_step, {"comment": "%ITER%%BATCH%", "log": log})],
            tracker_inputs=[],
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
            interval=BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            input_files=[],
            preprocessing=[(dummy_step, {"comment": "%ITER%", "log": log})],
            tracker_inputs=[],
            config=config,
        )

        assert log == [0, 1]
        assert isinstance(log[0], int)

    def test_registry_single_field(self, config) -> None:
        """Test a field can be loaded from memory with "store" and "use"."""
        tracker = DummyTracker()
        log: list[str] = []

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[],
            preprocessing=[
                (make_field, {"name": "p"}, {"store": "field1"}),
                (load_field, {"log": log}, {"use": "field1"}),
            ],
            tracker_inputs=[],
            config=config,
        )

        assert log == ["loaded p"]

    def test_registry_name_reused(self, config) -> None:
        """Test storing with the same name replaces the previous field."""
        tracker = DummyTracker()
        log: list[str] = []

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[],
            preprocessing=[
                (make_field, {"name": "p1"}, {"store": "field1"}),
                (make_field, {"name": "p2"}, {"store": "field1"}),
                (load_field, {"log": log}, {"use": "field1"}),
            ],
            tracker_inputs=[],
            config=config,
        )

        assert log == ["loaded p2"]

    def test_registry_missing(self, config) -> None:
        """Test an error is raised when trying to load a field that is missing."""
        tracker = DummyTracker()

        with pytest.raises(KeyError, match=r"fields are not available.*: field1"):
            batching(
                tracker,
                interval=BATCH_INTERVAL,
                time_range=BATCH_TIME_RANGE,
                input_files=[],
                preprocessing=[(load_field, {}, {"use": "field1"})],
                config=config,
            )

    @pytest.mark.parametrize(
        "store, use, load_log, error",
        [
            pytest.param(["a", "b", "c"], "b", "u", None, id="store all, use one"),
            pytest.param(["a", "b", "c"], ["b", "a"], "u p", None, id="use multiple"),
            pytest.param("u", "u", "u", None, id="single var name"),
            pytest.param(["v", "u"], ["u", "v"], "u v", None, id="multiple var names"),
            pytest.param(
                ["a", "b", "c", "d"],
                None,
                None,
                (ValueError, "Number of fields to store exceeds the number provided"),
                id="too many to store",
            ),
            pytest.param(
                ["a"],
                None,
                None,
                (ValueError, "Fields with the following names are not provided .*: a"),
                id="invalid var name",
            ),
        ],
    )
    def test_registry_multiple_fields(
        self, config, store, use, load_log, error
    ) -> None:
        """Test the registry when using multiple fields."""
        tracker = DummyTracker()
        log: list[str] = []

        # Check for failure if expected
        with (
            nullcontext() if error is None else pytest.raises(error[0], match=error[1])
        ):
            batching(
                tracker,
                interval=BATCH_INTERVAL,
                time_range=BATCH_TIME_RANGE,
                input_files=[],
                preprocessing=[
                    (make_fields, {"names": ["p", "u", "v"]}, {"store": store}),
                    (load_field, {"log": log}, {"use": use}),
                ],
                tracker_inputs=[],
                config=config,
            )

            # Check the fields were stored / loaded as expected
            assert log == ["loaded " + load_log]


class TestBatching:
    """Tests for the batching utility."""

    def test_expand_batch_time_range(self) -> None:
        """Test batch buffer periods are applied and clamped to the overall range."""
        expanded_range = _expand_batch_time_range(
            (cf.dt("2000-02-01", calendar=None), cf.dt("2000-03-01", calendar=None)),
            (cf.dt("2000-01-01", calendar=None), cf.dt("2000-03-01", calendar=None)),
            {
                "buffer_period": timedelta(days=10),
                "start_buffer_period": timedelta(days=1),
            },
        )

        assert expanded_range == (
            cf.dt("2000-01-31", calendar=None),
            cf.dt("2000-03-01", calendar=None),
        )

    def test_calendar_from_input(self, config) -> None:
        """Test batch ranges use the calendar of the input file."""
        input_file = config["output_dir"] / "input.nc"
        field = make_field("input", calendar="360_day")
        cf.write(field, str(input_file))  # type: ignore[operator]

        calendar = _get_calendar([str(input_file)])
        time_range = (
            cf.dt("2000-01-01", calendar=calendar),
            cf.dt("2000-03-01", calendar=calendar),
        )
        ranges = _batch_time_ranges(time_range, "month")

        assert calendar == "360_day"
        assert ranges == [
            (
                cf.dt("2000-01-01", calendar=calendar),
                cf.dt("2000-02-01", calendar=calendar),
            ),
            (
                cf.dt("2000-02-01", calendar=calendar),
                cf.dt("2000-03-01", calendar=calendar),
            ),
        ]

    def test_batch_directories(self, config) -> None:
        """Test batching creates the expected batch directories."""
        tracker = DummyTracker()
        n_iter = 2

        batching(
            tracker,
            [],
            BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            tracker_inputs=[],
            config=config,
        )

        # Check the directories have been created
        batch_dirs = [config["output_dir"] / f"batch_{i}" for i in range(n_iter)]
        for batch_dir in batch_dirs:
            assert batch_dir.exists()

    def test_batch_directories_deletion(self, config) -> None:
        """Test batching deletes the batch directories when delete_batch_dirs=True."""
        tracker = DummyTracker()
        n_iter = 2

        del config["delete_batch_dirs"]  # Default is True

        batching(
            tracker,
            input_files=[],
            interval=BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            retrieve_data=lambda _, batch_dir: Path.touch(batch_dir / "file.txt"),
            tracker_inputs=[],
            config=config,
        )

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
            input_files=[],
            interval=BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            retrieve_data=retrieve_data,
            tracker_inputs=[],
            config=config,
        )

        assert log == [f"retrieve_data {i}" for i in range(n_iter)]

    def test_input_time_selection(self, config, mocker) -> None:
        """Test each batch selects its corresponding time range."""
        tracker = DummyTracker()

        # Create the input and mock the time selection function
        input_file1 = str(config["output_dir"] / "input1.nc")
        input_file2 = str(config["output_dir"] / "input2.nc")
        field1 = make_field("input1", calendar="360_day")
        cf.write(field1, input_file1)  # type: ignore[operator]
        cf.write(make_field("input2", calendar="360_day"), input_file2)  # type: ignore[operator]
        select_time_range = mocker.patch(
            "tctrack.utils.batching.select_time_range", return_value=field1
        )

        batching(
            tracker,
            [
                (input_file1, {"batch_file": None}),
                (input_file2, {"batch_file": None, "time_varying": False}),
            ],
            interval="month",
            time_range=("2000-01-01", "2000-03-01"),
            tracker_inputs=[],
            config=config,
        )

        # Check the time-varying input was subspaced in time for each batch
        assert select_time_range.call_args_list == [
            call(
                [input_file1],
                (
                    cf.dt("2000-01-01", calendar="360_day"),
                    cf.dt("2000-02-01", calendar="360_day"),
                ),
            ),
            call(
                [input_file1],
                (
                    cf.dt("2000-02-01", calendar="360_day"),
                    cf.dt("2000-03-01", calendar="360_day"),
                ),
            ),
        ]

    def test_input_time_selection_with_buffers(self, config, mocker) -> None:
        """Test batch input selection includes the requested time buffers."""
        tracker = DummyTracker()
        input_file = str(config["output_dir"] / "input.nc")
        field = make_field("input")
        cf.write(field, input_file)  # type: ignore[operator]
        select_time_range = mocker.patch(
            "tctrack.utils.batching.select_time_range", return_value=field
        )

        batching(
            tracker,
            [(input_file, {"batch_file": None})],
            interval="month",
            time_range=("2000-01-01", "2000-03-01"),
            tracker_inputs=[],
            config={
                **config,
                "buffer_period": timedelta(days=10),
            },
        )

        time_range1 = (
            cf.dt("2000-01-01", calendar=None),
            cf.dt("2000-02-11", calendar=None),
        )
        time_range2 = (
            cf.dt("2000-01-31", calendar=None),
            cf.dt("2000-03-01", calendar=None),
        )
        assert select_time_range.call_args_list == [
            call([input_file], time_range1),
            call([input_file], time_range2),
        ]

    def test_input_file_copied_to_batch(self, config) -> None:
        """Test an input file is copied to the batch directory."""
        tracker = DummyTracker()
        input_file = config["output_dir"] / "input.nc"
        cf.write(make_field("input"), str(input_file))  # type: ignore[operator]

        batching(
            tracker,
            input_files=str(input_file),
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            config=config,
        )

        batch_file = config["output_dir"] / "batch_0" / input_file.name
        assert batch_file.is_file()
        assert tracker.inputs == [[str(batch_file)]]

    def test_input_file_batch_name(self, config) -> None:
        """Test an input file can be given a filename in the batch directory."""
        tracker = DummyTracker()
        input_file = str(config["output_dir"] / "input.nc")
        cf.write(make_field("input"), input_file)  # type: ignore[operator]

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[(input_file, {"batch_file": "renamed.nc"})],
            config=config,
        )

        batch_file = config["output_dir"] / "batch_0" / "renamed.nc"
        assert batch_file.is_file()
        assert tracker.inputs == [[str(batch_file)]]

    def test_input_file_batch_name_none(self, config) -> None:
        """Test an input file can be excluded from the batch directory."""
        tracker = DummyTracker()
        input_file = str(config["output_dir"] / "input.nc")
        cf.write(make_field("input"), input_file)  # type: ignore[operator]

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[(input_file, {"batch_file": None})],
            tracker_inputs=[],
            config=config,
        )

        batch_dir = config["output_dir"] / "batch_0"
        assert not any(batch_dir.iterdir())
        assert tracker.inputs == [[]]

    def test_input_fields_are_stored_for_preprocessing(self, config) -> None:
        """Test input fields can be stored in the preprocessing registry."""
        tracker = DummyTracker()
        input_file = str(config["output_dir"] / "input.nc")
        cf.write(make_field("input"), input_file)  # type: ignore[operator]
        log: list[str] = []

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=[(input_file, {"store": "input", "batch_file": None})],
            preprocessing=[(load_field, {"log": log}, {"use": "input"})],
            tracker_inputs=[],
            config=config,
        )

        assert log == ["loaded input"]

    def test_input_file_wildcards(self, config) -> None:
        """Test wildcard inputs are correctly written to the batch directory."""
        tracker = DummyTracker()
        cf.write(make_field("input_1"), str(config["output_dir"] / "input_1.nc"))  # type: ignore[operator]
        cf.write(make_field("input_2"), str(config["output_dir"] / "input_2.nc"))  # type: ignore[operator]

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=BATCH_TIME_RANGE,
            input_files=str(config["output_dir"] / "input_*.nc"),
            config=config,
        )

        batch_file = config["output_dir"] / "batch_0" / "input_*.nc"
        assert batch_file.is_file()
        assert tracker.inputs == [[str(batch_file)]]

    def test_input_file_wildcards_missing(self, config) -> None:
        """Test batching fails if there are no input files that match a wildcard."""
        tracker = DummyTracker()

        with pytest.raises(FileNotFoundError, match="No files matched"):
            batching(
                tracker,
                input_files="input_*.nc",
                interval=BATCH_INTERVAL,
                time_range=BATCH_TIME_RANGE,
                config=config,
            )

    def test_tracker_input_output_files(self, config) -> None:
        """Test batching correctly sets the per-batch input and output files."""
        tracker = DummyTracker()
        n_iter = 2
        tracker_inputs = ["file1", "file2"]

        batching(
            tracker,
            [],
            BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            tracker_inputs=tracker_inputs,
            config=config,
        )

        # Check the input and output files were set correctly
        batch_dirs = [config["output_dir"] / f"batch_{i}" for i in range(n_iter)]
        assert tracker.inputs == [
            [str(batch_dir / file) for file in tracker_inputs]
            for batch_dir in batch_dirs
        ]
        assert tracker.outputs == [
            str(config["output_dir"] / f"tracks_{i}.nc") for i in range(n_iter)
        ]

    def test_tracker_input_wildcards(self, config) -> None:
        """Test batching expands input file wildcards in each batch directory."""
        tracker = DummyTracker()

        def retrieve_data(_: int, batch_dir: Path) -> None:
            (batch_dir / "input_2.nc").touch()
            (batch_dir / "input_1.nc").touch()

        batching(
            tracker,
            interval=BATCH_INTERVAL,
            time_range=TWO_BATCH_TIME_RANGE,
            input_files=[],
            retrieve_data=retrieve_data,
            tracker_inputs="input_*.nc",
            config=config,
        )

        assert tracker.inputs == [
            [
                str(config["output_dir"] / f"batch_{i}" / "input_1.nc"),
                str(config["output_dir"] / f"batch_{i}" / "input_2.nc"),
            ]
            for i in range(2)
        ]
