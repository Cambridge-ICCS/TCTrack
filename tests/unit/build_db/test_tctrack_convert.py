"""Tests for tctrack.build_db.tctrack_convert."""

# ruff: noqa: SLF001
from typing import Any
from unittest.mock import MagicMock

import pytest

from tctrack.build_db.tctrack_convert import Converter, get_argument_parser, main


@pytest.fixture
def create_class():
    """Create a mock class."""

    def _create_class(tracker):
        class MockClass:
            tctrack_tracker = tracker

            def __init__(self, **kwargs: Any) -> None:
                pass

            def convert(self) -> None:
                print(f"MockClass: {tracker}")

        return MockClass

    return _create_class


class TestGetArgumentParser:
    """Tests for get_argument_parser function."""

    def test_parser_description(self):
        """Test the description provided for the parser."""
        parser = get_argument_parser()

        assert parser.description == "CLI tool to convert datasets into TCTrack format."

    @pytest.mark.parametrize(
        "is_required,expected_args",
        [
            (True, ["tracker", "output_path", "input_path"]),
            (False, ["help", "start_date", "end_date"]),
        ],
    )
    def test_arguments(self, is_required, expected_args):
        """Test the required and optional arguments."""
        parser = get_argument_parser()

        actual_args = []
        for action in parser._actions:
            if action.required == is_required:
                actual_args.append(action.dest)

        assert actual_args == expected_args

    def test_argument_help_text(self):
        """Test the help text for each argument."""
        expected_help_text = {
            "help": "show this help message and exit",
            "tracker": "The tracker/source for the input file. Ex. IBTrACS.",
            "start_date": "The date (in YYYY-MM-DD format) from which the data "
            "needs to be filtered from the input file.",
            "end_date": "The date (in YYYY-MM-DD format) to which the data "
            "needs to be filtered from the input file.",
            "output_path": "Path to the output NetCDF file.",
            "input_path": "Path to the input NetCDF file.",
        }
        parser = get_argument_parser()

        actual_help_text = {}
        for action in parser._actions:
            actual_help_text[action.dest] = action.help

        assert actual_help_text == expected_help_text


class TestMain:
    """Tests for main function."""

    @pytest.mark.parametrize("tracker", ["A", "C"])
    def test_valid_tracker(self, create_class, capsys, tracker):
        """Test by passing a valid tracker argument."""
        mock_subclasses = MagicMock()
        mock_subclasses.return_value = [
            create_class("A"),
            create_class("B"),
            create_class("C"),
        ]
        Converter.__subclasses__ = mock_subclasses

        main(["-t", tracker, "-o", "/path/to/output.nc", "/path/to/input.nc"])

        output = capsys.readouterr()
        assert output.out == f"MockClass: {tracker}\n"

    def test_invalid_tracker(self, create_class):
        """Test by passing an invalid tracker argument."""
        mock_subclasses = MagicMock()
        mock_subclasses.return_value = [
            create_class("A"),
            create_class("B"),
            create_class("C"),
        ]
        Converter.__subclasses__ = mock_subclasses

        with pytest.raises(SystemExit) as excinfo:
            main(["-t", "D", "-o", "/path/to/output.nc", "/path/to/input.nc"])

        assert str(excinfo.value) == "Unsupported tracker..."
