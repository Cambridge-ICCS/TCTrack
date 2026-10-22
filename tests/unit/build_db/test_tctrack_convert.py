"""Tests for tctrack.build_db.tctrack_convert."""

from typing import Any
from unittest.mock import MagicMock

import pytest

from tctrack.build_db.tctrack_convert import Converter, main


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


class TestMain:
    """Tests for main function."""

    @pytest.mark.parametrize("source", ["A", "C"])
    def test_valid_source(self, create_class, capsys, source):
        """Test by passing a valid source argument."""
        mock_subclasses = MagicMock()
        mock_subclasses.return_value = [
            create_class("A"),
            create_class("B"),
            create_class("C"),
        ]
        Converter.__subclasses__ = mock_subclasses

        main(["-s", source, "-o", "/path/to/output.nc", "/path/to/input.nc"])

        output = capsys.readouterr()
        assert output.out == f"MockClass: {source}\n"

    def test_invalid_source(self, create_class):
        """Test by passing an invalid source argument."""
        mock_subclasses = MagicMock()
        mock_subclasses.return_value = [
            create_class("A"),
            create_class("B"),
            create_class("C"),
        ]
        Converter.__subclasses__ = mock_subclasses

        with pytest.raises(SystemExit) as excinfo:
            main(["-s", "D", "-o", "/path/to/output.nc", "/path/to/input.nc"])

        assert str(excinfo.value) == "Unsupported source..."
