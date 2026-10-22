"""CLI tool to convert datasets into TCTrack format."""

import argparse
import logging
import sys

from tctrack.build_db.converters import Converter


def get_argument_parser() -> argparse.ArgumentParser:
    """Create argument parser for CLI arguments."""
    parser = argparse.ArgumentParser(
        description="CLI tool to convert datasets into TCTrack format."
    )
    parser.add_argument(
        "--source",
        "-s",
        required=True,
        help="The source for the input file. E.g. IBTrACS.",
    )
    parser.add_argument(
        "--start-date",
        "-d",
        help="The date (in YYYY-MM-DD format) from which the"
        " data needs to be filtered from the input file.",
    )
    parser.add_argument(
        "--end-date",
        "-e",
        help="The date (in YYYY-MM-DD format) to which the"
        " data needs to be filtered from the input file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        dest="output_path",
        help="Path to the output NetCDF file.",
    )
    parser.add_argument(
        "input_path", metavar="INPUT_PATH", help="Path to the input NetCDF file."
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s", force=True
    )

    parser = get_argument_parser()
    args = parser.parse_args(argv)

    Converter_ = None
    for cls in Converter.__subclasses__():
        if cls.tctrack_tracker.lower() == args.source.lower():
            Converter_ = cls
            break

    if not Converter_:
        msg = "Unsupported source..."
        raise SystemExit(msg)

    Converter_(**vars(args)).convert()  # type: ignore
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
