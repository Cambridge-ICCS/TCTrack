"""Unit tests for preprocessing helper utilities."""

from pathlib import Path

import cf
import numpy as np
import pytest

from tctrack.utils.preprocessing import (
    FieldSelect,
    _load_field,
    calculate_vorticity,
    collapse_field,
    combine_time,
    gaussian_grid,
    read_files,
    regrid_to_field,
    regrid_to_gaussian,
    replace_fill_value,
    separate_variables,
    set_netcdf_variable_name,
    subsample_field,
)


def make_field(var_name: str, time: str | None = None) -> cf.Field:
    """Create a small example field with an optional time value."""
    standard_names = {
        "mslp": "air_pressure_at_mean_sea_level",
        "u": "eastward_wind",
        "v": "northward_wind",
    }

    field = cf.example_field(0).copy()
    field.nc_set_variable(var_name)
    field.set_property("standard_name", standard_names[var_name])
    if time is not None:
        field.coordinate("T").set_data([cf.dt(time)], inplace=True)
    return field


def write_fields(fields: list[cf.Field], path: Path) -> str:
    """Write one or more fields to a file and return the file path."""
    cf.write(fields, str(path))  # type: ignore[operator]
    return str(path)


class TestPreprocessing:
    """Tests for preprocessing functions."""

    def test_read_files_combines_fields(self, tmp_path):
        """Test read_files accepts files with multiple fields."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )

        fields = read_files(input_file)

        assert [field.nc_get_variable() for field in fields] == ["mslp", "u"]

    def test_read_files_combines_time(self, tmp_path):
        """Test read_files accepts and combines fields split temporally over files."""
        input_files = [
            write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc"),
            write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc"),
        ]

        fields = read_files(input_files)

        assert len(fields) == 1
        assert fields[0].coordinate("T").size == 2

    def test_read_files_wildcard(self, tmp_path):
        """Test read_files correctly expands wildcard filepaths."""
        write_fields(make_field("mslp"), tmp_path / "a.nc")
        write_fields(make_field("u"), tmp_path / "b.nc")

        fields = read_files(str(tmp_path / "*.nc"))

        assert len(fields) == 2

    def test_read_files_wildcard_no_matches(self, tmp_path):
        """Test read_files fails for wildcard paths with no matches."""
        with pytest.raises(FileNotFoundError, match="No files matched input pattern"):
            read_files(str(tmp_path / "*.nc"))

    def test_read_files_output_file(self, tmp_path):
        """Test read_files writes the combined output when requested."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")
        output_file = tmp_path / "output.nc"

        fields = read_files(input_file, str(output_file))

        assert len(fields) == 1
        assert output_file.exists()
        assert cf.read(str(output_file))[0].nc_get_variable() == "mslp"

    def test_combine_time_bounds(self, tmp_path):
        """Test combine_time correctly selects data in time bounds."""
        input_files = [
            write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc"),
            write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc"),
        ]

        fields = combine_time(input_files, ("2000-01-01", "2000-01-02"))

        assert len(fields) == 1  # Same field (mslp)
        assert fields[0].coordinate("T").size == 1  # Upper bound is excluded

    def test_separate_varibles(self, tmp_path):
        """Test separate_variables correctly splits variables across multiple files."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )
        output_files = {
            "mslp": str(tmp_path / "mslp.nc"),
            "u": str(tmp_path / "u.nc"),
        }

        fields = separate_variables(input_file, output_files)

        assert [field.nc_get_variable() for field in fields] == ["mslp", "u"]
        assert read_files(output_files["mslp"])[0] == fields[0]
        assert read_files(output_files["u"])[0] == fields[1]

    def test_separate_varibles_invalid(self, tmp_path):
        """Test separate_variables fails if an invalid variable name is given."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")

        with pytest.raises(ValueError, match=r"A variable to save \(invalid\)"):
            separate_variables(input_file, {"invalid": str(tmp_path / "output.nc")})

    def test_load_field_accepts_fields(self):
        """Test _load_field accepts in-memory fields."""
        field = make_field("mslp")

        assert _load_field(field) is field

    def test_load_field_accepts_files(self, tmp_path):
        """Test _load_field accepts a filename / list of files."""
        file_a = write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc")
        file_b = write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc")

        field = _load_field([file_a, file_b])

        assert field.nc_get_variable() == "mslp"
        assert field.coordinate("T").size == 2

    def test_load_field_rejects_multifield_files(self, tmp_path):
        """Test _load_field rejects files with multiple fields."""
        input_file = write_fields(
            [make_field("u"), make_field("v")],
            tmp_path / "input.nc",
        )

        with pytest.raises(
            ValueError,
            match=r"Use FieldSelect\(files, variable_name\) to select a field",
        ):
            _load_field(input_file)

    def test_load_field_accepts_field_select(self, tmp_path):
        """Test _load_field selects a field when given a FieldSelect object."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )

        field = _load_field(FieldSelect(input_file, "u"))

        assert field.nc_get_variable() == "u"

    def test_load_field_rejects_missing_field_select(self, tmp_path):
        """Test _load_field fails when a selected variable is missing."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")

        with pytest.raises(
            ValueError,
            match="No field with NetCDF variable name 'invalid' was found",
        ):
            _load_field(FieldSelect(input_file, "invalid"))

    def test_subsample_field(self):
        """Test subsample_field works correctly."""
        field = make_field("mslp")

        subset = subsample_field(field, {"X": slice(2)})

        assert subset.shape == (field.shape[0], 2)

    def test_subsample_field_rejects_empty_subspace_kwargs(self):
        """Test subsample_field fails when no selectors are provided."""
        with pytest.raises(ValueError, match="At least one subspace selector"):
            subsample_field(make_field("mslp"), {})

    def test_collapse_field(self):
        """Test collapse_field works correctly."""
        field = make_field("mslp")

        collapsed = collapse_field(field, "mean", "X")

        assert collapsed.shape == (field.axis_size("latitude"),)
        assert np.allclose(
            collapsed.array, field.collapse("mean", axes="X", squeeze=True).array
        )

    def test_calculate_vorticity(self):
        """Test calculate_vorticity works correctly."""
        field_u = make_field("u")
        field_v = make_field("v")

        vorticity = calculate_vorticity(field_u, field_v)

        assert vorticity.nc_get_variable() == "vorticity"
        assert (
            vorticity.get_property("standard_name")
            == "atmosphere_upward_absolute_vorticity"
        )
        assert vorticity.get_property("units") == "s-1"

    def test_replace_fill_value(self):
        """Test replace_fill_value works correctly."""
        field = make_field("mslp")
        field[0, 0] = cf.masked

        filled = replace_fill_value(field, -1.0)

        assert filled.array[0, 0] == pytest.approx(-1.0)

    def test_set_netcdf_variable_name(self):
        """Test set_netcdf_variable_name works correctly."""
        field = make_field("mslp")

        renamed = set_netcdf_variable_name(
            field,
            "pressure",
            coord_names={"X": "longitude", "Y": "latitude"},
        )

        assert renamed.nc_get_variable() == "pressure"
        assert renamed.coordinate("X").nc_get_variable() == "longitude"
        assert renamed.coordinate("Y").nc_get_variable() == "latitude"

    def test_regrid_to_field(self):
        """Test regrid_to_field works correctly."""
        target = make_field("u")

        regridded = regrid_to_field(make_field("v"), target)

        assert regridded.shape == target.shape

    def test_gaussian_grid(self):
        """Test Gaussian grid helper returns the expected coordinate sizes."""
        latitude, longitude = gaussian_grid(4)

        assert len(latitude) == 8
        assert len(longitude) == 16
        assert longitude[0] == pytest.approx(0.0)
        assert longitude[-1] == pytest.approx(337.5)
        assert latitude[0] == pytest.approx(-latitude[-1])

    def test_regrid_to_gaussian(self):
        """Test regrid_to_gaussian works correctly."""
        field = make_field("mslp")

        regridded = regrid_to_gaussian(field, 4)

        assert regridded.shape == (8, 16)
