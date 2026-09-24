"""Unit tests for preprocessing functions."""

from pathlib import Path

import cf
import numpy as np
import pytest

import tctrack.preprocessing
from tctrack.preprocessing import (
    _load_field,
    calculate_vorticity,
    calculate_wind_speed,
    collapse_field,
    flip_axis,
    gaussian_grid,
    multiply_field,
    read_files,
    regrid_to_field,
    regrid_to_gaussian,
    replace_fill_value,
    select_time_range,
    separate_variables,
    set_netcdf_info,
    set_time_units,
    squeeze_field,
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


def write_fields(fields: cf.Field | list[cf.Field], path: Path) -> str:
    """Write one or more fields to a file and return the file path."""
    cf.write(fields, str(path))  # type: ignore[operator]
    return str(path)


class TestPreprocessing:
    """Tests for preprocessing functions."""

    def test_read_files_combines_fields(self, tmp_path: Path):
        """Test read_files accepts files with multiple fields."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )

        fields = read_files(input_file)

        assert [field.nc_get_variable() for field in fields] == ["mslp", "u"]

    def test_read_files_combines_time(self, tmp_path: Path):
        """Test read_files accepts and combines fields split temporally over files."""
        input_files = [
            write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc"),
            write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc"),
        ]

        fields = read_files(input_files)

        assert len(fields) == 1
        assert fields[0].coordinate("T").size == 2

    def test_read_files_wildcard(self, tmp_path: Path):
        """Test read_files correctly expands wildcard filepaths."""
        write_fields(make_field("mslp"), tmp_path / "a.nc")
        write_fields(make_field("u"), tmp_path / "b.nc")

        fields = read_files(str(tmp_path / "*.nc"))

        assert len(fields) == 2

    def test_read_files_wildcard_no_matches(self, tmp_path: Path):
        """Test read_files fails for wildcard paths with no matches."""
        with pytest.raises(FileNotFoundError, match="No files matched input pattern"):
            read_files(str(tmp_path / "*.nc"))

    def test_read_files_output_file(self, tmp_path: Path):
        """Test read_files writes the combined output when requested."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")
        output_file = tmp_path / "output.nc"

        fields = read_files(input_file, output_file=str(output_file))

        assert len(fields) == 1
        assert output_file.exists()
        assert cf.read(str(output_file))[0].nc_get_variable() == "mslp"  # type: ignore[operator]

    def test_select_time_range_bounds(self, tmp_path: Path):
        """Test select_time_range correctly selects data in time bounds."""
        input_files = [
            write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc"),
            write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc"),
            write_fields(make_field("u", "2000-01-01"), tmp_path / "c.nc"),
            write_fields(make_field("u", "2000-01-02"), tmp_path / "d.nc"),
        ]

        fields = select_time_range(input_files, ("2000-01-01", "2000-01-02"))

        # Same fields (mslp, u)
        assert len(fields) == 2
        assert fields[0].nc_get_variable() == "mslp"
        assert fields[1].nc_get_variable() == "u"
        # Upper bound is excluded
        assert fields[0].coordinate("T").size == 1
        assert fields[1].coordinate("T").size == 1

    def test_select_time_range_squeeze(self, tmp_path: Path):
        """Test select_time_range squeezes size-1 list outputs."""
        input_files = [
            write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc"),
            write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc"),
        ]

        output = select_time_range(input_files, ("2000-01-01", "2000-01-02"))

        # Check it returns just the mslp field, not a list
        assert isinstance(output, cf.Field)
        assert output.nc_get_variable() == "mslp"
        assert output.coordinate("T").size == 1

    def test_separate_varibles(self, tmp_path: Path):
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

    def test_separate_varibles_invalid(self, tmp_path: Path):
        """Test separate_variables fails if an invalid variable name is given."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")

        with pytest.raises(ValueError, match=r"A variable to save \(invalid\)"):
            separate_variables(input_file, {"invalid": str(tmp_path / "output.nc")})

    def test_separate_varibles_return_order(self, tmp_path: Path):
        """Test separate_variables returns fields in the requested order."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )

        fields = separate_variables(input_file, {}, return_order=["u", "mslp"])

        assert [field.nc_get_variable() for field in fields] == ["u", "mslp"]

    def test_separate_varibles_invalid_return_order(self, tmp_path: Path):
        """Test separate_variables fails if an invalid return order name is given."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")

        with pytest.raises(ValueError, match=r"A variable to return \(invalid\)"):
            separate_variables(input_file, {}, return_order=["invalid"])

    def test_load_field_accepts_fields(self):
        """Test _load_field accepts in-memory fields."""
        field = make_field("mslp")

        assert _load_field(field) is field

    def test_load_field_accepts_files(self, tmp_path: Path):
        """Test _load_field accepts a filename / list of files."""
        file_a = write_fields(make_field("mslp", "2000-01-01"), tmp_path / "a.nc")
        file_b = write_fields(make_field("mslp", "2000-01-02"), tmp_path / "b.nc")

        field = _load_field([file_a, file_b])

        assert field.nc_get_variable() == "mslp"
        assert field.coordinate("T").size == 2

    def test_load_field_rejects_multifield_files(self, tmp_path: Path):
        """Test _load_field rejects files with multiple fields."""
        input_file = write_fields(
            [make_field("u"), make_field("v")],
            tmp_path / "input.nc",
        )

        with pytest.raises(ValueError, match=r"Use {.*} to select a field"):
            _load_field(input_file)

    def test_load_field_accepts_field_select(self, tmp_path: Path):
        """Test _load_field selects a field when given a FieldSelect dictionary."""
        input_file = write_fields(
            [make_field("mslp"), make_field("u")],
            tmp_path / "input.nc",
        )

        field = _load_field({"files": input_file, "var_name": "u"})

        assert field.nc_get_variable() == "u"

    def test_load_field_rejects_missing_field_select(self, tmp_path: Path):
        """Test _load_field fails when a selected variable is missing."""
        input_file = write_fields(make_field("mslp"), tmp_path / "input.nc")

        with pytest.raises(
            ValueError,
            match="No field with NetCDF variable name 'invalid' was found",
        ):
            _load_field({"files": input_file, "var_name": "invalid"})

    def test_squeeze_field(self):
        """Test squeeze_field removes size-1 dimensions."""
        field = make_field("mslp")
        field.insert_dimension("T", 1, inplace=True)

        squeezed = squeeze_field(field)

        assert squeezed.shape == field.squeeze().shape
        assert "T" not in squeezed.domain_axes()

    def test_flip_axis(self):
        """Test flip_axis reverses the axis and its coordinate values."""
        field = make_field("mslp")

        before = field.copy()
        flipped = flip_axis(field, "Y")

        assert flipped.coordinate("Y").array[0] == before.coordinate("Y").array[-1]
        assert np.allclose(
            flipped.coordinate("Y").array, before.coordinate("Y").array[::-1]
        )

    def test_flip_axis_stored_direction(self):
        """Test flip_axis updates the stored_direction property if present."""
        field = make_field("mslp")
        field.coordinate("Y").set_property("stored_direction", "increasing")
        # Make sure it is increasing
        assert field.coordinate("Y").array[0] < field.coordinate("Y").array[-1]

        flipped = flip_axis(field, "Y")

        assert flipped.coordinate("Y").get_property("stored_direction") == "decreasing"

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
            == "atmosphere_upward_relative_vorticity"
        )
        assert vorticity.get_property("units") == "s-1"

    def test_calculate_wind_speed(self):
        """Test calculate_wind_speed works correctly."""
        field_u = make_field("u")
        field_v = make_field("v")

        wind_speed = calculate_wind_speed(field_u, field_v)

        assert wind_speed.nc_get_variable() == "wind_speed"
        assert wind_speed.get_property("standard_name") == "wind_speed"
        expected = np.hypot(field_u.array, field_v.array)
        assert np.allclose(wind_speed.array, expected)

    def test_multiply_field(self):
        """Test multiply_field works correctly."""
        field = make_field("mslp")

        scaled = multiply_field(field, 0.5)

        assert np.allclose(scaled.array, field.array * 0.5)

    def test_set_time_units(self):
        """Test set_time_units converts values without changing datetimes."""
        field = make_field("mslp", "2000-01-01")

        updated = set_time_units(field, "days since 1999-12-31")

        assert str(updated.coordinate("T").Units) == "days since 1999-12-31"
        assert updated.coordinate("T").array[0] == 1
        assert updated.coordinate("T").datetime_array[0] == cf.dt(
            2000, 1, 1, calendar="gregorian"
        )

    def test_replace_fill_value(self):
        """Test replace_fill_value works correctly."""
        field = make_field("mslp")
        field[0, 0] = cf.masked

        filled = replace_fill_value(field, -1.0)

        assert filled.array[0, 0] == pytest.approx(-1.0)

    def test_set_netcdf_info(self):
        """Test set_netcdf_info sets names and properties correctly."""
        field = make_field("mslp")

        updated = set_netcdf_info(
            field,
            nc_name="pressure",
            properties={"standard_name": "air_pressure"},
            coord_nc_names={"X": "longitude", "Y": "latitude"},
        )

        assert updated.nc_get_variable() == "pressure"
        assert updated.get_property("standard_name") == "air_pressure"
        assert updated.coordinate("X").nc_get_variable() == "longitude"
        assert updated.coordinate("Y").nc_get_variable() == "latitude"

    def test_set_netcdf_info_optional_name(self):
        """Test set_netcdf_info leaves the name unchanged when nc_name is None."""
        field = make_field("mslp")

        updated = set_netcdf_info(field)

        assert updated.nc_get_variable() == "mslp"

    def test_set_netcdf_info_unlimited_axis(self, tmp_path: Path):
        """Test set_netcdf_info sets and removes unlimited axis status."""
        field = make_field("mslp")

        updated = set_netcdf_info(field, axis_unlimited="T")
        assert updated.domain_axis("T").nc_is_unlimited() is True

        output = str(tmp_path / "unlimited.nc")
        updated = set_netcdf_info(
            updated, axis_unlimited=("T", False), output_file=output
        )
        assert updated.domain_axis("T").nc_is_unlimited() is False

        # Check the status is applied in the written file
        fields = cf.read(output)  # type: ignore[operator]
        assert fields[0].domain_axis("T").nc_is_unlimited() is False

    def test_regrid_esmpy_guard(self, monkeypatch):
        """Test regridding fails clearly when esmpy is unavailable."""
        monkeypatch.setattr(tctrack.preprocessing, "ESMPY_AVAILABLE", False)

        with pytest.raises(ImportError, match="Regridding requires esmpy"):
            regrid_to_field(make_field("v"), make_field("u"))

    @pytest.mark.skipif(
        not tctrack.preprocessing.ESMPY_AVAILABLE,
        reason="esmpy is not pip-installable so this currently fails in CI",
    )
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

    @pytest.mark.skipif(
        not tctrack.preprocessing.ESMPY_AVAILABLE,
        reason="esmpy is not pip-installable so this currently fails in CI",
    )
    def test_regrid_to_gaussian(self):
        """Test regrid_to_gaussian works correctly."""
        field = make_field("mslp")

        regridded = regrid_to_gaussian(field, 4)

        assert regridded.shape == (8, 16)
