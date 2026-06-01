"""Module with preprocessing functions that can be used with batching."""

import glob
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias, TypeVar

import cf
import numpy as np

FilePaths: TypeAlias = str | Sequence[str]


@dataclass(frozen=True)
class FieldSelect:
    """File name(s) plus the NetCDF variable name to select."""

    files: FilePaths
    var_name: str


FieldSource: TypeAlias = FilePaths | FieldSelect | cf.Field


def _expand_input_paths(paths: FilePaths) -> list[str]:
    """Expand wildcard input paths into a list of file paths."""
    path_list = [paths] if isinstance(paths, str) else list(paths)
    expanded_paths: list[str] = []

    for path in path_list:
        matches = sorted(glob.glob(path))
        if matches:
            expanded_paths.extend(matches)
            continue
        if glob.has_magic(path):
            msg = f"No files matched input pattern '{path}'."
            raise FileNotFoundError(msg)
        expanded_paths.append(path)

    return expanded_paths


T = TypeVar("T", cf.Field, list[cf.Field])


def _write_output(result: T, output_file: str | None) -> T:
    """Optionally write output before returning."""
    if output_file is not None:
        cf.write(result, output_file)  # type: ignore[operator]
    return result


def read_files(
    input_files: FilePaths,
    output_file: str | None = None,
    *,
    select: str | None = None,
) -> list[cf.Field]:
    """Read fields from separate files. Optionally write to a single output file."""
    fields = list(cf.read(_expand_input_paths(input_files), select=select))  # type: ignore[operator]
    return _write_output(fields, output_file)


def combine_time(
    input_files: FilePaths,
    time_bounds: tuple[str, str] | None = None,
    output_file: str | None = None,
) -> list[cf.Field]:
    """Combine files in time, with optional time bounds of [start, end)."""
    fields = read_files(input_files)

    if time_bounds is not None:
        time_interval = cf.wi(
            cf.dt(time_bounds[0]), cf.dt(time_bounds[1]), open_upper=True
        )
        fields = [field.subspace(T=time_interval) for field in fields]

    return _write_output(fields, output_file)


def separate_variables(
    input_files: FilePaths,
    output_files: dict[str, str],
) -> list[cf.Field]:
    """Split variables into separate files. The output keys are nc variable names."""
    fields = {field.nc_get_variable(): field for field in read_files(input_files)}

    for var_name, output_file in output_files.items():
        if var_name not in fields:
            msg = f"A variable to save ({var_name}) is not provided in the inputs."
            raise ValueError(msg)
        cf.write(fields[var_name], output_file)  # type: ignore[operator]

    return list(fields.values())


def _load_field(source: FieldSource) -> cf.Field:
    """Load a field from NetCDF file(s) ."""
    if isinstance(source, cf.Field):
        return source

    if isinstance(source, FieldSelect):
        fields = read_files(source.files, select=f"ncvar%{source.var_name}")
        if not fields:
            msg = f"No field with NetCDF variable name '{source.var_name}' was found."
            raise ValueError(msg)
        return fields[0]

    if isinstance(source, (str, Sequence)):
        fields = read_files(source)
        if len(fields) != 1:
            msg = (
                f"Expected one field from '{source}', but found {len(fields)}. "
                "Use FieldSelect(files, variable_name) to select a field."
            )
            raise ValueError(msg)
        return fields[0]

    msg = (
        "Invalid input type for the field source. "
        "Allowed types are cf.Field, FieldSelect, or string filepath(s)."
    )
    raise ValueError(msg)


def subsample_field(
    input_: FieldSource,
    subspace_kwargs: dict[str, Any],
    output_file: str | None = None,
    *,
    squeeze: bool = True,
) -> cf.Field:
    """Subsample a field using ``cf.Field.subspace``."""
    if not subspace_kwargs:
        msg = "At least one subspace selector must be provided to 'subspace_kwargs'."
        raise ValueError(msg)

    field = _load_field(input_)
    subset = field.subspace(**subspace_kwargs)
    if squeeze:
        subset.squeeze(inplace=True)
    return _write_output(subset, output_file)


def collapse_field(
    input_: FieldSource,
    method: str,
    axes: str | Sequence[str],
    output_file: str | None = None,
    *,
    squeeze: bool = True,
) -> cf.Field:
    """Collapse a field over one or more axes."""
    field = _load_field(input_)
    collapsed = field.collapse(method, axes=axes)
    if squeeze:
        collapsed.squeeze(inplace=True)
    return _write_output(collapsed, output_file)


def calculate_curl_xy(
    input_x: FieldSource,
    input_y: FieldSource,
    variable_name: str,
    variable_info: dict[str, str],
    output_file: str | None = None,
) -> cf.Field:
    """Calculate the curl of two fields which are x and y components of a vector."""
    field_x = _load_field(input_x)
    field_y = _load_field(input_y)

    curl = cf.curl_xy(field_x, field_y, radius="earth")
    curl.nc_set_variable(variable_name)
    for name, value in variable_info.items():
        if name == "units":
            curl.override_units(value, inplace=True)
        else:
            curl.set_property(name, value)
    return _write_output(curl, output_file)


def calculate_vorticity(
    input_u: FieldSource,
    input_v: FieldSource,
    output_file: str | None = None,
) -> cf.Field:
    """Calculate vorticity from colocated velocity fields."""
    return calculate_curl_xy(
        input_u,
        input_v,
        variable_name="vorticity",
        variable_info={
            "standard_name": "atmosphere_upward_absolute_vorticity",
            "units": "s-1",
        },
        output_file=output_file,
    )


def replace_fill_value(
    input_: FieldSource,
    fill_value: float,
    output_file: str | None = None,
) -> cf.Field:
    """Replace masked values in a field using ``cf.Field.filled``."""
    field = _load_field(input_)
    field.filled(fill_value=fill_value, inplace=True)
    return _write_output(field, output_file)


def set_netcdf_variable_name(
    input_: FieldSource,
    field_name: str,
    output_file: str | None = None,
    *,
    coord_names: dict[str, str] | None = None,
) -> cf.Field:
    """Set NetCDF variable names for a field and (optionally) its coordinates."""
    field = _load_field(input_)
    field.nc_set_variable(field_name)
    for coordinate, variable_name in (coord_names or {}).items():
        field.coordinate(coordinate).nc_set_variable(variable_name)
    return _write_output(field, output_file)


def regrid_to_field(
    input_: FieldSource,
    target: FieldSource | cf.Domain,
    output_file: str | None = None,
    *,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto the grid of another field / domain."""
    field = _load_field(input_)
    if not isinstance(target, cf.Domain):
        target = _load_field(target)

    regridded = field.regrids(target, method=method)
    regridded.nc_clear_dataset_chunksizes()  # Avoids a possible error when writing
    return _write_output(regridded, output_file)


def regrid_to_lat_lon(
    input_: FieldSource,
    latitude: np.ndarray,
    longitude: np.ndarray,
    output_file: str | None = None,
    *,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto a latitude-longitude grid."""
    field = _load_field(input_)

    domain = field.domain.copy()
    lat_coord = domain.dimension_coordinate("latitude")
    lat_coord.set_data(latitude, inplace=True)
    lat_coord.del_bounds()
    lon_coord = domain.dimension_coordinate("longitude")
    lon_coord.set_data(longitude, inplace=True)
    lon_coord.del_bounds()

    regridded = field.regrids((lat_coord, lon_coord), method=method)
    regridded.nc_clear_dataset_chunksizes()  # Avoids a possible error when writing
    return _write_output(regridded, output_file)


def gaussian_grid(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Create regular Gaussian latitude and longitude coordinates."""
    latitude = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(2 * n)[0]))
    longitude = np.arange(0.0, 360.0, 360.0 / (4 * n))
    return latitude, longitude


def regrid_to_gaussian(
    input_: FieldSource,
    n: int,
    output_file: str | None = None,
    *,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto a regular Gaussian grid with n lat points per hemisphere."""
    lat, lon = gaussian_grid(n)
    return regrid_to_lat_lon(input_, lat, lon, output_file=output_file, method=method)


__all__ = [
    "FieldSelect",
    "calculate_curl_xy",
    "calculate_vorticity",
    "collapse_field",
    "combine_time",
    "gaussian_grid",
    "read_files",
    "regrid_to_field",
    "regrid_to_gaussian",
    "regrid_to_lat_lon",
    "replace_fill_value",
    "separate_variables",
    "set_netcdf_variable_name",
    "subsample_field",
]
