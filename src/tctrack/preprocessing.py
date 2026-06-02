"""Module with preprocessing functions that can be used with batching."""

import glob
import importlib.util
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias, TypeVar

import cf
import numpy as np

ESMPY_AVAILABLE = importlib.util.find_spec("esmpy") is not None


def _require_esmpy() -> None:
    """Guard to ensure that esmpy is installed for regridding."""
    if not ESMPY_AVAILABLE:
        msg = (
            "Regridding requires esmpy to be installed. "
            "See the dependency installation documentation:\n"
            "https://tctrack.readthedocs.io/en/latest/getting-started/index.html#esmpy"
        )
        raise ImportError(msg)


FilePaths: TypeAlias = str | Sequence[str]


@dataclass(frozen=True)
class FieldSelect:
    """Class containing the file name(s) plus the NetCDF variable name to select.

    Necessary for choosing a variable from files which contain multiple.

    Parameters
    ----------
    files : str | Sequence[str]
        Input file path(s) to read from. ``glob`` pattern matching allowed.
    var_name : str
        NetCDF variable name to select from the input files.
    """

    files: str | Sequence[str]
    var_name: str


FieldSource: TypeAlias = str | Sequence[str] | FieldSelect | cf.Field
"""Type alias for the allowed sources for ``cf.Field`` arguments.

The ``cf.Field`` can be passed directly or using the path(s) CF-NetCDF file(s).
If the file(s) contain multiple fields then :class:`FieldSelect` should be used to
specify which to use.
"""


def _expand_input_paths(paths: str | Sequence[str]) -> list[str]:
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
    input_files: str | Sequence[str],
    *,
    output_file: str | None = None,
    select: str | None = None,
) -> list[cf.Field]:
    """Read fields from one or more files.

    Parameters
    ----------
    input_files : str | Sequence[str]
        Input file path(s) to read. ``glob`` pattern matching allowed.
    output_file : str | None, optional
        Output file to write the loaded fields to.
    select : str | None, optional
        Optional field selection for ``cf.read``.

    Returns
    -------
    list[cf.Field]
        The list of fields read from the input files.
    """
    fields = list(
        cf.read(  # type: ignore[operator]
            _expand_input_paths(input_files), select=select, netcdf_backend="netCDF4"
        )
    )
    return _write_output(fields, output_file)


def select_time_range(
    input_files: str | Sequence[str],
    time_bounds: tuple[str, str],
    *,
    output_file: str | None = None,
) -> list[cf.Field]:
    """Combine files in time and select a time range.

    Parameters
    ----------
    input_files : str | Sequence[str]
        Input file path(s) to combine. ``glob`` pattern matching allowed.
    time_bounds : tuple[str, str]
        Start and end datetime strings in format ``"YYYY-MM-DD[ HH:MM]"``.
        The end bound is open / exclusive.
    output_file : str | None, optional
        Output file to write the result to.

    Returns
    -------
    list[cf.Field]
        The list of combined fields.
    """
    fields = read_files(input_files)

    time_interval = cf.wi(cf.dt(time_bounds[0]), cf.dt(time_bounds[1]), open_upper=True)
    fields = [field.subspace(T=time_interval) for field in fields]

    return _write_output(fields, output_file)


def separate_variables(
    input_files: str | Sequence[str],
    output_files: dict[str, str],
) -> list[cf.Field]:
    """Split variables into separate files.

    Parameters
    ----------
    input_files : str | Sequence[str]
        Input file path(s) to read. ``glob`` pattern matching allowed.
    output_files : dict[str, str]
        Mapping from NetCDF variable name to output file path.

    Returns
    -------
    list[cf.Field]
        The list of fields read from the input files.
    """
    fields = {field.nc_get_variable(): field for field in read_files(input_files)}

    for var_name, output_file in output_files.items():
        if var_name not in fields:
            msg = f"A variable to save ({var_name}) is not provided in the inputs."
            raise ValueError(msg)
        cf.write(fields[var_name], output_file)  # type: ignore[operator]

    return list(fields.values())


def _load_field(source: FieldSource) -> cf.Field:
    """Load a single field from an in-memory field or file input."""
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
    *,
    output_file: str | None = None,
    squeeze: bool = False,
) -> cf.Field:
    """Subsample a field using ``cf.Field.subspace``.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing which field to load.
    subspace_kwargs : dict[str, Any]
        Keyword arguments passed to ``cf.Field.subspace``.
    output_file : str | None, optional
        Output file to write the result to.
    squeeze : bool, optional
        Whether to squeeze size-1 dimensions after subspacing. Default: ``False``.

    Returns
    -------
    cf.Field
        The subsampled field.
    """
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
    *,
    output_file: str | None = None,
    squeeze: bool = True,
) -> cf.Field:
    """Collapse a field over one or more axes.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing which field to load.
    method : str
        Collapse method passed to ``cf.Field.collapse``. E.g. ``"mean"``, ``"minimum"``.
    axes : str | Sequence[str]
        Axis or axes to collapse over.
    output_file : str | None, optional
        Output file to write the collapsed field to.
    squeeze : bool, optional
        Whether to squeeze size-1 dimensions after collapsing.

    Returns
    -------
    cf.Field
        The collapsed field.
    """
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
    *,
    output_file: str | None = None,
) -> cf.Field:
    """Calculate the curl of x and y vector components.

    Parameters
    ----------
    input_x : FieldSource
        Field for the x component.
    input_y : FieldSource
        Field for the y component.
    variable_name : str
        NetCDF variable name for the output field.
    variable_info : dict[str, str]
        Field properties to set on the output.
    output_file : str | None, optional
        Output file to write the curl field to.

    Returns
    -------
    cf.Field
        Curl field derived from the two inputs.
    """
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
    *,
    output_file: str | None = None,
) -> cf.Field:
    """Calculate vorticity from colocated velocity fields.

    Parameters
    ----------
    input_u : FieldSource
        Field for the eastward velocity component.
    input_v : FieldSource
        Field for the northward velocity component.
    output_file : str | None, optional
        Output file to write the vorticity field to.

    Returns
    -------
    cf.Field
        Vorticity field.
    """
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
    *,
    output_file: str | None = None,
) -> cf.Field:
    """Replace masked values in a field using ``cf.Field.filled``.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing which field to load.
    fill_value : float
        Value for missing data.
    output_file : str | None, optional
        Output file to write the updated field to.

    Returns
    -------
    cf.Field
        Field with fill value replaced.
    """
    field = _load_field(input_)
    field.filled(fill_value=fill_value, inplace=True)
    return _write_output(field, output_file)


def set_netcdf_variable_name(
    input_: FieldSource,
    field_name: str,
    *,
    output_file: str | None = None,
    coord_names: dict[str, str] | None = None,
) -> cf.Field:
    """Set NetCDF variable names for a field and, optionally, its coordinates.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing which field to load.
    field_name : str
        NetCDF variable name for the field.
    output_file : str | None, optional
        Output file to write the updated field to.
    coord_names : dict[str, str] | None, optional
        Optional updated NetCDF variable names for coordinates. Keys are the standard
        names.

    Returns
    -------
    cf.Field
        Field with updated NetCDF variable names.
    """
    field = _load_field(input_)
    field.nc_set_variable(field_name)
    for coordinate, variable_name in (coord_names or {}).items():
        field.coordinate(coordinate).nc_set_variable(variable_name)
    return _write_output(field, output_file)


def regrid_to_field(
    input_: FieldSource,
    target: FieldSource | cf.Domain,
    *,
    output_file: str | None = None,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto the grid of another field or domain.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing the field to regrid.
    target : FieldSource | cf.Domain
        Target field or domain that supplies the destination grid.
    output_file : str | None, optional
        Output file to write the regridded field to.
    method : str, optional
        Regridding method passed to ``cf.Field.regrids``.

    Returns
    -------
    cf.Field
        Regridded field.
    """
    _require_esmpy()
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
    *,
    output_file: str | None = None,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto a latitude-longitude grid.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing the field to regrid.
    latitude : np.ndarray
        Latitude coordinate values for the target grid.
    longitude : np.ndarray
        Longitude coordinate values for the target grid.
    output_file : str | None, optional
        Output file to write the regridded field to.
    method : str, optional
        Regridding method passed to ``cf.Field.regrids``.

    Returns
    -------
    cf.Field
        Regridded field on the requested latitude-longitude grid.
    """
    _require_esmpy()
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
    """Create regular Gaussian latitude and longitude coordinates.

    Parameters
    ----------
    n : int
        Number of latitude points per hemisphere.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Latitude and longitude coordinate arrays.
    """
    latitude = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(2 * n)[0]))
    longitude = np.arange(0.0, 360.0, 360.0 / (4 * n))
    return latitude, longitude


def regrid_to_gaussian(
    input_: FieldSource,
    n: int,
    *,
    output_file: str | None = None,
    method: str = "linear",
) -> cf.Field:
    """Regrid a field onto a regular Gaussian grid.

    Parameters
    ----------
    input_ : FieldSource
        A field, file path(s), or :class:`FieldSelect` describing the field to regrid.
    n : int
        Number of latitude points per hemisphere for the target gaussian grid.
    output_file : str | None, optional
        Output file to write the regridded field to.
    method : str, optional
        Regridding method passed to ``cf.Field.regrids``.

    Returns
    -------
    cf.Field
        Regridded field on the Gaussian grid.
    """
    lat, lon = gaussian_grid(n)
    return regrid_to_lat_lon(input_, lat, lon, output_file=output_file, method=method)


__all__ = [  # noqa: RUF022  # Prevent reorder for a more logical order in the api docs
    "read_files",
    "select_time_range",
    "separate_variables",
    "subsample_field",
    "collapse_field",
    "calculate_curl_xy",
    "calculate_vorticity",
    "replace_fill_value",
    "set_netcdf_variable_name",
    "regrid_to_field",
    "regrid_to_lat_lon",
    "gaussian_grid",
    "regrid_to_gaussian",
    "FieldSource",
    "FieldSelect",
]
