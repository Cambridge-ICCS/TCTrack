"""Module containing parameter classes for the TSTORMS tracker."""

import warnings
from dataclasses import dataclass

from tctrack.core import TCTrackerParameters


@dataclass(repr=False)
class TSTORMSBaseParameters(TCTrackerParameters):
    """
    Dataclass containing values used in configuring the TSTORMS install.

    This class is intended to configure the install so that executables can be run
    correctly and data stored in a desired location.
    Configuration of the actual cyclone detection algorithm is done using
    :class:`TSTORMSDetectParameters` and :class:`TSTORMSStitchParameters`.
    """

    tstorms_dir: str
    """
    Full path to the TSTORMS installation directory. This will be likely be a directory
    named `tropical_storms_pub/` and should contain `tstorms_driver/` and
    `trajectory_analysis/` subdirectories.
    """

    output_dir: str
    """
    Full path to the directory where TSTORMS outputs should be deposited.
    """

    input_dir: str = ""
    """
    Full path to the directory where TSTORMS input files can be found.
    Defaults to empty string which will load from current directory.
    """


@dataclass(repr=False)
class TSTORMSDetectParameters(TCTrackerParameters):
    """
    Dataclass of values used by the Driver operation of TSTORMS for candidate detection.

    Default values are set to match those in the TSTORMS tstorms_driver source-code.

    Raises
    ------
    ValueError
        - If northern latitude bound is less than the southern latitude bound.
        - If any of the input file filenames are left as ``None``.
    UserWarning
        - If do_thickness is set to True as this has no effect.
    """

    u_in_file: str | None = None
    """
    Filename of the u (zonal) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the u field named
    ``u_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``u850``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    If unset, it will be automatically determined from the input file list.
    """

    v_in_file: str | None = None
    """
    Filename of the v (meridional) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the v field named
    ``v_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``v850``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    If unset, it will be automatically determined from the input file list.
    """

    vort_in_file: str | None = None
    """
    Filename of the vorticity input file.
    This should be a NetCDF file containing a single lat-lon slice of the vorticity
    field named ``vort850`` at the 850 hPa level.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    If unset, it will be automatically determined from the input file list.
    """

    tm_in_file: str | None = None
    """
    Filename of the temperature input file.
    This should be a NetCDF file containing a single lat-lon slice of the mean
    temperature of the warm-core layer named ``tm``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    If unset, it will be automatically determined from the input file list.
    """

    slp_in_file: str | None = None
    """
    Filename of the sea-level pressure input file.
    This should be a NetCDF file containing a single lat-lon slice of sea-level pressure
    named ``slp``.
    This should be the full path unless the location of input files was specified in
    :class:`TSTORMSBaseParameters`' ``input_dir``.
    If unset, it will be automatically determined from the input file list.
    """

    use_sfc_wind: bool = True
    """
    Whether to use surface winds (``True``), or winds at 850 hPa level.
    """

    vort_crit: float = 3.5e-5
    """
    Critical vorticity threshold [s-1] to be exceeded to qualify as a candidate storm.
    """

    tm_crit: float = 0.5
    """
    Critical warm core local maximum to be exceeded to qualify as a candidate storm.
    """

    thick_crit: float = 50.0
    """
    Critical thickness threshold to be exceeded to qualify as a candidate storm.
    Note that thickness calculations are not yet implemented in TSTORMS.
    """

    dist_crit: float = 4.0
    """
    Critical radius [degrees] within which vorticity, sea-level pressure, and other
    maxima/minima must lie within each other to qualify as a candidate storm.
    """

    lat_bound_n: float = 90.0
    """
    Northern latitude bound [degrees] below which storm detection should occur.
    """

    lat_bound_s: float = -90.0
    """
    Southern latitude bound [degrees] above which storm detection should occur.
    """

    do_spline: bool = False
    """
    Whether to use splines for detecting minma instead of gridpointwise search.
    """

    do_thickness: bool = False
    """
    Whether to use thickness of the 200-1000 hPa layer as a variable for detecting
    candidate storms. Note that this functionality is not yet implemented in TSTORMS.
    """

    def __post_init__(self):
        """Validate parameters."""
        if self.lat_bound_n < self.lat_bound_s:
            msg = (
                f"Northern latitude bound ({self.lat_bound_n}) is less than "
                f"Southern latitude bound ({self.lat_bound_s})."
            )
            raise ValueError(msg)
        if self.do_thickness:
            msg = (
                "`do_thickness` is set, but will have no effect as this feature is not "
                "implemented in TSTORMS."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)


@dataclass(repr=False)
class TSTORMSStitchParameters(TCTrackerParameters):
    """
    Dataclass containing values used by the stitching trajectory operation of TSTORMS.

    Default values are set to match those in the TSTORMS trajectory_analysis
    source-code.

    Raises
    ------
    ValueError
        If northern latitude bound is less than the southern latitude bound.
    UserWarning
        If do_thickness is set to True as this has no effect.
    """

    r_crit: float = 900.0
    """Maximum daily track length [km] between succcessive points in a trajectory."""

    wind_crit: float = 17.0
    """Critical wind speed [m/s] for trajectory calculations."""

    vort_crit: float = 3.5e-5
    """Critical vorticity threshold [s-1] for trajectory calculations."""

    tm_crit: float = 0.5
    """Critical warm core threshold for trajectory calculations."""

    thick_crit: float = 50.0
    """Critical thickness threshold for trajectory calculations."""

    n_day_crit: int = 2
    """Minimum number of days a trajectory must last to be valid."""

    do_filter: bool = True
    r"""
    Whether to apply filtering of trajectories.
    Filtering is based on landmask and lat-lon bounds to generate \*_filt output files.
    """

    lat_bound_n: float = 40.0
    """Northern latitude bound [degrees] for trajectory filtering."""

    lat_bound_s: float = -40.0
    """Southern latitude bound [degrees] for trajectory filtering."""

    do_spline: bool = False
    """
    Whether to use splines for trajectory calculations.
    Should match the value used in the Driver routine.
    If splines used then :attr:`twc_crit` and :attr:`thick_crit` will not be used in
    comparisons.
    """

    do_thickness: bool = False
    """
    Whether to use thickness of the 200-1000 hPa layer as a variable for detecting
    candidate storms. Note that this functionality is not yet implemented in TSTORMS.
    """

    def __post_init__(self):
        """Validate parameters."""
        if self.lat_bound_n < self.lat_bound_s:
            msg = (
                f"Northern latitude bound ({self.lat_bound_n}) is less than "
                f"Southern latitude bound ({self.lat_bound_s})."
            )
            raise ValueError(msg)
        if self.do_thickness:
            msg = (
                "`do_thickness` is set, but will have no effect as this feature is not "
                "implemented in TSTORMS."
            )
            warnings.warn(msg, category=UserWarning, stacklevel=3)


def parameter_set(
    name: str,
    tstorms_dir: str,
    output_dir: str,
) -> tuple[
    TSTORMSBaseParameters,
    TSTORMSDetectParameters,
    TSTORMSStitchParameters,
]:
    """Return a default TSTORMS parameter set.

    Parameters
    ----------
    name : str
        The parameter set to return. Valid values are:

        - ``"default"`` or ``"Vitart2001"`` for the original TSTORMS algorithm presented
          in [Vitart2001]_.
    tstorms_dir : str
        Path to the TSTORMS installation.
    output_dir : str
        Path where TSTORMS output is written.

    Returns
    -------
    tuple[TSTORMSBaseParameters, TSTORMSDetectParameters, TSTORMSStitchParameters]
        Base, detection, and stitching parameters.

    Raises
    ------
    ValueError
        If ``name`` is not a supported parameter set.
    """
    parameter_sets = [
        (
            ["default", "Vitart2001"],
            (
                TSTORMSBaseParameters(tstorms_dir=tstorms_dir, output_dir=output_dir),
                TSTORMSDetectParameters(
                    vort_crit=3.5e-5,
                    tm_crit=0.0,
                    thick_crit=50.0,
                    dist_crit=4.0,
                    lat_bound_n=70.0,
                    lat_bound_s=-70.0,
                    do_spline=False,
                    do_thickness=False,
                    use_sfc_wind=True,
                ),
                TSTORMSStitchParameters(
                    r_crit=900.0,
                    wind_crit=17.0,
                    vort_crit=3.5e-5,
                    tm_crit=0.0,
                    n_day_crit=2,
                    do_filter=True,
                    lat_bound_n=70.0,
                    lat_bound_s=-70.0,
                ),
            ),
        ),
    ]
    for aliases, parameter_set in parameter_sets:
        if name in aliases:
            return parameter_set

    available = ", ".join("/".join(aliases) for aliases, _ in parameter_sets)
    msg = f"Unknown parameter set: {name}. Available sets: {available}."
    raise ValueError(msg)
