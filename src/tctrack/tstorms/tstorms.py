"""Module providing classes that bind the TSTORMS code.

References
----------
- TSTORMS NOAA page: https://www.gfdl.noaa.gov/tstorms/
"""

import warnings
from dataclasses import dataclass

from tctrack.core import TCTrackerParameters


@dataclass(repr=False)
class TSTORMSParameters(TCTrackerParameters):
    """
    Dataclass containing values used in configuring the TSTORMS install.

    This class is intended to configure the install to that executables can be run
    correctly and data stored in a desired location.
    Configuration of the actual cyclone detection algorithm is done using
    :class:`DriverParameters` and :class:`TrajectoryParameters`.
    """

    tstorms_dir: str
    """
    Full path to the TSTORMS installation directory. This will be likely be a directory
    named `tropical_storms_pub/` and should contain `tstorms_driver/` and
    `trajectory_analysis/` subdirectories.
    """

    output_dir: str
    """
    Full path to trhe directory where TSTORMS outputs should be deposited.
    """


@dataclass(repr=False)
class DriverParameters(TCTrackerParameters):
    """
    Dataclass containing values used by the Driver operation of TSTORMS.

    Default values are set to match those in the TSTORMS source-code.

    Raises
    ------
    ValueError
        - If northern latitude bound is less than the southern latitude bound.
        - If any of the input file filenames are left as ``None``.
    UserWarning
        - If do_thickness is set to True as this has no effect.
    """

    u_in_file: str
    """
    Filename of the u (zonal) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the u field named
    ``u_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``u850``.
    """

    v_in_file: str
    """
    Filename of the v (meridional) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the v field named
    ``v_ref`` (if :attr:`use_sfc_wind` is ``True``) or ``v850``.
    """

    vort_in_file: str
    """
    Filename of the vorticity input file.
    This should be a NetCDF file containing a single lat-lon slice of the vorticity
    field named ``vort850`` at the 850 hPa level.
    """

    tm_in_file: str
    """
    Filename of the temperature input file.
    This should be a NetCDF file containing a single lat-lon slice of the mean
    temperature of the warm-core layer named ``tm``.
    """

    slp_in_file: str
    """
    Filename of the sea-level pressure input file.
    This should be a NetCDF file containing a single lat-lon slice of sea-level pressure
    named ``slp``.
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
    Critical warm core threshold to be exceeded to qualify as a candidate storm.
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
class TrajectoryParameters(TCTrackerParameters):
    """
    Dataclass containing values used by the Trajectory operation of TSTORMS.

    Default values are set to match those in the TSTORMS source-code.

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
    """
    Whether to apply filtering of trajectories.
    Filtering is based on landmask and lat-lon bounds to generate *_filt output files.
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
