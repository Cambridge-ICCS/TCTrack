"""Module providing classes that bind the TSTORMS code.

References
----------
- TSTORMS NOAA page: https://www.gfdl.noaa.gov/tstorms/
"""

import warnings
from dataclasses import dataclass

from tctrack.core import TCTrackerParameters


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

    u_in_file: str | None = None
    """
    Filename of the u (zonal) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the u field named
    `u_ref` (if `use_sfc_wind=True`) or `u850`.
    """

    v_in_file: str | None = None
    """
    Filename of the v (meridional) velocity input file.
    This should be a NetCDF file containing a single lat-lon slice of the v field named
    `v_ref` (if `use_sfc_wind=True`) or `v850`.
    """

    vort_in_file: str | None = None
    """
    Filename of the vorticity input file.
    This should be a NetCDF file containing a single lat-lon slice of the vorticity
    field named `vort850` at the 850 hPa level.
    """

    tm_in_file: str | None = None
    """
    Filename of the temperature input file.
    This should be a NetCDF file containing a single lat-lon slice of the mean
    temperature of the warm-core layer named `tm`.
    """

    slp_in_file: str | None = None
    """
    Filename of the sea-level pressure input file.
    This should be a NetCDF file containing a single lat-lon slice of sea-level pressure
    named `slp`.
    """

    use_sfc_wind: bool = True
    """
    Whether to use surface winds (`True`), or winds at 850 hPa level.
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
        if any(
            f is None
            for f in [
                self.u_in_file,
                self.v_in_file,
                self.vort_in_file,
                self.tm_in_file,
                self.slp_in_file,
            ]
        ):
            msg = "Input file not provided for one of u, v, vort, tm, slp."
            raise ValueError(msg)
