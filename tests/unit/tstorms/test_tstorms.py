"""
Unit tests for the TSTORMS Python bindings.

Note that these do not require a TSTORMS installation to run and make use of
pytest-mock to mock the results of subprocess calls to the system.
"""

import pytest

from tctrack.tstorms import DriverParameters, TrajectoryParameters


@pytest.fixture
def tstorms_filenames() -> dict[str, str]:
    """Provide filenames for DriverParameters as non-optional/no default."""
    return {
        "u_in_file": "u.nc",
        "v_in_file": "v.nc",
        "vort_in_file": "vort.nc",
        "tm_in_file": "tm.nc",
        "slp_in_file": "slp.nc",
    }


class TestTSTORMSTypes:
    """Tests for the different Classes and Types defined for TSTORMS."""

    def test_driver_parameters_defaults(self, tstorms_filenames: dict) -> None:
        """Check the default values for DriverParameters."""
        params = DriverParameters(**tstorms_filenames)
        # Check values of all defaults except filenames
        assert params.use_sfc_wind is True
        assert params.vort_crit == 3.5e-5
        assert params.tm_crit == 0.5
        assert params.thick_crit == 50.0
        assert params.dist_crit == 4.0
        assert params.lat_bound_n == 90.0
        assert params.lat_bound_s == -90.0
        assert params.do_spline is False
        assert params.do_thickness is False

    def test_driver_parameters_lat_bounds_error(self, tstorms_filenames: dict) -> None:
        """Check ValueError when northern lat bound is less than southern lat bound."""
        with pytest.raises(
            ValueError,
            match=r"Northern latitude bound.*is less than.*Southern latitude bound",
        ):
            DriverParameters(**tstorms_filenames, lat_bound_n=-10.0, lat_bound_s=10.0)

    def test_driver_parameters_do_thickness_warning(
        self, tstorms_filenames: dict
    ) -> None:
        """Check UserWarning when do_thickness is set to True."""
        with pytest.warns(
            UserWarning, match="`do_thickness` is set, but will have no effect.*"
        ):
            DriverParameters(**tstorms_filenames, do_thickness=True)

    def test_trajectory_parameters_defaults(self) -> None:
        """Check the default values for TrajectoryParameters."""
        params = TrajectoryParameters()
        # Check values of all defaults
        assert params.r_crit == 900.0
        assert params.wind_crit == 17.0
        assert params.vort_crit == 3.5e-5
        assert params.tm_crit == 0.5
        assert params.thick_crit == 50.0
        assert params.n_day_crit == 2
        assert params.do_filter is True
        assert params.lat_bound_n == 40.0
        assert params.lat_bound_s == -40.0
        assert params.do_spline is False
        assert params.do_thickness is False

    def test_trajectory_parameters_lat_bounds_error(self) -> None:
        """Check ValueError when northern lat bound is less than southern lat bound."""
        with pytest.raises(
            ValueError,
            match=r"Northern latitude bound.*is less than.*Southern latitude bound",
        ):
            TrajectoryParameters(lat_bound_n=-10.0, lat_bound_s=10.0)

    def test_trajectory_parameters_do_thickness_warning(self) -> None:
        """Check UserWarning when do_thickness is set to True."""
        with pytest.warns(
            UserWarning, match="`do_thickness` is set, but will have no effect.*"
        ):
            TrajectoryParameters(do_thickness=True)
