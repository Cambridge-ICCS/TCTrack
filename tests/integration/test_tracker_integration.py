
"""Real integration tests using tutorial scripts with actual data."""

import tempfile
from pathlib import Path

import pytest
import xarray as xr

import tctrack.tempest_extremes as te


@pytest.fixture
def subset_data_files():
    """Create subset of real data with only first 5 time steps."""
    data_dir = Path("/home/sg2147/rds/rds-inspire-tc-TqEGHMWTn8A/sg2147/data_processed")

    if not data_dir.exists():
        pytest.skip("Real data directory not found.")

    original_files = {
        "zg7h": data_dir
        / "zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
        "psl": data_dir
        / "psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
        "sfcWind": data_dir
        / "sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc",
        "orog": data_dir / "orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
    }

    if not all(f.exists() for f in original_files.values()):
        pytest.skip("Real data files not found.")

    temp_dir = tempfile.TemporaryDirectory()
    temp_path = Path(temp_dir.name)

    # Create subset files
    for _var_name, original_file in original_files.items():
        ds = xr.open_dataset(original_file)
        if "time" not in ds.dims:
            subset_ds = ds
        else:
            subset_ds = ds.isel(time=slice(0, 5))

        subset_file = temp_path / original_file.name
        subset_ds.to_netcdf(subset_file)
        ds.close()

    # Return the file paths directly
    input_files = [
        str(
            temp_path
            / "zg7h_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc"
        ),
        str(
            temp_path
            / "psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc"
        ),
        str(
            temp_path
            / "sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500801-19501030.nc"
        ),
        str(temp_path / "orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc"),
    ]

    yield input_files  # ← Return the paths

    temp_dir.cleanup()


class TestRunTrackerSubprocessIntegration:
    """Integration tests with real tracker and subset data."""

    def test_tempest_extremes_with_subset_data(self, subset_data_files):  # ← indented
        """Test full TempestExtremes pipeline with small subset of real data."""
        input_files = subset_data_files

        dn_params = te.TEDetectParameters(
            in_data=input_files,
            search_by_min="psl",
            merge_dist=6.0,
            closed_contours=[
                te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
            ],
            out_header=True,
            output_commands=[
                te.TEOutputCommand(var="psl", operator="min", dist=0.0),
                te.TEOutputCommand(var="sfcWind", operator="max", dist=2.0),
            ],
            output_dir="te_outputs_test",
        )

        sn_params = te.TEStitchParameters(
            caltype="360_day",
            max_sep=8.0,
            min_time="2",
        )

        tracker = te.TETracker(dn_params, sn_params)
        output_file = "tracks_test_subset.nc"

        tracker.run_tracker(output_file)

        assert Path(output_file).exists(), "Output file not created."
        assert Path(output_file).stat().st_size > 0, "Output file is empty."

        Path(output_file).unlink()
