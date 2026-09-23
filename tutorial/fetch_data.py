"""Script to download the ERA5 data used in the TCTrack tutorial."""

import os

import cdsapi


def fetch_data_cds(data_dir: str = "data/") -> None:
    """Download the ERA5 tutorial data from CDS.

    This requires registering for a CDS API key:
    https://cds.climate.copernicus.eu/how-to-api

    Then accept the licences for the data:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download#manage-licences

    Parameters
    ----------
    data_dir : str, optional
        Directory to download the data to. Created if it does not exist.
    """
    print("Downloading data. This may take several minutes.")
    os.makedirs(data_dir, exist_ok=True)
    client = cdsapi.Client()

    TIME_PERIOD = {
        "year": ["1950"],
        "month": ["09"],
        "day": ["01", "02", "03", "04", "05", "06", "07"],
        "time": ["00:00", "06:00", "12:00", "18:00"],
    }

    BASE_REQUEST = {
        "product_type": ["reanalysis"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [90, -180, -90, 180],
    }

    def retrieve(
        dataset: str, variables: list[str], target: str, request: dict
    ) -> None:
        """Download a dataset to a named file in the data directory."""
        full_request = {"variable": variables, **BASE_REQUEST, **TIME_PERIOD, **request}
        client.retrieve(dataset, full_request).download(
            target=os.path.join(data_dir, target)
        )

    # Pressure-level data:
    # - geopotential
    # - temperature
    # - vorticity
    pressure_levels_dataset = "reanalysis-era5-pressure-levels"

    retrieve(
        pressure_levels_dataset,
        ["geopotential"],
        "era5_z.nc",
        {"pressure_level": ["300", "500"]},
    )
    retrieve(
        pressure_levels_dataset,
        ["temperature"],
        "era5_t.nc",
        {"pressure_level": ["200", "250", "300", "400", "500"]},
    )
    retrieve(
        pressure_levels_dataset,
        ["vorticity"],
        "era5_vo.nc",
        {"pressure_level": ["850"]},
    )

    # Single-level data
    # - surface wind components
    # - mean sea-level pressure
    # - surface geopotential (time invariant so use 1 time)
    single_level_dataset = "reanalysis-era5-single-levels"

    # Note: wind speed is calculated from the wind components in
    # preprocess_data.py as the "10m_wind_speed" product is not available
    # from the MARS archive for this period.
    retrieve(
        single_level_dataset,
        [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "mean_sea_level_pressure",
        ],
        "era5_sfc.nc",
        {},
    )
    retrieve(
        single_level_dataset,
        ["geopotential"],
        "era5_sfc_z.nc",
        {"year": ["1940"], "month": ["01"], "day": ["01"], "time": ["00:00"]},
    )


if __name__ == "__main__":
    fetch_data_cds()
