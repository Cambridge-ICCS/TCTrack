"""Script to download the ERA5 data used in the TCTrack tutorial.

By default this downloads a pre-prepared bundle from the TCTrack GitHub
releases. Alternatively the data can be fetched directly from the Copernicus
Climate Data Store (CDS) using the ``fetch_data_cds`` function, which requires
registering for a CDS API key:
https://cds.climate.copernicus.eu/how-to-api
"""

import hashlib
import os
import tarfile
import urllib.request


def fetch_data() -> None:
    """Download the ERA5 tutorial data from the TCTrack releases."""
    data_url = (
        "https://github.com/Cambridge-ICCS/TCTrack/releases/download/"
        "data/tutorial-v1/tctrack-tutorial-data.tar.gz"
    )
    expected_sha256 = "d256b48d904be0fa31093453eab36ed1641ca645c1a78a1824b912b1318fef3c"

    print("Downloading data...")
    bundle = "tctrack-tutorial-data.tar.gz"
    urllib.request.urlretrieve(data_url, bundle)  # noqa: S310 - URL audit not required

    # Verify checksum
    digest = hashlib.sha256()
    with open(bundle, "rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != expected_sha256:
        msg = (
            "Checksum of the downloaded data does not match the expected value. "
            "Please report this at https://github.com/Cambridge-ICCS/TCTrack/issues"
        )
        raise RuntimeError(msg)

    # Extract data
    with tarfile.open(bundle) as tar:
        tar.extractall(filter="data")
    os.remove(bundle)

    print("Done.")


def fetch_data_cds() -> None:
    """Download the ERA5 tutorial data from CDS.

    This requires registering for a CDS API key:
    https://cds.climate.copernicus.eu/how-to-api

    Then accept the licences for the data:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download#manage-licences
    """
    import cdsapi  # noqa: PLC0415

    print("Downloading data. This may take several minutes.")
    data_dir = "data/"
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
    fetch_data()
