# Sample data for the machine-learning tracker

`era5_dikeledi_2025-01-10.nc` is a small ERA5 subset used by the machine-learning
tracker's tests and tutorial. It contains the input variables of
`tctrack.machine_learning.MLTracker` for the period in which Cyclone Dikeledi crossed
the northern Mozambique Channel.

| | |
|---|---|
| Time | 2025-01-10 00:00 to 2025-01-12 06:00 UTC, 6-hourly (10 timesteps) |
| Region | 22.25°S to 2.5°S, 41.0°E to 60.75°E (80 × 80 points, 0.25° spacing) |
| Pressure levels | 1000, 750 and 500 hPa |
| Pressure-level variables | relative humidity, air temperature, eastward wind, northward wind, relative vorticity |
| Surface variables | sea surface temperature (`sst`), 2 m temperature (`t2m`) |

Run on this file with the default model, the tracker finds Dikeledi at a threshold of
0.25 within about 1 degree of the IBTrACS positions, along with some spurious tracks.

## Source and licence

Generated from ERA5 reanalysis data by cropping to the region and period above. No other
processing was applied.

> Hersbach, H., et al. (2023): ERA5 hourly data on pressure levels from 1940 to
> present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS).
> DOI: [10.24381/cds.bd0915c6](https://doi.org/10.24381/cds.bd0915c6)
>
> Hersbach, H., et al. (2023): ERA5 hourly data on single levels from 1940 to
> present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS).
> DOI: [10.24381/cds.adbb2d47](https://doi.org/10.24381/cds.adbb2d47)

Contains modified Copernicus Climate Change Service information 2025. Neither the
European Commission nor ECMWF is responsible for any use that may be made of the
information it contains. ERA5 is distributed under the
[CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) licence.
