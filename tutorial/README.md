# TCTrack Tutorial

This directory contains Jupyter notebooks used in the TCTrack tutorial. For full
details please see the [online tutorial
documentation](https://tctrack.readthedocs.io/en/latest/getting-started/tutorial.html).

Contained in this directory are:

- `tempest_extremes.ipynb` - a notebook running the full TCTrack workflow with
  Tempest Extremes
- `tstorms.ipynb` - a notebook running the full TCTrack workflow with TSTORMS
- `install_tempest_extremes.sh` - to install Tempest Extremes locally
- `install_tstorms.sh` - to install TSTORMS locally
- `run_ml_tracker.py` - to run the machine-learning tracker on a small ERA5 sample (in
  `data/machine_learning`) to generate tracks and detections

The notebooks are self-contained and will download the example data they
require, but the tracking algorithms must be installed first by sourcing the
installation scripts above.
