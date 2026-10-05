"""Script to run the machine-learning tracker as part of the TCTrack tutorial."""

from tctrack import machine_learning

# Small ERA5 sample (10-12 January 2025, 6-hourly) covering the northern Mozambique
# Channel, where Cyclone Dikeledi crossed. See data/machine_learning/README.md.
# The model is downloaded from the HuggingFace Hub on first use. If the repository
# requires a token, set the HF_TOKEN environment variable.
ml_params = machine_learning.MLParameters(
    input_file="../data/machine_learning/era5_dikeledi_2025-01-10.nc",
    # Lower than the default of 0.5, so that the model follows the whole of the storm
    # track in this sample, at the cost of some spurious detections.
    threshold=0.25,
)

ml_tracker = machine_learning.MLTracker(ml_params)
ml_tracker.run_tracker("tracks_ml.nc")
ml_tracker.detections_to_netcdf("detections_ml.nc")
