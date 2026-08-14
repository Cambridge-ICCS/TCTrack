"""Module containing parameter classes for the machine learning tracker."""

from collections.abc import Sequence
from dataclasses import dataclass

from tctrack.core import TCTrackerParameters
from tctrack.core.ml_tracker import TCMLParameters


@dataclass(repr=False)
class MLParameters(TCMLParameters):
    """Dataclass containing values for parameters used by MLTracker.

    Extends :class:`~tctrack.core.ml_tracker.TCMLParameters` with the input
    file path and sets the default HuggingFace repository for this model.

    See Also
    --------
    MLTracker : The tracker class that uses these parameters.
    """

    hf_repo_id: str = "surbhigoel456/cyclone-TC-ML"
    """HuggingFace Hub repository ID for the cyclone detection model."""

    pressure_variables: Sequence[str] = (
        "relative_humidity",
        "air_temperature",
        "eastward_wind",
        "northward_wind",
        "atmosphere_relative_vorticity",
    )
    """CF identities of the pressure-level input variables, in model channel order.

    Correspond to ``pressure_var_codes`` in the reference training pipeline
    (`cyclone-track-ml <https://github.com/robert-edwin-rouse/cyclone-track-ml>`__),
    but given as the actual CF ``standard_name`` each variable is converted to
    (via cfgrib) in ERA5 output, which differs from the CDS API request name
    used there for several variables (e.g. ``"temperature"`` in the reference
    pipeline is ``"air_temperature"`` here).
    """

    pressure_levels: Sequence[float] = (1000, 750, 500)
    """Pressure levels (hPa) to extract for each of :attr:`pressure_variables`."""

    sst_variable: str = "ncvar%sst"
    """CF identity of the sea surface temperature variable.

    ERA5 sea surface temperature has no ``standard_name`` (cfgrib leaves it
    ``"unknown"``), so it must be selected by NetCDF variable name instead.
    """

    t2m_variable: str = "ncvar%t2m"
    """CF identity of the 2-metre temperature, used to gap-fill SST over land."""

    normalisation_stats_path: str | None = None
    """Path to a ``.nc`` file with per-channel ``mean`` and ``range`` variables
    (each shape ``(channel,)``, ordered as in :meth:`MLTracker.preprocess`),
    used to normalise input as ``(x - mean) / range``.

    These must be the statistics computed over the model's *training* set;
    they cannot be recomputed from a single inference file. Matches the
    ``data/normalisation_parameters.nc`` file saved by ``compiler.py`` in the
    reference training pipeline. If ``None``, the statistics bundled with
    this package - computed over the default model's training set - are
    used instead.
    """

    merge_distance_deg: float = 2.0
    """Distance (degrees lat/lon) within which candidates in the same timestep
    are merged, keeping only the highest-scoring one.

    The model often splits a single storm into several adjacent clusters -
    either side of a lifecycle-class boundary, or where confidence dips below
    :attr:`~TCMLParameters.threshold` mid-storm - which would otherwise each
    become a separate track. Mirrors ``merge_dist`` in TempestExtremes'
    ``DetectNodes``, which resolves the same problem for its own candidates.
    Set to ``0`` to disable merging.
    """


@dataclass(repr=False)
class MLStitchParameters(TCTrackerParameters):
    """Dataclass containing values used by :meth:`MLTracker.stitch`.

    See Also
    --------
    MLTracker : The tracker class that uses these parameters.
    """

    max_distance_deg: float = 3.0
    """Maximum distance (degrees lat/lon) between a track's last point and a
    candidate in the next timestep for them to be linked as the same storm.

    Not part of the reference training pipeline, which only trains
    per-pixel classification and has no timestep-to-timestep linking step.
    The default assumes 6-hourly input (as used for training) and typical
    cyclone translation speeds; reduce for higher-frequency input.
    """

    max_gap: int = 1
    """Number of consecutive timesteps a track may go unmatched before it is
    closed, tolerating brief drops below :attr:`~TCMLParameters.threshold`.
    """

    min_length: int = 2
    """Minimum number of observations a track must have to be kept, filtering
    out single-timestep detections likely to be noise.
    """


def parameter_set(name: str) -> tuple[MLParameters, MLStitchParameters]:
    """Return a default machine learning parameter set.

    Parameters
    ----------
    name : str
        The parameter set to return. Valid values are:

        - ``"default"`` for the default parameter values used by the model.

    Returns
    -------
    tuple[MLParameters, MLStitchParameters]
        Detection and stitching parameters.

    Raises
    ------
    ValueError
        If ``name`` is not a supported parameter set.
    """
    parameter_sets = [
        (["default"], (MLParameters(), MLStitchParameters())),
    ]
    for aliases, parameter_set in parameter_sets:
        if name in aliases:
            return parameter_set

    available = ", ".join("/".join(aliases) for aliases, _ in parameter_sets)
    msg = f"Unknown parameter set: {name}. Available sets: {available}."
    raise ValueError(msg)
