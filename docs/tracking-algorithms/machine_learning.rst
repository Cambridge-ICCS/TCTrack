Machine Learning
================

The machine learning (ML) tracker detects tropical cyclones using a U-Net
segmentation model trained to classify each grid point of an input weather
field into one of five tropical cyclone lifecycle stages (or background).
The model and training pipeline are described in the reference
`cyclone-track-ml <https://github.com/robert-edwin-rouse/cyclone-track-ml>`_
repository.

Unlike the other trackers in TCTrack, the ML tracker does not wrap an
external tracking executable. Detection is performed directly in Python
using a `PyTorch <https://pytorch.org/>`_ TorchScript model, with weights
downloaded from a `HuggingFace Hub <https://huggingface.co/>`_ repository
(or supplied locally).

For full details of the Machine Learning API in TCTrack see the
:doc:`TCTrack Machine Learning API documentation <../api/machine_learning_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the machine_learning module to use references throughout this page.
.. py:module:: tctrack.machine_learning
   :no-index:

Model
-----

Using the ML tracker requires no external software to be installed; the
model weights are fetched automatically.

By default, model weights are downloaded from the HuggingFace Hub repository
given by :attr:`MLParameters.hf_repo_id` (which defaults to the pretrained
model at
`surbhigoel456/cyclone-TC-ML <https://huggingface.co/surbhigoel456/cyclone-TC-ML>`_)
the first time a :class:`MLTracker` is constructed, and cached locally by
``huggingface_hub`` for subsequent use. If the repository is private or
gated, an access token must be supplied either via
:attr:`~tctrack.core.ml_tracker.TCMLParameters.hf_token` or the ``HF_TOKEN``
environment variable.

Alternatively, a model file already available on disk can be used directly
by setting :attr:`~tctrack.core.ml_tracker.TCMLParameters.model_path`, in
which case no download takes place.

Usage
-----

Cyclone tracking with the ML tracker consists of two phases: detection,
which runs the model over each timestep of the input grid to locate
candidate storms, and stitching, which links candidates across timesteps
into trajectories - mirroring the detect/stitch split used by the other
trackers in TCTrack, even though the underlying algorithm is unrelated.

Usage of the ML tracker in TCTrack is through the ``machine_learning`` module.

This provides the :class:`MLTracker` class that stores algorithm parameters
and provides access to the detection and stitching methods. These are
configured through the :class:`MLParameters` and :class:`MLStitchParameters`
dataclasses.

Detection
^^^^^^^^^

In the following example we demonstrate the approach for detecting tropical
cyclones using the ML tracker. First, we set up the :meth:`~MLTracker.detect`
functionality to run on an input file. Full details of each input are given
in the :class:`MLParameters` documentation:

.. code-block:: python

    from tctrack.machine_learning import MLParameters, MLTracker

    ml_params = MLParameters(
        input_file="era5_2020.nc",
        threshold=0.5,
        merge_distance_deg=2.0,
    )

    # Initialize the tracker - this downloads (or loads) the model.
    tracker = MLTracker(ml_params)

    tracker.detect()

For each timestep, the model produces a per-pixel class probability map over
the input grid. Pixels whose highest-probability class is not background,
and whose confidence exceeds :attr:`~tctrack.core.ml_tracker.TCMLParameters.threshold`,
are grouped into discrete storm candidates and, where multiple candidates lie
within :attr:`~MLParameters.merge_distance_deg` of one another, merged into
the single strongest candidate. The resulting candidate locations are stored
internally, ready for stitching.

Stitching
^^^^^^^^^

These candidates may be stitched into trajectories using the
:meth:`~MLTracker.stitch` functionality. Stitching is based on locating,
at each timestep, the nearest unmatched candidate to each active track's
last known position. Tracks that go unmatched for too long are closed, and
short-lived tracks are discarded. Full details of each input are given in
the :class:`MLStitchParameters` documentation.

.. code-block:: python

    from tctrack.machine_learning import MLStitchParameters

    stitch_params = MLStitchParameters(
        stitch_max_distance_deg=3.0,
        stitch_max_gap=1,
        stitch_min_length=2,
    )
    tracker = MLTracker(ml_params, stitch_params)
    tracker.detect()
    trajectories = tracker.stitch()

Output
^^^^^^

After running detect and stitch to generate storm trajectories, the tracks
can be written to a NetCDF file fully compliant with the
`CF-Conventions <https://cfconventions.org/>`_ (specifically the
`trajectory data format <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#trajectory-data>`_)
using :meth:`~MLTracker.to_netcdf`:

.. code-block:: python

    tracker.to_netcdf("my_ml_cf_trajectories.nc")

This can be read using any NetCDF reading utility, though
`cf-python <https://ncas-cms.github.io/cf-python/>`_ will load it following
the `CF data model <https://ncas-cms.github.io/cf-python/#cf-data-model>`_.

The underlying detections - the individual storm locations found by
:meth:`~MLTracker.detect`, before stitching connects them into tracks - can
also be written out directly, as independent CF point data, using
:meth:`~MLTracker.detections_to_netcdf`:

.. code-block:: python

    tracker.detections_to_netcdf("my_ml_detections.nc")

This is useful for inspecting or evaluating detections independently of the
stitching step, since it retains every candidate regardless of whether it
was ultimately linked into a trajectory of the required minimum length.

Combined run
^^^^^^^^^^^^

The above examples demonstrate running :meth:`~MLTracker.detect`,
:meth:`~MLTracker.stitch`, and :meth:`~MLTracker.to_netcdf` separately.
However, it is likely that users will want to run these in succession, which
can be done using the :meth:`~MLTracker.run_tracker` method after defining
an :class:`MLTracker` object with appropriate :class:`MLParameters` and
:class:`MLStitchParameters`:

.. code-block:: python

    from tctrack.machine_learning import (
        MLParameters,
        MLStitchParameters,
        MLTracker,
    )

    ml_params = MLParameters(input_file="era5_2020.nc")
    stitch_params = MLStitchParameters()
    tracker = MLTracker(ml_params, stitch_params)

    tracker.run_tracker("my_ml_cf_trajectories.nc")

Input data
----------

The ML tracker requires a single CF-NetCDF input file containing the
variables the model was trained on, all on the same latitude/longitude grid:

* Pressure-level variables:

  * Each of :attr:`~MLParameters.pressure_variables` (by default relative
    humidity, air temperature, eastward wind, northward wind, and relative
    vorticity) at each of :attr:`~MLParameters.pressure_levels`
    (by default 1000, 750, and 500 hPa).
  * Selected from the input file by CF identity (``standard_name``), so the
    file must expose these variables with a ``Z`` (pressure) dimension
    covering the configured levels.

* Sea surface temperature:

  * Selected via :attr:`~MLParameters.sst_variable`. ERA5 sea surface
    temperature has no ``standard_name``, so by default this is selected by
    NetCDF variable name (``ncvar%sst``) instead.
  * Used together with :attr:`~MLParameters.t2m_variable` (2-metre
    temperature) to gap-fill sea surface temperature over land, and to
    derive a static land/sea mask from where sea surface temperature is
    undefined.

The file must also have a ``time`` coordinate, used both to order the
timesteps that :meth:`~MLTracker.detect` iterates over and to timestamp
output.

Channels are normalised before being passed to the model using per-channel
statistics computed over the training set. By default the statistics
bundled with TCTrack (matching the pretrained model) are used; a different
set of statistics may be supplied via
:attr:`~MLParameters.normalisation_stats_path` if using a model trained
elsewhere.

To extract these variables from an input dataset and write to individual
files, combine multiple files over times, or regrid variables to a
consistent grid see :doc:`../data/preprocessing_data`.
