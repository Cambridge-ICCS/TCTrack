Machine Learning
================

The machine learning (ML) tracker detects tropical cyclones using a U-Net
segmentation model, which classifies each grid point of an input weather field into
one of five tropical cyclone lifecycle stages or background. The detections found at
each timestep are then linked into trajectories.

The model and its training pipeline are described in the reference `cyclone-track-ml
<https://github.com/robert-edwin-rouse/cyclone-track-ml>`_ repository, and the
pretrained model is available from the `HuggingFace Hub
<https://huggingface.co/surbhigoel456/cyclone-TC-ML>`_.

Unlike the other trackers in TCTrack, the ML tracker does not wrap an external
executable. The model is a `PyTorch <https://pytorch.org/>`_ TorchScript file that is
run directly from Python.

For full details of the Machine Learning API in TCTrack see the
:doc:`TCTrack Machine Learning API documentation <../api/machine_learning_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the machine_learning module to use references throughout this page.
.. py:module:: tctrack.machine_learning
   :no-index:

Installation
------------

Using the ML tracker does not require any external software to be installed. The
dependencies, including PyTorch and ``huggingface_hub``, are installed with TCTrack.

The model weights are downloaded from the HuggingFace Hub repository given by
:attr:`MLParameters.hf_repo_id` (which defaults to the pretrained model at
`surbhigoel456/cyclone-TC-ML <https://huggingface.co/surbhigoel456/cyclone-TC-ML>`_)
the first time a :class:`MLTracker` is constructed, and are cached locally by
``huggingface_hub`` for subsequent use. If the repository is private or gated, an
access token must be supplied either through
:attr:`~tctrack.core.ml_tracker.TCMLParameters.hf_token` or the ``HF_TOKEN``
environment variable. The token is removed from the parameters once the tracker is
constructed, so that it is not written to any output file.

Alternatively, a model file that is already available on disk can be used by setting
:attr:`~tctrack.core.ml_tracker.TCMLParameters.model_path`, in which case nothing is
downloaded. The model is run on the CPU by default; set
:attr:`~tctrack.core.ml_tracker.TCMLParameters.device` to use another PyTorch device,
such as ``"cuda"``.

Algorithm
---------

The model takes 17 input channels at each timestep: five pressure-level variables
(relative humidity, air temperature, eastward and northward wind, and relative
vorticity) at three pressure levels (1000, 750 and 500 hPa), a land/sea mask, and sea
surface temperature. Each channel is normalised using statistics computed over the
training set. The model outputs a probability for each of the five classes (background
and four lifecycle stages: non-developing storm, cyclolysis, cyclogenesis and active
cyclone) at every grid point.

Tropical cyclone tracking then proceeds in two steps, which mirror the detection and
stitching steps of the other trackers in TCTrack.

*Detection*: at each timestep, grid points whose most probable class is not background
and whose probability is at least :attr:`~tctrack.core.ml_tracker.TCMLParameters.threshold`
are identified as storm points. Neighbouring storm points are grouped into a single
candidate, regardless of their class, because a single storm can span two lifecycle
stages. The location of a candidate is the centroid of its points weighted by their
probability, and its class and score are those of its most probable point. Candidates
that lie within :attr:`~MLParameters.merge_distance_deg` of each other are merged,
keeping only the one with the highest score, because the model often splits a single
storm into several adjacent candidates.

*Stitching*: candidates are linked into trajectories by working through the timesteps
in order. At each timestep, every active track is extended with the nearest unmatched
candidate within :attr:`~MLStitchParameters.max_distance_deg` of its last position.
A track that goes unmatched for more than :attr:`~MLStitchParameters.max_gap`
timesteps is closed, and tracks with fewer than
:attr:`~MLStitchParameters.min_length` points are discarded.

Distances between candidates are calculated in degrees. The longitude difference is
wrapped across the date line and scaled by the cosine of the latitude, which is an
approximation that is accurate enough for the small distances involved.

The stitching step is not part of the reference training pipeline, which only trains
the per-point classification. The model was trained on full global ERA5 fields, and
the detections it makes on small regional subsets can be less reliable.

Usage
-----

Usage of the ML tracker is performed using the :mod:`tctrack.machine_learning` module.
This contains the :class:`MLTracker` class, which is used to run the algorithm, and the
:class:`MLParameters` and :class:`MLStitchParameters` dataclasses for specifying the
parameters of the detection and stitching steps.

The example below illustrates how to use these to run the tracker with the default
parameters. The output trajectories are written to the ``"trajectories.nc"`` file in a
CF-compliant NetCDF format.

.. code-block:: python

    from tctrack.machine_learning import MLParameters, MLStitchParameters, MLTracker

    ml_params = MLParameters(input_file="era5_2020.nc")
    stitch_params = MLStitchParameters()

    tracker = MLTracker(ml_params, stitch_params)
    tracker.run_tracker("trajectories.nc")

Creating the :class:`MLTracker` downloads or loads the model. The stitching
parameters are optional, and the defaults are used if they are not provided. Further
parameters can be set in :class:`MLParameters`, the most notable of which is
:attr:`~tctrack.core.ml_tracker.TCMLParameters.threshold`, the minimum probability for
a point to be treated as part of a storm.

The :meth:`~MLTracker.run_tracker` method performs several steps in succession. These
are:

* :meth:`~MLTracker.detect`, to run the model over each timestep and locate the
  candidate storms. This reads and normalises the input data itself, using
  :meth:`~MLTracker.preprocess`.
* :meth:`~MLTracker.stitch`, to link the candidates into trajectories.
* :meth:`~MLTracker.to_netcdf`, to write the trajectories to a CF-compliant NetCDF
  file.

It is possible to call some or all of these to perform the tracking manually, for
example to inspect the detections before stitching them:

.. code-block:: python

    tracker = MLTracker(ml_params, stitch_params)
    tracker.detect()
    trajectories = tracker.stitch()
    tracker.to_netcdf("trajectories.nc")

The trajectories are written in the `CF-Conventions
<https://cfconventions.org/>`_ `trajectory data format
<https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#trajectory-data>`_.
In addition to the time, latitude and longitude, each point of a trajectory has the
lifecycle stage (``class_index``), the model confidence (``score``) and the values of
the input variables at the storm location. The file can be read using any NetCDF
reading utility, though `cf-python <https://ncas-cms.github.io/cf-python/>`_ will load
it following the `CF data model <https://ncas-cms.github.io/cf-python/#cf-data-model>`_.

The individual storm locations found by :meth:`~MLTracker.detect`, before stitching
links them into tracks, can also be written as independent CF point data using
:meth:`~MLTracker.detections_to_netcdf`:

.. code-block:: python

    tracker.detections_to_netcdf("detections.nc")

This is useful for inspecting or evaluating the detections independently of the
stitching step, as it contains every candidate whether or not it was linked into a
trajectory of the required length.

Input Data
----------

The input data to the ML tracker must be contained in a single CF-NetCDF file, with all
of the variables on the same latitude/longitude grid:

* Pressure-level variables:

  * Each of :attr:`~MLParameters.pressure_variables` (by default relative humidity,
    air temperature, eastward wind, northward wind and relative vorticity) at each of
    :attr:`~MLParameters.pressure_levels` (by default 1000, 750 and 500 hPa).
  * These are selected from the input file by CF identity (``standard_name``), so the
    file must contain these variables with a ``Z`` (pressure) dimension covering the
    configured levels.

* Sea surface temperature:

  * Selected using :attr:`~MLParameters.sst_variable`. ERA5 sea surface temperature has
    no ``standard_name``, so by default it is selected using its NetCDF variable name
    (``ncvar%sst``).
  * This is used together with :attr:`~MLParameters.t2m_variable` (2-metre
    temperature) to fill the sea surface temperature over land, and to derive a
    land/sea mask from where the sea surface temperature is undefined.

The file must also have a ``time`` coordinate, which is used to order the timesteps
that :meth:`~MLTracker.detect` iterates over and to timestamp the output.

The statistics used to normalise the channels are those computed over the training set
of the pretrained model, which are included with TCTrack. If you are using a model
trained elsewhere, the statistics can be provided in a NetCDF file containing ``mean``
and ``range`` for each channel, using :attr:`~MLParameters.normalisation_stats_path`.

Some preprocessing of the input data is likely to be required in order to extract these
variables, combine multiple files over time, and regrid onto a consistent grid. See
:doc:`../data/preprocessing_data` for examples.

Example Data
------------

A small ERA5 sample, covering the northern Mozambique Channel from 10 to 12 January
2025 (6-hourly), when Cyclone Dikeledi crossed it, is included in the TCTrack
repository in ``data/machine_learning/``. It contains all of the variables above, and
the script ``tutorial/run_ml_tracker.py`` runs the ML tracker on it.

Modifying Parameters
--------------------

The following parameters can be changed to suit the data:

* :attr:`~tctrack.core.ml_tracker.TCMLParameters.threshold` is the minimum probability
  for a point to be considered part of a storm. Lowering it finds more of a storm, and
  can allow tracks to be followed for longer, but it also gives more spurious
  detections. The tutorial script uses ``0.25`` for the example data, as opposed to
  the default of ``0.5``.
* :attr:`~MLParameters.merge_distance_deg` is the distance within which candidates at
  the same timestep are merged. Setting it to ``0`` disables merging.
* :attr:`~MLStitchParameters.max_distance_deg` is the furthest a storm can move between
  timesteps and still be linked. The default assumes 6-hourly input data, as used for
  training, so it should be reduced for higher-frequency data.
* :attr:`~MLStitchParameters.max_gap` and :attr:`~MLStitchParameters.min_length`
  control how many missed timesteps a track can survive, and how many points it needs
  to be kept.
