TSTORMS
=======

TSTORMS is an early tracking software for tropical cyclones originally developed
by NOAA GFDL.
It is made available under a GPL-2.0 license.

TCTrack provides bindings for both the ``Driver`` routine for detecting candidate
storms, also the ``Trajectories`` functionality for stitching nodes to generate trajectories.

For full details of the TSTORMS API in TCTrack see the
:doc:`TCTrack TSTORMS API documentation <../api/tstorms_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the tempest_extremes module to use references throughout this page.
.. py:module:: tctrack.tstorms
   :no-index:

Installation
------------

TODO

Usage
-----

Cyclone tracking in TSTORMS consists of two phases: detection, which performs
node detection of candidate storms for each snapshot in time, and stitching, which
combines suitable nodes across timesteps to generate trajectories. [*]_

Usage of TSTORMS in TCTrack is through the ``tstorms`` module.

This provides the :class:`TSTORMSTracker` class that stores algorithm parameters and provides
access to the methods. The detection and stitching algorithms can be configured through
the various parameters in the :class:`TSTORMSBaseParameters`, :class:`TSRTORMSDetectParameters`, and
:class:`TSTORMSStitchParameters` dataclasses.

In the following example we demonstrate the approach for detecting tropical cyclones
using TSTORMS. First, we set up the :meth:`~TSTORMSTracker.detect`
functionality to run on a series of input files to generate output.
Detection of candidate storms is based on locating points that satisfy certain minima or
maxima, with full details of each input given in the :class:`TSTORMSDetectParameters`
documentation:

.. code-block:: python

    from tctrack.tstorms import (
        TSTORMSBaseParameters,
        TSTORMSDetectParameters,
        TSTORMSTracker,
    )

    tstorms_params = TSTORMSBaseParameters(
        tstorms_dir="path/to/tstorms/installation/",
        output_dir="path/to/place/outputs/",
        input_dir="path/to/input/files/",
    )
    detect_params = TSTORMSDetectParameters(
        u_in_file="u_input.nc",
        v_in_file="v_input.nc",
        vort_in_file="vort_input.nc",
        tm_in_file="tm_input.nc",
        slp_in_file="slp_input.nc",
        vort_crit=3.5e-5,
        tm_crit=0.0,
        thick_crit=50.0,
        dist_crit=4.0,
        lat_bound_n=70.0,
        lat_bound_s=-70.0,
        do_spline=False,
        do_thickness=False,
        use_sfc_wind=True,
    )

    # Initialize the tracker
    tracker = TSTORMSTracker(tstorms_params, detect_params)

    detect_call = tracker.detect(verbose=True)

The result of this will be a ``cyclones`` file in the prescribed output directory with
candidate storms organised by date in the format described in :meth:`~TSTORMSTracker.detect`.

These candidates may be stitched into plausible trajectories using the
:meth:`~TSTORMSTracker.stitch` functionality.
Stitching is based on locating subsequent candidate points within a maximum distance of
the previous point. Valid trajectories must last for a set time period and can be
filtered within latitude bounds and to satisfy variable thresholds.
Full details of each input are given in the :class:`TSTORMSStitchParameters` documentation.

.. code-block:: python

    from tctrack.tstorms import TSTORMSStitchParameters

    stitch_params = TSTORMSStitchParameters(
        r_crit=900.0,
        wind_crit=17.0,
        vort_crit=3.5e-5,
        tm_crit=0.0,
        n_day_crit=2,
        do_filter=True,
        lat_bound_n=70.0,
        lat_bound_s=-70.0,
    )
    tracker = TSTORMSTracker(tstorms_params, detect_params, stitch_params)
    stitching_call = tracker.stitch(verbose=True)

The result of this will be the following set of files at the output location
provided in :class:`TSTORMSStitchParameters`:

- ``ori`` containing the origins of each storm in the form
  ``lon, lat, YY, MM, DD, HH``,
- ``traj`` containing trajectory point data in the form
  ``lon, lat, wind, psl, YY, MM, DD, HH`` for each trajectory,
- ``trav`` containing trajectory point data and vorticity in the form
  ``lon, lat, wind, psl, vort_max, YY, MM, DD, HH`` for each trajectory,
- ``stats`` containing information about the number of storms per month/year
  in each basin.

If filtering is applied (:class:`TSTORMSStitchParameters` ``do_filter``) there will also be
additional files ``ori_filt``, ``traj_filt``, and ``trav_filt`` from after this
takes place.

Finally, after running stitch to generate storm trajectories, the tracks can be written to
a NetCDF file fully compliant with the `CF-Conventions <https://cfconventions.org/>`_
(specifically the `trajectory data format <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#trajectory-data>`_)
using :meth:`~TSTORMSTracker.to_netcdf`:

.. code-block:: python

    tracker.to_netcdf("my_tstorms_cf_trajectories.nc")

This can be read using any NetCDF reading utility, though
`cf-python <https://ncas-cms.github.io/cf-python/>`_ will load it following the
`CF data model <https://ncas-cms.github.io/cf-python/#cf-data-model>`_.

The above examples demonstrate running :meth:`~TSTORMSTracker.detect`,
:meth:`~TSTORMSTracker.stitch`, and :meth:`~TSTORMSTracker.to_netcdf` separately.
However, it is likely that users will want to these in succession which can be done
using the :meth:`~TSTORMSTracker.run_tracker` method after defining a :class:`TSTORMSTracker`
object with appropriate :class:`TSTORMSDetectParameters` and :class:`TSTORMSStitchParameters`:

.. code-block:: python

    from tctrack.tstorms import (
        TSTORMSBaseParameters,
        TSTORMSDetectParameters,
        TSTORMSStitchParameters,
        TSTORMSTracker,
    )

    tstorms_params = TSTORMSBaseParameters(...)
    detect_params = TSTORMSDetectParameters(...)
    stitch_params = TSTORMSStitchParameters(...)
    tracker = TSTORMSTracker(tstorms_params, detect_params, stitch_params)

    tracker.run_tracker("my_tstorms_cf_trajectories.nc")



.. [*] These are named detect and stitch for consistency across the TCTrack package.
   However, in the TSTORMS source they are referred to as tstorms_driver and
   trajectory_analysis respectively.

