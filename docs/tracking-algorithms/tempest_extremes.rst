Tempest Extremes
================

Tempest Extremes is a popular software for detecting and tracking tropical cyclones
and other weather systems in climate datasets.
It is written by Paul Ullrich and released under a BSD 2-clause license.

The code can be found `on GitHub <https://github.com/ClimateGlobalChange/tempestextremes>`_,
with `online documentation available <https://climate.ucdavis.edu/tempestextremes.php>`_.

TCTrack provides bindings for tropical cyclone tracking sections of Tempest Extremes,
namely the ``DetectNodes`` and ``StitchNodes`` functionalities.
Other components may be added in future is required by users.

For full details of the Tempest Extremes API in TCTrack see the
:doc:`TCTrack Tempest Extremes API documentation <../api/tempest_extremes_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the tempest_extremes module to use references throughout this page.
.. py:module:: tctrack.tempest_extremes
   :no-index:

Installation
------------

Using the Tempest Extremes module in TCTrack requires Tempest Extremes to be installed
on a user's system.

The code is built using CMake.
It requires a C++ compiler and an installation of NetCDF with C++ bindings.
Full instructions for this can be found in the
`Tempest Extremes Documentation <https://github.com/ClimateGlobalChange/tempestextremes/tree/master?tab=readme-ov-file#installation-via-cmake-recommended>`_
but are summarised here::

    git clone https://github.com/ClimateGlobalChange/tempestextremes.git
    cd tempestextremes/
    git checkout 5feb3a04d29fd62a1f13fa9c0b85daeefcbecd6f
    mkdir build
    cmake -B build/ -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF .
    cmake --build build/

This will clone the Tempest Extremes code and check out the most recent commit that
TCTrack has been built against (5feb3a0 - a more recent version may work if you require
the latest features).
It will then build using CMake
Once this is complete the Tempest Extremes executeables can then be found in
``tempestextremes/build/bin/``.

To use these from TCTrack you will need to add the directory to your ``PATH`` so that
the executeables can be found at runtime::

    export PATH=$PATH:~/tempestextremes/build/bin

Note: you will need to modify this as appropriate if you cloned Tempest Extremes
somewhere other than ``~/``.

Usage
-----

Cyclone tracking in Tempest Extremes consists of two phases: node detection of candidate
storms for each snapshot in time, and stitching of nodes across timesteps to generate
tracks.

Usage of Tempest Extremes in TCTrack is done through the ``tempest_extremes`` module.

This provides the :class:`TETracker` class that stores algorithm parameters and provides
access to the methods. The detection and stitching algorithm can be configured through
the various parameters in the :class:`DetectNodesParameters` and
:class:`StitchNodesParameters` dataclasses.

Detection
^^^^^^^^^

In the following example we demonstrate the approach for detecting tropical cyclones as
detailed in [Ullrich2021]_. First, we set up the :meth:`~TETracker.detect_nodes`
functionality to run on a series of input files to generate output. We configure
detection to be done based on minima in sea-level pressure (psl) every 6 hours, with
filters based on closed contours of psl and the geopotential height difference [*]_ [*]_, and
merging of candidates within 6 degrees of one another. Additional output fields are also
added for psl, surface elevation (orog), and surface windspeed (sfcWind) so that they
may be used with :meth:`~TETracker.stitch_nodes`:

.. code-block:: python

    import tctrack.tempest_extremes as te

    input_files = [
        "psl_E3hr_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001010300-195006302100.nc",
        "zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001-195003.nc",
        "orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
        "sfcWind_Prim3hr_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001010130-195006302230.nc",
    ]

    closed_contours = [
        te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
        te.TEContour(
            var="_DIFF(zg7h(250hPa),zg7h(500hPa))", delta=-6.0, dist=6.5, minmaxdist=1.0
        ),
    ]

    output_commands = [
        te.TEOutputCommand(var="psl", operator="min", dist=0.0),
        te.TEOutputCommand(var="orog", operator="max", dist=0.0),
        te.TEOutputCommand(var="sfcWind", operator="max", dist=2.0),
    ]

    dn_params = te.DetectNodesParameters(
        in_data=input_files,
        search_by_min="psl",
        time_filter="6hr",
        merge_dist=6.0,
        closed_contours=closed_contours,
        out_header=True,
        output_commands=output_commands,
        output_file="nodes_out.dat",
    )

    te_tracker = te.TETracker(dn_params)

    run_info = te_tracker.detect_nodes()

Stitching
^^^^^^^^^

This can then followed by StitchNodes which is set up to combine nodes into a track that
are less than 8 degrees from one another, with a track length of at least 54 hours and 8
degrees end-to-end, with up to 24 hours of missing data between adjacent nodes. This is
then filtered based upon the latitude, windspeed, and surface altitude. The format of
the ``"tracks_out.txt"`` output file is described in the documentation for
:meth:`~TETracker.stitch_nodes`:

.. code-block:: python

    threshold_filters = [
        te.TEThreshold(var="lat", op="<=", value=50, count=10),
        te.TEThreshold(var="lat", op=">=", value=-50, count=10),
        te.TEThreshold(var="orog", op="<=", value=150, count=10),
        te.TEThreshold(var="sfcWind", op=">=", value=10, count=10),
    ]

    sn_params = te.StitchNodesParameters(
        output_file="tracks_out.txt",
        caltype="360_day",
        max_sep=8.0,
        min_time="54h",
        max_gap="24h",
        min_endpoint_dist=8.0,
        threshold_filters=threshold_filters,
    )

    te_tracker = te.TETracker(dn_params, sn_params)

    run_info = te_tracker.stitch_nodes()

Output
^^^^^^

Finally, after running stitch nodes to generate tracks, one can write the tracks to
a NetCDF file fully compliant with the `CF-Conventions <https://cfconventions.org/>`_
(specifically the `trajectory data format <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#trajectory-data>`_)
using :meth:`~TETracker.to_netcdf`:

.. code-block:: python

    te_tracker.to_netcdf("my_cf_tracks.nc")

This can be read using any NetCDF reading utility, though
`cf-python <https://ncas-cms.github.io/cf-python/>`_ will load it following the
`CF data model <https://ncas-cms.github.io/cf-python/#cf-data-model>`_.

Combined run
^^^^^^^^^^^^

The above examples demonstrate running :meth:`~TETracker.detect_nodes`,
:meth:`~TETracker.stitch_nodes`, and :meth:`~TETracker.to_netcdf` separately. However,
it is likely that users will want to these in succession which can be done using the
:meth:`~TETracker.run_tracker` method after defining a :class:`TETracker` object with
appropriate :class:`DetectNodesParameters` and :class:`StitchNodesParameters`:

.. code-block:: python

    dn_params = te.DetectNodesParameters(...)
    sn_params = te.StitchNodesParameters(...)
    te_tracker = te.TETracker(dn_params, sn_params)

    te_tracker.run_tracker("my_cf_tracks.nc")

When running in this way there are certain parameters that *should not* be
set. In :class:`StitchNodesParameters`, :attr:`~StitchNodesParameters.in_fmt`
and :attr:`~StitchNodesParameters.in_file` will be determined from
:class:`DetectNodesParameters`. Also,
:attr:`~StitchNodesParameters.out_seconds` should be left ``False`` to enable
proper conversion to the netCDF output file.

If the intermediate files generated by Tempest Extremes are not of interest then the
values of :attr:`DetectNodesParameters.output_file` and
:attr:`StitchNodesParameters.output_file` can be left as ``None`` to avoid saving them.
They will instead be stored in a temporary directory that lasts for the lifetime of the
:class:`TETracker` instance.

Input data
----------

The input data must be stored in NetCDF files with dimensions named ``lon``, ``lat``,
and ``time`` (or otherwise specified using the :attr:`~DetectNodesParameters.lon_name`
and :attr:`~DetectNodesParameters.lat_name` parameters). There can be multiple input
files and multiple variables per file. However, each variable must only appear in one
file, i.e. variables split over time into multiple files must first be combined (see
:ref:`combine_time`).

.. [*] In the example code a functional operation (`_DIFF`) is used for the geopotential
   height difference. View the `TempestExtremes online documentation
   <https://climate.ucdavis.edu/tempestextremes.php#VariableExpressions>`_ for details
   about this and other functional operations. As is done here, functional operations
   can be used for the closed contours, but they should not be used elsewhere as the
   reading of the `detect_nodes` output file and setting the output metadata will not
   work.

.. [*] The geopotential height difference in [Ullrich2021]_ is calculated between 300
   hPa and 500 hPa. However, if the input data does not have a 300 hPa level then using
   a 250 hPa level is a reasonable substitute.

.. rubric:: References

.. [Ullrich2021] Ullrich, Paul A., et al. “TempestExtremes v2.1: A Community Framework for Feature Detection, Tracking, and Analysis in Large Datasets.” Geoscientific Model Development 14, no. 8 (2021): 5023–48. https://doi.org/10.5194/gmd-14-5023-2021.

