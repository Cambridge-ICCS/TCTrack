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
``tempestextremes/build_serial/bin/``.


Usage
-----

Cyclone tracking in Tempest Extremes consists of two phases: node detection of candidate
storms for each snapshot in time, and stitching of nodes across timesteps to generate
tracks.

Usage of Tempest Extremes in TCTrack is done through the ``tempest_extremes`` module.

This provides the ``TETracker`` class that stores algorithm parameters and provides access
to the methods.
The detection and stitching algorithm can be configured through the various parameters
in the ``DetectNodesParameters`` and ``StitchNodesParameters`` dataclasses.

In the following example we set up the DetectNodes functionality to run on a series of
input files to generate output. We configure detection to be done based on minima in
psl, with closed contours of psl and zgdiff, and merging of candidates within 6 degrees
of one another. Additional output fields are also added for psl and orog so that they
may be used with StitchNodes:

.. code-block:: python

    import tctrack.tempest_extremes as te

    input_files = [
        "psl_E3hr_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001010300-195006302100.nc",
        "zgdiff_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195001-195003.nc",
        "orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc",
    ]

    closed_contours = [
        te.TEContour(var="psl", delta=200.0, dist=5.5, minmaxdist=0.0),
        te.TEContour(var="zgdiff", delta=-6.0, dist=6.5, minmaxdist=1.0),
    ]

    output_commands = [
        te.TEOutputCommand(var="psl", operator="min", dist=0.0),
        te.TEOutputCommand(var="orog", operator="max", dist=0.0),
    ]

    dn_params = te.DetectNodesParameters(
        in_data=input_files,
        search_by_min="psl",
        merge_dist=6.0,
        closed_contours=closed_contours,
        out_header=True,
        output_commands=output_commands,
        output_file="nodes_out.dat",
    )

    te_tracker = te.TETracker(dn_params)

    run_info = te_tracker.detect_nodes()

This can then followed by StitchNodes which is set up to combine nodes into a track that
are less than 8 degrees from one another, with a track length of at least 10 nodes and 8
degrees end-to-end, with a maximum of 3 times missing between each pair of nodes. This
is then filtered based upon the lattitude and surface altitude. The format of the
``"tracks_out.txt"`` output file is described in
:meth:`~tctrack.tempest_extremes.TETracker.stitch_nodes`:

.. code-block:: python

    threshold_filters = [
        te.TEThreshold(var="lat", op="<=", value=40, count=10),
        te.TEThreshold(var="lat", op=">=", value=-40, count=10),
        te.TEThreshold(var="orog", op="<=", value=1500, count=10),
        te.TEThreshold(var="orog", op="<=", value=10, count=4),
    ]

    sn_params = te.StitchNodesParameters(
        output_file="tracks_out.txt",
        caltype="360_day",
        max_sep=8.0,
        min_time=10,
        max_gap=3,
        min_endpoint_dist=8.0,
        threshold_filters=threshold_filters,
    )

    te_tracker = te.TETracker(dn_params, sn_params)

    run_info = te_tracker.stitch_nodes()

However, it is likely preferable to run both
:meth:`~tctrack.tempest_extremes.TETracker.detect_nodes` and
:meth:`~tctrack.tempest_extremes.TETracker.stitch_nodes` together. Which can be done
with a single :class:`~tctrack.tempest_extremes.TETracker` object:

.. code-block:: python

    te_tracker = te.TETracker(dn_params, sn_params)
    te_tracker.detect_nodes()
    te_tracker.stitch_nodes()
