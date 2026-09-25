TCTrack Tutorial
================

This page outlines a tutorial for using TCTrack to generate cyclone tracks using
Tempest Extremes and TSTORMS.

We work though all steps of the process from obtaining and preparing data,
installing and using TCTrack and the wrapped algorithms, to plotting some simple
outputs.

Details are provided for installing and running multiple tracking algorithms, though
you may choose to follow just one.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import tctrack to use references throughout this page.
.. py:module:: tctrack
   :no-index:


Installation
------------

First, follow the :ref:`installation instructions <getting-started/index:installation>`.
Make sure to use a conda environment and install TCTrack from source (not from PyPI) so
that the tutorial scripts are downloaded.

The next step is to install the trackers we want to call from TCTrack, in this case
:doc:`Tempest Extremes <../tracking-algorithms/tempest_extremes>` and
:doc:`TSTORMS <../tracking-algorithms/tstorms>`.

To do this, first ensure that the following dependencies have also been installed:

* NetCDF (with C++ bindings and Tempest Extremes and Fortran bindings for TSTORMS)
* A C++ Compiler (for Tempest Extremes)
* A Fortran Compiler (for TSTORMS) -- ifort is assumed,
  though :ref:`others are available <tracking-algorithms/tstorms:installation>`

From your cloned version of TCTrack navigate to the ``tutorial/`` directory at the top
level::

    cd tutorial

This directory contains a number of pre-prepared scripts demonstrating a simple
TCTrack workflow to provide an example of usage and take you through the process of
detecting cyclone tracks from climate data.


To install Tempest Extremes source the installation script which will
clone and build Tempest Extremes locally under the tutorial directory, as
described on the :doc:`Tempest Extremes pages <../tracking-algorithms/tempest_extremes>`,
and add the executables to the ``PATH``.
Read carefully to understand and check that you are happy with what it will do before
running::

    source install_tempest_extremes.sh


To install TSTORMS source the installation script to clone and build TSTORMS locally
under the tutorial directory, as described on the :doc:`TSTORMS pages <../tracking-algorithms/tstorms>`.
Read carefully to understand and check that you are happy with what it will do before
running::

    source install_tstorms.sh


Obtaining Data
--------------

We will use example data from the ERA5 reanalysis dataset from 1950-09-01 to 1950-09-07.

This can be downloaded from the TCTrack GitHub releases using the included fetch data
script::

    python fetch_data.py

This will download and place the NetCDF files in a ``data/`` directory.


Pre-processing of Data
----------------------

From inside the conda environment run the script to pre-process the data::

    python preprocess_data.py

This will pre-process the downloaded data as required for our codes and place it in
``data_processed/``.
This includes the following processes:

* 10m wind speed is calculated from the wind components
* Surface geopotential is converted to orography in metres
* Single variables are extracted for TSTORMS inputs with the required dimension
  names, orientation, and time units
* A mean is taken over pressure levels of temperature

Note the use of the Python ``del`` command where appropriate as the data can consume a
large amount of memory which we want to free when possible.


Running the code
----------------

To run Tempest Extremes over the data we use the enclosed Python script which will
execute using the parameters from the [Ullrich2021]_ paper::

    python run_tempest_extremes.py

Intermediate files will be placed in ``te_outputs/`` whilst the resulting tracks file
will be output as ``tracks_tempest_extremes.nc``.
These can be inspected using::

    ncdump -h tracks_tempest_extremes.nc


To run TSTORMS over the data we use the enclosed Python script which will
execute similarly to the [Vitart2001]_ paper::

    python run_tstorms.py

Intermediate files will be placed in ``tstorms_outputs/`` whilst the resulting tracks file
will be output as ``tracks_tstorms.nc``.
These can be inspected using::

    ncdump -h tracks_tstorms.nc


Visualising Results
-------------------

Finally we can visualise the results by plotting the reacks from the output data files.
This can be done by running the included plotting script::

    python plot_tracks.py

which will generate a png figure of tracks plotted on a map using windspeed as a measure
of intensity.

By default this will read from ``tracks_tempest_extremes.nc``, but this can be changed in
the file.

Note that the plotting script requires the following Python packages to be installed in
the local environment: ``numpy``, ``NetCDF4``, ``matplotlib``, and ``cartopy``.
