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

Install TCTrack from GitHub.
It is recommended that a Python virtual environment is used.
For more information see the :ref:`installation instructions <getting-started/index:installation>`.::

    git clone git@github.com:Cambridge-ICCS/TCTrack.git
    cd TCTrack
    pip install .


The next step is to install the trackers we want to call from TCTrack, in this case
:ref:`Tempest Extremes <../tracking-algorithms/tempest_extremes>` and
:ref:`TSTORMS <../tracking-algorithms/tstorms>`.

Ensure that the following dependencies have been installed

* NetCDF (with C++ bindings and Tempest Extremes and Fortran bindings for TSTORMS)
* A C++ Compiler (for Tempest Extremes)
* A Fortran Compiler (for TSTORMS) -- ifort is assumed,
  though :ref:`others are available <../tracking-algorithms/tstorms:installation>`


From your cloned version of TCTrack navigate to the ``tutorial/`` directory at the top
level::

    cd tutorial

This directory contains a number of pre-prepared scripts demonstrating a simple
TCTrack workflow to provide an example of usage and take you through the process of
detecting cyclone tracks from climate data.


To install Tempest Extremes source the installation script which will
clone and build Tempest Extremes locally under the tutorial directory, as
described on the :ref:`Tempest Extremes pages <../tracking-algorithms/tempest_extremes>`,
and add the executables to the ``PATH``.
Read carefully to understand and check that you are happy with what it will do before
running::

    source install_tempest_extremes.sh


To install TSTORMS source the installation script to clone and build TSTORMS locally
under the tutorial directory, as described on the :ref:`TSTORMS pages <../tracking-algorithms/tstorms>`.
Read carefully to understand and check that you are happy with what it will do before
running::

    source install_tstorms.sh


Obtaining Data
--------------

We will use example data from CMIP6, specifically from the HadGEM model and the
1950-historical experiment. We will use just a subset of ASO (August, September, October)
from 1951 to demonstrate the code and some preprocessing techniques.

This data can be obtained from the `ESGF CEDA archive <https://esgf.ceda.ac.uk/thredds/catalog/catalog.html>`_
or direct from `CEDA <https://data.ceda.ac.uk/>`_ using the included fetch data script::

    ./fetch_data.sh

This will fetch the NetCDF data files using ``wget`` and place them in a ``data/``
directory.
Again, read carefully to understand and check that you are happy with what it will do
before running.


Pre-processing of Data
----------------------

Some preprocessing of the data is required.
We will do this using cf-python and esmf as described on the
:doc:`../data/preprocessing_data` pages.

The recommended method for installing esmf is via Conda for which we recommend using
`miniforge <https://github.com/conda-forge/miniforge>`_.
From a conda base environment create a new environment with cf-python and esmpy installed::

    conda create -n regrid -c conda-forge cf-python cf-plot udunits2 esmpy
    conda activate regrid

From inside the conda environment run the regridding script to pre-process the data::

    python regrid.py

This will pre-process the downloaded data as required for our codes and place it in
``data_processed/``.
This includes the following processes:

* Target months (ASO) are extracted from the yearly files
* The 3hr data is mapped to the same daily time values as the daily data
* Single variables at single levels are extracted for TSTORMS inputs
* Vorticity is calculated from velocity data (following regridding to surface grid)
* A mean is taken over pressure levels of temperature

Note that since this requires significant IO it may take a little time to complete.
Note also the use of the Python ``del`` command where appropriate as the data consumes a
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

By default this will read from tracks_tempest_extremes.nc, but this can be changed in
the file.

Note that the plotting script requires the following Pythion packages to be installed in
the local environment: ``numpy``, ``NetCDF4``, ``matplotlib``, and ``cartopy``.



.. rubric:: References

.. [Ullrich2021] Ullrich, Paul A., et al. “TempestExtremes v2.1: A Community Framework for Feature Detection, Tracking, and Analysis in Large Datasets.” Geoscientific Model Development 14, no. 8 (2021): 5023–48. https://doi.org/10.5194/gmd-14-5023-2021.

.. [Vitart2001] Vitart, Frédéric, and Timothy N. Stockdale.
   "Seasonal Forecasting of Tropical Storms Using Coupled GCM Integrations",
   Monthly Weather Review 129, 10 (2001): 2521-2537.
   `https://doi.org/10.1175/1520-0493(2001)129\<2521:SFOTSU\>2.0.CO;2 <https://doi.org/10.1175/1520-0493(2001)129\<2521:SFOTSU\>2.0.CO;2>`_

