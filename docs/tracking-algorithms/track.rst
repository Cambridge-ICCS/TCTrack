TRACK
=====

TRACK is a software package and algorithm for tracking and analysing tropical
cyclones and other weather systems in meteorological and oceanic datasets. It
is developed by Kevin Hodges at the University of Reading, UK, and released
under an AGPLv3 license.

The TRACK code can be found `on the University of Reading GitLab
<https://gitlab.act.reading.ac.uk/track/track>`_ and the algorithm is described in
[Hodges2017]_.

TCTrack provides bindings to perform tropical cyclone tracking using the basic TRACK
algorithm. Other components may be added in future as required by users.

For full details of the TRACK API in TCTrack see the
:doc:`TCTrack TRACK API documentation <../api/track_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the track module to use references throughout this page.
.. py:module:: tctrack.track
   :no-index:

Installation
------------

Using the TRACK module in TCTrack requires TRACK to be installed on a user's system.

TRACK requires NetCDF as well as a C compiler, a FORTRAN 77 compiler, Make, and
`makedepend`. It can then be installed with the following commands:

.. code-block:: bash

  # Download the TRACK code
  git clone https://gitlab.act.reading.ac.uk/track/track
  cd track
  git checkout 6ded301a5f5183d73e5b49c16019024b9a53eff7

  # Set necessary environment variables
  export PATH=${PATH}:.
  export CC=gcc
  export FC=gfortran
  export NETCDF=$(nc-config --prefix)
  export ARFLAGS=

  # Update the source file dependencies
  make -C src -f Makefile.linux depend

  # Compile the code
  ./master -build -fext=run -inpf=inputs.nc -upath=$(cd .. && pwd)

This will clone TRACK and checkout the most recent commit that TCTrack has been built
against before compiling the code using `make`. The compiled executable will then be
``track/bin/track.run``.

For further details about installation see the instructions `here
<https://gitlab.act.reading.ac.uk/track/track/-/blob/master/INSTALL>`_.

Usage
-----

.. TODO

.. rubric:: References

.. [Hodges2017] Hodges, Kevin, Alison Cobb, and Pier Luigi Vidale. "How Well Are Tropical Cyclones Represented in Reanalysis Datasets?", Journal of Climate 30, 14 (2017): 5243-5264, doi: https://doi.org/10.1175/JCLI-D-16-0557.1
