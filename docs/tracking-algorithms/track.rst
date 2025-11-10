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

TRACK requires `NetCDF <https://www.unidata.ucar.edu/software/netcdf>`_ (it can
be installed through system package managers, or see the `detailed installation
instructions
<https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html>`_)
as well as a C compiler, a FORTRAN 77 compiler, Make, and `makedepend`. It can
then be installed with the following commands:

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

This will clone TRACK and checkout the most recent commit that TCTrack has been tested
against, however there may be newer versions available. The code is then compiled using
`make`, creating the executable at ``track/bin/track.run``.

For further details about installation see the instructions `here
<https://gitlab.act.reading.ac.uk/track/track/-/blob/master/INSTALL>`_.

Algorithm
---------

Tropical cyclone tracking in TRACK uses relative vorticity at 850 hPa. This is
preprocessed using spectral filtering to retain only wavenumbers 6 - 63. Tropical
cyclone candidates are then identified at each time step where there is a local maximum
in the relative vorticity exceeding a specific threshold. These candidates are then
stitched into trajectories by minimising a cost function. Finally, the trajectories are
filtered to remove any that do not travel a sufficient distance or last for a sufficient
duration.

Tracking can also be performed using different fields, including mean sea level
pressure, but this is not currently implemented in TCTrack. In addition, the TRACK
algorithm describes performing checks for a warm core based upon the difference in
vorticity between 850 and 200 hPa. However, this is also not currently implemented.

Further details can be found in [Hodges2017]_.

Usage
-----

TRACK takes a single NetCDF file as input which must contain the northward and eastward
windspeeds on a single `Gaussian grid <https://en.wikipedia.org/wiki/Gaussian_grid>`_
and in units of ``m s-1``. Some preprocessing of the input data may be required. For
example, using the following `NCO <https://nco.sourceforge.net/nco.html>`_ and `CDO
<https://code.mpimet.mpg.de/projects/cdo/wiki/tutorial>`_ commands:

.. code-block:: bash

   # Combine the windspeeds onto a single grid / file (grid given by ua_file.nc)
   ncremap -i va_file.nc -d ua_file.nc -o ua_va_file.nc

   # Map to a Gaussian grid with 512 (2x256) latitude points
   cdo remapbil,F256 ua_va_file.nc ua_va_file_F256.nc

Usage of TRACK is performed using the :mod:`tctrack.track` module. This contains the
:class:`TRACKTracker` class which is used to run the algorithm, and the
:class:`TRACKParameters` dataclass for specifying the various parameters.

The example below illustrates how to use these to run TRACK using the default
parameters. The output trajectories are written to the ``"trajectories.nc"`` file in a
CF-compliant NetCDF format. Note that the location of the TRACK source directory must be
provided using the :attr:`~TRACKParameters.base_dir` parameter and the NetCDF variable
names for the windspeed should be specified with
:attr:`~TRACKParameters.wind_var_names`.

.. code-block:: python

    from tctrack.track import TRACKParameters, TRACKTracker

    params = TRACKParameters(
        base_dir="/path/to/track",
        input_file="ua_va_file_F256.nc",
        wind_var_names=("ua", "va"),
    )

    tracker = TRACKTracker(params)
    tracker.run_tracker("trajectories.nc")

In addition to :attr:`~TRACKParameters.base_dir` and
:attr:`~TRACKParameters.wind_var_names` there are a number of other parameters
that can be used in :class:`TRACKParameters`. Most notable is
:attr:`~TRACKParameters.filter_distance` which sets the minimum distance (in degrees) of
trajectories.

The :meth:`~TRACKTracker.run_tracker` method performs several steps in sucession. These
are:

* :meth:`~TRACKTracker.calculate_vorticity`, to calculate the relative vorticity from
  the windspeed.
* :meth:`~TRACKTracker.spectral_filtering`, to perform the spectral filtering on the
  vorticity.
* :meth:`~TRACKTracker.tracking`, to perform the tracking to identify trajectories.
* :meth:`~TRACKTracker.filter_trajectories`, to filter trajectories based on distance
  and duration.
* :meth:`~TRACKTracker.to_netcdf`, to convert the output to a CF-compliant NetCDF file.

It is possible to call some or all of these to manually perform the tracking or repeat
individual steps without rerunning the whole process. For each step, intermediate input
and output files are stored in the 'indat' and 'outdat' folders in the TRACK source
directory - the same location as when calling TRACK manually. Therefore, any generated
files that you wish to keep should be moved to prevent them from being overwritten.

Modifying Parameters
--------------------

Each step of the algorithm makes a call out to the TRACK software with a list of inputs
that are built up using :class:`TRACKParameters`. If you wish to experiment with
different inputs that are not currently implemented in TCTrack then you can export these
input lists, modify the resulting files, and then run TRACK using these files as input.

To export the inputs you can set the :attr:`~TRACKParameters.export_inputs` parameter to
``True``. The result will be a ``.in`` file for each step of the algorithm with one
input per line. It is a bit opaque what each input corresponds to, so refer to the
comments in the ``_get_<step>_inputs`` methods in :class:`TRACKTracker`, or look at the
output of the TRACK calls.

To read in the files, simply set the :attr:`~TRACKParameters.read_inputs` parameter to
``True``. Be aware that the files will then supercede any values set in
:class:`TRACKParameters`. If any of the inputs files do not exist it will fall back to
the generated inputs for that step. So it is possible to only read inputs for a single
step by deleting the exported files for the other steps. The location of the input files
can also be changed by setting the value of :attr:`~TRACKParameters.inputs_directory`.

If you identify any useful inputs it would be appreciated if you can contribute by
adding these to :class:`TRACKParameters` and the relevant input generation method(s).

.. rubric:: References

.. [Hodges2017] Hodges, Kevin, Alison Cobb, and Pier Luigi Vidale. "How Well Are Tropical Cyclones Represented in Reanalysis Datasets?", Journal of Climate 30, 14 (2017): 5243-5264, doi: https://doi.org/10.1175/JCLI-D-16-0557.1
