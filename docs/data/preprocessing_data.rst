Preprocessing Data
==================

Many of the tracking algorithms require the input data to be in specific formats. Here
we detail how to perform some of the typically required preprocessing steps using the
cf-python library. Other tools can be used for the same tasks, however we focus on
cf-python since it provides a uniform interface and it is a dependency of TCTrack.

.. _preprocessing_deps:

Dependencies
------------

Any regridding with cf-python requires `esmpy <https://earthsystemmodeling.org/esmpy/>`_
and `ESMF <https://earthsystemmodeling.org/>`_ as dependencies. These are not
pip-installable but can be installed in a conda environment.

.. _combine_time:

Combining in Time
-----------------

Typically, data is separated into different files in time, but must often be combined
into a single file. The below code example illustrates how this can be done with
cf-python:

.. code-block:: python

    # Read the list input files. This automatically concatenates in time.
    input_files = [...]
    field = cf.read(input_files)[0]

    # (Optionally) Select a time interval. This uses the first three months of 1950
    time_interval = cf.wi(cf.dt("1950-01-01"), cf.dt("1950-04-01"), open_upper=True)
    field = field.subspace(T=time_interval)

    # Write the combined data to a single file
    cf.write(field, "combined-output.nc")

Combine Variables
-----------------

Variables will often be stored in separate files. To combine them with cf-python simply
read them in separately and then write them together:

.. code-block:: python

    # Read the separate input files
    field1 = cf.read("var1_file.nc")[0]
    field2 = cf.read("var2_file.nc")[0]

    # Write the combined fields to a single file
    cf.write([field1, field2], "combined_file.nc")

Regridding
----------

.. note::
    To regrid using cf-python requires esmpy and ESMF to be installed as dependencies
    (:ref:`see above <preprocessing_deps>`). There are also other tools available
    including xarray, `NCO <https://nco.sourceforge.net/nco.html>`_ (ncremap), and `CDO
    <https://code.mpimet.mpg.de/projects/cdo/wiki/tutorial>`_ (cdo remap...).

Regridding variables will involve either using the grid of an existing variable or
creating a new grid. Each of which is shown below. The interpolation method can be
specified using the ``method`` argument, with options such as ``"linear"``,
``"conservative"``, and nearest neighbour search (`see here for details
<https://ncas-cms.github.io/cf-python/method/cf.Field.regrids.html>`_).

To use an existing variable:

.. code-block:: python

    # Get the fields for the two variables
    field1 = cf.read("var1_file.nc")[0]
    field2 = cf.read("var2_file.nc")[0]

    # Regrid field1 onto the grid of field2
    field1 = field1.regrids(field2, method="linear")
    field1.nc_clear_dataset_chunksizes()  # Avoids a possible error when writing

To regrid onto a new grid:

.. code-block:: python

    field = cf.read("var1_file.nc")[0]

    # Create a new grid at regular longitude and latitude coordinates
    domain = cf.Domain.create_regular((-180, 180, 1), (-90, 90, 1))

    # Regrid
    field = field.regrid(domain, method="linear")
    field.nc_clear_dataset_chunksizes()  # Avoids a possible error when writing

Gaussian Grid
^^^^^^^^^^^^^

If, as in :doc:`TRACK <../tracking-algorithms/track>`, a regular `Gaussian grid
<https://en.wikipedia.org/wiki/Gaussian_grid>`_ is required (i.e. the latitude points
satisfy the arcsin of the roots of a Legendre polynomial), the new longitudes and
latitudes need to be defined. These are used to define new ``cf.DimensionCoordinate``
objects to be used for the regridding.

.. code-block:: python

    field = cf.read("var1_file.nc")[0]

    # Define a regular Gaussian grid with 'n' points per hemisphere
    n = 256
    lon = np.arange(0, 360, 360 / (4 * n))
    lat = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(2 * n)[0]))

    # Copy and modify the latitude and longitude DimensionCoordinates
    domain = field.domain.copy()
    lat_coord = domain.dimension_coordinate("latitude")
    lat_coord.set_data(lat, inplace=True)
    lat_coord.del_bounds()
    lon_coord = domain.dimension_coordinate("longitude")
    lon_coord.set_data(lon, inplace=True)
    lon_coord.del_bounds()

    # Regrid
    field = field.regrids((lat_coord, lon_coord), method="linear")
    field.nc_clear_dataset_chunksizes()  # Avoids a possible error when writing
