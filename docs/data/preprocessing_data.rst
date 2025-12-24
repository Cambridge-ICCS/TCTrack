Preprocessing Data
==================

Many of the tracking algorithms require the input data to be in specific formats. Here
we detail how to perform some of the typically required preprocessing steps using the
cf-python library. Other tools can be used for the same tasks, however we focus on
cf-python since it provides a uniform interface and it is a dependency of TCTrack.

For full documentation of the routines described on these pages and more see the
`cf python documentation <https://ncas-cms.github.io/cf-python/>`_.

.. _preprocessing_deps:

Dependencies
------------

Any regridding with cf-python requires `esmpy <https://earthsystemmodeling.org/esmpy/>`_
and `ESMF <https://earthsystemmodeling.org/>`_ as dependencies. These are not
pip-installable but can be installed in a conda environment:

.. code-block:: bash

    conda create -n tctrack_env
    conda activate tctrack_env
    conda install -c conda-forge esmpy h5py
    pip install tctrack

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

Separating Variables
--------------------

If variables instead need to be separated into multiple files, such as in :doc:`TSTORMS
<../tracking-algorithms/tstorms>`, the opposite proceedure is followed:

.. code-block:: python

    # Read in the combined file
    field1, field2 = cf.read("combined_file.nc")

    # Write to separate files
    cf.write(field1, "var1_file.nc")
    cf.write(field2, "var2_file.nc")

Subsampling
-----------

Sometimes we wish to subsample, e.g. to move from hourly data to daily.
This can be done again using cf-python's ``subspace`` command, this time providing a
``slice`` or indices to extract the values of interest:

.. code-block:: python

    # Read the separate input files
    field1 = cf.read("var1_file.nc")[0]

    # Generate subspaces as required
    # Take the 5th element of the 'Z' coordinate
    field2 = field1.subspace(Z=[5])
    # Take the zeroth and fifth elements of the 'X' coordinate
    field3 = field1.subspace(X=[0, 5])
    # Every second elements of the 'Y' coordinate between 3 and -3
    field4 = field1.subspace(Y=slice(3, -3, 2))

Note that if only a single element is taken (e.g. a slice of a single pressure level)
then the field will retain this as a coordinate dimension.
To remove the single-valued coordinate from the field use cf-python's
``squeeze`` before writing to file:

.. code-block:: python

    # Read the separate input files
    field1 = cf.read("var1_file.nc")[0]

    # Slice the 5th pressure level ('Z' coordinate)
    field2 = field1.subspace(Z=[5])

    # Squeeze to remove the single-valued Z from field dimensions
    field_ta.squeeze(inplace=True)
    # or, for a new field
    field3 = field2.squeeze()

Operations
----------

cf-python provides various operations to calculate new fields.
These include both matehmatical operations and statistical collapses.

For example, to calculate vorticity from coincident velocity date we can use ``curl_xy``:

.. code-block:: python

    # Read the separate input files
    u_field = cf.read("u_file.nc")[0]
    v_field = cf.read("v_file.nc")[0]

    # calculate vorticity
    w_field = cf.curl_xy(u_field, v_field, radius="earth")
    w_field.nc_set_variable("vorticity")
    w_field.set_property("standard_name", "atmosphere_upward_absolute_vorticity")
    w_field.set_property("units", "s-1")

    # Save the new variable to NetCDF
    cf.write(w_field, "vorticity_file.nc")

Or to take a mean over a coordinate:

.. code-block:: python

    # Read the separate input files
    field = cf.read("file.nc")[0]

    # Take the mean in the zonal 'X' coordinate and squeeze to remove 'X' dimension
    field_zonal_mean = field.collapse("mean", axes="X")
    field_zonal_mean.squeeze(inplace=True)

    # Save the new variable to NetCDF
    cf.write(field_zonal_mean, "zonal_mean_file.nc")

Setting Fill Values
^^^^^^^^^^^^^^^^^^^

Sometimes it us useful to replace fill values after an operation before writing to file.
This can be done using cf-python's ``filled`` routine.
For example, after to set any null or masked values to ``0.0`` after calculating
vorticity above use:

.. code-block:: python

    w_field.filled(fill_value=0.0, inplace=True)

before writing to file.

Set NetCDF Variable Name
------------------------

To set specfic NetCDF variable names for the fields and coordinates you can use the
``nc_set_variable`` methods:

.. code-block:: python

    field = cf.read("var1_file.nc")[0]

    # Set the new netcdf variable names for the field and coordinates
    field.nc_set_variable("slp")
    field.coordinate("latitude").nc_set_variable("lat")

    # Save with the new netcdf variable names
    cf.write(field, "slp_file.nc")

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

Note that regridding can be performed inplace using ``inplace=True``.

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
