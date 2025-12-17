Preprocessing Data
==================

Many of the tracking algorithms require the input data to be in specific formats. Here
we detail how to perform some of the typically required preprocessing steps using the
cf-python library.

Any regridding will require `esmpy <https://earthsystemmodeling.org/esmpy/>`_ and `ESMF
<https://earthsystemmodeling.org/>`_ as dependencies. These are not pip-installable but
can be installed in a conda environment.

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
