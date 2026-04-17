Getting Started
===============

This getting started guide is intended for new users of TCTrack to get them up
and running as quickly as possible whilst introducing the main concepts.
More advanced users should consult the full :doc:`api documentation <../api/index>`.

After reading the installation information here new users may wish to work through the
:doc:`tutorial <tutorial>` to see examples of how TCTrack is used and familiarise
themselves with the workflow before using their own data.

.. toctree::
   :maxdepth: 3
   :hidden:


Installation
------------

Environment
~~~~~~~~~~~

TCTrack is a Python package and requires Python version 3.10 or later.

Before installing TCTrack ensure that any external :ref:`dependencies
<getting-started/index:dependencies>` are installed. This can be done manually as
:ref:`described below <getting-started/index:dependencies>`. However, the recommended
approach is to install these non-python dependencies using a conda virtual environment.

To install conda we recommend using `miniforge
<https://github.com/conda-forge/miniforge>`_. Once miniforge is installed conda can be
initialised with the below command. This can also be added to your shell rc file. ::

    source /path/to/miniforge3/etc/profile.d/conda.sh

Then create and activate the conda virtual environment with the relevant dependencies
installed::

    conda create -n tctrack-env -c conda-forge cf-python cf-plot udunits2 esmpy
    conda activate tctrack-env

When finished using TCTrack this environment can be turned off with ``conda deactivate``.

The individual tracking algorithms will also need to be :ref:`installed separately
<getting-started/index:tracking algorithms>`.


Installation from PyPI
~~~~~~~~~~~~~~~~~~~~~~

Once the environment has been set up TCTrack can be installed. It is easiest to do this
from PyPI using pip::

     pip install tctrack


Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, it can be installed from source by cloning the repository from GitHub::

    git clone https://github.com/Cambridge-ICCS/TCTrack
    cd TCTrack

pip can then be run to install the code from the source into the environment::

    pip install .

If you are developing TCTrack you should install as an editable package with the
additional developer dependencies::

    pip install --editable .[dev]

The `dev` optional dependencies include the `test`, `lint`, and `doc` subgroups.


Dependencies
------------

UDUNITS
~~~~~~~

TCTrack makes use of the `cf-python <https://ncas-cms.github.io/cf-python/>`_ package to
generate CF-compliant NetCDF files.
This brings with it the dependency of `UDUNITS <https://docs.unidata.ucar.edu/udunits/current/>`_
which needs to be installed by the user.

Binaries are available on most systems, for example on Ubuntu::

    apt-get install -y libudunits2-dev

or on mac one can use homebrew::

    brew install udunits

noting that we may need to add the library to the dynamic path e.g.::

    export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/udunits/2.2.28/lib


esmpy
~~~~~

Any regridding of data with cf-python requires `esmpy
<https://earthsystemmodeling.org/esmpy/>`_ and `ESMF
<https://earthsystemmodeling.org/>`_ as dependencies. This is not needed directly in the
TCTrack package but may be needed for initial pre-processing of data, such as in the
tutorial and described in the :doc:`../data/preprocessing_data` page.

These are not pip-installable but can be installed in a conda environment::

    conda install -c conda-forge esmpy


Tracking algorithms
~~~~~~~~~~~~~~~~~~~

To use some tracking algorithms requires additional dependencies to be installed, namely
the libraries that they are wrapping.

For detailed information see the specific
:doc:`tracking algorithms documentation<../tracking-algorithms/index>`.
