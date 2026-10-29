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

TCTrack is a Python package and requires Python version 3.11 or later.

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

Docker container as environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative approach for setting up the development environment is via Docker container.
The TCTrack repository contains a `Dockerfile` to build the required Docker image.
Installing `Docker Desktop <https://docs.docker.com/desktop/>`_ is a prerequisite for
this approach.

Once Docker Desktop is setup, the Docker image can be built by running the following
command from the root directory of the TCTrack repository::

    docker build -t tctrack:latest .

The Docker image has the TCTrack dependencies installed. The optional dependencies can
be installed later from within the Docker container as required. 

The Docker image has `Tempest Extremes` tracking algorithm already installed in it. The 
other tracking algorithms need to be installed manually as required. See
:doc:`tracking algorithms documentation<../tracking-algorithms/index>`.

To create a Docker container using the built image, run::

    docker run --name iccs_tctrack -it -v REPO_PATH_ON_HOST:REPO_PATH_ON_CONTAINER tctrack:latest

The above command creates a container named `iccs_tctrack` and provides a shell to access it.
The command also mounts the TCTrack repository on the host to the `REPO_PATH_ON_CONTAINER`
path within the container. As a result, any code changes made from the host will be
visible within the container as well.

Refer to the `Docker Docs <https://docs.docker.com/reference/cli/docker/container/>`_
to stop a running container and start it again later.


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
additional developer dependencies and the pre-commit hooks installed::

    pip install --editable .[dev]
    pre-commit install

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

Any preprocessing that involves regridding of data using :mod:`tctrack.preprocessing` or
cf-python requires `esmpy <https://earthsystemmodeling.org/esmpy/>`_ and `ESMF
<https://earthsystemmodeling.org/>`_ as dependencies. This is not needed directly in the
tracking algorithms.

These are not pip-installable but can be installed in a conda environment::

    conda install -c conda-forge esmpy


Tracking algorithms
~~~~~~~~~~~~~~~~~~~

To use some tracking algorithms requires additional dependencies to be installed, namely
the libraries that they are wrapping.

For detailed information see the specific
:doc:`tracking algorithms documentation<../tracking-algorithms/index>`.
