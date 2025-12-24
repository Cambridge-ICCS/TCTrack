Getting Started
===============

This getting started guide is intended for new users of TCTrack to get them up
and running as quickly as possible whilst introducing the main concepts.
More advanced users should consult the full :doc:`api documentation`<../api/index>`.

After reading the installation information here new users may wish to work through the
:doc:`tutorial `<tutorial>` to see examples of how TCTrack is used and familiarise
themselves with the workflow before using their own data.

.. toctree::
   :maxdepth: 3
   :hidden:


Installation
------------

Installation Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

Before using TCTrack ensure that any external :ref:`getting-started/index:dependencies` as required are
installed as :ref:`described below <getting-started/index:dependencies>`.

TCTrack is a Python package (Python 3.10+).
It can be obtained by cloning the repository from GitHub::

    git clone https://github.com/Cambridge-ICCS/TCTrack
    cd TCTrack

Installation is performed using pip, run from within the cloned source directory.
It is recommended that this is done from within a virtual environment::

    python3 -m venv .venv
    source .venv/bin/activate

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


Tracking algorithms
~~~~~~~~~~~~~~~~~~~

To use some tracking algorithms requires additional dependencies to be installed, namely
the libraries that they are wrapping.

For detailed information see the specific
:doc:`tracking algorithms documentation<../tracking-algorithms/index>`.
