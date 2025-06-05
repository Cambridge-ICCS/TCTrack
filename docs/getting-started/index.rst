Getting Started
===============

This getting started guide is intended for new users of TCTrack to get them up
and running as quickly as possible whilst introducing the main concepts.
More advanced users should consult the full :doc:`API documentation`<../api/index>`.

.. toctree::
   :maxdepth: 3
   :hidden:


Installation
------------

Installation Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

TCTrack is a Python package.
It can be installed by cloning the repository from GitHub and using pip::

    git clone https://github.com/Cambridge-ICCS/TCTrack
    cd TCTrack
    pip install .

It is recommended that this is done from within a virtual environment::

    python3 -m venv .venv
    source .venv/bin/activate

If you are developing TCTrack you should install as an editable package with the
additional developer dependencies::

    pip install --editable .[dev]

The `dev` optional dependencies include the `test`, `lint`, and `doc` subgroups.


Dependencies
------------

To use some tracking algorithms requires additional dependencies to be installed, namely
the libraries that they are wrapping.

For detailed information see the specific
:doc:`tracking algorithms documentation<../tracking-algorithms/index>`.
