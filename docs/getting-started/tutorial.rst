TCTrack Tutorial
================

We have included some tutorials for using TCTrack to generate cyclone tracks using
Tempest Extremes and TSTORMS.

These consists of Jupyter notebooks in the ``tutorial/`` directory of
the repository. These go through the steps of preprocessing data and running the
tracking algorithms for a week of data. These are shown on the pages below:

* :doc:`Tempest Extremes tutorial <tutorial_tempest_extremes>`
* :doc:`TSTORMS tutorial <tutorial_tstorms>`

.. Adds the tutorials to the contents menu
.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial_tempest_extremes
   tutorial_tstorms


Running manually
----------------

To run the tutorial manually, first follow the :ref:`installation instructions
<getting-started/index:installation>`. Make sure to use a conda environment and install
TCTrack from source (not from PyPI) so that the tutorial notebooks are downloaded.

Use the additional tutorial dependencies when installing tctrack::

    pip install .[tutorial]

From your cloned version of TCTrack navigate to the ``tutorial/`` directory and run
jupyter::

    cd tutorial
    jupyter notebook

Next, open the desired notebook and follow the instructions.
