TSTORMS
=======

TSTORMS is an early tracking software for tropical cyclones originally developed
by NOAA GFDL.
It is made available under a GPL-2.0 license.

TCTrack provides bindings for both the ``Driver`` routine for detecting candidate
storms, also the ``Trajectories`` functionality for generating trajectories.

For full details of the TSTORMS API in TCTrack see the
:doc:`TCTrack TSTORMS API documentation <../api/tstorms_api>`.

.. toctree::
   :maxdepth: 2
   :hidden:

.. Import the tempest_extremes module to use references throughout this page.
.. py:module:: tctrack.tstorms
   :no-index:

Installation
------------

TODO

Usage
-----

Cyclone tracking in TSTORMS consists of two phases: the "driver" which performs
node detection of candidate storms for each snapshot in time, and "trajectories" which
stitches together suitable nodes across timesteps to generate trajectories.

Usage of TSTORMS in TCTrack is through the ``tstorms`` module.

This provides the :class:`TSTORMSTracker` class that stores algorithm parameters and provides
access to the methods. The driver and trajectory algorithms can be configured through
the various parameters in the :class:`DriverParameters` and
:class:`TrajectoryParameters` dataclasses.

.. TODO
.. Example of driver
.. Example of trajectories
.. Example of run_tracker
.. Notes about duplicated parameters or other quirks
.. Simple description of the algorithm
