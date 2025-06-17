Tempest Extremes
================

Tempest Extremes is a popular software for detecting and tracking tropical cyclones
and other weather systems in climate datasets.
It is written by Paul Ullrich and released under a BSD 2-clause license.

For an overview of the functionalities, installation, and usage see the
:doc:`Tempest Extremes section of the documentation <../tracking-algorithms/tempest_extremes>`.

.. automodule:: tctrack.tempest_extremes

.. autoclass:: TETracker
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: DetectNodesParameters
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: StitchNodesParameters
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TEContour
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TEOutputCommand
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TEThreshold
   :members:
   :undoc-members:
   :show-inheritance:


.. _stitch-nodes-output-format:
``stitch_nodes()`` Output Format
--------------------------------

The default ``"gfdl"`` output of :meth:`TETracker.stitch_nodes` is a plain-text "nodefile" format which
contains a number of tracks, each of which in the form below.

.. code-block:: text

   start <N> <year> <month> <day> <hour>
   	 <i> <j> <lon> <lat> <var1> <var2> ... <year> <month> <day> <hour>
   	 ...
   	 <i> <j> <lon> <lat> <var1> <var2> ... <year> <month> <day> <hour>

- ``N`` is the number of nodes in the track (and number of lines below header).
- ``i``, ``j`` are grid indices.
- ``var1``, ``var2``, etc., are scalar variables as defined by :attr:`StitchNodesParameters.in_fmt` (typically, psl, orog).
- ``hour`` may instead be seconds if :attr:`StitchNodesParameters.out_seconds` is ``True``.
