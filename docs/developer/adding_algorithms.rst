Adding a Tracking Algorithm to TCTrack
======================================

If you want to extend TCTrack to add a new tracking algorithm, the documentation on
this page provides guidance about how to do this to fit with the existing software.

.. contents:: Contents:
   :local:
   :depth: 2


Code
----

Packaging and Location
~~~~~~~~~~~~~~~~~~~~~~

Each new tracking algorithm should be added as a subpackage to the ``TCTrack``
package.

New code should go into a directory ``src/tctrack/my_tracking_algorithm/``,
include a ``__init__.py`` file defining the user-facing namespace, and be
added as an import to the main package namespace in ``src/tctrack/__init__.py``.
Further information about packaging in Python can be found in the
`Python Packaging User Guide <https://packaging.python.org/en/latest/guides/section-build-and-publish/>`_.

Structure
~~~~~~~~~

New algorithms should follow a similar structure to the existing implementations in
order to provide a familiar and standardised interface.

The code is written in an object-oriented approach with the main item being a class
``MyAlgorithmTracker``.
Initialisation of this class will take as input(s) the various specific parameters
required by the algorithm, setting them as instance attributes.
This will include input data, algorithm-specific parameters, and any information
controlling outputs.
Default values for parameters (instance attributes) should be set where possible such
that they can be used to generate sensible results without modification.

Any functionalities from the algorithm/code should be wrapped as class methods,
callable from an instance of that class, with behaviour governed by the instance
attributes set during initialisation.

The class should contain a ``run_tracker`` method that will perform the end-to-end
function of reading inputs and generating tracks as cf-compliant output.
This is intended to provide easy top-level usage and match the interface of other
algorithms in the ``TCTrack`` package.
It may be assembled from smaller, more specific methods that allow experienced users
increased control over the functionality.

The coding style described in the :ref:`developer guidelines <linting>` should be
followed for all new code.

Data
~~~~

Outputs from ``TCTrack`` algorithms are expected to be in a format that is compliant
with the `cf-conventions <https://cfconventions.org/>`_, specifically in the
`Trajectory format <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#_multidimensional_array_representation_of_trajectories>`_.
This can be achieved using the `cf-python package <https://ncas-cms.github.io/cf-python/>`_
which is already a dependency of ``TCTrack``.

Testing
~~~~~~~

New algorithms should come with a corresponding `pytest <https://docs.pytest.org>`_
unit test suite added at ``tests/unit/test_my_tracking_algorithm.py``.
See the :ref:`testing <testing>` section of the developer guidelines for more
information.

Integration tests for small dummy datasets can be added to
``tests/integration/test_my_tracking_algorithm.py``


Documentation
-------------

There are two pieces of documentation required for new algorithms.
Developers should follow the :ref:`general documentation guidelines <documentation>`
for ``TCTrack``, with further guidance below.

The first documentation is of the code itself.
All Python code should be annotated using
`docstrings <https://peps.python.org/pep-0257/>`_ following the
`numpy style <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
These should then be used to auto-generate API documentation by adding a file
at ``docs/api/my_tracking_algorithm_api.rst`` and including it in the index.

The second documentation is user information about your tracking algorithm.
This should be added as a page at
``docs/tracking-algorithms/my_tracking_algorithm.rst``, again adding to the main index.
This should include a brief description of the algorithm, comprehensive instructions for
how to obtain and install any external dependencies to ``TCTrack``, examples for
usage and getting started, and any key references or other information as required.
