Developer/Contributor Information
=================================

TCTrack is an open-source project developed by the Instutute of Computing for
Climate Science (ICCS) at the University of Cambridge.

Contributions are welcome, be they bug reports, bug fixes,
features/enhancements, documentation improvements etc.
For details of how to get involved with any of these see below.

.. contents:: Contents:
   :local:
   :depth: 1

.. note::

  Parts of this document are based on the 
  `xarray contributing guide <https://docs.xarray.dev/en/stable/contribute/contributing.html>`_
  which itself is heavily based on the 
  `Pandas Contributing Guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_.


.. _issues:

Bug Reports and Feature Requests
--------------------------------

Bug reports are an important part of making *TCTrack* more stable.
Having a complete bug report allows others to reproduce the bug and provides insight
into fixing.

Trying out the bug-producing code on the main branch is often a worthwhile exercise
to confirm that the bug still exists.
It is also worth searching existing bug reports and pull requests to see if the issue
has already been reported and/or fixed.


Submitting a bug report
~~~~~~~~~~~~~~~~~~~~~~~

If you find a bug in the code or documentation, submit a ticket to the
`Issue Tracker <https://github.com/Cambridge-ICCS/TCTrack/issues>`_.

When reporting a bug, please include the following:

#. A short, self-contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github-flavored-markdown/>`_::

      ```python
      import tctrack

      """Your code here."""
      ...
      ```

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the *TCTrack* community and be open to
comments/ideas from others.

See this `stackoverflow article <https://stackoverflow.com/help/mcve>`_ 
for more detailed tips on writing a good bug report.

If you have fixed the issue yourself in your own version of *TCTrack* please note
this on the issue and follow up by opening a :ref:`pull request <pull_request>`.


Submitting a feature request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there is a feature you think would be useful to add to *TCTrack* please make
a request using the
`Issue Tracker <https://github.com/Cambridge-ICCS/TCTrack/issues>`_.
When doing so please try to include the following details:

#. A clear description of the functionality you would like to see added and why you feel
   it would be useful.

#. Any external references that are helpful/required to understand and implement
   your request in detail.

#. A short code snippet indicating what the desired API might look like.

If you have already implemented the feature yourself in your own version of
*TCTrack* please note this on the issue and follow up by opening a
:ref:`pull request <pull_request>`.


.. _pull_request:

Pull Requests
-------------

If you have something to contribute to the *TCTrack* codebase this is done by
opening a pull request to the
`main repository on GitHub <https://github.com/Cambridge-ICCS/TCTrack>`_.

Here is a summary of the expected workflow when contributing:

#. Make sure there is an open issue on the
   `Issue Tracker <https://github.com/Cambridge-ICCS/TCTrack/issues>`_ as
   :ref:`described above <issues>` detailing the bug/feature that you are addressing.

#. `Fork the main repository <https://github.com/Cambridge-ICCS/TCTrack/fork>`_
   into your own personal GitHub space, and then clone and work on this fork.
   You should work on a branch within this fork with a sensible name that reflects
   what you are working on.

#. As you work on the code, commit your changes in sensible-sized chunks with clear
   commit messages.
   A commit should detail any changes made to perform a particular action en route
   to the overall goal. When writing commit messages remember that it needs to be
   clearly understandable to other developers as to what they contribute.
   See previous commits in the project for examples.

   As you work keep the following aspects in mind:

   a. Do not place large changes to multiple files in a single commit.

   b. Try and remember to apply the :ref:`stylistic checks and balances <linting>`
      to your code before committing.

   c. Make sure that you include appropriate tests alongside your
      code contributions. Code without tests will not be merged.

   d. Make sure that you include/update any docstrings in the code, and that they
      conform to the `numpy style <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
      See the rest of the code for examples.

   e. Make sure that you :ref:`update the documentation <documentation>` where
      necessary to reflect the additions you have made. If adding a significant
      top-level feature to the code you may want to update the
      relevant pages to showcase your additions.

#. Once you push code back to your GitHub fork you can open a pull request.
   For small bug-fixes and features you may wait until you feel things are complete
   before opening the pull request.
   However, if you wish for feedback/intermediate review then please open the pull
   request in draft mode during development.

#. When opening a pull request ensure that it contains:

   * A sensible title summarising its contribution.
   * A `reference <https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/autolinked-references-and-urls>`_
     to the issue number(s) that it is addressing.
   * A description of what has been done making it easy for the maintainers to review.

Once a pull request is opened it will be reviewed by the project maintainers and any
requests for changes/improvement fed back to the author.
Once the maintainers are happy, your code will be approved and the pull request merged!


.. _linting:

Code quality and formatting
---------------------------

Writing good code is not just about what you write.
It is also about *how* you write it.
During continuous integration several tools will be run to check your code
for stylistic errors.
Generating any warnings will cause these tests to fail.
Thus, good style is a requirement for submitting code to *TCTrack*.

*TCTrack* uses tools to ensure consistent and quality code formatting throughout:

- `ruff <https://docs.astral.sh/ruff/>`_ for:

  - standardized code formatting
  - code quality checks
  - checking docstrings against the numpy conventions

- `mypy <http://mypy-lang.org/>`_ for static type checking of
  `type hints <https://docs.python.org/3/library/typing.html>`_.

These will be checked on all pull requests and commits to main, so it is suggested you
run them on your code before committing.

This can be done with a development install by running the following bash commands from
the root directory:

.. code-block:: shell

    ruff format src/
    ruff check src/
    mypy src/
    blackdoc docs/

Sometimes it makes sense to
`disable a ruff warning <https://docs.astral.sh/ruff/linter/#error-suppression>`_.
We generally prefer that this is done on a case-by-case basis in the code.
If you have justification for turning off any warnings in your contribution please
document them in your pull request.

The full *ruff* configuration for the project is contained in the 
`pyproject.toml <https://github.com/Cambridge-ICCS/TCTrack/blob/main/pyproject.toml>`_
file.

Code quality is enforced for pull requests through a code-quality continuous integration
workflow run using GitHub actions.


.. _testing:

Testing
-------

All code contributions should have accompanying unit and integration tests to ensure
that all parts of the code are functioning properly.

TCTrack uses the `pytest <https://docs.pytest.org>`_ framework for testing,
with subprocess calls mocked with `pytest-mock <https://pytest-mock.readthedocs.io>`_ in
unit tests.

Tests are stored separately from the main code in the ``tests/`` directory at the root
of the package. There are separate subdirectories for ``unit`` (testing the TCTrack
Python code) and ``integration`` (testing the interaction with other libraries) tests.

To run the tests from a development install use, from the command-line:

.. code-block:: shell

    pytest tests/unit
    pytest tests/integration

to run both the unit and integration tests sequentially.
Note that the integration tests require installation of third-party libraries that
TCTrack wraps.

Code quality of the tests is maintained using ruff and mypy (see :ref:`linting`).
Check these from a development install by running:

.. code-block:: shell

    ruff format tests/
    ruff check tests/
    mypy src tests

Note that mypy is run over both src and tests to pick up the TCTrack type hints.

Testing standards are enforced for pull requests through continuous integration
workflows run using GitHub actions.
These are run for a number of Python versions and operating systems.


.. _documentation:

Documentation
-------------

The documentation is written in
`reStructuredText <https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html>`_
and built using Sphinx.
The `Sphinx Documentation <https://www.sphinx-doc.org/en/master/contents.html>`_
has an excellent
`introduction to reST <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
in addition to other aspects of Sphinx.

Sphinx allows much of the API documentation to be be generated automatically
from the docstrings in the source code.
This is why it is important to put time into these.

The rest of the documentation, such as the installation and getting started pages, and
the contribution guidelines that you are reading right now, are written out and stored
in the ``docs/`` directory of the code.

To build the documentation on a development install run::

    cd docs/
    make clean
    make html

This will generate HTML output files in the folder ``docs/_build/html/`` that can be
viewed in a browser.
