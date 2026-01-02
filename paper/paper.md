---
title: 'TCTrack: FAIR Tropical Cyclone Tracking'
tags:
  - Python
  - climate science
  - weather
  - metadata
  - FAIR
  - CF conventions
  - environment
authors:
  - name: Jack W. Atkinson
    orcid: 0000-0001-5001-4812
    affiliation: "1"
    corresponding: true
  - name: Sam Avis
    affiliation: "1"
    orcid: 0000-0001-9637-9489
affiliations:
 - name: Institute of Computing for Climate Science, University of Cambridge, UK
   index: 1
date: 02 January 2026
bibliography: paper.bib

---

# Summary

TCTrack is a Python package providing a unified, user-friendly interface to multiple
established tropical cyclone tracking algorithms: TRACK [@Hodges2017],
TempestExtremes [@Ullrich2021], and TSTORMS [@Vitart1997].
Whilst these algorithms are widely used in the climate science community, each has its
own unique interface, input requirements, and output format, making intercomparison
studies and routine application challenging. \mbox{TCTrack} addresses these usability issues by
wrapping each algorithm with a consistent Python API based on abstract base classes.
Crucially, TCTrack transforms the output from each algorithm into a standardised NetCDF format that is fully
CF-compliant [@EarthScienceDataSystems2024] using the discrete sampling
geometry trajectory format, ensuring that all metadata from input data and processing
parameters are preserved.
This approach makes TCTrack particularly valuable for climate model evaluation studies,
such as those in CMIP, and for tracking algorithm intercomparison work.
By providing standardised, metadata-rich outputs, TCTrack aligns with the FAIR
principles for research data and software [@Wilkinson2016; @Barker2022], ensuring outputs are
Findable, Accessible, Interoperable, and Reusable.


# Statement of need

Tropical cyclone tracking is an essential component of climate model evaluation and
climate change impact assessment.
However, researchers face significant barriers when attempting to use existing tracking codes.
Each established algorithm—TRACK, TempestExtremes, and TSTORMS—requires separate installation,
has distinct input data requirements and parameter specifications, and produces output
in custom formats with limited or no metadata.
Users who wish to compare results across multiple tracking algorithms must invest
substantial effort learning each tool independently, converting between different data
formats, and manually adding provenance metadata to outputs.

TCTrack streamlines this workflow by providing a uniform interface across all three algorithms.
Users need only:

1) configure tracking parameters using algorithm-specific parameter dataclasses,
2) instantiate a tracker object and call the `run_tracker()` method,
3) use the CF-compliant netCDF output data.

All output files follow the CF conventions' discrete sampling geometry trajectory format
(sections 9 and H4 of CF-1.13), with coordinate variables, auxiliary coordinates, and
complete metadata including algorithm name, TCTrack version, and all parameters required
for reproducibility.
Cell methods are correctly specified when algorithms compute area means, maxima, or minima.

This standardised approach is particularly valuable for two key use cases.
First, climate model evaluation studies—such as the assessment of tropical cyclones in
Neural GCM [@Kochkov2024]—require consistent application of tracking algorithms to model
outputs.
Second, intercomparison studies are essential given documented variations in cyclone
statistics between different CMIP models [@Bourdin2024], tracking methods [@Roberts2020],
and reanalysis datasets [@Hodges2017].
TCTrack facilitates such studies by providing a common framework whilst preserving the
distinct characteristics of each underlying algorithm.


# Software description

TCTrack is implemented as a Python package built around abstract base classes for
`Tracker` and `TrackerParameters`.
Each supported tracking algorithm has concrete implementations of these classes that
handle algorithm-specific details whilst presenting a uniform interface.
Users configure tracking by instantiating the appropriate parameter class, then pass
this to a tracker class which provides standard methods including `run()` and
`to_netcdf()`.
The tracker classes manage all interaction with the underlying Fortran and C++ codes
through Python's `subprocess` module, eliminating the need for users to understand the
native interfaces of each algorithm.

A key design principle of TCTrack is the separation of concerns: the package focuses
specifically on the tracking step, whilst comprehensive documentation guides users in
preparing input data using cf-python [@Hassell2020]—a core dependency that provides
robust CF-compliant data handling.
The documentation includes tutorials demonstrating complete workflows and clearly
specifies the input requirements for each algorithm.

The output transformation is TCTrack's most significant contribution to reproducibility.
Each tracker class captures metadata from input files and augments it with processing
provenance, including the tracking algorithm used, TCTrack version, and complete
parameter specifications.
Variable metadata is preserved from inputs, and appropriate CF cell methods are added
when algorithms perform spatial aggregations.
The resulting trajectory-format netCDF files are fully standards-compliant and can be
readily consumed by downstream analysis tools without additional processing.

TCTrack is packaged using modern Python standards with a `pyproject.toml` configuration
and supports Python 3.10 and above.
The codebase includes comprehensive unit tests using pytest to verify the interface
implementations, and continuous integration via GitHub Actions ensures compatibility
across multiple Python versions.
Code quality is maintained through automated checks with ruff for linting and mypy for
static type checking, ensuring robustness and maintainability.

The package is open source, extensively documented at https://tctrack.readthedocs.io/,
and available via GitHub at https://github.com/Cambridge-ICCS/TCTrack.
The documentation includes installation guidance for the external tracking codes, API
reference, and worked examples.


# Comparison to other approaches

The three tracking algorithms wrapped by TCTrack—TRACK, TempestExtremes, and
TSTORMS—represent established, peer-reviewed methods widely used in the climate science
community.
TRACK uses a feature-tracking approach based on vorticity maxima and has been applied
to assess tropical cyclones in reanalysis datasets [@Hodges2017],
TempestExtremes employs geometric criteria on multiple variables [@Ullrich2021],
and TSTORMS searches for co-located vorticity and temperature maxima and sea-level
pressure minima [@Vitart1997].
Each has been validated against observations and used in numerous climate studies.
However, their native interfaces are command-line based with algorithm-specific input
formats and minimal output metadata, creating barriers to adoption and intercomparison.

Alternative Python-based tracking tools exist, such as tobac [@Heikenfeld2019], which
provides general feature tracking capabilities primarily for meteorological applications.
However, tobac implements new tracking methodologies rather than providing access to the
established algorithms that have been extensively used and validated in the climate community.
TCTrack's approach differs fundamentally: it enables researchers to use proven methods
whilst addressing the practical challenges of usability, interoperability, and metadata
standards.
The package architecture is also extensible—additional tracking algorithms can be integrated
following the contributor guidelines, as the abstract base class design accommodates
diverse implementations.

The emphasis on CF-compliant output distinguishes TCTrack from wrapper scripts or
informal tools that researchers may develop for personal use.
By leveraging cf-python [@Hassell2020] for data handling and strictly adhering to
CF conventions for trajectory data, TCTrack ensures that outputs meet community standards
for FAIR (Findable, Accessible, Interoperable, and Reusable) data
[@Wilkinson2016; @Barker2022].
This is particularly important for CMIP-related workflows and institutional data
repositories.

# Examples of Use

TODO

# Future development

TODO

# Acknowledgments

This project is supported by a philantropic donation from INIGO Insurance as part of
the InSPIRe project.
We also thank the Institute of Computing for Climate Science for their support.
We thank Alison Ming and Charles Powell for their testing and feedback during
development.


# References
