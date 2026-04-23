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
  - name: Samuel J. Avis
    affiliation: "1"
    orcid: 0000-0001-9637-9489
    corresponding: true
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
compliant with the Climate and Forecasting (CF) conventions [@EarthScienceDataSystems2024]
using the discrete sampling geometry trajectory format, ensuring that all metadata from
input data and processing parameters are preserved.
This approach makes TCTrack particularly valuable for climate model evaluation studies,
such as those in CMIP (the Coupled Model Intercomparison Project), and for tracking algorithm intercomparison work.
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
(sections 9 and H4 of CF-1.13 [@EarthScienceDataSystems2024]), with coordinate variables,
auxiliary coordinates, and complete metadata including algorithm name, TCTrack version,
and all parameters required for reproducibility.

This standardised framework that TCTrack provides will support researchers, forecasters,
and the natural hazard industry (e.g. insurers) for using different algorithms to project
the impact of tropical cyclones.
It is particularly valuable for two key use cases.
First, climate model evaluation studies—such as the assessment of tropical cyclones in
Neural GCM [@Kochkov2024]—require consistent application of tracking algorithms to model
outputs.
Second, intercomparison studies are essential given documented variations in cyclone
statistics between different CMIP models [@Bourdin2024], tracking methods [@Roberts2020],
and reanalysis datasets [@Hodges2017].
TCTrack facilitates such studies by providing a common framework whilst preserving the
distinct characteristics of each underlying algorithm.


# Software Design

TCTrack is implemented as a Python package built around abstract base classes
`TCTrackerParameters` and `TCTracker`.
The first captures the unique aspects of each algorithm in a dataclass as
standardised inputs.
The second is the interface to a tracking algorithm itself, which makes use of
`subprocess` calls to the underlying Fortran and C++ codes. This abstracts away the diverse
native interfaces of each algorithm from the end user.
Each supported tracking algorithm has concrete implementations of these classes that
handle algorithm-specific details whilst presenting a standardised Python API.
Users configure tracking by instantiating the appropriate parameter class, then pass
this to a tracker class which provides standard methods including `run_tracker()` and
`to_netcdf()`.

The use of abstract base classes introduces modularity to make it straightforward to add
additional tracking algorithms, with clear interfaces to the rest of the package.
This is designed so that new algorithms can either wrap existing, proven software—as
with those currently implemented—or be implemented directly in TCTrack. For example,
we are in currently in discussions with researchers to add a new machine-learning based
detection algorithm to the package.

A particularly important aspect of TCTrack is its standardised, CF-compliant output
format into which we transform output tracks from all algorithms before presenting
them to the end-user.
This ensures that outputs are FAIR and immediately usable in downstream analysis without
additional processing.
This reflects the growing importance of FAIR principles in climate science, and is
achieved by building on the cf-python software package [@Hassell2017, @Hassell2020].
Each tracker class preserves the metadata from input files and augments it to reflect any
processing of variables, such as spatial aggregation. It also provides the tracking
algorithm used, TCTrack version, and complete parameter specifications for
reproducibility.

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

# Research Impact Statement

TCTrack addresses a need within the community for simple, accessible tooling for exploring
the differences in tropical cyclone tracks and statistics produced by different tracking
algorithms and conducting intercomparison studies.
It has been used in climate science research as part of CMIP6 model comparison studies
by researchers at the University of Cambridge, with results presented to
the scientific community in @Atkinson2026. This work also aims to reproduce and
update the findings of @Knutson2020 to project the response of tropical cyclone activity
to future climate change using more recent CMIP6 data.
TCTrack is also being used in industry by INIGO Insurance to generate tracks in a
robust, standardised format for downstream use in data analysis and presentation.


# State of the Field

The three tracking algorithms provided by TCTrack (TRACK, TempestExtremes, and
TSTORMS) are established, peer-reviewed methods used in the climate science
community.
TRACK uses a feature-tracking approach based on vorticity maxima and has been applied
to assess tropical cyclones in reanalysis datasets [@Hodges2017],
TempestExtremes employs geometric criteria on multiple variables [@Ullrich2021],
and TSTORMS searches for co-located vorticity and temperature maxima and sea-level
pressure minima [@Vitart1997].
However, their native interfaces are command-line based with algorithm-specific input
formats and minimal output metadata, creating barriers to adoption and intercomparison.

Alternative Python-based tracking tools exist, such as tobac [@Heikenfeld2019], which
provides general feature tracking capabilities primarily for meteorological applications.
However, tobac implements new tracking methodologies rather than providing access to the
established algorithms that have been extensively used and validated in the climate community.
TCTrack's approach differs fundamentally: it enables researchers to use proven methods
whilst addressing the practical challenges of usability, interoperability, and metadata
standards.
TCTrack is also extensible in that additional tracking algorithms can be integrated
following the contributor guidelines and building on the abstract base class design.

The decision to build TCTrack rather than change the data format/running procedures of
multiple existing projects comes the desire to also provide a unified interface
for researchers as well as improved data representation.
TCTrack's focus on usability and FAIR data builds on what the existing tools offer,
making it a valuable part of the ecosystem.


# Acknowledgments

This project is supported by a philantropic donation from INIGO Insurance as part of
the InSPIRe project.
We also thank the Institute of Computing for Climate Science for their support.
We thank Alison Ming and Charles Powell for their testing and feedback during
development.

# AI Usage Disclosure

During the development of this software AI tools (ChatGPT, Devstral) were used to understand usage
of functions, review written code to suggest issues and improvements, and for generating
boilerplate and parameterisation of tests. Any code suggestions by AI were thoroughly
checked and double checked during code review by a second person.
During the preparation of this paper, AI was used to refactor an initial draft to align
with the changes to the JOSS scope. This was thoroughly checked and built upon
thereafter.

# References
