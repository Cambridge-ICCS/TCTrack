<div align="center">
  <img src="docs/images/TCTrack_Logo.svg" alt="TCTrack Logo">
</div>

# TCTrack

![GitHub License](https://img.shields.io/github/license/Cambridge-ICCS/TCTrack)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/Cambridge-ICCS/TCTrack/unit-tests.yaml?label=unit-tests)
[![Documentation Status](https://readthedocs.org/projects/TCTrack/badge/?version=latest)](https://tctrack.readthedocs.io/en/latest/?badge=latest)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)


TCTrack is a Python library providing bindings to tracking algorithms for tropical
cyclones in an accessible manner to generate high-quality and
[FAIR](https://en.wikipedia.org/wiki/FAIR_data) output data.

It can be used for tracking cyclones in simulations and observations, and to compare
the output of different algorithms for a variety of data sources.

> [!WARNING]  
> The software is currently in pre-release and may therefore be subject to changes.


## Installation

### Dependencies
The package requires Python 3 (>=3.10).

### Package Installation
First clone the repository using `git`:
```sh
git clone https://github.com/Cambridge-ICCS/TCTrack
```

Set up and activate a virtual environment:
```sh
python3 -m venv .venv
source .venv/bin/activate
```
When finished using TCTrack this can be turned off with `deactivate`.

Then install the package using `pip`:
```sh
pip install --editable .
```


## Using TCTrack

To use TCTrack from within Python after install import it as you would any other
library:

```python
import tctrack
#[Example]
```

Full details can be found in the
[getting-started documentation online](https://tctrack.readthedocs.io/developer/index.html).
For a complete description of the library API see 
[API documentation](https://tctrack.readthedocs.io/developer/index.html).


## Contributing

Contributions and collaborations are welcome.

For bugs, feature requests, and clear suggestions for improvement please
[open an issue](https://github.com/Cambridge-ICCS/TCTrack/issues).

If you have added something to _TCTrack_ that would be useful to others, or can
address an [open issue](https://github.com/Cambridge-ICCS/TCTrack/issues), please
[fork the repository](https://github.com/Cambridge-ICCS/TCTrack/fork) and open a
pull request.

Additional dependencies for deleopment can be installed as follows:
```sh
pip install --editable .[dev]
```

Full details for contribution and developers can be found in the
[online documentation](https://tctrack.readthedocs.io/developer/index.html).

### Code of Conduct

Everyone participating in the _TCTrack_ project, and in particular in the
issue tracker, pull requests, and social media activity, is expected to treat other
people with respect and, more generally, to follow the guidelines articulated in the
[Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/).


## License

Copyright &copy; ICCS

*TCTrack* is distributed under the [GPL 3](https://github.com/Cambridge-ICCS/TCTrack/blob/main/LICENSE).


## Acknowledgments

This work was funded by a philantropic donation to the
[University of Cambridge](https://www.cam.ac.uk/) from
[INIGO Insurance](https://inigoinsurance.com/) as part of the InSPIRe project.

<p align="left">
  <img src="docs/images/inigo_inspire.png" width="355">
</p>

The TCTrack logo was designed by [Jack Atkinson](https://jackatkinson.net/) - [@jatkinson1000](https://github.com/jatkinson1000).
