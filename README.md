# TCTrack

The TCTrack package aims to provide functionality for running tropical cyclone (TC)
tracking algorithms and analysing the results. This can be used to compare the output
of different algorithms for a variety of data sources.

> [!WARNING]  
> The software is currently in early development and therefore subject to significant
> changes.


## Installation

### Dependencies
The package requires Python 3 (>=3.9).

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
```python
import tctrack
#[Example]
```


## Contributing

Contributions and collaborations are welcome.

For bugs, feature requests, and clear suggestions for improvement please
[open an issue](https://github.com/Cambridge-ICCS/TCTrack/issues).

If you have added something to _TCTrack_ that would be useful to others, or can
address an [open issue](https://github.com/Cambridge-ICCS/TCTrack/issues), please
[fork the repository](https://github.com/Cambridge-ICCS/TCTrack/fork) and open a
pull request.

### Additional dependencies for development

- [pytest](https://docs.pytest.org/en/stable/) is used for unit testing.
- [ruff](https://docs.astral.sh/ruff/) is used for formatting and linting.
- [mypy](https://www.mypy-lang.org/) is used for static checking of type hints.

These optional developer dependencies can be installed as follows:
```sh
pip install --editable .[dev]
```

### Testing

All code contributions should have accompanying unit tests to ensure that all parts of
the code are functioning properly.

When added, the testing framework will utilise `pytest`.

### Code quality

All code should be annotated using [type hints](https://peps.python.org/pep-0484/) and
formatted using `ruff format`:
```sh
ruff format
```
In addition, the code should be linted using `ruff check` to check for errors:
```sh
ruff check
```
Details of the ruff configuration and specific rules can be found in the
`pyproject.toml` file.

Code quality is enforced on all pull requests and merges through a continuous
integration workflow using GitHub actions.

### Code of Conduct
Everyone participating in the _TCTrack_ project, and in particular in the
issue tracker, pull requests, and social media activity, is expected to treat other
people with respect and, more generally, to follow the guidelines articulated in the
[Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/).


## License

Copyright &copy; ICCS

*TCTrack* is distributed under the [MIT Licence](https://github.com/Cambridge-ICCS/TCTrack/blob/main/LICENSE).
