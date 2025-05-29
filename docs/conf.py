# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "TCTrack"
copyright = "2025, Jack Atkinson, Sam Avis"
author = "Jack Atkinson, Sam Avis"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",  # numpy docstring support
    "sphinx.ext.viewcode",  # links to highlighted source code
    "sphinx.ext.autodoc",  # pull in documentation from docstrings
    "sphinx_toolbox.more_autodoc.typehints",  # use type hints from code/docstrings
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "*tests*"]

pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
autodoc_member_order = "bysource"
napoleon_preprocess_types = True
