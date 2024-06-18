# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from docutils import nodes


# -- Project information -----------------------------------------------------

project = "pyranges_plot"
copyright = "2024, Ester Muñoz del Campo, Marco Mariotti"
author = "Ester Muñoz del Campo, Marco Mariotti"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
]
autosummary_generate = True  # Enable summary table generation

autodoc_default_options = {
    "members": True,
    "imported-members": True,
}

templates_path = ["_templates"]
html_static_path = ["_static"]

# exclude_patterns = ['_generated_hidden*']

master_doc = "index"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"


def monospaced_link(name, rawtext, text, lineno, inliner, options={}, content=[]):
    url = text.split(" ")[-1].strip("<>")
    clickable_text = " ".join(text.split(" ")[:-1])
    # Create a reference node, which is the docutils node for hyperlinks
    # unescaped_text = utils.unescape(text)

    node = nodes.reference(rawtext, clickable_text, refuri=url, **options)

    # Add a special class to this node
    node["classes"].append("monospaced-link")
    return [node], []


def setup(app):
    app.add_role("mslink", monospaced_link)
    app.add_css_file("custom.css")
