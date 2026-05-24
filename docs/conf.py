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
import doctest
import sphinx_rtd_theme
import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------

project = "pyrangeyes"
copyright = "2024, Ester Muñoz del Campo, Marco Mariotti"
author = "Ester Muñoz del Campo, Marco Mariotti"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
]
autosummary_generate = True  # Enable summary table generation

doctest_default_flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE

doctest_global_setup = """
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyrangeyes as pe
import pyrangeyes.plot_main as _pe_plot_main
import pyrangeyes.pr_register_plot as _pe_pr_register_plot

if not hasattr(pe, '_doctest_real_plot'):
    pe._doctest_real_plot = _pe_plot_main.plot
_real_plot = pe._doctest_real_plot
pe.reset_options()
pe.set_id_col(None)

def _doctest_plot(*args, **kwargs):
    explicit_return = kwargs.get('return_plot') is not None or kwargs.get('to_file') is not None
    kwargs.setdefault('warnings', False)
    if not explicit_return:
        kwargs['return_plot'] = 'fig'
    result = _real_plot(*args, **kwargs)
    if not explicit_return:
        if hasattr(result, 'close'):
            result.close()
        else:
            plt.close('all')
        return None
    return result

pe.plot = _doctest_plot
_pe_plot_main.plot = _doctest_plot
_pe_pr_register_plot.plot = _doctest_plot
"""

doctest_global_cleanup = """
import matplotlib.pyplot as plt
plt.close('all')
"""


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
    clickable_label = " ".join(text.split(" ")[:-1])
    # Create a reference node, which is the docutils node for hyperlinks
    # unescaped_label = utils.unescape(text)

    node = nodes.reference(rawtext, clickable_text, refuri=url, **options)

    # Add a special class to this node
    node["classes"].append("monospaced-link")
    return [node], []


def setup(app):
    app.add_role("mslink", monospaced_link)
    app.add_css_file("custom.css")
