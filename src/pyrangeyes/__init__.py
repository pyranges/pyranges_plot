from ._version import __version__  # noqa: F401
from .core import (
    set_engine,  # noqa: F401
    get_engine,  # noqa: F401
    set_id_col,  # noqa: F401
    get_id_col,  # noqa: F401
    set_warnings,  # noqa: F401
    get_warnings,  # noqa: F401
    set_theme,  # noqa: F401
    get_theme,  # noqa: F401
    print_options,  # noqa: F401
    get_options,  # noqa: F401
    set_options,  # noqa: F401
    reset_options,  # noqa: F401
)
from .plot_main import plot  # noqa: F401
from .track import Track  # noqa: F401
from . import adapters  # noqa: F401
from .pr_register_plot import register_methods  # noqa: F401
from . import example_data  # noqa: F401
from . import vcf  # noqa: F401
from .make_subsets import make_scatter  # noqa: F401


def browse(*args, **kwargs):
    """Experimental interactive browser.

    Plotly is an optional dependency; importing pyrangeyes must not require it.
    Import the browser implementation lazily so environments installed without
    the plotly extra can still use the rest of pyrangeyes.
    """
    try:
        from .browser import browse as _browse
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("plotly"):
            raise ModuleNotFoundError(
                "pyrangeyes.browse is experimental and requires the optional "
                "plotly dependencies. Install with `pip install pyrangeyes[plotly]` "
                "or `pip install pyrangeyes[all]` to use it."
            ) from exc
        raise
    return _browse(*args, **kwargs)
