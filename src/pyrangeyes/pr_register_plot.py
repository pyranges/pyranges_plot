from pyranges1 import PyRanges
from .core import set_engine
from .plot_main import plot
from .track import Track


def _track(self, adapter=None, **options):
    return Track(self, adapter, **options)


def register_methods(engine=None):
    """
    Register pyrangeyes convenience methods on PyRanges.

    Adds ``PyRanges.plot(...)`` and ``PyRanges.track(...)``. Optionally, set the
    plotting engine.

    Parameters
    ----------
    engine: {str}, default None

        Optional string to set the engine for plotting: Matplotlib ('plt', 'matplotlib') or Plotly ('ply', 'plotly').

    Examples
    --------
    >>> import pyrangeyes as pe

    >>> pe.register_methods()

    >>> pe.register_methods("matplotlib")

    """

    if engine is not None:
        set_engine(engine)

    PyRanges.plot = plot
    PyRanges.track = _track
