from pyranges import PyRanges
from .core import set_engine
from .plot_main import plot


def register_plot(engine=None):
    """
    Register the plot function as a method to PyRanges.

    Allows to use the plot function as a method of PyRanges, as PyRanges.plot().
    Optionally, set the plotting engine.

    Parameters
    ----------
    engine: {str}, default None

        Optional string to set the engine for plotting: Matplotlib ('plt', 'matplotlib') or Plotly ('ply', 'plotly').

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.register_plot()

    >>> prp.register_plot("matplotlib")

    """

    if engine is not None:
        set_engine(engine)

    # Attach the wrapper as a method to PyRanges
    PyRanges.plot = plot
