Installation
~~~~~~~~~~~~

Pyranges Plot requires Python ≥ 3.12 and can be installed using pip.

As the plot production can based on either Plotly or Matplotlib, the graphic library-specific
installations are enabled. This way if the user wants to install exclusively the dependencies
correspondent to one of those libraries it can be done by running just one of the following
commands: ::

    pip install pyranges-plot[plotly]

    pip install pyranges-plot[matplotlib]

To install all dependencies in order to be able to use both engines' functionalities,
the [all] option must be specified: ::

    pip install pyranges-plot[all]

Note that the minimal installation by :code:`pip install pyranges-plot` is not able to produce
plots since the graphical dependencies are not installed.
