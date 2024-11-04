Installation
~~~~~~~~~~~~

Pyranges Plot requires Python ≥ 3.12 and can be installed using pip.

Pyranges Plot supports two alternative graphical libraries ("engines"): plotly and matplotlib.
At least one must be installed. Use these commands to install Pyranges Plot together with
your engine of choice::

    pip install pyranges-plot[plotly]

    pip install pyranges-plot[matplotlib]

To install both engines, use instead::

    pip install pyranges-plot[all]

Note that the minimal installation by :code:`pip install pyranges-plot` is not able to produce
plots since the graphical dependencies are not installed.
