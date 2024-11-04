Tutorial
~~~~~~~~

This tutorial assumes some familiarity with pyranges v1.
If necessary, go through its tutorial first: https://pyranges1.readthedocs.io/

.. contents:: Contents of Tutorial
   :depth: 3


Getting started
---------------

The first compulsory step to obtain a plot is setting the **engine**, using function
:func:`set_engine <pyranges_plot.set_engine>` after importing. We also **register** the plot function
using :func:`register_plot <pyranges_plot.register_plot>`, which is optional but convenient:
it allows to use the plot function directly from PyRanges objects (further explained later).

    >>> import pyranges_plot as prp
    >>> prp.set_engine("plotly")  # possible engines: "plotly" and "matplotlib"
    >>> prp.register_plot()


Pyranges Plot centralizes the interface to producing graphics in
the :func:`plot <pyranges_plot.plot>` function. It offers plenty of options to
customize the appearance of the plot, showcased in this tutorial.
To that end, we will use some example data included in the Pyranges Plot package.
Yet, any PyRanges object can be used, e.g. loaded from gff, gtf, bam files.

    >>> p = prp.example_data.p1
    >>> print(p)
      index  |      Chromosome  Strand      Start      End  transcript_id    feature1    feature2
      int64  |           int64  object      int64    int64  object           object      object
    -------  ---  ------------  --------  -------  -------  ---------------  ----------  ----------
          0  |               1  +               1       11  t1               a           A
          1  |               1  +              40       60  t1               a           A
          2  |               2  -              10       25  t2               b           B
          3  |               2  -              70       80  t2               b           B
          4  |               2  +              85      100  t3               c           C
          5  |               2  +             110      115  t3               c           C
          6  |               2  +             150      180  t3               c           C
          7  |               3  +             140      152  t4               d           D
    PyRanges with 8 rows, 7 columns, and 1 index columns.
    Contains 3 chromosomes and 2 strands.

By default, :func:`plot <pyranges_plot.plot>` produces an interactive plot. If the Matplotlib engine is selected,
a window appears. If the Plotly engine is selected, a server is automatically opened, and
an address is printed in the console. The plot can be accessed by opening this address in a browser.

    >>> prp.plot(p)

.. image:: images/prp_rtd_01.png

Interactive navigation is intuitive:

* Hover over intervals to see their details in a **tooltip**
* Click and drag to zoom in on a region.
* Double-click to reset the zoom level.
* Inspect the rest of buttons on the top-right to see other available actions.

To create a pdf or png image file instead of opening an interactive plot,
use the ``to_file`` parameter of :func:`plot <pyranges_plot.plot>`.

    >>> prp.plot(p, to_file="my_plot.png")

Because we **registered** the plot function, we can also invoke it like a method of the PyRanges object, as
``PyRanges.plot(...)``. This is equivalent to the previous code:

    >>> p.plot(to_file="my_plot.png")

In the figure above, intervals are displayed individually, i.e. each PyRanges row is treated as a separate entity.
To link the intervals instead, as to represent a transcript composed of exons, use the ``id_column`` parameter,
indicating the column name that defines the groups of intervals.

    >>> prp.plot(p, id_col="transcript_id")

.. image:: images/prp_rtd_02.png

Because the ``id_col`` parameter is used frequently, it can be set as default for all plots using function
:func:`set_id_col <pyranges_plot.set_id_col>`. The following code is equivalent to the previous one:

    >>> prp.set_id_col("transcript_id")
    >>> prp.plot(p)


Selecting what to plot
----------------------
The data above has only 4 interval groups (hereafter, "transcripts") so all of them were included in the plot.
By default, a **maximum of 25 transcripts** are plotted, customizable with the ``max_shown`` parameter of
:func:`plot <pyranges_plot.plot>`.
Below, we can set the maximum number of transcripts show as 2. Note the warning shown:

    >>> prp.plot(p, max_shown=2)

.. image:: images/prp_rtd_03.png

To plot only a subset of the data, use the Pandas/PyRanges object's slicing capabilities.
For example, this plots the intervals on chromosome 2, positive strand, between positions 100 and 200:

    >>> (p.loci[2, '+', 100:200]).plot()

By default, the **limits of plot coordinates** are set to show all the data, and leave some margin at the edges.
This is customizable with the ``limits`` parameter.
The user can decide to change all or some of the coordinate limits leaving the rest as default if desired.
The ``limits`` parameter accepts different input types:

* Dictionary with chromosome names as keys, and a tuple of two integer numbers indicating the limits (or ``None`` to leave as default).

* Tuple of two integer numbers, which sets the same limits for all plotted chromosomes.

* PyRanges object, wherein Start and End columns define the limits for the corresponding Chromosome.

    >>> prp.plot(p, limits={1: (None, 100), 2: (60, 200), 3: None})

.. image:: images/prp_rtd_04.png

    >>> prp.plot(p, limits=(0,300))

.. image:: images/prp_rtd_05.png

Coloring
--------
By default, the intervals are **colored** according to the ID column
(``transcript_id`` in this case,  previously set as default with :func:`set_id_col <pyranges_plot.set_id_col>`).

We can select any other column to color the intervals by using the ``color_col`` parameter
of :func:`plot <pyranges_plot.plot>`.
For example, let's color by the Strand column:

    >>> prp.plot(p, color_col="Strand")

.. image:: images/prp_rtd_06.png

Now the "+" strand transcripts are displayed in one color and the ones on the "-" strand in another color.
Note that pyranges_plot used its default color scheme, and mapped each value in the  ``color_col`` column to a color.

The  **colormap** parameter of :func:`plot <pyranges_plot.plot>` centralizes coloring customization.
It is a versatile parameter, accepting many different types of input.
Using a dictionary allows to exert full control over the coloring, explicitly setting each value-color pair:

    >>> prp.plot(p, color_col="Strand",
    ...          colormap={"+": "green", "-": "red"})

.. image:: images/prp_rtd_07.png

Alternatively, the user may just define the sequence of colors used
(letting pyranges_plot pick which color to assign to each value).
One can provide a list of colors in hex or rgb; or a string recognized as the name of an available
Matplotlib or Plotly colormap;
or an actual Matplotlib or Plotly colormap object. Below, we invoke the "Dark2" Matplotlib colormap:

    >>> prp.plot(p, colormap="Dark2")

.. image:: images/prp_rtd_08.png

.. @maxtico: please add a plot showcasing the legend=True option. Add some short text before it

In this section, we have seen how to color intervals based on their attributes.
Next, we will see how to customize the appearance of the plot itself.


Appearance customization options: cheatsheet
--------------------------------------------

A wide range of **options** are available to customize appearance, as summarized below:

.. image:: images/options_fig_wm.png

These options can be provided as parameters to the :func:`plot <pyranges_plot.plot>` function, or
set as default beforehand. Let's see an example of providing them as parameters:

    >>> prp.plot(p, plot_bkg="rgb(173, 216, 230)", plot_border="#808080", title_color="magenta")

.. image:: images/prp_rtd_15.png

To instead set these options as default, use the :func:`set_options <pyranges_plot.set_options>` function:

    >>> prp.set_options('plot_bkg', 'rgb(173, 216, 230)')
    >>> prp.set_options('plot_border', '#808080')
    >>> prp.set_options('title_color', 'magenta')
    >>> prp.plot(p)  # this will now open a plot identical to the previous one

To inspect the current default options, use the
:func:`print_options <pyranges_plot.print_options>` function.
Note that any modified values from the built-in defaults will be marked with an asterisk (*):

    >>> prp.print_options()
    +------------------+--------------------+---------+--------------------------------------------------------------+
    |     Feature      |       Value        | Edited? |                         Description                          |
    +------------------+--------------------+---------+--------------------------------------------------------------+
    |     colormap     |       popart       |         | Sequence of colors to assign to every group of intervals     |
    |                  |                    |         | sharing the same “color_col” value. It can be provided as a  |
    |                  |                    |         | Matplotlib colormap, a Plotly color sequence (built as       |
    |                  |                    |         | lists), a string naming the previously mentioned color       |
    |                  |                    |         | objects from Matplotlib and Plotly, or a dictionary with     |
    |                  |                    |         | the following structure {color_column_value1: color1,        |
    |                  |                    |         | color_column_value2: color2, ...}. When a specific           |
    |                  |                    |         | color_col value is not specified in the dictionary it will   |
    |                  |                    |         | be colored in black.                                         |
    |   exon_border    |        None        |         | Color of the interval's rectangle border.                    |
    |     fig_bkg      |       white        |         | Bakground color of the whole figure.                         |
    |    grid_color    |     lightgrey      |         | Color of x coordinates grid lines.                           |
    |     plot_bkg     | rgb(173, 216, 230) |    *    | Background color of the plots.                               |
    |   plot_border    |      #808080       |    *    | Color of the line delimiting the plots.                      |
    |    shrunk_bkg    |    lightyellow     |         | Color of the shrunk region background.                       |
    |     tag_bkg      |        grey        |         | Background color of the tooltip annotation for the gene in   |
    |                  |                    |         | Matplotlib.                                                  |
    |   title_color    |      magenta       |    *    | Color of the plots' titles.                                  |
    |    title_size    |         18         |         | Size of the plots' titles.                                   |
    |     x_ticks      |        None        |         | Int, list or dict defining the x_ticks to be displayed.      |
    |                  |                    |         | When int, number of ticks to be placed on each plot. When    |
    |                  |                    |         | list, it corresponds to de values used as ticks. When dict,  |
    |                  |                    |         | the keys must match the Chromosome values of the data,       |
    |                  |                    |         | while the values can be either int or list of int; when int  |
    |                  |                    |         | it corresponds to the number of ticks to be placed; when     |
    |                  |                    |         | list of int it corresponds to de values used as ticks. Note  |
    |                  |                    |         | that when the tick falls within a shrunk region it will not  |
    |                  |                    |         | be diplayed.                                                 |
    +------------------+--------------------+---------+--------------------------------------------------------------+
    |   arrow_color    |        grey        |         | Color of the arrow indicating strand.                        |
    | arrow_line_width |         1          |         | Line width of the arrow lines                                |
    |    arrow_size    |       0.006        |         | Float corresponding to the fraction of the plot or int       |
    |                  |                    |         | corresponding to the number of positions occupied by a       |
    |                  |                    |         | direction arrow.                                             |
    |   exon_height    |        0.6         |         | Height of the exon rectangle in the plot.                    |
    |   intron_color   |        None        |         | Color of the intron lines. When None, the color of the       |
    |                  |                    |         | first interval will be used.                                 |
    |     text_pad     |       0.005        |         | Space where the id annotation is placed beside the           |
    |                  |                    |         | interval. When text_pad is float, it represents the          |
    |                  |                    |         | percentage of the plot space, while an int pad represents    |
    |                  |                    |         | number of positions or base pairs.                           |
    |    text_size     |         10         |         | Fontsize of the text annotation beside the intervals.        |
    |     v_spacer     |        0.5         |         | Vertical distance between the intervals and plot border.     |
    +------------------+--------------------+---------+--------------------------------------------------------------+
    |   plotly_port    |        8050        |         | Port to run plotly app.                                      |
    | shrink_threshold |        0.01        |         | Minimum length of an intron or intergenic region in order    |
    |                  |                    |         | for it to be shrunk while using the “shrink” feature. When   |
    |                  |                    |         | threshold is float, it represents the fraction of the plot   |
    |                  |                    |         | space, while an int threshold represents number of           |
    |                  |                    |         | positions or base pairs.                                     |
    +------------------+--------------------+---------+--------------------------------------------------------------+

To reset options to built-in defaults,  use :func:`reset_options <pyranges_plot.reset_options>`.
By default, it will reset all options. Providing arguments, you can select which options to reset:

    >>> prp.reset_options('plot_background')  # reset one feature
    >>> prp.reset_options(['plot_border', 'title_color'])  # reset a few features
    >>> prp.reset_options()  # reset all features


Built-in and custom themes
--------------------------

A pyranges_plot **theme** is a collection of options for appearance customization (those displayed above
with :func:`print_options <pyranges_plot.print_options>`) each with a set value.
Themes are implemented as dictionaries, that are passed to the :func:`set_theme <pyranges_plot.set_theme>` function.
In practice, setting a theme is equivalent to setting options like we did above
with :func:`set_options <pyranges_plot.set_options>`, but with a single command.

For example, below we create a theme corresponding to the appearance of our last plot:

    >>> my_theme = {
    ...     "plot_bkg": "rgb(173, 216, 230)",
    ...     "plot_border": "#808080",
    ...     "title_color": "magenta"
    ... }
    >>> prp.set_theme(my_theme)
    >>> prp.plot(p)  # this will now open a plot identical to the previous one

Pyranges_plot comes with a few built-in themes, listed in the :func:`set_theme <pyranges_plot.set_theme>` function's
documentation. For example, here's the "dark" theme:

    >>> prp.set_theme('dark')
    >>> prp.plot(p)

.. @maxtico: please add this plot

To reset the theme, you can resort again to :func:`reset_options <pyranges_plot.reset_options>`.



Managing space: packed/unpacked, shrink
---------------------------------------

By default, pyranges_plot tries to save as much vertical space as possible,
so the transcripts are placed one beside the other, in a "packed" disposition.
To instead display one transcript per row, set the ``packed`` parameter as ``False``:

.. code-block::

    prp.plot(p, packed=False, legend = True)

.. image:: images/prp_rtd_09.png

.. @maxtico: please remove legend from text and replace the plot accordingly


Pyranges_plot offers the option to reduce horizontal space, occupied by introns or intergenic regions,
by activating the ``shrink`` parameter.
The  ``shrink_threshold`` determines the minimum length of a region without visible intervals to be shrunk.
When a float is provided, it will be interpreted as a fraction of the visible coordinate limits,
while when an int is given it will be interpreted as number of base pairs.

.. code-block::

    ppp = prp.example_data.p3
    print(ppp)


.. code-block::

    index    |    Chromosome    Strand    Start    End      transcript_id
    int64    |    object        object    int64    int64    object
    -------  ---  ------------  --------  -------  -------  ---------------
    0        |    1             +         90       92       t1
    1        |    1             +         61       64       t1
    2        |    1             +         104      113      t1
    3        |    1             +         228      229      t1
    ...      |    ...           ...       ...      ...      ...
    16       |    2             -         42       46       t5
    17       |    2             -         37       40       t5
    18       |    2             +         60       70       t6
    19       |    2             +         80       90       t6
    PyRanges with 20 rows, 5 columns, and 1 index columns.
    Contains 2 chromosomes and 2 strands.


.. code-block::

    prp.plot(ppp, shrink=True)

.. image:: images/prp_rtd_13.png

.. code-block::

    prp.plot(ppp, shrink=True, shrink_threshold=0.2)

.. image:: images/prp_rtd_14.png


Showing mRNA structure
----------------------

A familiar visualization to many bioinformaticians involves showing the mRNA structure with coding sequences (CDS)
displayed thicker than UTR (untranslated) regions. This is achieved by setting the ``thick_cds`` parameter to ``True``.
Note that data must be coded like standard GFF/GTF files,
with different rows for exons and for CDS, wherein CDS are subsets of exons. A "Feature" column must be present
and contain "exon" or "CDS" values:

.. code-block::

    pp = prp.example_data.p2
    print(pp)


.. code-block::

    index    |    Chromosome    Strand    Start    End      transcript_id    feature1    feature2    Feature
    int64    |    int64         object    int64    int64    object           object      object      object
    -------  ---  ------------  --------  -------  -------  ---------------  ----------  ----------  ---------
    0        |    1             +         1        11       t1               1           A           exon
    1        |    1             +         40       60       t1               1           A           exon
    2        |    2             -         10       25       t2               1           B           CDS
    3        |    2             -         70       80       t2               1           B           CDS
    ...      |    ...           ...       ...      ...      ...              ...         ...         ...
    10       |    4             -         30500    30700    t5               2           E           CDS
    11       |    4             -         30647    30700    t5               2           E           exon
    12       |    4             +         29850    29900    t6               2           F           CDS
    13       |    4             +         29970    30000    t6               2           F           CDS
    PyRanges with 14 rows, 8 columns, and 1 index columns.
    Contains 4 chromosomes and 2 strands.


.. code-block::

    prp.plot(pp, thick_cds=True)

.. image:: images/prp_rtd_12.png



Displaying multiple PyRanges objects
------------------------------------

In some cases, the data intervals might overlap. An example could be when some intervals in
the PyRanges object correspond to exons and others correspond to "GCA" appearances. For such
cases, the ``thickness_col`` and ``depth_col`` parameters are implemented.

The :func:`plot <pyranges_plot.plot>` function can accept more than one PyRanges object, provided as a list.
In this case, pyranges_plot will display them in the same plot, one on top of the other, for each common chromosome.
The intervals of different PyRanges object are separated by a vertical spacer.

Let's see an example with two PyRanges objects, mapping the occurrences of two amino acids, alanine and cysteine:

.. code-block::

    p_ala = prp.example_data.p_ala
    p_cys = prp.example_data.p_cys

    print(p_ala)
    print(p_cys)



.. code-block::

      index  |      Start      End    Chromosome  id        trait1    trait2      depth
      int64  |      int64    int64         int64  object    object    object      int64
    -------  ---  -------  -------  ------------  --------  --------  --------  -------
          0  |         10       20             1  gene1     exon      gene_1          0
          1  |         50       75             1  gene1     exon      gene_1          0
          2  |         90      130             1  gene1     exon      gene_1          0
          3  |         13       16             1  gene1     aa        Ala             1
          4  |         60       63             1  gene1     aa        Ala             1
          5  |         72       75             1  gene1     aa        Ala             1
          6  |        120      123             1  gene1     aa        Ala             1
    PyRanges with 7 rows, 7 columns, and 1 index columns.
    Contains 1 chromosomes.

      index  |      Start      End    Chromosome  id        trait1    trait2      depth
      int64  |      int64    int64         int64  object    object    object      int64
    -------  ---  -------  -------  ------------  --------  --------  --------  -------
          0  |         10       20             1  gene1     exon      gene_1          0
          1  |         50       75             1  gene1     exon      gene_1          0
          2  |         90      130             1  gene1     exon      gene_1          0
          3  |         15       18             1  gene1     aa        Cys             1
          4  |         55       58             1  gene1     aa        Cys             1
          5  |         62       65             1  gene1     aa        Cys             1
          6  |        100      103             1  gene1     aa        Cys             1
          7  |        110      113             1  gene1     aa        Cys             1
    PyRanges with 8 rows, 7 columns, and 1 index columns.
    Contains 1 chromosomes.



.. code-block::

    prp.plot([p_ala, p_cys])

.. @maxtico: please make this plot

When providing multiple PyRanges objects, it is useful to differentiate them in the plot. The ``y_labels`` parameter
allows to provide a list of strings, one for each PyRanges object, to be displayed on the left side of the plot:

.. code-block::

    prp.plot(
        [p_ala, p_cys],
        y_labels=["pr Alanine", "pr Cysteine"]
    )

.. @maxtico: make this plot

Customizing depth and thickness
-------------------------------

When dealing with overlapping intervals (e.g. see data above), the default visualization may fail to show
relevant information, because some intervals are hidden behind others. To address this, the
``depth_col`` parameter can be used to highlight overlapping intervals. This parameter accepts a
column name from the PyRanges object, which must contain integer values. The higher the value, the
closer the interval will be to the top of the plot, ensuring its visibility:

.. code-block::

    prp.plot(
        [p_ala, p_cys],
        id_col="id",
        y_labels=["pr Alanine", "pr Cysteine"],
        depth_col="depth"
    )

.. @maxtico: make this plot

Another way to highlight overlapping regions is by playing with the height (or thickness) of the blocks representing
intervals. This is achieved by using the ``thickness_col`` parameter, which defines a data column name whose values
determine thickness of the corresponding intervals:

.. code-block::
    prp.plot(
        [p_ala, p_cys],
        id_col="id",
        color_col="trait1",
        y_labels=["pr Alanine", "pr Cysteine"],
        thickness_col="trait1",
    )


.. image:: images/prp_rtd_11.png

.. @maxtico: replace this last plot (I changed the code but didn't update the plot)


Additional information: tooltips and titles
-------------------------------------------

In interactive plots there is the option of showing information about the gene when the
mouse is placed over its structure. This information always shows the gene's strand if
it exists, the start and end coordinates and the ID. To add information contained in other
dataframe columns to the tooltip, a string should be given to the ``tooltip`` parameter. This
string must contain the desired column names within curly brackets as shown below.

Similarly, the title of the chromosome plots can be customized giving the desired string to
the ``title_chr`` parameter, where the correspondent chromosome value of the data is referred
to as {chrom}. An example could be the following:

.. code-block::

    prp.plot(
        p,
        tooltip="first feature: {feature1}\nsecond feature: {feature2}",
        title_chr='Chr: {chrom}'
        )

.. image:: images/prp_rtd_10.png