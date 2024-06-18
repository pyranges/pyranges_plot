Tutorial
~~~~~~~~

Getting started
---------------

The first step to obtain a plot is always setting the **engine**. The way to do it is using
the ``set_engine`` function after importing.::
    import pyranges_plot as prp

    prp.set_engine("plotly")

Similarly, some other variables can be set prior to the plot call, like ``id_col``,
``warnings`` and ``theme``; though unlike engine, they can be given as parameters to
the :code:`plot`` function.

Pyranges Plot evolves around the :code:`plot` function, which admits output definition
through its parameters and appearance customization options through ``kargs``. To showcase
its functionalities we will load some example data using a dictionary, however Pyranges
includes a series of data loading options like gff, gtf, bam... (take a look at `Pyranges
documentation<https://pyranges1.readthedocs.io/en/latest/>`_ to know more!)::
    import pyranges as pr

    p = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 2, 2, 2, 3],
                    "Strand": ["+", "+", "-", "-", "+", "+", "+", "+"],
                    "Start": [1, 40, 10, 70, 85, 110, 150, 140],
                    "End": [11, 60, 25, 80, 100, 115, 180, 152],
                    "transcript_id":["t1", "t1", "t2", "t2", "t3", "t3", "t3", "t4"],
                    "feature1": ["a", "a", "b", "b", "c", "c", "c", "d"],
                    "feature2": ["A", "A", "B", "B", "C", "C", "C", "D"]})
    print(p)

.. code-block::
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

Once the set up is ready, a minimal plot can be obtained with just:
.. code-block::
    prp.plot(p)

.. image:: images/prp_rtd_01.png