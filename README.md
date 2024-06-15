# pyranges_plot
Gene visualization package for dataframe objects generated with [PyRanges](https://pyranges.readthedocs.io/en/latest/index.html).




## Overview
The goal is getting a plot displaying a series of genes, transcripts, or any kind
of ranges contained in a PyRanges object. It displays the genes' intron-exon structure 
in its corresponding chromosome, enabling easy visualization of your PyRanges data. The 
Pyranges version compatible with Pyranges Plot is >= 1.0.0.

To obtain the plot there are some features to be defined by the user, one is the 
**engine** since it can be based on Matplotlib or Plotly, the other is the name 
of the **gene ID** column in your data. The rest of features can either be left 
as default or be customized. In example, the plot shows the first 25 genes of the 
dataframe by default, but this can be modified. It is worth noting that the order 
of the genes will be conserved when performing the subset.

In the case of coloring, Pyranges Plot offers a wide versatility. The data feature 
(column) according to which the genes will be colored is by default the gene ID, but 
this "color column" can be selected manually. Color specifications can be left as the 
default colormap (``plotly.colors.qualitative.Alphabet``) or be provided as dictionaries, 
lists or color objects from either Matplotlib or Plotly regardless of the chosen engine. 
When a colormap or list of colors is specified, the colors assigned to the genes will 
iterate over the provided ones following the color column pattern. In the case of 
concrete color instructions such as dictionary, the genes will be colored according 
to it while the non-specified ones will be colored in black.

<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/general_ex.png">
</p>




## Installation
PyRanges-Plot can be installed using pip. To install all dependencies in order to be able to 
use all the functionalities of the package and both engines the `[all]` option must be 
specified:

```
pip install pyranges-plot[all]
```

If the user wishes to use only one of the engines, the installation af all dependencies 
can be avoided by using the engine-specific installation option:
```
# For matplotlib
pip install pyranges-plot[plt]

# For plotly
pip install pyranges-plot[plotly]
```

Note that the minimal installation by `pip install pyranges-plot` is not able to produce plots 
since the plot dependencies are not installed.



## Examples
### Minimal version
Next we will test pyranges_plot visualization options using the ``plot`` function. 
For that we will be using a PyRanges object generated from a dictionary.

```python
import pyranges as pr
import pyranges_plot as prp

p = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 2, 2, 2, 3],
                  "Strand": ["+", "+", "-", "-", "+", "+", "+", "+"],
                  "Start": [1, 40, 10, 70, 85, 110, 150, 140],
                  "End": [11, 60, 25, 80, 100, 115, 180, 152], 
                  "transcript_id":["t1", "t1", "t2", "t2", "t3", "t3", "t3", "t4"], 
                  "feature1": ["a", "a", "b", "b", "c", "c", "c", "d"], 
                  "feature2": ["A", "A", "B", "B", "C", "C", "C", "D"]})
print(p)

```
```
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
```


The generated data is a stranded PyRanges object containing 4 genes in 3 chromosomes 
as shown above. Having this example data stored in the variable ``p``, we are able to 
start exploring Pyranges Plot options. We can get a plot in a single line:



```python
prp.plot(p, engine="plt", id_col="transcript_id")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex01.png">
</p>



The output is an interactive Matplotlib plot. To obtain it we just need to provide the data, 
the engine and the name of the ID column. However, the engine and the ID column can be set 
previously so there is no need to specify them anymore while plotting:

```python
# As engine use 'plotly' or 'ply' for Plotly plots and 'matplotlib' or 'plt' for Matplotlib plots
prp.set_engine('plotly')
prp.set_id_col('transcript_id')
```

Now the plots will be based on Plotly because we set it as the engine, though they will look 
the same as the Matplotlib ones. Also, both libraries offer interactive zoom options. For 
Matplotlib…
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_fixex03.png">
</p>

and for Plotly.
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_fixex04.png">
</p>


### :left_right_arrow: Playing with limits

Since the data has only 4 genes all of them are plotted, but the function has a default limit 
of 25, so in a case where the data contains more genes it will only show the top 25, unless 
the ``max_ngenes`` parameter is specified. For example, we can set the maximum number of genes 
as 2. Note that in the case of plotting more than 25 a warning about the plot's integrity 
will appear.



```python
prp.plot(p, max_shown=2)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex02.png">
</p>



Another pyranges_plot functionality is allowing to define the plots' coordinate limits through 
the ``limits`` parameter. The default limits show some space between the first and last 
plotted exons of each chromosome, but these can be customized. The user can decide to change 
all or some of the coordinate limits leaving the rest as default if desired. The limits can 
be provided as a dictionary, tuple or PyRanges object:
 - Dictionary where the keys should be the data's chromosome names in string format and the 
 values can be either ``None`` or a tuple indicating the limits. When a chromosome is not 
 specified in the dictionary, or it is assigned ``None`` the coordinates will appear as default. 
 - Tuple option sets the limits of all plotted chromosomes as specified.
 - PyRanges object can also be used to define limits, allowing the visualization of one object's 
 genes in another object's range window.



```python
prp.plot(p, limits={1: (None, 100), 2: (60, 200), 3: None})
prp.plot(p, limits=(0,300))
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex05.png">
</p>
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex06.png">
</p>


### :rainbow: Coloring

We can try to color the genes according to the strand column instead of the ID (default). For 
that the ``color_col`` parameter should be used.



```python
prp.plot(p, color_col="Strand")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex07.png">
</p>



This way we see the "+" strand genes in one color and the "-" in another color. Additionally, 
these colors can be customized through the ``colormap`` parameter. For 
this case we can specify it as a dictionary in the following way:



```python
prp.plot(
    p,
    color_col="Strand",
    colormap={"+": "green", "-": "red"}
)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex08.png">
</p>



The parameter ``colormap`` is very versatile because it accepts dictionaries for specific coloring, 
but also Matplotlib and Plotly color objects such as colormaps (or even just the string name of 
these objects) as well as lists of colors. For example, we can use the Dark2 Matplotlib colormap, 
even if the plot is based on Plotly:



```python
prp.plot(p, colormap="Dark2")
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex09.png">
</p>

### :eyes: Display options

The disposition of the genes is by default a packed disposition, so the genes are 
preferentially placed one beside the other. But this disposition can be displayed 
as 'full' if the user wants to display one gene under the other by setting the 
``packed`` parameter as ``False``. Also, a legend can be added by setting the `legend`
parameter to `True`.



```python
prp.plot(p, packed=False, legend = True)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex10.png">
</p>



In interactive plots there is the option of showing information about the gene when the mouse is 
placed over its structure. This information always shows the gene's strand if it exists, the start and 
end coordinates and the ID. To add information contained in other dataframe columns to the tooltip, 
a string should be given to the ``tooltip`` parameter. This string must contain the desired column 
names within curly brackets as shown in the example. Similarly, the title of the chromosome plots can be customized giving the desired string to 
the `title_chr` parameter, where the correspondent chromosome value of the data is referred 
to as {chrom}. An example could be the following:

```python
prp.plot(
 p,
 tooltip="first feature: {feature1}\nsecond feature: {feature2}",
 title_chr='Chr: {chrom}'
)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex11.png">
</p>



#### :dizzy: Show transcript structure

Another interesting feature is showing the transcript structure, so the exons appear as 
wider rectangles than UTR regions. For that the proper information should be stored in the 
`"Feature"` column of the data. A usage example is:

```python
pp = pr.from_dict({
 "Chromosome": [1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4],
 "Strand": ["+", "+", "-", "-", "+", "+", "+", "+", "-", "-", "-", "-", "+", "+"],
 "Start": [1, 40, 10, 70, 85, 110, 150, 140, 30100, 30150, 30500, 30647, 29850, 29970],
 "End": [11, 60, 25, 80, 100, 115, 180, 152, 30300, 30300, 30700, 30700, 29900, 30000],
 "transcript_id": ["t1", "t1", "t2", "t2", "t3", "t3", "t3", "t4", "t5", "t5", "t5", "t5", "t6", "t6"],
 "feature1": ["1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2"],
 "feature2": ["A", "A", "B", "B", "C", "C", "C", "D", "E", "E", "E", "E", "F", "F"],
 "Feature": ["exon", "exon", "CDS", "CDS", "CDS", "CDS", "CDS", "exon", "exon", "CDS", "CDS", "exon", "CDS", "CDS"]

})

prp.plot(pp, thick_cds=True)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex12.png">
</p>

#### :dizzy: Reduce intron size

In order to facilitate visualization, pyranges_plot offers the option to reduce the introns 
which exceed a given threshold size. For that the `shrink` parameter should be used. 
Additionally, the threshold can be defined by the user through kargs or setting the default options 
as explained in the next section, using `shrink_threshold`, when a float is provided as 
shrink_threshold it will be interpreted as a fraction of the original coordinate range, 
while when an int is given it will be interpreted as number of base pairs.

```python
ppp = pr.from_dict({'Chromosome': ['1'] * 10 + ['2'] * 10,
                    'Strand': ['+', '+', '+', '+', '-', '-', '-', '-', '+', '+'] + ["+", "+", "+", "+", "-", "-", "-",
                                                                                    "-", "+", "+"],
                    'Start': [90, 61, 104, 228, 9, 142, 52, 149, 218, 151] + [5, 27, 37, 47, 1, 7, 42, 37, 60, 80],
                    'End': [92, 64, 113, 229, 12, 147, 57, 155, 224, 153] + [8, 32, 40, 50, 5, 10, 46, 40, 70, 90],
                    'transcript_id': ['t1', 't1', 't1', 't1', 't2', 't2', 't2', 't2', 't3', 't3'] + ["t4", "t4", "t4",
                                                                                                     "t4", "t5", "t5",
                                                                                                     "t5", "t5", "t6",
                                                                                                     "t6"],
                    'Feature': ["exon"] * 20
                    }
                   )

prp.plot(ppp, shrink=True)
prp.plot(ppp, shrink=True, shrink_threshold=0.2)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex13.png">
</p>
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex14.png">
</p>


### :ribbon: Appearance customizations

There are some features of the plot appearance which can also be customized, like the 
background, plot border or titles. To check these customizable features and its default 
options values, the `print_options` function should be used. These values con be modified for all 
the following plots through the `set_options` function; However, for a single plot, these 
features can be given as kargs to the `plot` function (see `shrink_threshold` in the example 
above).

```python
# Check the default options values
prp.print_options()
```
```
+------------------------+-------+---------+--------------------------------------------------------------+
|        Feature         | Value | Edited? |                         Description                          |
+------------------------+-------+---------+--------------------------------------------------------------+
|     tag_background     | grey  |         | Background color of the tooltip annotation for the gene in   |
|                        |       |         | Matplotlib.                                                  |
|    plot_background     | white |         | Background color for the chromosomes plots.                  |
|      plot_border       | black |         | Color of the line defining the chromosome plots.             |
|       title_size       |  18   |         | Size of the plots' titles.                                   |
|      title_color       | black |         | Color of the plots' titles.                                  |
+------------------------+-------+---------+--------------------------------------------------------------+
|       exon_width       |  0.4  |         | Height of the exon rectangle in the plot.                    |
|    arrow_line_width    |   1   |         | Line width of the arrow lines (for stranded PyRanges).       |
|      arrow_color       | grey  |         | Direction arrow color (for stranded PyRanges).               |
|       arrow_size       | 0.006 |         | Fraction or percentage of the plot occupied by a direction   |
|                        |       |         | arrow.                                                       |
| arrow_intron_threshold | 0.04  |         | Minimum size of the intron to plot a direction arrow in it.  |
|                        |       |         | Provided as a float correspondig to the plot fraction or     |
|                        |       |         | percentage.                                                  |
+------------------------+-------+---------+--------------------------------------------------------------+
|    shrink_threshold    | 0.05  |         | Minimum lenght of an intron in order for it to be shrunk   |
|                        |       |         | while using the shrink feature. When threshold is       |
|                        |       |         | float, it represents the percentage of the plot space,       |
|                        |       |         | while an int threshold represents number of positions or     |
|                        |       |         | base pairs.                                                  |
|      plotly_port       | 8050  |         | Port to run plotly app.                                      |
+------------------------+-------+---------+--------------------------------------------------------------+


```

Once you found the feature you would like to customize, it can be modified:

```python

# Change the default options values
prp.set_options('plot_background', 'rgb(173, 216, 230)')
prp.set_options('plot_border', '#808080')
prp.set_options('title_color', 'magenta')

# Make the customized plot
prp.plot(p)
```
<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/prplot_ex15.png">
</p>


Now the modified values will be marked when checking the options values:

```python
prp.print_options()
```
```
+------------------------+--------------------+---------+--------------------------------------------------------------+
|        Feature         |       Value        | Edited? |                         Description                          |
+------------------------+--------------------+---------+--------------------------------------------------------------+
|     tag_background     |        grey        |         | Background color of the tooltip annotation for the gene in   |
|                        |                    |         | Matplotlib.                                                  |
|    plot_background     | rgb(173, 216, 230) |    *    | Background color for the chromosomes plots.                  |
|      plot_border       |      #808080       |    *    | Color of the line defining the chromosome plots.             |
|       title_size       |         18         |         | Size of the plots' titles.                                   |
|      title_color       |      magenta       |    *    | Color of the plots' titles.                                  |
+------------------------+--------------------+---------+--------------------------------------------------------------+
|       exon_width       |        0.4         |         | Height of the exon rectangle in the plot.                    |
|    arrow_line_width    |         1          |         | Line width of the arrow lines (for stranded PyRanges).       |
|      arrow_color       |        grey        |         | Direction arrow color (for stranded PyRanges).               |
|       arrow_size       |       0.006        |         | Fraction or percentage of the plot occupied by a direction   |
|                        |                    |         | arrow.                                                       |
| arrow_intron_threshold |        0.04        |         | Minimum size of the intron to plot a direction arrow in it.  |
|                        |                    |         | Provided as a float correspondig to the plot fraction or     |
|                        |                    |         | percentage.                                                  |
+------------------------+--------------------+---------+--------------------------------------------------------------+
|    shrink_threshold    |        0.05        |         | Minimum lenght of an intron in order for it to be shrunk   |
|                        |                    |         | while using the shrink feature. When threshold is       |
|                        |                    |         | float, it represents the percentage of the plot space,       |
|                        |                    |         | while an int threshold represents number of positions or     |
|                        |                    |         | base pairs.                                                  |
|      plotly_port       |        8050        |         | Port to run plotly app.                                      |
+------------------------+--------------------+---------+--------------------------------------------------------------+

```

To return to the original appearance of the plot, the `reset_options` function can restore 
all or some parameters. By default, it will reset all the features, but it also accepts a 
string for resetting a single feature or a list of strings to reset a few.

```python
prp.reset_options()  # reset all
prp.reset_options('plot_background')  # reset one feature
prp.reset_options(['plot_border', 'title_color'])  # reset a few features
```



Once we are able to get the plot we want, it can be exported to pdf or png format using the 
``to_file`` parameter. This parameter takes a string with the name or path of the file including
its extension. Additionally, the size can be customized through the ``file_size`` parameter by 
providing a tuple containing the height and width values.



```python
# Build the plot and save it in pdf or png
prp.plot(p, to_file='my_plot.pdf', file_size=(1300, 600))

# An example of some pyranges adjustments and save
p_subset = p[p.transcript_id.isin(['t3', 't4'])]
prp.plot(p_subset, colormap='Set3', to_file='t3_t4_plot.png')
```

<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/t3_t4_plot.png">
</p>



## Coming soon
* Accept different PyRanges objects or DataFrames as input for the same plot.
* Bases will be displayed along coordinates.
* Colorblind friendly.
