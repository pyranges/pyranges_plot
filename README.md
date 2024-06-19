# pyranges_plot
Gene visualization package for dataframe objects generated with [PyRanges](https://pyranges.readthedocs.io/en/latest/index.html).




## Overview
The goal is getting a plot displaying a series of genes, transcripts, or any kind
of ranges contained in a PyRanges object. It displays the genes' intron-exon structure 
in its corresponding chromosome, enabling easy visualization of your PyRanges data. The 
Pyranges version compatible with Pyranges Plot is >= 1.0.0 (find it at https://github.com/pyranges/pyranges_1.x.git).

To obtain the plot there are some features to be defined by the user, one is the 
**engine** since it can be based on Matplotlib or Plotly, the other is optional and 
refers to the name of the **gene ID** column in your data. The rest of features can 
either be left as default or be customized. In example, the plot shows the first 25 
genes of the dataframe by default, but this can be modified. 

In the case of coloring, Pyranges Plot offers a wide versatility. The data feature 
(column) according to which the genes will be colored is by default the gene ID, but 
this "color column" can be selected manually. Color specifications can be left as the 
default colormap or be provided as dictionaries, lists or color objects from either 
Matplotlib or Plotly regardless of the chosen engine. When a colormap or list of colors 
is specified, the colors assigned to the genes will iterate over the provided ones 
following the color column pattern. In the case of concrete color instructions such as 
dictionary, the genes will be colored according to it while the non-specified ones will 
be colored in black.

<p align="center">
    <img src="https://github.com/emunozdc/pyranges_plot/raw/main/images/general_ex.png">
</p>




## Installation
PyRanges-Plot can be installed using pip. To install all dependencies in order to be able to 
use all the functionalities of the package and both engines, the `[all]` option must be 
specified:

```
pip install pyranges-plot[all]
```

If the user wishes to use only one of the engines, the installation of all dependencies 
can be avoided by using the engine-specific installation options:
```
# For matplotlib
pip install pyranges-plot[plt]

# For plotly
pip install pyranges-plot[plotly]
```

Note that the minimal installation by `pip install pyranges-plot` is not able to produce plots 
since the graphical dependencies are not installed.


## Documentation
Pyranges Plot documentation and tutorial can be found at [readthedocs](https://pyranges-plot.readthedocs.io/en/latest/).


## Coming soon
* Bases will be displayed along coordinates.
