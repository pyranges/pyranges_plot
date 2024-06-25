# pyranges_plot
Gene visualization package for dataframe objects generated with [PyRanges](https://pyranges1.readthedocs.io/).

PyRanges Plot significantly facilitates genomic data
interpretation by providing powerful, customizable, and user-friendly visualizations. It
effectively covers the gap between data manipulation and visualization, thereby accelerating
the analysis workflow in genomic research.

## Overview
The goal is getting a plot displaying a series of genes, transcripts, or any kind
of ranges contained in a PyRanges object. It displays the genes' intron-exon structure 
in its corresponding chromosome, enabling easy visualization of your PyRanges data. The 
Pyranges version compatible with Pyranges Plot is >= 1.0.0 (find it at https://github.com/pyranges/pyranges_1.x.git).

To obtain a plot, the variable `engine` must be specified by the user first. This variable 
defines the graphic library on which the plots will be based: the valid `engine` options 
are "matplotlib" or "plt" for Matplotlib and "plotly" or "ply" for Plotly. 

Every other functionality can be defined during the `plot` function call. These 
functionalities include the ID column to group the intervals belonging to the same item 
(transcript, gene, protein...), items disposition, coloring criteria and palette, labels 
and output form among others. The input for the `plot` fucntion is 1 or more PyRanges 
objects, and the output is by default an interactive plot with zooming options and tooltip 
information, but if desired the plot can be directly exported to a png or pdf file.



<p align="center">
    <img src="https://github.com/pyranges/pyranges_plot/blob/main/images/general_ex.png">
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
pip install pyranges-plot[matplotlib]

# For plotly
pip install pyranges-plot[plotly]
```

Note that the minimal installation by `pip install pyranges-plot` is not able to produce plots 
since the graphical dependencies are not installed.


## Documentation
Pyranges Plot documentation and tutorial can be found at https://pyranges-plot.readthedocs.io/.


## Coming soon
* Bases will be displayed along coordinates.
