import copy


plot_features_dict = {
    "arrow_color": ("grey", "Color of the arrow indicating strand.", " "),
    "arrow_line_width": (
        1,
        "Line width of the arrow lines",
        " ",
    ),
    "arrow_size": (
        0.006,
        "Float corresponding to the fraction of the plot or int corresponding to the number of positions occupied by a direction arrow.",
        " ",
    ),
    "colormap": (
        "Alphabet",
        "Sequence of colors to assign to every group of intervals sharing the same “color_col” value. It can be provided as a Matplotlib colormap, a Plotly color sequence (built as lists), a string naming the previously mentioned color objects from Matplotlib and Plotly, or a dictionary with the following structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific color_col value is not specified in the dictionary it will be colored in black.",
        " ",
    ),
    "exon_border": (None, "Color of the interval's rectangle border.", " "),
    "exon_height": (0.6, "Height of the exon rectangle in the plot.", " "),
    "fig_bkg": ("white", "Bakground color of the whole figure.", " "),
    "grid_color": ("lightgrey", "Color of x coordinates grid lines.", " "),
    "intron_color": (
        None,
        "Color of the intron lines. When None, the color of the first interval will be used.",
        " ",
    ),
    "plot_bkg": ("white", "Background color of the plots.", " "),
    "plot_border": ("black", "Color of the line delimiting the plots.", " "),
    "plotly_port": (8050, "Port to run plotly app.", " "),
    "shrink_threshold": (
        0.01,
        "Minimum length of an intron or intergenic region in order for it to be shrunk while using the “shrink” feature. When threshold is float, it represents the fraction of the plot space, while an int threshold represents number of positions or base pairs.",
        " ",
    ),
    "shrunk_bkg": (
        "lightyellow",
        "Color of the shrunk region background.",
        " ",
    ),
    "tag_bkg": (
        "grey",
        "Background color of the tooltip annotation for the gene in Matplotlib.",
        " ",
    ),
    "text_pad": (
        0.005,
        "Space where the id annotation is placed beside the interval. When text_pad is float, it represents the percentage of the plot space, while an int pad represents number of positions or base pairs.",
        " ",
    ),
    "text_size": (10, "Fontsize of the text annotation beside the intervals.", " "),
    "title_color": ("black", "Color of the plots' titles.", " "),
    "title_size": (18, "Size of the plots' titles.", " "),
    "v_spacer": (0.5, "Vertical distance between the intervals and plot border.", " "),
}

# Normal (light theme)
plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
plot_features_dict_vals = {}
for key, val in plot_features_dict.items():
    plot_features_dict_vals[key] = val[0]

# Dark theme
theme_dark = {
    "colormap": "G10",
    "fig_bkg": "#1f1f1f",
    "plot_border": "white",
    "title_color": "goldenrod",
    "plot_bkg": "grey",
    "grid_color": "darkgrey",
    "arrow_color": "lightgrey",
    "shrunk_bkg": "lightblue",
}

# Store themes
builtin_themes = {
    "light": plot_features_dict_vals,
    "dark": theme_dark,
}
