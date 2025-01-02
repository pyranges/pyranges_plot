import copy

ori_l = ["#f05f89", "#f0db36", "#7bc45f", "#5e4699", "#f7943a", "#537ebf", "#ee3a36"]
darker_l = ["#9b3c59", "#a28b22", "#4d6e3a", "#3c285f", "#a46327", "#345a7d", "#9c2523"]
lighter_l = [
    "#ffadc9",
    "#ffee76",
    "#a8e89a",
    "#816bb9",
    "#ffc56b",
    "#82b3ff",
    "#ff7a74",
]
prp_cmap = ori_l + lighter_l + darker_l

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
        "popart",
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
    "return_plot": (None, "Whether the plot is returned or not.", " "),
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
    "title_font": ("Arial", "Font of the plots' titles.", " "),
    "v_spacer": (0.5, "Vertical distance between the intervals and plot border.", " "),
    "x_ticks": (
        None,
        "Int, list or dict defining the x_ticks to be displayed. When int, number of ticks to be placed on each plot. When list, it corresponds to de values used as ticks. When dict, the keys must match the Chromosome values of the data, while the values can be either int or list of int; when int it corresponds to the number of ticks to be placed; when list of int it corresponds to de values used as ticks. Note that when the tick falls within a shrunk region it will not be diplayed.",
        " ",
    ),
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

# Mariotti theme
ori_l = ["#f05f89", "#f0db36", "#7bc45f", "#5e4699", "#f7943a", "#537ebf", "#ee3a36"]
darker_l = ["#9b3c59", "#a28b22", "#4d6e3a", "#3c285f", "#a46327", "#345a7d", "#9c2523"]
lighter_l = [
    "#ffadc9",
    "#ffee76",
    "#a8e89a",
    "#816bb9",
    "#ffc56b",
    "#82b3ff",
    "#ff7a74",
]

theme_mariotti = {
    "colormap": ori_l + lighter_l + darker_l,
    "shrunk_bkg": "#e7e0f5",
    "fig_bkg": "#fff5ee",
    "title_color": "#590000",
    "plot_border": "#4c644c",
}

# Swimming pool theme
theme_sp = {
    "fig_bkg": "#696969",
    "plot_bkg": "#71E2E8",
    "colormap": ["#0D61AF", "#B82C10", "white"],
    "shrunk_bkg": "#c6e6c6",
    "plot_border": "#011334",
    "title_color": "#011334",
}

# Store themes
builtin_themes = {
    "light": plot_features_dict_vals,
    "dark": theme_dark,
    "Mariotti_lab": theme_mariotti,
    "swimming_pool": theme_sp,
}
