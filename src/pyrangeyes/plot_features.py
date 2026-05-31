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
        "Colors to assign to interval fills. Use 'direct' when fill_col/outline_col/label_color_col already contain literal colors. "
        "A dict channel mapping must have 'fill' and may also have 'outline' and 'label'; 'fill'/'outline' aliases reuse channels. "
        "Values may be Matplotlib/Plotly colormap names, color lists, value-to-color mappings, or quantitative specs.",
        " ",
    ),
    "outline_color": (
        None,
        "Fixed color for interval outlines. When None, outlines use the resolved interval fill colors.",
        " ",
    ),
    "interval_height": (
        0.6,
        "Default (and maximum) height of rendered interval blocks.",
        " ",
    ),
    "auto_height_px_per_unit": (
        80,
        "Pixels assigned to one vertical layout unit when pyrangeyes infers figure height automatically.",
        " ",
    ),
    "figure_bg": ("white", "Bakground color of the whole figure.", " "),
    "grid_color": ("lightgrey", "Color of x coordinates grid lines.", " "),
    "intron_color": (
        None,
        "Color of the intron lines. When None, the color of the first interval will be used.",
        " ",
    ),
    "track_bg": ("white", "Background color of the plots.", " "),
    "plot_border": ("black", "Color of the line delimiting the plots.", " "),
    "plotly_port": (8050, "Port to run plotly app.", " "),
    "return_plot": (None, "Whether the plot is returned or not.", " "),
    "shrink_threshold": (
        0.01,
        "Minimum length of an intron or intergenic region in order for it to be shrunk while using the “shrink” feature. When threshold is float, it represents the fraction of the plot space, while an int threshold represents number of positions or base pairs.",
        " ",
    ),
    "shrunk_bg": (
        "lightyellow",
        "Color of the shrunk region background.",
        " ",
    ),
    "squish_factor": (
        0.3,
        "Factor applied to rendered interval height and stacked-row spacing for tracks with squish=True.",
        " ",
    ),
    "tag_bkg": (
        "grey",
        "Background color of the tooltip annotation for the gene in Matplotlib.",
        " ",
    ),
    "label_pad": (
        1,
        "Space, in percent of the visible plot span, between interval labels and intervals. For example, label_pad=1 means 1%.",
        " ",
    ),
    "label_size": (12, "Fontsize of the text annotation beside the intervals.", " "),
    "label_color": (
        "black",
        "Fixed color of interval labels unless label_color_col or colormap['label'] maps them.",
        " ",
    ),
    "label_angle": (0, "Rotation angle of interval labels, in degrees.", " "),
    "label_position": (
        "above",
        "Position of interval labels: 'left', 'right', 'center', 'top'/'above', or 'bottom'/'below'.",
        " ",
    ),
    "label_fit": (
        True,
        "Whether text labels reserve space during pack layout to reduce overlaps.",
        " ",
    ),
    "title_color": ("black", "Color of panel titles.", " "),
    "title_size": (18, "Font size of panel titles.", " "),
    "title_font": ("Arial", "Font family of panel titles.", " "),
    "v_spacer": (0.25, "Vertical distance between the intervals and plot border.", " "),
    "x_ticks": (
        None,
        "Int, list or dict defining the x_ticks to be displayed. When int, number of ticks to be placed on each plot. When list, it corresponds to de values used as ticks. When dict, the keys must match the Chromosome values of the data, while the values can be either int or list of int; when int it corresponds to the number of ticks to be placed; when list of int it corresponds to de values used as ticks. Note that when the tick falls within a shrunk region it will not be diplayed.",
        " ",
    ),
}

print_option_categories = [
    (
        "General and plot appearance",
        [
            "figure_bg",
            "plot_border",
            "plotly_port",
            "return_plot",
        ],
    ),
    (
        "Panels and axes",
        [
            "grid_color",
            "shrunk_bg",
            "shrink_threshold",
            "v_spacer",
            "x_ticks",
            "auto_height_px_per_unit",
            "title_color",
            "title_size",
            "title_font",
        ],
    ),
    (
        "Tracks and interval appearance",
        [
            "track_bg",
            "interval_height",
            "squish_factor",
            "colormap",
            "outline_color",
            "intron_color",
            "arrow_color",
            "arrow_line_width",
            "arrow_size",
        ],
    ),
    (
        "Text options per-track",
        [
            "label_pad",
            "label_size",
            "label_color",
            "label_angle",
            "label_position",
            "label_fit",
            "tag_bkg",
        ],
    ),
]

# Normal (light theme)
plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
plot_features_dict_vals = {}
for key, val in plot_features_dict.items():
    plot_features_dict_vals[key] = val[0]

# Dark theme
theme_dark = {
    "colormap": "G10",
    "figure_bg": "#1f1f1f",
    "plot_border": "white",
    "label_color": "white",
    "title_color": "goldenrod",
    "track_bg": "grey",
    "grid_color": "darkgrey",
    "arrow_color": "lightgrey",
    "shrunk_bg": "lightblue",
}

# pastel default theme
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

theme_pastel = {
    "colormap": ori_l + lighter_l + darker_l,
    "shrunk_bg": "#e7e0f5",
    "figure_bg": "#fff5ee",
    "title_color": "#590000",
    "plot_border": "#4c644c",
}

# Swimming pool theme
theme_sp = {
    "figure_bg": "#696969",
    "track_bg": "#71E2E8",
    "colormap": ["#0D61AF", "#B82C10", "white"],
    "shrunk_bg": "#c6e6c6",
    "plot_border": "#011334",
    "title_color": "#011334",
}

# Store themes
builtin_themes = {
    "light": plot_features_dict_vals,
    "dark": theme_dark,
    "pastel": theme_pastel,
    "swimming_pool": theme_sp,
}
