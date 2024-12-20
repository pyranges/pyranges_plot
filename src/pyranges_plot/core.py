import pandas as pd
import importlib
from pyranges.core.names import END_COL

from .names import CUM_DELTA_COL
from .plot_features import (
    plot_features_dict,
    plot_features_dict_in_use,
    builtin_themes,
)


def check4dependency(name):
    """Check if a module is installed"""
    try:
        importlib.import_module(name)
        return True
    except ImportError:
        return False


# CORE FUNCTIONS
# id_col
ID_COL = None


def set_id_col(name):
    """
    Defines the ID column for the data.

    Parameters
    ----------
    name: str

         Indicates the name of the ID column to be used when dealing with data.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_id_col('gene_id')

    """

    global ID_COL
    ID_COL = name


def get_id_col():
    """Returns the currently defined ID column (id_col)."""

    return ID_COL


# engine
ENGINE = None


def set_engine(name):
    """
    Defines the engine for the plots.

    Parameters
    ----------
    name: str
        Indicates if Matplotlib ('plt', 'matplotlib') or Plotly ('ply', 'plotly') should be used.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_engine('plt')

    """

    global ENGINE
    ENGINE = name


def get_engine():
    """Returns the currently defined engine."""

    return ENGINE


# warnings
WARNINGS = True


def set_warnings(option):
    """
    Defines if the warnings should be shown.

    Parameters
    ----------
    option: bool

         True for showing the warnings and false to turn them off.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_warnings(False)

    """

    global WARNINGS
    WARNINGS = option


def get_warnings():
    """Returns the current warnings state."""

    return WARNINGS


theme = None


def set_theme(name):
    """
    Defines the color theme for the plots.

    Parameters
    ----------
    name: {str, dict, None}
        Name of the predefined theme or dictionary with defined options to be set as new default.
        Currently available themes are "dark", "light", "Mariotti_lab" and "swimming_pool".

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_theme("dark")
    >>> prp.set_theme({"title_color": "goldenrod", "exon_height": 0.8})

    """

    global theme
    theme = name

    if name is None:
        return

    if isinstance(theme, str):
        if theme not in builtin_themes.keys():
            raise Exception(
                f'The name "{theme}" is not a valid theme name. Accepted themes are: {builtin_themes.keys()}'
            )
        else:
            name = builtin_themes[theme]

    if isinstance(name, dict):
        for key, value in name.items():
            # is it different from default?
            mod_tag = " "
            if name[key] != plot_features_dict[key][0]:
                mod_tag = "*"

            plot_features_dict_in_use[key] = (
                value,
                plot_features_dict[key][1],
                mod_tag,
            )  # (value, description, modified tag)


def get_theme():
    """Returns the currently defined color theme."""

    return theme


# Related to default features (options)


def set_options(varname, value=None):
    """
    Define some feature options of the plot layout.

    Parameters
    ----------
    varname: {str, dict}

        Name of the variable to change, or dictionary with the {variable: value} pairs.
        To check the available customizable options use prp.print_default().

    value:

        New value of the variable to be assigned to varname if needed.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_options('plot_background', 'magenta')

    >>> prp.set_options('title_size', 20)

    >>> prp.set_options({'plot_background': 'magenta', 'title_size': 20})

    """

    if isinstance(varname, str):
        varname = {varname: value}

    for key, val in varname.items():
        mod_tag = " "
        if varname[key] != plot_features_dict[key][0]:
            mod_tag = "*"
        plot_features_dict_in_use[key] = (
            val,
            plot_features_dict[key][1],
            mod_tag,
        )  # (value, description, modified tag)


def get_options(varname="all"):
    """
    Obtain the deafault value for a plot layout variable/s and its description.

    Parameters
    ----------
    varname: {str, list}, default 'all'

        Name of the variable/s to get the value and description. Use "values" to get
        only the variables values excluding the description and modified tag.

    """

    # list of variables
    if isinstance(varname, list):
        vars_list = []
        for var in varname:
            vars_list.append(plot_features_dict_in_use[var][0])
        return vars_list

    # all variables
    elif varname == "all":
        return plot_features_dict_in_use
    elif varname == "values":
        val_features_dict_in_use = {}
        for key, val in plot_features_dict_in_use.items():
            val_features_dict_in_use[key] = val[0]
        return val_features_dict_in_use

    # one variable
    else:
        if varname in plot_features_dict_in_use:
            return plot_features_dict_in_use[varname][0]
        else:
            raise Exception(
                f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}"
            )


def get_original_options():
    """Returns the dictionary with the original plot features."""

    return plot_features_dict


def reset_options(varname="all"):
    """
    Reset the deafault value for one, some or all plot layout variables to their original vlaue.

    Parameters
    ----------
    varname: {str, list}, default 'all'

        Name of the variable/ to reset the value.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.reset_options()

    >>> prp.reset_options('all')

    >>> prp.reset_options('tag_bkg')

    >>> prp.reset_options(['title_size', 'tag_background'])

    """

    plot_features_dict_in_use = get_options()
    plot_features_dict = get_original_options()

    # list of variables
    if isinstance(varname, list):
        for var in varname:
            plot_features_dict_in_use[var] = plot_features_dict[var]

    # all variables
    elif varname == "all":
        for var in plot_features_dict_in_use.keys():
            plot_features_dict_in_use[var] = plot_features_dict[var]

    # one variable
    else:
        try:
            if varname in plot_features_dict_in_use.keys():
                plot_features_dict_in_use[varname] = plot_features_dict[varname]
            else:
                raise Exception(
                    f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}"
                )
        except SystemExit as e:
            print("An error occured:", e)


def divide_desc(desc, cutoff):
    """Divide long feature description in lines."""

    lines_l = []
    while len(desc) > cutoff:
        for i in range(59, -1, -1):
            if desc[i] == " ":
                lines_l.append(desc[:i])
                desc = desc[i + 1 :]
                break
    lines_l.append(desc)

    return lines_l


def print_options(return_keys=False):
    """Prints the customizable features default values and description."""

    # store data
    plot_features_dict_in_use = get_options()

    # prepare data to print
    if not return_keys:
        feat_df = pd.DataFrame.from_dict(
            plot_features_dict_in_use,
            orient="index",
            columns=["Value", "Description", "Modified"],
        )

        # Calculate column sizes
        name_sz = max([len(val) for val in plot_features_dict_in_use])
        value_sz = max([len(str(val)) for val in feat_df["Value"]])
        if value_sz < 5:  # value has a minimum of 5
            value_sz = 5
        mod_sz = 7  # according to "Edited?" length
        desc_sz = 60

        # Function to format row
        def format_row(key, value):
            if len(value.iloc[1]) <= 60:
                return f"| {key:^{name_sz}} | {str(value.iloc[0]):^{value_sz}} | {value.iloc[2]:^{mod_sz}} | {value.iloc[1]:<{desc_sz}} |"

            else:
                lines_l = divide_desc(value.iloc[1], cutoff=desc_sz)
                fstr = f"| {key:^{name_sz}} | {str(value.iloc[0]):^{value_sz}} | {value.iloc[2]:^{mod_sz}} | {lines_l[0]:<{desc_sz}} |"
                empty = " "
                for i in range(1, len(lines_l)):
                    fstr += f"\n| {empty:^{name_sz}} | {empty:^{value_sz}} | {empty:^{mod_sz}} | {lines_l[i]:<{desc_sz}} |"

                return fstr

        # Create table header
        header = f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+\n"
        header += f"| {'Feature':^{name_sz}} | {'Value':^{value_sz}} | {'Edited?':^{mod_sz}} | {'Description':^{desc_sz}} |\n"
        header += f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"

        # Divide features
        extragen_feat_df = feat_df[
            feat_df.index.isin(
                [
                    "colormap",
                    "tag_bkg",
                    "fig_bkg",
                    "plot_bkg",
                    "plot_border",
                    "title_size",
                    "title_color",
                    "title_font",
                    "grid_color",
                    "exon_border",
                    "shrunk_bkg",
                    "x_ticks",
                ]
            )
        ].copy()

        intragen_feat_df = feat_df[
            feat_df.index.isin(
                [
                    "exon_height",
                    "v_spacer",
                    "text_size",
                    "text_pad",
                    "arrow_line_width",
                    "arrow_color",
                    "arrow_size",
                    "intron_color",
                ]
            )
        ].copy()

        other_feat_df = feat_df[
            feat_df.index.isin(["shrink_threshold", "plotly_port"])
        ].copy()

        # Create table rows
        rows_eg = "\n".join(
            [format_row(key, value) for key, value in extragen_feat_df.iterrows()]
        )
        rows_ig = "\n".join(
            [format_row(key, value) for key, value in intragen_feat_df.iterrows()]
        )
        rows_o = "\n".join(
            [format_row(key, value) for key, value in other_feat_df.iterrows()]
        )

        # Print table
        print(header)
        print(rows_eg)
        print(
            f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"
        )
        print(rows_ig)
        print(
            f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
        )
        print(rows_o)
        print(
            f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
        )

    if return_keys:
        return set(plot_features_dict_in_use.keys())


def cumdelting(num_l, ts_data, chrom):
    """Update a list of coordinates according to cumdelta."""

    for i in range(len(num_l)):
        cdel = 0
        # get proper cumdelta
        for ix, row in ts_data[chrom].iterrows():
            if row[END_COL] <= num_l[i]:
                cdel = row[CUM_DELTA_COL]
            else:
                break
        num_l[i] -= cdel

    return num_l
