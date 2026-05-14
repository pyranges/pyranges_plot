import pandas as pd
import importlib
from pyranges1.core.names import END_COL

from .names import CUM_DELTA_COL
from .plot_features import (
    plot_features_dict,
    plot_features_dict_in_use,
    builtin_themes,
)
from . import adapters


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
    >>> import pyrangeyes as pre

    >>> pre.set_id_col('gene_id')

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
    >>> import pyrangeyes as pre

    >>> pre.set_engine('plt')

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
    >>> import pyrangeyes as pre

    >>> pre.set_warnings(False)

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
        Currently available themes are "dark", "light", "pastel" and "swimming_pool".

    Examples
    --------
    >>> import pyrangeyes as pre

    >>> pre.set_theme("dark")
    >>> pre.set_theme({"title_color": "goldenrod", "interval_height": 0.8})

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
            if value != plot_features_dict[key][0]:
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


def set_options(varname=None, value=None, *, adapter=None, variable=None):
    """
    Define plot-layout options, or adapter options when ``adapter`` is given.

    Parameters
    ----------
    varname : str or dict, optional
        Plot option name to change, or a dictionary with ``{option: value}``
        pairs. Use :func:`print_options` to inspect available plot options.

        When ``adapter`` is provided, this is the adapter option name unless
        ``variable`` is also provided.

    value : object, optional
        New value assigned to ``varname`` or ``variable``.

    adapter : str, optional keyword-only
        Adapter whose options should be changed, for example ``"mRNA"``.

    variable : str or dict, optional keyword-only
        Alias for ``varname`` when setting adapter options. This makes calls
        such as ``set_options(adapter="mRNA", variable="utr_height", value=0.5)``
        explicit while preserving the original positional API.

    Examples
    --------
    >>> import pyrangeyes as pre

    >>> pre.set_options('plot_background', 'magenta')

    >>> pre.set_options('title_size', 20)

    >>> pre.set_options({'plot_background': 'magenta', 'title_size': 20})

    >>> pre.set_options(adapter='mRNA', variable='utr_height', value=0.5)

    """

    if adapter is not None:
        if variable is not None:
            varname = variable
        if varname is None:
            raise ValueError("Please provide variable when setting adapter options.")
        adapters.set_options(adapter, varname, value)
        return

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


def get_options(varname="all", *, adapter=None):
    """
    Obtain plot-layout options, or adapter options when ``adapter`` is given.

    Parameters
    ----------
    varname : str or list, default 'all'
        Option name(s) to retrieve.

        - ``"all"`` returns the full ``{option: (value, description, modified)}``
          mapping.
        - ``"values"`` returns only current values as ``{option: value}``.
        - A list returns option values in the same order.
        - A single option name returns one current value.

    adapter : str, optional keyword-only
        Adapter whose options should be read, for example ``"mRNA"``. Adapter
        functions use their ``DEFAULT`` sentinel to pull current values from
        this option store at runtime.

    """

    if adapter is not None:
        return adapters.get_options(adapter, varname)

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


def reset_options(varname="all", *, adapter=None):
    """
    Reset one, some, or all plot-layout options to their original value.

    When ``adapter`` is provided, reset options for that adapter instead.

    Parameters
    ----------
    varname : str or list, default 'all'
        Option name, list of names, or ``"all"``.

    adapter : str, optional keyword-only
        Adapter whose options should be reset, for example ``"mRNA"``.

    Examples
    --------
    >>> import pyrangeyes as pre

    >>> pre.reset_options()

    >>> pre.reset_options('all')

    >>> pre.reset_options('tag_bkg')

    >>> pre.reset_options(['title_size', 'tag_background'])

    >>> pre.reset_options(adapter='mRNA')

    """

    if adapter is not None:
        adapters.reset_options(adapter, varname)
        return

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


def _format_options_table(options_dict, *, feature_label="Feature"):
    feat_df = pd.DataFrame.from_dict(
        options_dict,
        orient="index",
        columns=["Value", "Description", "Modified"],
    )
    feat_df["Value"] = feat_df["Value"].map(
        lambda value: "<infer>" if value is None else value
    )

    name_sz = max([len(val) for val in options_dict] + [len(feature_label)])
    value_sz = max([len(str(val)) for val in feat_df["Value"]])
    if value_sz < 5:
        value_sz = 5
    mod_sz = 7
    desc_sz = 60

    def format_row(key, value):
        if len(value.iloc[1]) <= 60:
            return f"| {key:^{name_sz}} | {str(value.iloc[0]):^{value_sz}} | {value.iloc[2]:^{mod_sz}} | {value.iloc[1]:<{desc_sz}} |"

        lines_l = divide_desc(value.iloc[1], cutoff=desc_sz)
        fstr = f"| {key:^{name_sz}} | {str(value.iloc[0]):^{value_sz}} | {value.iloc[2]:^{mod_sz}} | {lines_l[0]:<{desc_sz}} |"
        empty = " "
        for i in range(1, len(lines_l)):
            fstr += f"\n| {empty:^{name_sz}} | {empty:^{value_sz}} | {empty:^{mod_sz}} | {lines_l[i]:<{desc_sz}} |"
        return fstr

    header = f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+\n"
    header += f"| {feature_label:^{name_sz}} | {'Value':^{value_sz}} | {'Edited?':^{mod_sz}} | {'Description':^{desc_sz}} |\n"
    header += f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
    rows = "\n".join([format_row(key, value) for key, value in feat_df.iterrows()])
    footer = f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
    print(header)
    print(rows)
    print(footer)


def print_options(return_keys=False, *, adapter=None):
    """Print customizable plot or adapter options.

    Parameters
    ----------
    return_keys : bool, default False
        If True, return the option-name set instead of printing a table.
    adapter : str, optional keyword-only
        Adapter whose options should be printed, for example ``"mRNA"``.
    """

    if adapter is not None:
        adapter_options = get_options(adapter=adapter)
        if not return_keys:
            _format_options_table(adapter_options, feature_label="Adapter option")
        else:
            return set(adapter_options.keys())
        return

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
        header = f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+\n"
        header += f"| {'Feature':^{name_sz}} | {'Value':^{value_sz}} | {'Edited?':^{mod_sz}} | {'Description':^{desc_sz}} |\n"
        header += f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"

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
                    "outline_color",
                    "shrunk_bkg",
                    "x_ticks",
                ]
            )
        ].copy()

        intragen_feat_df = feat_df[
            feat_df.index.isin(
                [
                    "interval_height",
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
            f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
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
