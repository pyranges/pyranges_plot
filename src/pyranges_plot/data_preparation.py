import numpy as np
from intervaltree import IntervalTree
import pyranges as pr
from pyranges.core.names import CHROM_COL, START_COL, END_COL

# Check for matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    missing_plt_flag = 0
except ImportError:
    missing_plt_flag = 1
# Check for plotly
try:
    import plotly.colors as pc

    missing_ply_flag = 0
except ImportError:
    missing_ply_flag = 1


from .names import (
    PR_INDEX_COL,
    SHRTHRES_COL,
    TEXT_PAD_COL,
    COLOR_INFO,
    COLOR_TAG_COL,
    BORDER_COLOR_COL,
)
from .core import cumdelting, get_engine, get_warnings, check4dependency
from .matplotlib_base.core import plt_popup_warning


############ COMPUTE INTRONS OFF THRESHOLD
def compute_thresh(df, chrmd_df_grouped):
    """Get shrink threshold from limits"""

    chrom = df[CHROM_COL].iloc[0]
    chrmd = chrmd_df_grouped.loc[chrom]
    limit_range = chrmd["max"] - chrmd["min"]
    df[SHRTHRES_COL] = [int(df[SHRTHRES_COL].iloc[0] * limit_range)] * len(df)

    return df


############ COMPUTE TEXT PAD SIZE
def compute_tpad(df, chrmd_df_grouped):
    """Get text pad size from limits"""

    chrom = df[CHROM_COL].iloc[0]
    chrmd = chrmd_df_grouped.loc[chrom]
    limit_range = chrmd["max"] - chrmd["min"]
    df[TEXT_PAD_COL] = [int(df[TEXT_PAD_COL].iloc[0] * limit_range)] * len(df)

    return df


############ SUBSET
def make_subset(df, id_col, max_shown):
    """Reduce the number of genes to work with."""

    # create a column indexing all the genes in the df
    df["gene_index"] = df.groupby(id_col, group_keys=False).ngroup()
    tot_ngenes = max(df["gene_index"])

    # select maximum number of genes
    if max(df.gene_index) + 1 <= max_shown:
        subdf = df
    else:
        subdf = df[df.gene_index < max_shown]

    # remove the gene_index column from the original df
    df.drop("gene_index", axis=1, inplace=True)

    return subdf, tot_ngenes


############ GENESMD_DF


###packed
def genesmd_packed(genesmd_df):
    """xxx"""

    # Initialize IntervalTree and used y-coordinates list
    trees = [IntervalTree()]

    def find_tree(row):
        for tree in trees:
            if not tree.overlaps(row[START_COL], row[END_COL]):
                return tree
        trees.append(IntervalTree())
        return trees[-1]

    # Assign y-coordinates
    for idx, row in genesmd_df.iterrows():
        tree = find_tree(row)
        tree.addi(row[START_COL], row[END_COL], idx)
        genesmd_df.at[idx, "ycoord"] = trees.index(tree)

    return genesmd_df


def update_y(genesmd_df, exon_height, v_spacer):
    """Update y coords according to previous prs"""

    # Consider pr dividing lines spot and the height of the previous pr to update y coords
    y_prev_df = (
        genesmd_df.groupby(PR_INDEX_COL)["ycoord"]
        .max()
        .shift(-1, fill_value=-(exon_height + v_spacer * 2))
        .apply(lambda x: x + (exon_height + v_spacer * 2))
        .loc[::-1]
        .cumsum()[::-1]
    )
    y_prev_df.name = "update_y_prev"
    genesmd_df = genesmd_df.join(y_prev_df, on=PR_INDEX_COL)
    genesmd_df["ycoord"] += genesmd_df["update_y_prev"]

    return genesmd_df


###colors for genes
def is_pltcolormap(colormap_string):
    """Checks whether the string given is a valid plt colormap name."""

    if check4dependency("matplotlib"):
        try:
            colormap = plt.colormaps[colormap_string]
            if colormap is not None and isinstance(colormap, mcolors.Colormap):
                return True
            else:
                return False

        except KeyError:
            return False

    else:
        return False


def is_plycolormap(colormap_string):
    """Checks whether the string given is a valid plotly color object name."""

    if check4dependency("plotly"):
        if hasattr(pc.sequential, colormap_string):
            return True
        elif hasattr(pc.diverging, colormap_string):
            return True
        elif hasattr(pc.cyclical, colormap_string):
            return True
        elif hasattr(pc.qualitative, colormap_string):
            return True


def get_plycolormap(colormap_string):
    """Provides the plotly color object corresponding to the string given."""

    if hasattr(pc.sequential, colormap_string):
        return getattr(pc.sequential, colormap_string)
    elif hasattr(pc.diverging, colormap_string):
        return getattr(pc.diverging, colormap_string)
    elif hasattr(pc.cyclical, colormap_string):
        return getattr(pc.cyclical, colormap_string)
    elif hasattr(pc.qualitative, colormap_string):
        return getattr(pc.qualitative, colormap_string)


def subdf_assigncolor(subdf, colormap, color_col, exon_border):
    """Add color information to data."""

    # Create COLOR_COL column
    if len(color_col) > 1:
        subdf[COLOR_TAG_COL] = list(zip(*[subdf[c] for c in color_col]))
    else:
        subdf[COLOR_TAG_COL] = subdf[color_col[0]]

    # Assign colors to
    color_tags = subdf[COLOR_TAG_COL].drop_duplicates()
    n_color_tags = len(color_tags)

    # 0-string to colormap object if possible
    if isinstance(colormap, str):
        if is_pltcolormap(colormap):
            colormap = plt.get_cmap(colormap)
        elif is_plycolormap(colormap):
            colormap = get_plycolormap(colormap)
        else:
            raise Exception(
                "The provided string does not match any installed dependency colormap."
            )

    # 1-plt colormap to list
    if not missing_plt_flag and isinstance(colormap, mcolors.ListedColormap):
        colormap = list(colormap.colors)  # colors of plt object
        colormap = [
            "rgb({}, {}, {})".format(int(r * 255), int(g * 255), int(b * 255))
            for r, g, b in colormap
        ]  # compatible with plotly

    # 2-list to dict
    if isinstance(colormap, list):
        # adjust number of colors
        if n_color_tags > len(colormap):
            engine = get_engine()
            warnings = get_warnings()
            if engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "The genes are colored by iterating over the given color list."
                )
            elif engine in ["ply", "plotly"] and warnings:
                subdf["_iterwarning!"] = [1] * len(subdf)
        else:
            colormap = colormap[:n_color_tags]
        # make plotly rgb colors compatible with plt
        # if colormap[0][:3] == "rgb":
        #     numb_list = [
        #         rgb[rgb.find("(") + 1 : rgb.find(")")].split(",") for rgb in colormap
        #     ]
        #     colormap = [
        #         "#{:02x}{:02x}{:02x}".format(int(r), int(g), int(b))
        #         for r, g, b in numb_list
        #     ]
        # create dict of colors
        colormap = {
            str(color_tags.iloc[i]): colormap[i % len(colormap)]
            for i in range(n_color_tags)
        }

    # 3- Use dict to assign color to gene
    if isinstance(colormap, dict):
        subdf[COLOR_TAG_COL] = subdf[COLOR_TAG_COL].astype(str)
        subdf[COLOR_INFO] = subdf[COLOR_TAG_COL].map(colormap)

        # add black genes warning if needed
        if subdf[COLOR_INFO].isna().any():
            engine = get_engine()
            warnings = get_warnings()
            if engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "Some genes do not have a color assigned so they are colored in black."
                )
            elif engine in ["ply", "plotly"] and warnings:
                subdf["_blackwarning!"] = [1] * len(subdf)
            subdf[COLOR_INFO].fillna("black", inplace=True)  # black for not specified

    if exon_border:
        subdf[BORDER_COLOR_COL] = [exon_border] * len(subdf)
    else:
        subdf[BORDER_COLOR_COL] = subdf[COLOR_INFO]

    return subdf


def get_genes_metadata(df, id_col, color_col, packed, exon_height, v_spacer):
    """Create genes metadata df."""

    # Start df with chromosome and the column defining color
    # Define the aggregation functions for each column
    agg_funcs = {
        col: "first"
        for col in id_col + color_col
        # if col not in [START_COL, END_COL, PR_INDEX_COL]
    }
    agg_funcs[START_COL] = "min"
    agg_funcs[END_COL] = "max"
    # workaround for Chromosome in color_col list
    if CHROM_COL in color_col:
        genesmd_df = (
            df.groupby(
                [CHROM_COL, PR_INDEX_COL] + id_col, group_keys=False, observed=True
            ).agg(agg_funcs)
            # .reset_index(level=[PR_INDEX_COL, CHROM_COL])
        )
        genesmd_df["chromosome"] = genesmd_df[CHROM_COL]
        for i in range(len(color_col)):
            if color_col[i] == CHROM_COL:
                color_col[i] = "chromosome"

    else:
        genesmd_df = (
            df.groupby(
                [CHROM_COL, PR_INDEX_COL] + id_col, group_keys=False, observed=True
            ).agg(agg_funcs)
            # .reset_index(level=[PR_INDEX_COL, CHROM_COL])
        )

    genesmd_df["chrix"] = genesmd_df.groupby(
        CHROM_COL, group_keys=False, observed=True
    ).ngroup()
    # Sort by pr_ix and chromosome
    genesmd_df.sort_values(by=[PR_INDEX_COL, "chrix"], inplace=True)
    genesmd_df["gene_ix_xchrom"] = genesmd_df.groupby(
        ["chrix", PR_INDEX_COL], group_keys=False, observed=True
    ).cumcount()

    # Assign y-coordinate to genes
    if packed:
        genesmd_df["ycoord"] = -1
        genesmd_df = genesmd_df.groupby(
            [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
        ).apply(genesmd_packed)  # add packed ycoord column
        genesmd_df.reset_index(level=CHROM_COL, inplace=True)
        genesmd_df = genesmd_df.groupby(CHROM_COL, observed=True).apply(
            lambda x: update_y(x, exon_height, v_spacer)
        )
        genesmd_df.drop(CHROM_COL, axis=1, inplace=True)

    else:
        # one gene in each height
        # genesmd_df.set_index(PR_INDEX_COL, inplace=True, append=True)
        genesmd_df["ycoord"] = (
            genesmd_df.sort_values(by=PR_INDEX_COL, ascending=False)
            .groupby(CHROM_COL, group_keys=False, observed=True)
            .cumcount()
        )

        # now create col to update according to prev pr height if needed
        # only one chromosome (no matter how many pr)
        if len(genesmd_df.index.get_level_values(CHROM_COL).drop_duplicates()) == 1:
            genesmd_df = genesmd_df.assign(
                upd_yc=genesmd_df.groupby([PR_INDEX_COL], group_keys=False).ngroup(
                    ascending=False
                )
            )

            # increase proper spacing to place pr_line
            genesmd_df["ycoord"] += genesmd_df["upd_yc"] * (-1)
            genesmd_df["ycoord"] += genesmd_df["upd_yc"] * (exon_height + 2 * v_spacer)

        # +1 chromosome and +1 pr
        elif len(genesmd_df.index.get_level_values(PR_INDEX_COL).drop_duplicates()) > 1:
            genesmd_df = genesmd_df.assign(
                upd_yc=genesmd_df.groupby(CHROM_COL, group_keys=False).apply(
                    lambda x: x.groupby(PR_INDEX_COL).ngroup(ascending=False)
                )
            )
            # increase proper spacing to place pr_line
            genesmd_df["ycoord"] += genesmd_df["upd_yc"] * (-1)
            genesmd_df["ycoord"] += genesmd_df["upd_yc"] * (exon_height + 2 * v_spacer)

    return genesmd_df


############ CHRMD_DF


##limits
def chrmd_limits(chrmd_df, limits):
    """Compute 'min_max' column for chromosome metadata"""

    # 1- create min_max column containing (plot min, plot max)

    # no limits no info
    if limits is None:
        chrmd_df["min_max"] = [(np.nan, np.nan)] * len(chrmd_df)

    # one tuple for all chromosomes
    elif type(limits) is tuple:
        chrmd_df["min_max"] = [limits] * len(chrmd_df)

    # pyranges object
    elif type(limits) is pr.PyRanges:
        # create dict to map limits
        limits_chrmd_df = limits.groupby(
            CHROM_COL, group_keys=False, observed=True
        ).agg({START_COL: "min", END_COL: "max"})
        # limits_chrmd_dict = limits_chrmd_df.to_dict(orient="index")

        # function to get matching values from limits_chrmd_df
        def make_min_max(row):
            chromosome = row.name[0]
            if chromosome in limits_chrmd_df.index:
                limits = limits_chrmd_df.loc[chromosome]

                return (
                    limits[START_COL],
                    limits[END_COL],
                )  # chromosome in both sets of data
            else:
                return (np.nan, np.nan)  # chromosome does not match

        # create limits column in plotting data
        chrmd_df["min_max"] = chrmd_df.apply(make_min_max, axis=1)

    # dictionary as limits
    else:
        chrmd_df["min_max"] = [
            limits.get(index)
            for index in list(chrmd_df.index.get_level_values(CHROM_COL))
        ]  # fills with None the chromosomes not specified


def fill_min_max(row, ts_data):
    """Complete min_max column for chromosome metadata if needed."""

    minmax_t = row["min_max"]
    # deal with empty rows
    if minmax_t is None:
        minmax_t = (np.nan, np.nan)

    # check both items and put default if necessary
    minmax_l = list(minmax_t)

    # add default to lower limit
    if minmax_l[0] is None or np.isnan(minmax_l[0]):
        minmax_l[0] = row["min"]
    # add default to higher limit
    if minmax_l[1] is None or np.isnan(minmax_l[1]):
        minmax_l[1] = row["max"]
    # consider introns off for higher limit
    else:
        if len(row) == 5:
            new_upper_lim = cumdelting([minmax_l[1]], ts_data, row.name[0], row.name[1])
            minmax_l[1] = new_upper_lim[0]

    # put plot coordinates in min_max
    row["min_max"] = minmax_l
    return row


def get_chromosome_metadata(
    df, limits, genesmd_df, packed, v_spacer, exon_height, ts_data=None
):
    """Create chromosome metadata df."""

    # Start df
    agg_funcs = {
        START_COL: "min",
        END_COL: "max",
        "__id_col_2count__": "nunique",
    }

    chrmd_df = df.groupby([CHROM_COL, PR_INDEX_COL], observed=True).agg(agg_funcs)
    chrmd_df.rename(
        columns={START_COL: "min", END_COL: "max", "__id_col_2count__": "n_genes"},
        inplace=True,
    )

    # Adjust limits in case +1 pr
    if len(df[PR_INDEX_COL].drop_duplicates()) > 1:
        chrmd_df["min"] = chrmd_df.groupby(CHROM_COL, group_keys=False, observed=True)[
            "min"
        ].transform("min")
        chrmd_df["max"] = chrmd_df.groupby(CHROM_COL, group_keys=False, observed=True)[
            "max"
        ].transform("max")

    # Add limits
    chrmd_limits(chrmd_df, limits)  # unknown limits are nan
    chrmd_df = chrmd_df.apply(lambda x: fill_min_max(x, ts_data), axis=1)

    chrmd_df_grouped = (
        chrmd_df.reset_index(level=PR_INDEX_COL)
        .groupby(CHROM_COL, group_keys=False, observed=True)
        .agg(
            {
                "min": "first",
                "max": "first",
                "min_max": "first",
                PR_INDEX_COL: ["size", list],
            }
        )
    )
    chrmd_df_grouped.columns = ["min", "max", "min_max", "n_pr_ix", "present_pr"]

    # Store plot y height
    chrmd_df_grouped = chrmd_df_grouped.join(
        genesmd_df.groupby([CHROM_COL], group_keys=False, observed=True)["ycoord"].max()
    )
    chrmd_df_grouped.rename(columns={"ycoord": "y_height"}, inplace=True)
    chrmd_df_grouped["y_height"] += (
        0.5 + exon_height / 2
    )  # the middle of the rectangle is +.5 of ycoord

    # Obtain the positions of lines separating pr objects
    chrmd_df = chrmd_df.join(
        genesmd_df.groupby([CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True)[
            "ycoord"
        ].max()
    )
    chrmd_df.rename(columns={"ycoord": "pr_line"}, inplace=True)
    chrmd_df["pr_line"] = chrmd_df.groupby(CHROM_COL, observed=True)["pr_line"].shift(
        -1, fill_value=-(0.5 + exon_height / 2 + v_spacer)
    )

    chrmd_df["pr_line"] += (
        0.5 + exon_height / 2 + v_spacer
    )  # midle of rectangle is +.5 of ycoord

    # Set chrom_ix to get the right association to the plot index
    chrmd_df_grouped["chrom_ix"] = chrmd_df_grouped.groupby(
        CHROM_COL, group_keys=False, observed=True
    ).ngroup()

    return chrmd_df, chrmd_df_grouped
