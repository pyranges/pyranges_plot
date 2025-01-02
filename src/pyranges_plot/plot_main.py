import pandas as pd
import numpy as np


# import pyranges as pr
from .core import (
    get_id_col,
    get_engine,
    print_options,
    get_options,
    get_warnings,
    set_theme,
    get_theme,
    set_options,
)
from .plot_features import prp_cmap
from .data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
    compute_thresh,
    compute_tpad,
    subdf_assigncolor,
)
from .introns_off import introns_resize, recalc_axis
from pyranges.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL
from .names import (
    PR_INDEX_COL,
    ORISTART_COL,
    ORIEND_COL,
    SHRTHRES_COL,
    ADJSTART_COL,
    ADJEND_COL,
    CUM_DELTA_COL,
    EXON_IX_COL,
    TEXT_PAD_COL,
    THICK_COL,
)

# Check for matplotlib
try:
    from .matplotlib_base.plot_exons_plt import plot_exons_plt

    missing_plt_flag = 0
except ImportError:
    missing_plt_flag = 1

# Check for plotly
try:
    from .plotly_base.plot_exons_ply import plot_exons_ply

    missing_ply_flag = 0
except ImportError:
    missing_ply_flag = 1


def plot(
    data,
    *,
    id_col=None,
    warnings=None,
    max_shown=25,
    packed=True,
    return_plot=None,
    add_aligned_plots=None,
    color_col=None,
    thickness_col=None,
    depth_col=None,
    shrink=False,
    limits=None,
    thick_cds=False,
    text=True,
    legend=False,
    title_chr="Chromosome {chrom}",
    y_labels=None,
    tooltip=None,
    to_file=None,
    theme=None,
    **kargs,
):
    """
    Create genes plot from 1/+ PyRanges objects.

    Parameters
    ----------
    data: {pyranges.PyRanges or list of pyranges.PyRanges}
        Pyranges, derived dataframe or list of them with annotation data.

    id_col: str, default None
        Name of the column containing gene ID.

    warnings: bool, default True
        Whether the warnings should be shown or not.

    max_shown: int, default 20
        Maximum number of genes plotted in the dataframe order.

    packed: bool, default True
        Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).

    color_col: str, default None
        Name of the column used to color the genes. If not specified, id_col will be used.

    thickness_col: str, default None
        Name of the data column with max 2 different values to plot the intervals correspondig to one value to
        thicker than the others. The first value by alphabetical order will have the height specified
        as 'exon_height', and the second will be 0.3*'exon_height'. Note that this parameter will be
        overseen if the 'thick_cds' parameter is set to True.

    depth_col: str, default None
        Name of the data column to be used for setting the order to plot the intervals. The intervals with
        the lowest value in this column will be plotted first and the ones with higher values will plotted
        on top of them.

    shrink: bool, default False
        Whether to compress the intron ranges to facilitate visualization or not.

    limits: {None, dict, tuple, pyranges.pyranges_main.PyRanges}, default None
        Customization of coordinates for the chromosome plots.\n
        - None: minimum and maximum exon coordinate plotted plus a 5% of the range on each side.\n
        - dict: {chr_name1: (min_coord, max coord), chr_name2: (min_coord, max_coord), ...}. Not
        all the plotted chromosomes need to be specified in the dictionary and some coordinates
        can be indicated as None, both cases lead to the use of the default value.\n
        - tuple: the coordinate limits of all chromosomes will be defined as indicated.\n
        - pyranges.pyranges_main.PyRanges: for each matching chromosome between the plotted data
        and the limits data, the limits will be defined by the minimum and maximum coordinates
        in the pyranges object defined as limits. If some plotted chromosomes are not present they
        will be left as default.

    thick_cds: bool, default False
        Display differentially transcript regions belonging and not belonging to CDS. The CDS/exon information
        must be stored in the 'Feature' column of the PyRanges object or the dataframe. Note that any other
        Feature value other than exon and CDS will be discarded for plotting.

    text: {bool, '{string}'}, default True
        Whether an annotation should appear beside the gene in the plot. If True, the id/index will be used. To
        customize the annotation use the '{string}' option to choose another data column. Providing the text as
        a '{data_column_name}' allows slicing in the case of strings by using '{data_column_name[:4]}'.

    legend: bool, default False
        Whether the legend should appear in the plot.

    title_chr: str, default "Chromosome {chrom}"
        String providing the desired title for the chromosome plots. It should be given in a way where
        the chromosome value in the data is indicated as {chrom}.

    y_labels: list, default None
        Name to identify the PyRanges object/s in the plot.

    tooltip: str, default None
        Dataframe information to show in a tooltip when placing the mouse over a gene, the given
        information will be added to the default: strand, start-end coordinates and id. This must be
        provided as a string containing the column names of the values to be shown within curly brackets.
        For example if you want to show the value of the pointed gene for the column "col1" a valid tooltip
        string could be: "Value of col1: {col1}". Note that the values in the curly brackets are not
        strings. If you want to introduce a newline you can use a newline character "\" + "n".

    to_file: {str, tuple}, default None
        Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
        Optionally, a tuple can be privided where the file name is specified as a str in the first position and in the
        second position there is a tuple specifying the height and width of the figure in px.

    theme: str, default "light"
        General color appearance of the plot. Available modes: "light", "dark", "Mariotti_lab", "swimming_pool".

    **kargs
        Customizable plot features can be defined using kargs. Use print_options() function to check the variables'
        nomenclature, description and default values.



    Examples
    --------

    >>> import pyranges as pr, pyranges_plot as prp

    >>> prp.set_engine('plotly')

    >>> p = pr.PyRanges({"Chromosome": [1]*5, "Strand": ["+"]*3 + ["-"]*2, "Start": [10,20,30,25,40], "End": [15,25,35,30,50], "transcript_id": ["t1"]*3 + ["t2"]*2}, "feature1": ["A", "B", "C", "A", "B"])

    >>> plot(p, id_col="transcript_id",  max_shown=25, colormap='Set3', text=False)

    >>> plot(p, id_col="transcript_id", color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot(p, limits = {'1': (1000, 50000), '2': None, '3': (10000, None)}, title_chr="Chrom: {chrom}")

    >>> plot([p, p], id_col="transcript_id", shrink=True, tooltip = "Feature1: {feature1}")

    >>> plot([p, p], id_col="transcript_id", y_labels=["first_p", "second_p"], packed=False, to_file='my_plot.pdf')
    """

    # Treat input data as list
    if not isinstance(data, list):
        data = [data]

    # Ensure correct y_labels
    if y_labels:
        if len(y_labels) != len(data):
            raise Exception(
                f"The number of provided y_labels {y_labels} does not match the number of PyRanges objects ({len(data)})."
            )

    # Deal with export
    if to_file:
        # given str file name
        if isinstance(to_file, str):
            ext = to_file[-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = (1600, 800)
        # given tuple (name, size)
        else:
            ext = to_file[0][-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = to_file[1]
            to_file = to_file[0]
    # not given to_file, store default size
    else:
        file_size = (1600, 800)

    # Deal with id column
    if id_col is None:
        ID_COL = get_id_col()
        if not ID_COL:
            ID_COL = ["__interval_index__"]
    else:
        ID_COL = id_col
    # treat as list
    if isinstance(ID_COL, str):
        ID_COL = [ID_COL]

    for df_item in data:
        for id_str in ID_COL:
            # Ensure correct names
            if (
                id_str is not None
                and id_str not in df_item.columns
                and id_str != "__interval_index__"
            ):
                raise Exception(
                    "Please define a correct name of the ID column using either set_id_col() function or plot_generic parameter as plot_generic(..., id_col = 'your_id_col')"
                )
            # Avoid Nan in id column

    # Deal with transcript structure
    if thick_cds:
        for df_item in data:
            if "Feature" not in df_item.columns:
                raise Exception(
                    "The transcript structure information must be stored in 'Feature' column of the data."
                )

    # Deal with warnings
    if warnings is None:
        warnings = get_warnings()

    # Deal with engine
    engine = get_engine()

    # PREPARE DATA for plot
    # Deal with plot features as kargs
    wrong_keys = [k for k in kargs if k not in print_options(return_keys=True)]
    if wrong_keys:
        raise Exception(
            f"The following keys do not match any customizable features: {wrong_keys}.\nCheck the customizable variable names using the print_options function."
        )

    def getvalue(key):
        if key in kargs:
            value = kargs[key]
            return value  ## add invalid data type??
        else:
            return get_options(key)

    # Get default plot features
    # store old options to reset them after the plot
    oldtheme = get_theme()
    oldfeat_dict = get_options("values")

    # check option modifications in params
    if theme is None:  # not specified in params, check if it was set
        theme = get_theme()
    set_theme(theme)

    feat_dict = {
        "colormap": getvalue("colormap"),
        "intron_color": getvalue("intron_color"),
        "tag_bkg": getvalue("tag_bkg"),
        "fig_bkg": getvalue("fig_bkg"),
        "plot_bkg": getvalue("plot_bkg"),
        "plot_border": getvalue("plot_border"),
        "title_dict_plt": {
            "family": "sans-serif",
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")) - 5,
        },
        "title_dict_ply": {
            "family": getvalue("title_font"),
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")),
        },
        "grid_color": getvalue("grid_color"),
        "exon_border": getvalue("exon_border"),
        "exon_height": float(getvalue("exon_height")),
        "transcript_utr_width": 0.3 * float(getvalue("exon_height")),
        "v_spacer": getvalue("v_spacer"),
        "text_size": float(getvalue("text_size")),
        "text_pad": getvalue("text_pad"),
        "plotly_port": getvalue("plotly_port"),
        "arrow_line_width": float(getvalue("arrow_line_width")),
        "arrow_color": getvalue("arrow_color"),
        "arrow_size": getvalue("arrow_size"),
        "shrink_threshold": getvalue("shrink_threshold"),
        "shrunk_bkg": getvalue("shrunk_bkg"),
        "x_ticks": getvalue("x_ticks"),
    }
    shrink_threshold = feat_dict["shrink_threshold"]
    colormap = feat_dict["colormap"]
    if colormap == "popart":
        colormap = prp_cmap

    # restore options set before plot is called
    set_theme(oldtheme)
    set_options(oldfeat_dict)

    # Make DataFrame subset if needed
    df_d = {}
    tot_ngenes_l = []
    for pr_ix, df_item in enumerate(data):
        # deal with empty PyRanges
        if df_item.empty:
            continue
        df_item = df_item.copy()

        # consider not known id_col, plot each interval individually
        if ID_COL == ["__interval_index__"]:
            df_item["__interval_index__"] = [str(i) for i in range(len(df_item))]
            df_d[pr_ix], tot_ngenes = make_subset(
                df_item, "__interval_index__", max_shown
            )
            tot_ngenes_l.append(tot_ngenes)

        # known id_col
        else:
            df_d[pr_ix], tot_ngenes = make_subset(df_item, ID_COL, max_shown)
            tot_ngenes_l.append(tot_ngenes)

    for tot_ngenes in tot_ngenes_l:
        if tot_ngenes > max_shown:
            subset_warn = 1
            break
        else:
            subset_warn = 0

    # concat subset dataframes and create new column with input list index
    if not df_d:
        raise Exception("The provided PyRanges object/s are empty.")
    subdf = pd.concat(df_d, names=[PR_INDEX_COL]).reset_index(
        level=PR_INDEX_COL
    )  ### change to pr but doesn't work yet!!

    # group id_cols in one column to count genes in chrmd
    # if len(ID_COL) > 1:
    #   subdf["__id_col_2count__"] = list(zip(*[subdf[c] for c in ID_COL+[PR_INDEX_COL]+[CHROM_COL]]))
    # else:
    subdf["__id_col_2count__"] = list(
        zip(*[subdf[c] for c in [CHROM_COL] + [PR_INDEX_COL] + ID_COL])
    )

    # Deal with thickness_col
    # prioritize transcript structure
    if thick_cds:
        # keep only the "exon" and "CDS"
        subdf = subdf[subdf["Feature"].isin(["exon", "CDS"])]
        if subdf.empty:
            raise Exception(
                "The provided data does not contain any interval containing 'exon' or 'CDS' in the Feature column, no data wil be plotted using 'thick_cds'."
            )
        # set proper thickness column
        thickness_col = "Feature"

    # set proper height values
    if thickness_col:
        # Is it present in data?
        if thickness_col not in subdf.columns:
            raise Exception(
                f"The provided thickness_col {thickness_col} is not present in the given data."
            )

        # Does it have more than 2 values
        if len(subdf[thickness_col].drop_duplicates()) > 2:
            raise Exception("Thickness_col must have a max of 2 different values.")

        # add thickness_col
        thick_tags_l = sorted(
            list(subdf[thickness_col].drop_duplicates()), reverse=True
        )
        if len(thick_tags_l) == 1:
            thick_tags_l = 2 * thick_tags_l
        thick_tags_d = {
            thick_tags_l[0]: feat_dict["transcript_utr_width"],
            thick_tags_l[1]: feat_dict["exon_height"],
        }
        subdf[THICK_COL] = subdf[thickness_col].map(thick_tags_d)
    else:
        subdf[THICK_COL] = [feat_dict["exon_height"]] * len(subdf)

    # Store color information in data
    # color_col as list
    if color_col is None:
        color_col = ID_COL
    elif isinstance(color_col, str):
        color_col = [color_col]

    subdf = subdf_assigncolor(subdf, colormap, color_col, feat_dict["exon_border"])

    # Create genes metadata DataFrame
    genesmd_df = get_genes_metadata(
        subdf,
        ID_COL,
        color_col,
        packed,
        feat_dict["exon_height"],
        feat_dict["v_spacer"],
    )

    # Create chromosome metadata DataFrame
    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf,
        limits,
        genesmd_df,
        packed,
        feat_dict["v_spacer"],
        feat_dict["exon_height"],
    )

    # Deal with introns off
    # adapt coordinates to shrunk
    ts_data = {}
    subdf[ORISTART_COL] = subdf[START_COL]
    subdf[ORIEND_COL] = subdf[END_COL]
    tick_pos_d = {}
    ori_tick_pos_d = {}

    if shrink:
        # compute threshold
        if isinstance(shrink_threshold, int):
            subdf[SHRTHRES_COL] = [shrink_threshold] * len(subdf)
        elif isinstance(shrink_threshold, float):
            subdf[SHRTHRES_COL] = [shrink_threshold] * len(subdf)
            subdf = subdf.groupby(CHROM_COL, group_keys=False, observed=True).apply(
                lambda x: compute_thresh(x, chrmd_df_grouped) if not x.empty else None
            )

        subdf = subdf.groupby(CHROM_COL, group_keys=False, observed=True).apply(
            lambda x: introns_resize(x, ts_data, ID_COL)  # if not x.empty else None
        )  # empty rows when subset
        subdf[START_COL] = subdf[ADJSTART_COL]
        subdf[END_COL] = subdf[ADJEND_COL]

        # recompute limits
        chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
            subdf,
            limits,
            genesmd_df,
            packed,
            feat_dict["v_spacer"],
            feat_dict["exon_height"],
            ts_data=ts_data,
        )

        # compute new axis values and positions if needed
        if ts_data:
            tick_pos_d, ori_tick_pos_d = recalc_axis(
                ts_data, tick_pos_d, ori_tick_pos_d
            )

    else:
        subdf[CUM_DELTA_COL] = [0] * len(subdf)

    # Sort data to plot chromosomes and pr objects in order
    subdf.sort_values([CHROM_COL, PR_INDEX_COL] + ID_COL + [START_COL], inplace=True)
    chrmd_df.sort_values([CHROM_COL, PR_INDEX_COL], inplace=True)
    subdf[EXON_IX_COL] = subdf.groupby(
        [CHROM_COL, PR_INDEX_COL] + ID_COL, group_keys=False, observed=True
    ).cumcount()
    genesmd_df.sort_values([CHROM_COL, PR_INDEX_COL] + [START_COL], inplace=True)

    # Deal with text_pad
    text_pad = feat_dict["text_pad"]
    if isinstance(text_pad, int):
        subdf[TEXT_PAD_COL] = [text_pad] * len(subdf)
    elif isinstance(text_pad, float):
        subdf[TEXT_PAD_COL] = [text_pad] * len(subdf)
        subdf = subdf.groupby(CHROM_COL, group_keys=False, observed=True).apply(
            lambda x: compute_tpad(x, chrmd_df_grouped) if not x.empty else None
        )

    # Deal with added plots
    if (len(chrmd_df_grouped) > 1) and add_aligned_plots:
        raise Exception(
            f"The parameter add_aligned_plots accepts only one chromosome in the input data. The provided data contains {len(chrmd_df_grouped)}"
        )

    if "REF" in subdf.columns:
        subdf["REF"] = subdf["REF"].astype(str)
        subdf["REF"] = subdf["REF"].replace(["nan", "NaN", "None"], np.nan)

    if tooltip is None:
        # Create a list to store the updated tooltips
        updated_tooltips = []
        subdf["__tooltip__"] = ""
        for index, row in subdf.iterrows():
            if STRAND_COL in subdf.columns:
                strand = row.get(STRAND_COL)
            else:
                strand = ""
            if "REF" in subdf.columns:
                if pd.notna(row.get("REF", None)):
                    tool_str = row["REF"] + ">" + row["ALT"]
                    geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}<br>{tool_str}"
                else:
                    if strand:
                        geneinfo = f"[{strand}] ({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default with strand
                    else:
                        geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default without strand
            else:
                if strand:
                    geneinfo = f"[{strand}] ({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default with strand
                else:
                    geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default without strand

            updated_tooltips.append(geneinfo)
        # Assign the updated tooltips back to the DataFrame
        subdf["__tooltip__"] = updated_tooltips

    if tooltip is None:
        tooltip = "{__tooltip__}"

    # deal with engine and call proper plot
    if engine in ["plt", "matplotlib"]:
        if not missing_plt_flag:
            plot_exons_plt(
                subdf=subdf,
                depth_col=depth_col,
                tot_ngenes_l=tot_ngenes_l,
                feat_dict=feat_dict,
                genesmd_df=genesmd_df,
                chrmd_df=chrmd_df,
                chrmd_df_grouped=chrmd_df_grouped,
                ts_data=ts_data,
                max_shown=max_shown,
                id_col=ID_COL,
                transcript_str=thick_cds,
                tooltip=tooltip,
                legend=legend,
                y_labels=y_labels,
                text=text,
                title_chr=title_chr,
                packed=packed,
                to_file=to_file,
                file_size=file_size,
                warnings=warnings,
                tick_pos_d=tick_pos_d,
                ori_tick_pos_d=ori_tick_pos_d,
            )

        else:
            raise Exception(
                "Make sure to install matplotlib dependecies by running `pip install pyranges-plot[plt]`"
            )

    elif engine in ["ply", "plotly"]:
        if not missing_ply_flag:
            if return_plot is not None:
                return plot_exons_ply(
                    subdf=subdf,
                    depth_col=depth_col,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=thick_cds,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    y_labels=y_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                    subset_warn=subset_warn,
                )
            else:
                plot_exons_ply(
                    subdf=subdf,
                    depth_col=depth_col,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=thick_cds,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    y_labels=y_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                    subset_warn=subset_warn,
                )
        else:
            raise Exception(
                "Make sure to install plotly dependecies by running `pip install pyranges-plot[plotly]`"
            )

    else:
        raise Exception("Please define engine with set_engine().")
