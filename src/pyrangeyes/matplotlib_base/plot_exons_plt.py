import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from matplotlib.patches import Rectangle
from pyranges1.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL

from .core import plt_popup_warning, coord2percent, rgb_string_to_tuple
from .fig_axes import create_fig
from .data2plot import (
    apply_gene_bridge,
    plot_introns,
)
from ..names import (
    PR_INDEX_COL,
    BORDER_COLOR_COL,
    COLOR_INFO,
    COLOR_LEGEND_TITLE_COL,
    OUTLINE_LEGEND_TITLE_COL,
)
from ..legend import (
    categorical_fill_items,
    categorical_outline_items,
    quantitative_fill_info,
    quantitative_outline_info,
)

arrow_style = "round"


def _legend_title(subdf, title_col, fallback):
    if title_col in subdf.columns and not subdf.empty:
        title = subdf[title_col].iloc[0]
        if title:
            return str(title)
    return fallback


def _categorical_legend_entries(subdf):
    """Return ``[(label, artist), ...]`` for Matplotlib categorical legends."""
    fill_title = _legend_title(subdf, COLOR_LEGEND_TITLE_COL, "color")
    outline_title = _legend_title(subdf, OUTLINE_LEGEND_TITLE_COL, "outline")

    entries = [
        (
            f"{fill_title}: {label}",
            Rectangle((0, 0), 1, 1, facecolor=fill, edgecolor=outline),
        )
        for label, fill, outline in categorical_fill_items(subdf)
    ]
    entries.extend(
        (
            f"{outline_title}: {label}",
            Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=outline, linewidth=2),
        )
        for label, outline in categorical_outline_items(subdf)
    )
    return entries


def _legend_ncol(fig, n_items):
    if n_items <= 0:
        return 1
    # Keep labels in compact rows while avoiding very wide single-line legends.
    return max(1, min(n_items, int(fig.get_figwidth() // 1.5) or 1))


def _reserve_bottom_space(fig, categorical_entries, quantitative_infos):
    ncol = _legend_ncol(fig, len(categorical_entries))
    categorical_rows = 0
    if categorical_entries:
        categorical_rows = (len(categorical_entries) + ncol - 1) // ncol
    bottom = 0.035 + 0.045 * categorical_rows + 0.06 * len(quantitative_infos)
    bottom = min(0.45, bottom)
    fig.subplots_adjust(bottom=bottom)
    return ncol


def _add_bottom_legends(fig, categorical_entries, quantitative_infos):
    """Place all Matplotlib legends in a reserved bottom band.

    Figure-level legends and colorbars do not participate in Matplotlib's normal
    axes layout, so placing them inside the plotting area can overlap data,
    titles, or each other.  We reserve a bottom strip, then stack the categorical
    legend and any quantitative colorbars there in figure coordinates.
    """
    if not categorical_entries and not quantitative_infos:
        return

    ncol = _reserve_bottom_space(fig, categorical_entries, quantitative_infos)
    y = 0.015

    if categorical_entries:
        labels, handles = zip(*categorical_entries)
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, y),
            ncol=ncol,
            frameon=True,
        )
        categorical_rows = (len(categorical_entries) + ncol - 1) // ncol
        y += 0.045 * categorical_rows

    for qinfo in quantitative_infos:
        cmap = mcolors.LinearSegmentedColormap.from_list(
            f"pyrangeyes_{qinfo['title']}", qinfo["colors"]
        )
        norm = mcolors.Normalize(vmin=qinfo["vmin"], vmax=qinfo["vmax"])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.15, y, 0.7, 0.018])
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.set_label(qinfo["title"])
        y += 0.06


def plot_exons_plt(
    subdf,
    depth_col,
    tot_ngenes_l,
    feat_dict,
    genesmd_df,
    chrmd_df,
    chrmd_df_grouped,
    ts_data,
    id_col,
    max_shown=25,
    transcript_str=False,
    tooltip=None,
    legend=False,
    return_plot=None,
    add_aligned_plots=None,
    y_labels=False,
    text=True,
    title_chr=None,
    packed=True,
    to_file=None,
    file_size=None,
    warnings=None,
    tick_pos_d=None,
    ori_tick_pos_d=None,
):
    """Create Matplotlib plot."""

    # Get default plot features
    # deal with rgb color options
    for key, val in feat_dict.items():
        if isinstance(val, str) and val[:3] == "rgb":
            feat_dict[key] = rgb_string_to_tuple(val)

    # store values
    intron_color = feat_dict["intron_color"]
    tag_bkg = feat_dict["tag_bkg"]
    fig_bkg = feat_dict["fig_bkg"]
    plot_bkg = feat_dict["plot_bkg"]
    plot_border = feat_dict["plot_border"]
    title_dict_plt = feat_dict["title_dict_plt"]
    grid_color = feat_dict["grid_color"]
    exon_border = feat_dict["outline_color"]
    exon_height = feat_dict["exon_height"]
    v_spacer = feat_dict["v_spacer"]
    text_size = feat_dict["text_size"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size = feat_dict["arrow_size"]
    shrunk_bkg = feat_dict["shrunk_bkg"]
    x_ticks = feat_dict["x_ticks"]

    # Create figure and axes
    if legend:
        categorical_legend_entries = _categorical_legend_entries(subdf)
        quantitative_legend_infos = [
            qinfo
            for qinfo in [quantitative_fill_info(subdf), quantitative_outline_info(subdf)]
            if qinfo is not None
        ]
    else:
        categorical_legend_entries = []
        quantitative_legend_infos = []
    legend_item_d = {}
    # pixel in inches
    px = 1 / plt.rcParams["figure.dpi"]
    x = file_size[0] * px
    y = file_size[1] * px

    fig, axes = create_fig(
        x,
        y,
        chrmd_df,
        chrmd_df_grouped,
        genesmd_df,
        id_col,
        ts_data,
        legend_item_d,
        title_chr,
        title_dict_plt,
        plot_bkg,
        plot_border,
        grid_color,
        packed,
        legend,
        y_labels,
        x_ticks,
        tick_pos_d,
        ori_tick_pos_d,
        tag_bkg,
        fig_bkg,
        shrunk_bkg,
        v_spacer,
        exon_height,
        add_aligned_plots,
    )

    if legend:
        _add_bottom_legends(fig, categorical_legend_entries, quantitative_legend_infos)

    # Plot genes
    # pd.DataFrame.groupby(subdf,
    for _, grouped_subdf in subdf.groupby(
        id_col + [PR_INDEX_COL, CHROM_COL], group_keys=False, observed=True
    ):
        gby_plot_exons(
            grouped_subdf,
            axes,
            fig,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            tooltip,
            intron_color,
            tag_bkg,
            plot_border,
            transcript_str,
            text,
            text_size,
            exon_height,
            exon_border,
            arrow_line_width,
            arrow_color,
            arrow_size,
            depth_col,
        )

    # Prevent zoom in y axis
    # for ax in axes:
    #     initial_ylim = ax.get_ylim()
    #     # event handler function to fix y-axis range
    #     def on_xlims_change(axes):
    #         axes.set_ylim(initial_ylim)
    #     ax.callbacks.connect('xlim_changed', on_xlims_change)

    # Provide output
    if to_file is None:
        # evaluate warning
        one_warn = 0
        for tot_ngenes in tot_ngenes_l:
            if tot_ngenes > max_shown and warnings and not one_warn:
                one_warn = 1
                plt_popup_warning(
                    "The provided data contains more genes than the ones plotted."
                )
        if return_plot is None:
            plt.show()
        elif return_plot == "fig":
            return fig
        else:
            raise Exception("Matplotlib can not return a Dash app.")
    else:
        plt.savefig(to_file, format=to_file[-3:], dpi=100)


def gby_plot_exons(
    df,
    axes,
    fig,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    intron_color,
    tag_bkg,
    plot_border,
    transcript_str,
    text,
    text_size,
    exon_height,
    exon_border,
    arrow_line_width,
    arrow_color,
    arrow_size,
    depth_col,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df[CHROM_COL].iloc[0]
    chr_ix = chrmd_df_grouped.loc[chrom]["chrom_ix"]
    if isinstance(chr_ix, pd.Series):
        chr_ix = chr_ix.iloc[0]
    pr_ix = df[PR_INDEX_COL].iloc[0]
    genename = df["__id_col_2count__"].iloc[0]
    genemd = genesmd_df.loc[genename]  # store data for the gene
    geneid = list(df[id_col].iloc[0])
    if len(geneid) == 1:
        geneid = geneid[0]

    # in case same gene in +1 pr
    if not isinstance(genemd, pd.Series):
        genemd = genemd[genemd[PR_INDEX_COL] == pr_ix]  # in case same gene in +1 pr
        genemd = pd.Series(genemd.iloc[0])
    gene_ix = genemd["ycoord"] + 0.5
    # color of border of first interval will be used as intron color and utr color for simplicity
    if exon_border is None:
        exon_border = df[BORDER_COLOR_COL].iloc[0]

    ax = axes[chr_ix]
    if STRAND_COL in df.columns:
        strand = df[STRAND_COL].unique()[0]
    else:
        strand = ""

    # Make gene annotation
    # get the gene information to print on hover
    if strand:
        geneinfo = f"[{strand}] ({min(df.__oriStart__)}, {max(df.__oriEnd__)})\nID: {geneid}"  # default with strand
    else:
        geneinfo = f"({min(df.__oriStart__)}, {max(df.__oriEnd__)})\nID: {geneid}"  # default without strand

    # Plot INTRON lines
    # get introns
    df[START_COL] = df[START_COL].astype(int)
    df[END_COL] = df[END_COL].astype(int)
    introns = df.complement_ranges(id_col)
    introns["intron_dir_flag"] = [0] * len(introns)

    # consider shrunk
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    # deal with arrow and intron color
    if isinstance(arrow_size, int):
        arrow_size = coord2percent(ax, 0, arrow_size)
    if intron_color is None:
        intron_color = df[COLOR_INFO].iloc[0]

    dir_flag = plot_introns(
        introns,
        ts_chrom,
        fig,
        ax,
        geneinfo,
        tag_bkg,
        gene_ix,
        intron_color,
        strand,
        exon_height,
        arrow_color,
        arrow_style,
        arrow_line_width,
        arrow_size,
    )

    # Plot the gene rows as EXONS
    apply_gene_bridge(
        transcript_str,
        text,
        text_size,
        df,
        fig,
        ax,
        strand,
        gene_ix,
        tag_bkg,
        plot_border,
        geneid,
        showinfo,
        arrow_size,
        arrow_color,
        arrow_style,
        arrow_line_width,
        dir_flag,
        depth_col,
    )
