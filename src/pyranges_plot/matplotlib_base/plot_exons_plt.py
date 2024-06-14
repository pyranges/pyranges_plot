import matplotlib.pyplot as plt
import pandas as pd
from pyranges.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL

from .core import plt_popup_warning, coord2percent
from .fig_axes import create_fig
from .data2plot import (
    apply_gene_bridge,
    plot_introns,
)
from ..names import PR_INDEX_COL, BORDER_COLOR_COL

arrow_style = "round"


def plot_exons_plt(
    subdf,
    tot_ngenes_l,
    feat_dict,
    genesmd_df,
    chrmd_df,
    chrmd_df_grouped,
    ts_data,
    legend_item_d,
    id_col,
    max_shown=25,
    transcript_str=False,
    tooltip=None,
    legend=False,
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
    tag_bkg = feat_dict["tag_bkg"]
    fig_bkg = feat_dict["fig_bkg"]
    plot_bkg = feat_dict["plot_bkg"]
    plot_border = feat_dict["plot_border"]
    title_dict_plt = feat_dict["title_dict_plt"]
    grid_color = feat_dict["grid_color"]
    exon_border = feat_dict["exon_border"]
    exon_height = feat_dict["exon_height"]
    transcript_utr_width = feat_dict["transcript_utr_width"]
    v_spacer = feat_dict["v_spacer"]
    text_size = feat_dict["text_size"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size = feat_dict["arrow_size"]
    shrinked_bkg = feat_dict["shrinked_bkg"]

    # Create figure and axes
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
        tick_pos_d,
        ori_tick_pos_d,
        tag_bkg,
        fig_bkg,
        shrinked_bkg,
        v_spacer,
        exon_height,
    )

    # Plot genes
    subdf.groupby(id_col + [PR_INDEX_COL], group_keys=False, observed=True).apply(
        lambda subdf: gby_plot_exons(
            subdf,
            axes,
            fig,
            chrmd_df,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            tooltip,
            tag_bkg,
            plot_border,
            transcript_str,
            text,
            text_size,
            exon_height,
            exon_border,
            transcript_utr_width,
            arrow_line_width,
            arrow_color,
            arrow_size,
        )
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
        plt.show()
    else:
        plt.savefig(to_file, format=to_file[-3:], dpi=400)


def gby_plot_exons(
    df,
    axes,
    fig,
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    tag_bkg,
    plot_border,
    transcript_str,
    text,
    text_size,
    exon_height,
    exon_border,
    transcript_utr_width,
    arrow_line_width,
    arrow_color,
    arrow_size,
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
        geneinfo = f"[{strand}] ({min(df.__oriStart__)}, {max(df.__oriEnd__)})\nID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.__oriStart__)}, {max(df.__oriEnd__)})\nID: {genename}"  # default without strand

    # Plot INTRON lines
    sorted_exons = df[[START_COL, END_COL]].sort_values(by=START_COL)
    sorted_exons["intron_dir_flag"] = [0] * len(sorted_exons)
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    if isinstance(arrow_size, int):
        arrow_size = coord2percent(ax, 0, arrow_size)

    dir_flag = plot_introns(
        sorted_exons,
        ts_chrom,
        fig,
        ax,
        geneinfo,
        tag_bkg,
        gene_ix,
        exon_border,
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
        exon_border,  # this works as "color" used for utr (not interval)
        exon_border,
        tag_bkg,
        plot_border,
        genename,
        showinfo,
        exon_height,
        transcript_utr_width,
        arrow_size,
        arrow_color,
        arrow_style,
        arrow_line_width,
        dir_flag,
    )
