import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
from pyranges.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL

from .core import initialize_dash_app, coord2percent
from .fig_axes import create_fig
from .data2plot import plot_introns, apply_gene_bridge
from ..names import PR_INDEX_COL, BORDER_COLOR_COL, COLOR_INFO


def plot_exons_ply(
    subdf,
    depth_col,
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
    subset_warn=0,
):
    """Create Plotly plot."""

    # Get default plot features
    intron_color = feat_dict["intron_color"]
    fig_bkg = feat_dict["fig_bkg"]
    plot_bkg = feat_dict["plot_bkg"]
    plot_border = feat_dict["plot_border"]
    title_dict_ply = feat_dict["title_dict_ply"]
    grid_color = feat_dict["grid_color"]
    exon_border = feat_dict["exon_border"]
    exon_height = feat_dict["exon_height"]
    v_spacer = feat_dict["v_spacer"]
    text_size = feat_dict["text_size"]
    plotly_port = feat_dict["plotly_port"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size = feat_dict["arrow_size"]
    shrunk_bkg = feat_dict["shrunk_bkg"]
    x_ticks = feat_dict["x_ticks"]

    # Create figure and chromosome plots
    fig = create_fig(
        subdf,
        chrmd_df,
        chrmd_df_grouped,
        genesmd_df,
        id_col,
        ts_data,
        title_chr,
        title_dict_ply,
        grid_color,
        packed,
        y_labels,
        x_ticks,
        tick_pos_d,
        ori_tick_pos_d,
        shrunk_bkg,
        v_spacer,
        exon_height,
        plot_border,
        add_aligned_plots,
    )

    # Plot genes
    subdf.groupby(
        id_col + [PR_INDEX_COL, CHROM_COL], group_keys=False, observed=True
    ).apply(
        lambda subdf: gby_plot_exons(
            subdf,
            fig,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            tooltip,
            intron_color,
            legend,
            return_plot,
            transcript_str,
            text,
            text_size,
            exon_height,
            exon_border,
            plot_bkg,
            arrow_line_width,
            arrow_color,
            arrow_size,
            depth_col,
        )
    )  # .reset_index(level=PR_INDEX_COL)

    # Adjust plot display
    fig.update_layout(
        plot_bgcolor=plot_bkg,
        font_color=plot_border,
        showlegend=legend,
        paper_bgcolor=fig_bkg,
    )
    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor=plot_border,
        mirror=True,
        color=plot_border,
    )
    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor=plot_border,
        mirror=True,
        color=plot_border,
    )

    # Provide output
    # insert silent information for warnings
    if warnings:
        fig.data[0].customdata = np.array([0, 0, 0])  # [tot_ngenes_l, 0, 0])
        if "_blackwarning!" in subdf.columns and "_iterwarning!" in subdf.columns:
            fig.data[0].customdata = np.array(
                [subset_warn, 91124, 91321]
            )  # [subset_warn, 91124, 91321])
        elif "_blackwarning!" in subdf.columns and "_iterwarning!" not in subdf.columns:
            fig.data[0].customdata = np.array(
                [subset_warn, 91124, 0]
            )  # [subset_warn, 91124, 0])
        elif "_blackwarning!" not in subdf.columns and "_iterwarning!" in subdf.columns:
            fig.data[0].customdata = np.array(
                [subset_warn, 0, 91321]
            )  # [subset_warn, 0, 91321])
        elif (
            "_blackwarning!" not in subdf.columns
            and "_iterwarning!" not in subdf.columns
        ):
            fig.data[0].customdata = np.array(
                [subset_warn, 0, 0]
            )  # [subset_warn, 0, 91321])
    else:
        fig.data[0].customdata = np.array(["no warnings"])

    if to_file is None:
        if return_plot is None:
            app_instance = initialize_dash_app(fig, max_shown)
            app_instance.run(port=plotly_port)
        elif return_plot == "app":
            app_instance = initialize_dash_app(fig, max_shown)
            return app_instance
        elif return_plot == "fig":
            return fig
    else:
        fig.update_layout(width=file_size[0], height=file_size[1])
        pio.write_image(fig, to_file)


def gby_plot_exons(
    df,
    fig,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    intron_color,
    legend,
    return_plot,
    transcript_str,
    text,
    text_size,
    exon_height,
    exon_border,
    plot_background,
    arrow_line_width,
    arrow_color,
    arrow_size,
    depth_col,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df[CHROM_COL].iloc[0]
    pr_ix = df[PR_INDEX_COL].iloc[0]
    genename = df["__id_col_2count__"].iloc[0]
    genemd = genesmd_df.loc[genename]  # store data for the gene
    df["legend_tag"] = [1] + [0] * (
        len(df) - 1
    )  # only one legend entry/linked intervals
    geneid = list(df[id_col].iloc[0])
    if len(geneid) == 1:
        geneid = geneid[0]

    # in case same gene in +1 pr
    if not isinstance(genemd, pd.Series):
        genemd = genemd[genemd[PR_INDEX_COL] == pr_ix]  # in case same gene in +1 pr
        genemd = pd.Series(genemd.iloc[0])
    gene_ix = genemd["ycoord"] + 0.5
    # color of first interval will be used as intron color and utr color for simplicity
    if exon_border is None:
        exon_border = df[BORDER_COLOR_COL].iloc[0]

    chrom_ix = chrmd_df_grouped.loc[chrom]["chrom_ix"]

    if STRAND_COL in df.columns:
        strand = df[STRAND_COL].unique()[0]
    else:
        strand = ""

    # Get the gene information to print on hover
    # default
    if "vcf" in df.columns and df["vcf"].any():
        tooltip_col = df["Tooltip_col"]
        geneinfo = f"({min(df.__oriStart__)}, {max(df.__oriEnd__)})<br>ID: {genename}<br>{tooltip_col}"
    else:
        if strand:
            geneinfo = f"[{strand}] ({min(df.__oriStart__)}, {max(df.__oriEnd__)})<br>ID: {geneid}"  # default with strand
        else:
            geneinfo = f"({min(df.__oriStart__)}, {max(df.__oriEnd__)})<br>ID: {geneid}"  # default without strand

    # add annotation for introns to plot
    x0, x1 = min(df[START_COL]), max(df[END_COL])
    y0, y1 = gene_ix - exon_height / 160, gene_ix + exon_height / 160

    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=plot_background,
            mode="lines",
            line=dict(color=exon_border, width=0),
            hoverinfo="text",
            text=geneinfo,
            showlegend=False,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot INTRON lines
    # get introns
    df[START_COL] = df[START_COL].astype(int)
    df[END_COL] = df[END_COL].astype(int)
    introns = df.complement(transcript_id=id_col)
    introns["intron_dir_flag"] = [0] * len(introns)

    # consider shrunk
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    # deal with arrow and intron color
    if isinstance(arrow_size, int):
        arrow_size = coord2percent(fig, chrom_ix + 1, 0, arrow_size)
    if intron_color is None:
        intron_color = df[COLOR_INFO].iloc[0]

    dir_flag = plot_introns(
        introns,
        ts_chrom,
        fig,
        gene_ix,
        intron_color,
        chrom_ix,
        strand,
        geneid,
        exon_height,
        arrow_color,
        arrow_line_width,
        arrow_size,
    )

    # Plot the gene rows
    apply_gene_bridge(
        transcript_str,
        text,
        text_size,
        df,
        fig,
        strand,
        geneid,
        gene_ix,
        chrom_ix,
        showinfo,
        legend,
        return_plot,
        arrow_size,
        arrow_color,
        arrow_line_width,
        dir_flag,
        depth_col,
    )
