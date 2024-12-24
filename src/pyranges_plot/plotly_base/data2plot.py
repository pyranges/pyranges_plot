import pyranges as pr
from pyranges.core.names import START_COL, END_COL

from .core import coord2percent, percent2coord
import plotly.graph_objects as go
import pandas as pd

from ..names import (
    ADJSTART_COL,
    ADJEND_COL,
    EXON_IX_COL,
    TEXT_PAD_COL,
    COLOR_INFO,
    COLOR_TAG_COL,
    BORDER_COLOR_COL,
    THICK_COL,
)


def plot_direction(
    fig,
    strand,
    genename,
    item_size,
    item_threshold,
    start,
    stop,
    incl,
    gene_ix,
    chrom_ix,
    exon_height,
    arrow_color,
    arrow_line_width,
):
    """Plot the direction arrow in the given item if it proceeds."""

    dir_flag = 0

    if strand:
        # create and plot direction lines
        if item_size <= item_threshold:
            dir_flag = 1
            ##diagonal_line = OX arrow extension(item middle point +- incl), OY arrow extension (item middle point + half of exon width)
            top_plus = (
                [(start + stop) / 2 + incl, (start + stop) / 2 - incl],
                [gene_ix, gene_ix + exon_height / 2 - 0.01],
            )
            bot_plus = (
                [(start + stop) / 2 - incl, (start + stop) / 2 + incl],
                [gene_ix - exon_height / 2 + 0.01, gene_ix],
            )
            top_minus = (
                [(start + stop) / 2 + incl, (start + stop) / 2 - incl],
                [gene_ix - exon_height / 2 + 0.01, gene_ix],
            )
            bot_minus = (
                [(start + stop) / 2 - incl, (start + stop) / 2 + incl],
                [gene_ix, gene_ix + exon_height / 2 - 0.01],
            )

            if strand == "+":
                arrow_bot = go.Scatter(
                    x=bot_plus[0],
                    y=bot_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=str(genename),
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_plus[0],
                    y=top_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=str(genename),
                    hoverinfo="skip",
                )
                fig.add_trace(arrow_bot, row=chrom_ix + 1, col=1)
                fig.add_trace(arrow_top, row=chrom_ix + 1, col=1)

            elif strand == "-":
                arrow_bot = go.Scatter(
                    x=bot_minus[0],
                    y=bot_minus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=str(genename),
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_minus[0],
                    y=top_minus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=str(genename),
                    hoverinfo="skip",
                )
                fig.add_trace(arrow_bot, row=chrom_ix + 1, col=1)
                fig.add_trace(arrow_top, row=chrom_ix + 1, col=1)
    return dir_flag


def apply_gene_bridge(
    transcript_str,
    text,
    text_size,
    df,
    fig,
    strand,
    genename,
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
):
    """Evaluate data and provide plot_row with right parameters."""

    # If transcript structure subtract exons
    if transcript_str:
        cds = df[df["Feature"] == "CDS"]
        exons = df[df["Feature"] == "exon"]

        # if there are exons and cds, subtract
        if sum([cds.empty, exons.empty]) == 2:
            exons = exons.subtract_ranges(cds)
        df = pr.concat([cds, exons])

    # Define depth order
    if depth_col:
        df.sort_values(depth_col, inplace=True)

    df.apply(
        plot_row,
        args=(
            fig,
            strand,
            genename,
            gene_ix,
            chrom_ix,
            showinfo,
            legend,
            return_plot,
            arrow_size,
            arrow_color,
            arrow_line_width,
            dir_flag,
            text,
            text_size,
        ),
        axis=1,
    )


def plot_row(
    row,
    fig,
    strand,
    genename,
    gene_ix,
    chrom_ix,
    showinfo,
    legend,
    return_plot,
    arrow_size,
    arrow_color,
    arrow_line_width,
    dir_flag,
    text,
    text_size,
):
    """Plot elements corresponding to one row of one gene."""

    # Get the gene information to print on hover
    # default
    """
    if row.get("vcf", True):
        tooltip_col = row.get("Tooltip_col")
        geneinfo = f"({row.__oriStart__}, {row.__oriEnd__})<br>ID: {genename}<br>{tooltip_col}"
    else:
        if strand:
            geneinfo = f"[{strand}] ({row.__oriStart__}, {row.__oriEnd__})<br>ID: {genename}"  # default with strand
        else:
            geneinfo = f"({row.__oriStart__}, {row.__oriEnd__})<br>ID: {genename}"  # default without strand
    """
    geneinfo = ""
    # customized
    showinfo_dict = row.to_dict()  # first element of gene rows
    if showinfo:
        if showinfo.startswith("$") and showinfo[1:] in row.index:
            showinfo = showinfo_dict[showinfo[1:]]
        showinfo = showinfo.replace("\n", "<br>")
        geneinfo += "<br>" + showinfo.format(**showinfo_dict)

    # consider legend
    if legend:
        legend = bool(row["legend_tag"])

    # Exon start, stop and color
    start = int(row[START_COL])
    stop = int(row[END_COL])
    exon_color = row[COLOR_INFO]
    exon_border = row[BORDER_COLOR_COL]
    exon_height = row[THICK_COL]

    # convert to coordinates for rectangle
    x0, x1 = start, stop
    y0, y1 = (
        gene_ix - exon_height / 2,
        gene_ix + exon_height / 2,
    )  ##gene middle point -+ half of exon size

    # Plot EXON as rectangle
    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=exon_color,
            mode="lines",
            line=dict(color=exon_border),
            text=geneinfo,
            hoverinfo="text",
            name=str(row[COLOR_TAG_COL]),
            showlegend=legend,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Add ID annotation if it is the first exon
    if row[EXON_IX_COL] == 0 and text:
        text_pad = row[TEXT_PAD_COL]
        # text == True
        if isinstance(text, bool):
            ann = str(genename)
        # text == '{string}'
        else:
            row_dict = row.to_dict()
            ann = text.format(**row_dict)

        fig.add_annotation(
            dict(
                x=x0 - text_pad,
                y=(y0 + y1) / 2,
                showarrow=False,
                text=ann,
                textangle=0,
                xanchor="right",
            ),
            row=chrom_ix + 1,
            col=1,
            font={"size": text_size},
        )

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    # arrow_size = coord2percent(fig, chrom_ix + 1, 0.05 * start, 0.05 * stop)
    incl = percent2coord(fig, chrom_ix + 1, arrow_size / 2)  # how long in the plot (OX)

    # create and plot lines
    if not dir_flag:
        plot_direction(
            fig,
            strand,
            genename,
            incl * 2,
            stop - start,  # itself as threshold
            start,
            stop,
            incl,
            gene_ix,
            chrom_ix,
            exon_height,
            arrow_color,
            arrow_line_width,
        )


def plot_introns(
    introns,
    ts_chrom,
    fig,
    gene_ix,
    color,
    chrom_ix,
    strand,
    genename,
    exon_height,
    arrow_color,
    arrow_line_width,
    arrow_size,
):
    """Plot intron lines as needed."""

    dir_flag = []

    def apply_plot_intron(row):
        """Plot intron df as lines."""

        start = row[START_COL]
        stop = row[END_COL]

        # NOT introns off
        if ts_chrom.empty:
            ts_intron = pd.DataFrame()

        # INTRONS OFF
        else:
            ts_intron = ts_chrom[
                (ts_chrom[ADJSTART_COL] >= start) & (ts_chrom[ADJSTART_COL] < stop)
            ].reset_index()

        # Plot LINES binding exons
        # No to-shrink regions in intron
        if ts_intron.empty:
            # create continuous line
            x0, x1 = start, stop
            y0, y1 = gene_ix, gene_ix
            intron_line = go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(color=color, width=0.7, dash="solid"),
                hoverinfo="skip",
                showlegend=False,
            )
            fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

        # Intron has to-shrink regions
        else:
            # alternating fixed and to-shrink regions or vice versa
            prev_tsend = None

            # Iterate over to-shrink regions
            for ix, row in ts_intron.iterrows():
                # (1) Add previous fixed region if needed
                # consider intron starts with fixed region
                if not prev_tsend and row[ADJSTART_COL] != start:
                    prev_tsend = start

                # create continuous line
                x0, x1 = prev_tsend, row[ADJSTART_COL]
                y0, y1 = gene_ix, gene_ix
                intron_line = go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color=color, width=0.7, dash="solid"),
                    hoverinfo="skip",
                    showlegend=False,
                )
                fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # (2) Add to-shrink region
                x0, x1 = row[ADJSTART_COL], row[ADJEND_COL]
                y0, y1 = gene_ix, gene_ix
                intron_line = go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color=color, width=0.7, dash="dot"),
                    hoverinfo="skip",
                    showlegend=False,
                )
                fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # (3) Add final fixed region if needed
                if (ix == len(ts_intron) - 1) and (row[ADJEND_COL] != stop):
                    # add last fixed region
                    # create continuous line
                    x0, x1 = row[ADJEND_COL], stop
                    y0, y1 = gene_ix, gene_ix
                    intron_line = go.Scatter(
                        x=[x0, x1],
                        y=[y0, y1],
                        mode="lines",
                        line=dict(color=color, width=0.7, dash="solid"),
                        hoverinfo="skip",
                        showlegend=False,
                    )
                    fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # store interval end for next iteration
                prev_tsend = row[ADJEND_COL]

        # Plot DIRECTION ARROW in INTRONS if strand is known
        intron_size = coord2percent(fig, chrom_ix + 1, start, stop)
        incl = percent2coord(
            fig, chrom_ix + 1, arrow_size / 2
        )  # how long in the plot (OX)

        dir_flag.append(
            plot_direction(
                fig,
                strand,
                genename,
                arrow_size,  # size of arrow as threshold
                intron_size,
                start,
                stop,
                incl,
                gene_ix,
                chrom_ix,
                exon_height,
                arrow_color,
                arrow_line_width,
            )
        )

    introns.apply(apply_plot_intron, axis=1)

    if 1 in dir_flag:
        return 1
    else:
        return 0
