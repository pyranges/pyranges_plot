import pyranges as pr
from pyranges.core.names import START_COL, END_COL

from .core import coord2percent, percent2coord, make_annotation, rgb_string_to_tuple
from matplotlib.patches import Rectangle
import pandas as pd

from ..names import (
    ADJSTART_COL,
    ADJEND_COL,
    EXON_IX_COL,
    TEXT_PAD_COL,
    COLOR_INFO,
    BORDER_COLOR_COL,
    THICK_COL,
)


def plot_direction(
    ax,
    strand,
    item_size,
    item_threshold,
    start,
    stop,
    incl,
    gene_ix,
    exon_height,
    arrow_color,
    arrow_style,
    arrow_width,
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
                ax.plot(
                    bot_plus[0],
                    bot_plus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

                ax.plot(
                    top_plus[0],
                    top_plus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

            elif strand == "-":
                ax.plot(
                    bot_minus[0],
                    bot_minus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

                ax.plot(
                    top_minus[0],
                    top_minus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

    return dir_flag


def apply_gene_bridge(
    transcript_str,
    text,
    text_size,
    df,
    fig,
    ax,
    strand,
    gene_ix,
    tag_background,
    plot_border,
    genename,
    showinfo,
    arrow_size,
    arrow_color,
    arrow_style,
    arrow_width,
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

    # plot rows
    df.apply(
        plot_row,
        args=(
            fig,
            ax,
            strand,
            gene_ix,
            tag_background,
            plot_border,
            genename,
            showinfo,
            arrow_size,
            arrow_color,
            arrow_style,
            arrow_width,
            dir_flag,
            text,
            text_size,
        ),
        axis=1,
    )


def plot_row(
    row,
    fig,
    ax,
    strand,
    gene_ix,
    tag_background,
    plot_border,
    genename,
    showinfo,
    arrow_size,
    arrow_color,
    arrow_style,
    arrow_width,
    dir_flag,
    text,
    text_size,
):
    """Plot elements corresponding to one row of one gene."""

    # Make gene annotation
    # get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({row.__oriStart__}, {row.__oriEnd__})\nID: {genename}"  # default with strand
    else:
        geneinfo = f"({row.__oriStart__}, {row.__oriEnd__})\nID: {genename}"  # default without strand

    # customized
    showinfo_dict = row.to_dict()  # first element of gene rows
    if showinfo:
        geneinfo += "\n" + showinfo.format(**showinfo_dict)

    # Exon start, stop and color
    start = int(row[START_COL])
    stop = int(row[END_COL])
    exon_color = row[COLOR_INFO]
    if isinstance(exon_color, str) and exon_color[:3] == "rgb":
        exon_color = rgb_string_to_tuple(exon_color)
    exon_border = row[BORDER_COLOR_COL]
    if isinstance(exon_border, str) and exon_border[:3] == "rgb":
        exon_border = rgb_string_to_tuple(exon_border)
    exon_height = row[THICK_COL]
    if isinstance(arrow_color, str) and arrow_color[:3] == "rgb":
        arrow_color = rgb_string_to_tuple(arrow_color)

    # Plot EXON as rectangle
    exon_rect = Rectangle(
        (start, gene_ix - exon_height / 2),
        stop - start,
        exon_height,
        edgecolor=exon_border,
        facecolor=exon_color,
        fill=True,
    )
    ax.add_patch(exon_rect)

    # create annotation for exon
    make_annotation(exon_rect, fig, ax, geneinfo, tag_background)

    # Add ID annotation if it is the first exon
    if row[EXON_IX_COL] == 0 and text:
        text_pad = row[TEXT_PAD_COL]
        # text == True
        if isinstance(text, bool):
            ann = genename
        # text == '{string}'
        else:
            row_dict = row.to_dict()
            ann = text.format_map(row_dict)

        ax.annotate(
            ann,
            xy=(start - text_pad, gene_ix),
            horizontalalignment="right",
            verticalalignment="center",
            color=plot_border,
            fontsize=text_size,
        )

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    incl = percent2coord(ax, arrow_size / 2)  # how long is the arrow in the plot (OX)

    # create and plot lines
    if not dir_flag:
        plot_direction(
            ax,
            strand,
            incl * 2,
            stop - start,  # itself as threshold
            start,
            stop,
            incl,
            gene_ix,
            exon_height,
            arrow_color,
            arrow_style,
            arrow_width,
        )


def plot_introns(
    introns,
    ts_chrom,
    fig,
    ax,
    geneinfo,
    tag_background,
    gene_ix,
    exon_color,
    strand,
    exon_height,
    arrow_color,
    arrow_style,
    arrow_width,
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
            intron_line = ax.plot(
                [start, stop],
                [gene_ix, gene_ix],
                color=exon_color,
                linewidth=1,
                zorder=1,
            )
            # add to plot with annotation
            make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

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
                intron_line = ax.plot(
                    [prev_tsend, row[ADJSTART_COL]],
                    [gene_ix, gene_ix],
                    color=exon_color,
                    linewidth=1,
                    zorder=1,
                )
                # add to plot with annotation
                make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # (2) Add to-shrink region
                intron_line = ax.plot(
                    [row[ADJSTART_COL], row[ADJEND_COL]],
                    [gene_ix, gene_ix],
                    color=exon_color,
                    linewidth=0.5,
                    linestyle="--",
                    zorder=1,
                )
                # add to plot with annotation
                make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # (3) Add final fixed region if needed
                if (ix == len(ts_intron) - 1) and (row[ADJEND_COL] != stop):
                    # add last fixed region
                    # create continuous line
                    intron_line = ax.plot(
                        [row[ADJEND_COL], stop],
                        [gene_ix, gene_ix],
                        color=exon_color,
                        linewidth=1,
                        zorder=1,
                    )
                    # add to plot with annotation
                    make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # store interval end for next iteration
                prev_tsend = row[ADJEND_COL]

        intron_size = coord2percent(ax, start, stop)
        incl = percent2coord(
            ax, arrow_size / 2
        )  # how long is the arrow in the plot (OX)

        # Plot DIRECTION ARROW in INTRONS if strand is known
        dir_flag.append(
            plot_direction(
                ax,
                strand,
                arrow_size,
                intron_size,
                start,
                stop,
                incl,
                gene_ix,
                exon_height,
                arrow_color,
                arrow_style,
                arrow_width,
            )
        )

    introns.apply(apply_plot_intron, axis=1)

    if 1 in dir_flag:
        return 1
    else:
        return 0
