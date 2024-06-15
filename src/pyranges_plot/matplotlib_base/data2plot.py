from pyranges.core.names import START_COL, END_COL

from .core import coord2percent, percent2coord, make_annotation
from matplotlib.patches import Rectangle
import pandas as pd

from ..names import (
    ADJSTART_COL,
    ADJEND_COL,
    EXON_IX_COL,
    TEXT_PAD_COL,
    COLOR_INFO,
    BORDER_COLOR_COL,
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
        if item_size > item_threshold:
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
    exon_color,
    exon_border,
    tag_background,
    plot_border,
    genename,
    showinfo,
    exon_height,
    transcript_utr_width,
    arrow_size,
    arrow_color,
    arrow_style,
    arrow_width,
    dir_flag,
):
    """Evaluate data and provide plot_row with right parameters."""
    # NOT transcript strucutre
    if not transcript_str:
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
                exon_height,
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

    # WITH transcript structure
    else:
        ## add warning for no good transcript str here
        # transcript has only CDS and exon
        if (
            df.Feature.str.contains("CDS").any()
            and df.Feature.str.contains("exon").any()
        ):
            # get coordinates for utr and cds
            tr_start, cds_start = df.groupby(
                "Feature", group_keys=False, observed=True
            ).Start.apply(min)[["exon", "CDS"]]
            tr_end, cds_end = df.groupby(
                "Feature", group_keys=False, observed=True
            ).End.apply(max)[["exon", "CDS"]]

            # create utr
            start_utr = Rectangle(
                (tr_start, gene_ix - transcript_utr_width / 2),
                cds_start - tr_start,
                transcript_utr_width,
                edgecolor=exon_border,
                facecolor=exon_color,
                fill=True,
            )
            end_utr = Rectangle(
                (cds_end, gene_ix - transcript_utr_width / 2),
                tr_end - cds_end,
                transcript_utr_width,
                edgecolor=exon_border,
                facecolor=exon_color,
                fill=True,
            )
            ax.add_patch(start_utr)
            ax.add_patch(end_utr)
            # add ID annotation for utr
            if text:
                text_pad = df[TEXT_PAD_COL].iloc[0]
                # text == True
                if isinstance(text, bool):
                    ann = genename
                # text == '{string}'
                else:
                    row_dict = df.iloc[0].to_dict()  # use first row
                    ann = text.format_map(row_dict)
                ax.annotate(
                    ann,
                    xy=(tr_start - text_pad, gene_ix),
                    horizontalalignment="right",
                    verticalalignment="center",
                    color=plot_border,
                    fontsize=text_size,
                )

            # make annotation for utr
            # get the gene information to print on hover
            # default
            if strand:
                geneinfo_start = f"[{strand}] ({tr_start}, {cds_start})\nID: {genename}"  # default with strand
                geneinfo_end = f"[{strand}] ({cds_end}, {tr_end})\nID: {genename}"  # default with strand
            else:
                geneinfo_start = f"({tr_start}, {cds_start})\nID: {genename}"  # default without strand
                geneinfo_end = (
                    f"({cds_end}, {tr_end})\nID: {genename}"  # default without strand
                )

            make_annotation(start_utr, fig, ax, geneinfo_start, tag_background)
            make_annotation(end_utr, fig, ax, geneinfo_end, tag_background)

            # keep CDS data and plot it
            df = df.groupby("Feature", group_keys=False, observed=True).get_group("CDS")
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
                    exon_height,
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

        # transcript only has exon
        elif (
            not df.Feature.str.contains("CDS").any()
            and df.Feature.str.contains("exon").any()
        ):
            # plot just as utr
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
                    transcript_utr_width,
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

        # transcript only has CDS
        elif (
            df.Feature.str.contains("CDS").any()
            and not df.Feature.str.contains("exon").any()
        ):
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
                    exon_height,
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

        # transcript has neither, skip gene
        else:
            return


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
    exon_height,
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
    exon_border = row[BORDER_COLOR_COL]

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
    # arrow_size = coord2percent(ax, 0.05 * start, 0.05 * stop)
    incl = percent2coord(ax, arrow_size / 2)  # how long is the arrow in the plot (OX)

    # create and plot lines
    if not dir_flag:
        plot_direction(
            ax,
            strand,
            stop - start,
            arrow_size,  # itself as threshold
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
                intron_size,
                arrow_size,
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
