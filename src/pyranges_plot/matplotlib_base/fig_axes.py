import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from pyranges.core.names import CHROM_COL, START_COL, END_COL
from pyranges_plot.core import cumdelting
from .core import make_annotation
from ..names import PR_INDEX_COL, CUM_DELTA_COL


def ax_display(ax, title, chrom, t_dict, plot_back, plot_border):
    """Set plot features."""

    if title:
        ax.set_title(title.format(**locals()), fontdict=t_dict)

    ax.set_facecolor(plot_back)
    plt.setp(ax.spines.values(), color=plot_border)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=plot_border)
    ax.xaxis.set_tick_params(bottom=False)
    ax.yaxis.set_tick_params(left=False)


def ax_limits(ax, x_min, x_max, x_rang, grid_color):
    """Adapt plots coordinates."""

    ax.set_xlim(
        x_min - 0.05 * x_rang, x_max + 0.05 * x_rang
    )  # add 5% to limit coordinates range
    plt.ticklabel_format(style="plain")
    ax.grid(visible=True, axis="x", linestyle=":", color=grid_color)  # , zorder = -1)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)  # not scientific notation
    ax.xaxis.get_major_formatter().set_useOffset(False)  # not offset notation
    ax.xaxis.set_major_locator(
        MaxNLocator(integer=True)
    )  # only integer ticks for bases


def ax_shrink_rects(
    ax, fig, ts_data, chrom, y_min, y_max, shrinked_bkg, tag_background
):
    """Add shrinked regions rectangles to the plot."""

    rects_df = ts_data[chrom].copy()
    rects_df["cumdelta_end"] = rects_df[CUM_DELTA_COL]
    rects_df["cumdelta_start"] = rects_df[CUM_DELTA_COL].shift(periods=1, fill_value=0)
    rects_df[START_COL] -= rects_df["cumdelta_start"]
    rects_df[END_COL] -= rects_df["cumdelta_end"]

    for a, b, c, d in zip(
        rects_df[START_COL],
        rects_df[END_COL],
        rects_df["cumdelta_start"],
        rects_df["cumdelta_end"],
    ):
        ts_range = Rectangle(
            (a, y_min - 1),
            b - a,
            y_max + 3,
            edgecolor="grey",
            facecolor=shrinked_bkg,
            fill=True,
            linewidth=0,
        )
        ax.add_patch(ts_range)
        make_annotation(
            ts_range,
            fig,
            ax,
            f"Shrinked region:\n[{a + c} - {b + d}]",
            tag_background,
        )


def create_fig(
    x,
    y,
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    legend_item_d,
    title_chr,
    title_dict_plt,
    plot_background,
    plot_border,
    grid_color,
    packed,
    legend,
    y_labels,
    tick_pos_d,
    ori_tick_pos_d,
    tag_background,
    fig_bkg,
    shrinked_bkg,
    v_spacer,
    exon_height,
):
    """Generate the figure and axes fitting the data."""

    # Unify titles and start figure
    titles = [title_chr.format(**{"chrom": chrom}) for chrom in chrmd_df_grouped.index]
    fig = plt.figure(figsize=(x, y), facecolor=fig_bkg)

    gs = gridspec.GridSpec(
        len(titles),
        1,
        height_ratios=chrmd_df_grouped["y_height"].to_list(),
    )  # size of chromosome subplot according to number of gene rows

    # one plot per chromosome
    axes = []
    for i in range(len(titles)):
        chrom = chrmd_df_grouped.index[i]
        axes.append(plt.subplot(gs[i]))
        ax = axes[i]
        # Adjust plot display
        ax_display(ax, title_chr, chrom, title_dict_plt, plot_background, plot_border)

        # set x axis limits
        x_min, x_max = chrmd_df_grouped.loc[chrom]["min_max"]
        x_rang = x_max - x_min
        ax_limits(ax, x_min, x_max, x_rang, grid_color)

        # consider introns off
        if tick_pos_d:
            # get previous default ticks
            original_ticks = [
                int(tick.get_text().replace("−", "-")) for tick in ax.get_xticklabels()
            ][1:]

            # find previous ticks that should be conserved
            to_add_val = []
            # there is data to shrink
            if ori_tick_pos_d[chrom]:
                to_add_val += [
                    tick
                    for tick in original_ticks
                    if tick < min(ori_tick_pos_d[chrom])
                    or tick > max(ori_tick_pos_d[chrom])
                ]
                for ii in range(1, len(ori_tick_pos_d[chrom]) - 1, 2):
                    not_shr0 = ori_tick_pos_d[chrom][ii]
                    not_shr1 = ori_tick_pos_d[chrom][ii + 1]
                    to_add_val += [
                        i for i in original_ticks if not_shr0 < i <= not_shr1
                    ]

            # nothing to shrink
            else:
                to_add_val += original_ticks

            # compute new coordinates of conserved previous ticks
            to_add = to_add_val.copy()
            to_add = cumdelting(to_add, ts_data, chrom)

            # set new ticks
            x_ticks_val = sorted(to_add)
            # do not add ticks beyond adjusted limits
            x_ticks_val = [
                num
                for num in x_ticks_val
                if num <= chrmd_df_grouped.loc[chrom]["min_max"][1]
            ]
            x_ticks_name = sorted(to_add_val)[: len(x_ticks_val)]

            # adjust names
            ax.set_xticks(x_ticks_val)
            ax.set_xticklabels(x_ticks_name)

        # set y axis limits
        y_min = 0.5 - exon_height / 2
        y_max = chrmd_df_grouped.loc[chrom]["y_height"]
        ax.set_ylim(y_min - v_spacer, y_max + v_spacer)
        # gene name as y labels if not packed and not y_labels
        y_ticks_val = []
        y_ticks_name = []
        if not packed and not y_labels:
            y_ticks_val = (
                genesmd_df[genesmd_df["Chromosome"] == chrom]["ycoord"] + 0.5
            ).to_list()
            y_ticks_name_d = (
                genesmd_df[genesmd_df[CHROM_COL] == chrom]
                .groupby(PR_INDEX_COL, group_keys=False, observed=True)
                .groups
            )
            y_ticks_name_d = dict(sorted(y_ticks_name_d.items(), reverse=True))
            y_ticks_name = [list(id) for id in y_ticks_name_d.values()]
            y_ticks_name = [item for sublist in y_ticks_name for item in sublist]

        # Add shrink rectangles
        if ts_data:
            ax_shrink_rects(
                ax,
                fig,
                ts_data,
                chrom,
                y_min,
                y_max,
                shrinked_bkg,
                tag_background,
            )

        ax.tick_params(colors=plot_border, which="both")

        # Draw lines separating pr objects
        if chrmd_df["pr_line"].drop_duplicates().max() != 0:
            pr_line_y_l = chrmd_df.loc[chrom]["pr_line"].tolist()
            if isinstance(pr_line_y_l, int):
                pr_line_y_l = [pr_line_y_l]
            pr_line_y_l = [y_max + v_spacer] + pr_line_y_l
            present_pr_l = chrmd_df_grouped.loc[chrom]["present_pr"]

            # separate items with horizontal lines
            for j, pr_line_y in enumerate(pr_line_y_l):
                if pr_line_y != 0:
                    ax.plot(
                        [x_min - 0.1 * x_rang, x_max + 0.1 * x_rang],
                        [pr_line_y, pr_line_y],
                        color=plot_border,
                        linewidth=1,
                        zorder=1,
                    )

                    # add y_label in the middle of the subplot if needed
                    if y_labels:
                        if pr_line_y_l[j + 1] != 0:
                            y_ticks_val.append(((pr_line_y) + (pr_line_y_l[j + 1])) / 2)
                        else:
                            y_ticks_val.append((pr_line_y) / 2)
                        y_ticks_name.append(y_labels[int(present_pr_l[j])])

        ax.set_yticks(y_ticks_val)
        ax.set_yticklabels(list(y_ticks_name))

    plt.subplots_adjust(hspace=0.7)

    # Create legend
    if legend:
        labels, handles = zip(*legend_item_d.items())
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1))

    return fig, axes
