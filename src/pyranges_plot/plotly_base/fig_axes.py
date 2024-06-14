import plotly.subplots as sp
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from pyranges.core.names import CHROM_COL, START_COL, END_COL
from pyranges_plot.core import cumdelting
from pyranges_plot.names import PR_INDEX_COL, ORISTART_COL, ORIEND_COL, CUM_DELTA_COL


def calculate_ticks(subdf, num_ticks=10):
    """Calculate tick values for a given data range."""

    # Calculate range and initial tick interval
    data_min = subdf[ORISTART_COL].min()
    data_max = subdf[ORIEND_COL].max()
    data_range = data_max - data_min
    initial_interval = data_range / (num_ticks - 1)

    # Calculate a 'nice' interval between ticks
    # The exponent of the range
    exponent = np.floor(np.log10(initial_interval))
    # Fractional part of the range
    fractional_part = initial_interval / 10**exponent

    # Determine nice fractional part
    if fractional_part < 1.5:
        nice_fractional_part = 1
    elif fractional_part < 3:
        nice_fractional_part = 2
    elif fractional_part < 7:
        nice_fractional_part = 5
    else:
        nice_fractional_part = 10

    nice_interval = nice_fractional_part * 10**exponent

    # Calculate tick values
    tick_values = np.arange(
        start=data_min - (data_min % nice_interval), stop=data_max, step=nice_interval
    )
    # Adjust the last tick value if necessary
    if tick_values[-1] < data_max:
        tick_values = np.append(tick_values, tick_values[-1] + nice_interval)

    return tick_values


def create_fig(
    subdf,
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    title_chr,
    title_dict_ply,
    grid_color,
    packed,
    y_labels,
    tick_pos_d,
    ori_tick_pos_d,
    shrinked_bkg,
    v_spacer,
    exon_height,
    plot_border,
):
    """Generate the figure and axes fitting the data."""

    # Unify titles and start figure
    titles = [title_chr.format(**{"chrom": chrom}) for chrom in chrmd_df_grouped.index]
    titles = list(pd.Series(titles))
    fig = sp.make_subplots(
        rows=len(titles),
        cols=1,
        row_heights=chrmd_df_grouped["y_height"].to_list(),
        subplot_titles=titles,
    )

    # one subplot per chromosome
    for i in range(len(titles)):
        chrom = chrmd_df_grouped.index[i]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i + 1, col=1)

        # set title format if there are titles
        if fig.layout.annotations:
            fig.layout.annotations[i].update(font=title_dict_ply)

        # set x axis limits
        x_min, x_max = chrmd_df_grouped.loc[chrom]["min_max"]
        x_rang = x_max - x_min
        fig.update_xaxes(
            range=[x_min - 0.05 * x_rang, x_max + 0.05 * x_rang],
            tickformat="d",
            showgrid=True,
            gridcolor=grid_color,
            griddash="dot",
            zeroline=False,
            row=i + 1,
            col=1,
        )  # add 5% to limit coordinates range

        # consider introns off
        if tick_pos_d:
            # get previous default ticks
            chrom_subdf = subdf[subdf[CHROM_COL] == chrom]
            original_ticks = list(calculate_ticks(chrom_subdf))

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
                num for num in x_ticks_val if num <= chrmd_df_grouped.loc[chrom]["max"]
            ]
            x_ticks_name = sorted(to_add_val)[: len(x_ticks_val)]

            # set new ticks
            fig.update_xaxes(
                tickvals=x_ticks_val,
                ticktext=x_ticks_name,
                row=i + 1,
                col=1,
            )

        # set y axis limits
        y_min = 0.5 - exon_height / 2
        y_max = chrmd_df_grouped.loc[chrom]["y_height"]
        y_ticks_val = []
        y_ticks_name = []

        # gene names in y axis
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
            rects_df = ts_data[chrom]
            rects_df["cumdelta_end"] = rects_df[CUM_DELTA_COL]
            rects_df["cumdelta_start"] = rects_df[CUM_DELTA_COL].shift(
                periods=1, fill_value=0
            )
            rects_df[START_COL] -= rects_df["cumdelta_start"]
            rects_df[END_COL] -= rects_df["cumdelta_end"]

            for a, b, c, d in zip(
                rects_df[START_COL],
                rects_df[END_COL],
                rects_df["cumdelta_start"],
                rects_df["cumdelta_end"],
            ):
                x0, x1 = a, b
                y0, y1 = y_min - 1, y_max + 1
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, x1, x0, x0],
                        y=[y0, y0, y1, y1, y0],
                        fill="toself",
                        fillcolor=shrinked_bkg,
                        mode="lines",
                        line={"color": "lightyellow"},
                        text=f"Shrinked region:\n[{x0+c} - {x1+d}]",
                        hoverinfo="text",
                        line_width=0,
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )

            # Draw lines separating pr objects if +1
            if chrmd_df["pr_line"].drop_duplicates().max() != 0:
                pr_line_y_l = chrmd_df.loc[chrom]["pr_line"].tolist()
                if isinstance(pr_line_y_l, int):
                    pr_line_y_l = [pr_line_y_l]
                pr_line_y_l = [y_max + v_spacer] + pr_line_y_l
                present_pr_l = chrmd_df_grouped.loc[chrom]["present_pr"]

                # separate items with horizontal lines
                for j, pr_line_y in enumerate(pr_line_y_l):
                    if pr_line_y != 0:
                        # draw line
                        fig.add_hline(
                            y=pr_line_y,
                            line=dict(color=plot_border, width=1, dash="solid"),
                            row=i + 1,
                            col=1,
                        )

                        # add y_label in the middle of the subplot y axis if needed
                        if y_labels:
                            if pr_line_y_l[j + 1] != 0:
                                y_ticks_val.append(
                                    ((pr_line_y) + (pr_line_y_l[j + 1])) / 2
                                )
                            else:
                                y_ticks_val.append((pr_line_y) / 2)
                            y_ticks_name.append(y_labels[int(present_pr_l[j])])

            else:
                # pr names in y axis
                if y_labels:
                    y_ticks_val = [y_max / 2]
                    y_ticks_name = [str(y_labels)]

            fig.update_yaxes(
                range=[y_min - v_spacer, y_max + v_spacer],
                fixedrange=True,
                tickvals=y_ticks_val,
                ticktext=y_ticks_name,
                showgrid=False,
                zeroline=False,
                row=i + 1,
                col=1,
            )

    return fig
