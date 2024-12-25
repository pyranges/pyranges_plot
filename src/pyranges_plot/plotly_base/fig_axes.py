import plotly.subplots as sp
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from pyranges.core.names import START_COL, END_COL
from pyranges_plot.core import cumdelting
from pyranges_plot.names import PR_INDEX_COL, CUM_DELTA_COL


def calculate_ticks(chrom_md_grouped, chrom, num_ticks=10):
    """Calculate tick values for a given data range."""

    # Calculate range and initial tick interval
    data_min = chrom_md_grouped.loc[chrom]["min_max"][0]
    data_max = chrom_md_grouped.loc[chrom]["min_max"][1]
    data_range = data_max - data_min
    int_interval = int(data_range / (num_ticks - 1))
    if int_interval == 0:
        int_interval = 2

    # Calculate tick values
    tick_values = np.arange(start=data_min, stop=data_max, step=int_interval)
    # Adjust the last tick value if necessary
    if tick_values[-1] < data_max:
        tick_values = np.append(tick_values, tick_values[-1] + int_interval)

    return tick_values


def create_fig(
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
):
    """Generate the figure and axes fitting the data."""

    # Unify titles and start figure
    titles = [title_chr.format(**{"chrom": chrom}) for chrom in chrmd_df_grouped.index]
    titles = list(pd.Series(titles))

    # Additional plots configuration
    additional_plots = []
    custom_dict_list = []

    if add_aligned_plots is not None:
        # Check if add_aligned_plots is a list of go.Scatter or a list of tupples (go.Scatter,custom_dict)
        if isinstance(add_aligned_plots, list):
            if all(isinstance(item, go.Scatter) for item in add_aligned_plots):
                # Case: List of go.Scatter objects
                additional_plots = add_aligned_plots
            elif all(
                isinstance(item, tuple) and len(item) == 2 for item in add_aligned_plots
            ):
                # Case: List of tuples
                additional_plots = [
                    item[0] for item in add_aligned_plots
                ]  # First element of tuples
                custom_dict_list = [
                    item[1] for item in add_aligned_plots
                ]  # Second element of tuples
            else:
                raise ValueError(
                    "add_aligned_plots must be a list of go.Scatter objects or a list of tuples (go.Scatter, custom_dict)."
                )
        else:
            raise ValueError("add_aligned_plots must be a list.")

    # Defining vertical spacing between plots
    def_vertical_spacing = 0.25
    if custom_dict_list:
        for custom_dict in custom_dict_list:
            vertical_spacing = custom_dict.get("y_space", def_vertical_spacing)
    else:
        vertical_spacing = def_vertical_spacing

    # In case the user provides extra plots, these will be added to the main plot
    num_main_rows = len(titles)
    num_additional_rows = len(additional_plots) if add_aligned_plots else 0
    total_rows = num_main_rows + num_additional_rows
    row_heights_additional = [0.5] * num_additional_rows

    shared_axes = False
    if add_aligned_plots:
        shared_axes = True
        fig = sp.make_subplots(
            rows=total_rows,
            cols=1,
            row_heights=chrmd_df_grouped["y_height"].to_list() + row_heights_additional,
            subplot_titles=titles,
            shared_xaxes=shared_axes,
            vertical_spacing=vertical_spacing,
        )
    else:
        fig = sp.make_subplots(
            rows=total_rows,
            cols=1,
            row_heights=chrmd_df_grouped["y_height"].to_list() + row_heights_additional,
            subplot_titles=titles,
            shared_xaxes=shared_axes,
        )

    # one subplot per chromosome
    for i in range(total_rows):
        if i < num_main_rows:
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

            # Work with x labels
            x_ticks_val = list(calculate_ticks(chrmd_df_grouped, chrom))
            x_ticks_name = list(calculate_ticks(chrmd_df_grouped, chrom))

            # consider specified x_ticks
            if x_ticks:
                # Unpack if its dict
                if isinstance(x_ticks, dict):
                    if chrom in x_ticks.keys():
                        x_ticks_chrom = x_ticks[chrom]
                        if isinstance(x_ticks_chrom, int):
                            x_ticks_val = [
                                i
                                for i in np.linspace(
                                    int(x_ticks_val[0]),
                                    int(x_ticks_val[-1]),
                                    x_ticks_chrom,
                                )
                            ]
                            x_ticks_name = x_ticks_val

                        if isinstance(x_ticks_chrom, list):
                            x_ticks_val = x_ticks_chrom
                            x_ticks_name = x_ticks_val

                elif isinstance(x_ticks, int):
                    x_ticks_val = [
                        i
                        for i in np.linspace(
                            int(x_ticks_val[0]), int(x_ticks_val[-1]), x_ticks
                        )
                    ]
                    x_ticks_name = x_ticks_val

                elif isinstance(x_ticks, list):
                    x_ticks_val = x_ticks
                    x_ticks_name = x_ticks_val

            # adjust names, must fall within limits
            fig.update_xaxes(
                tickvals=[int(i) for i in x_ticks_val],
                ticktext=[int(i) for i in x_ticks_name],
                row=i + 1,
                col=1,
            )

            # consider introns off
            if tick_pos_d:
                # get previous default ticks
                # chrom_subdf = subdf[subdf[CHROM_COL] == chrom]
                original_ticks = x_ticks_val

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
                    if num <= chrmd_df_grouped.loc[chrom]["max"]
                ]
                x_ticks_name = sorted(to_add_val)[: len(x_ticks_val)]

                # set new ticks
                fig.update_xaxes(
                    tickvals=[int(i) for i in x_ticks_val],
                    ticktext=[int(i) for i in x_ticks_name],
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
                y_ticks_val = genesmd_df.loc[chrom]["ycoord"] + 0.5
                y_ticks_val.reset_index(PR_INDEX_COL, drop=True, inplace=True)
                y_ticks_name = y_ticks_val.index
                y_ticks_val = y_ticks_val.to_list()

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
                            fillcolor=shrunk_bkg,
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
    else:
        if add_aligned_plots is not None:
            # This corresponds to an additional aligned plot
            aligned_idx = i - num_main_rows
            fig.add_trace(additional_plots[aligned_idx], row=i + 1, col=1)

            # Adding customisation features
            if custom_dict_list and aligned_idx < len(custom_dict_list):
                custom_dict = custom_dict_list[aligned_idx]

            # Determine the y-axis domain of the current subplot
            yaxis_key = f"yaxis{i + 1}" if i > 0 else "yaxis"
            y_domain = fig.layout[yaxis_key].domain
            new_y_domain = [
                y_domain[0],
                y_domain[0] + custom_dict.get("height", 0.25),
            ]  # Increase the height by 25%
            fig.update_layout({yaxis_key: dict(domain=new_y_domain)})
            y_domain = fig.layout[yaxis_key].domain

            # Adding annotation for the subplot
            # Checking y axis length to adjust the title
            title_space = 0.025 if y_domain[1] >= 0.4 else 0
            fig.add_annotation(
                dict(
                    x=0.5,  # Center the title in the plot
                    y=y_domain[1] + title_space,
                    xref="paper",
                    yref="paper",
                    text=custom_dict["title"],
                    showarrow=False,
                    font=dict(
                        size=custom_dict.get("title_size", 18),
                        color=custom_dict.get("title_color", "black"),
                    ),
                )
            )

    return fig
