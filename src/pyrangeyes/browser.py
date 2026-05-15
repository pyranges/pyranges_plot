"""Interactive Plotly browser controls for pyrangeyes panels."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pyranges1.core.names import CHROM_COL, START_COL, END_COL

from .core import get_engine
from .plot_main import plot

BROWSER_MODES = ("squish", "packed", "list", "zip")


@dataclass(frozen=True)
class BrowserPanel:
    """A single browser track/panel."""

    index: int
    chrom: object
    data_index: int
    data: object
    adapter: object

    @property
    def label(self):
        if self.data_index == 0:
            return str(self.chrom)
        return f"{self.chrom} [{self.data_index + 1}]"


def _as_list(value):
    return value if isinstance(value, list) else [value]


def _adapter_for(adapter, ix, n):
    if isinstance(adapter, list):
        return adapter[ix]
    if n == 1:
        return adapter
    return adapter


def _split_panels(data, adapter=None):
    panels = []
    data_l = _as_list(data)
    for data_ix, pr_obj in enumerate(data_l):
        if CHROM_COL not in pr_obj.columns:
            panel = BrowserPanel(
                len(panels),
                None,
                data_ix,
                pr_obj,
                _adapter_for(adapter, data_ix, len(data_l)),
            )
            panels.append(panel)
            continue
        chroms = list(pr_obj[CHROM_COL].drop_duplicates())
        for chrom in chroms:
            panels.append(
                BrowserPanel(
                    len(panels),
                    chrom,
                    data_ix,
                    pr_obj[pr_obj[CHROM_COL] == chrom],
                    _adapter_for(adapter, data_ix, len(data_l)),
                )
            )
    return panels


def _mode_kwargs(mode, base_kwargs, base_text, base_interval_height):
    kwargs = deepcopy(base_kwargs)
    if mode == "squish":
        kwargs.update(
            packed=True,
            text=False,
            y_labels=False,
            interval_height=min(base_interval_height, 0.25),
            v_spacer=min(kwargs.get("v_spacer", 0.2), 0.08),
        )
    elif mode == "packed":
        kwargs.update(packed=True, text=base_text if base_text is not None else True)
    elif mode == "list":
        kwargs.update(packed=False, text=False)
    elif mode == "zip":
        return kwargs
    else:
        raise ValueError(
            f"Unknown browser mode {mode!r}; expected one of {BROWSER_MODES}."
        )
    return kwargs


def _trace_values(values):
    if values is None:
        return []
    try:
        return [v for v in values if v is not None]
    except TypeError:
        return []


def _add_source_figure(master, source, row, visible):
    trace_indices = []
    annotation_indices = []
    y_values = []
    x_values = []

    for trace in source.data:
        new_trace = deepcopy(trace)
        new_trace.visible = visible
        master.add_trace(new_trace, row=row, col=1)
        trace_indices.append(len(master.data) - 1)
        x_values.extend(_trace_values(getattr(trace, "x", None)))
        y_values.extend(_trace_values(getattr(trace, "y", None)))

    for ann in source.layout.annotations or []:
        ann_dict = ann.to_plotly_json()
        ann_dict.pop("xref", None)
        ann_dict.pop("yref", None)
        ann_dict["visible"] = visible
        master.add_annotation(ann_dict, row=row, col=1)
        annotation_indices.append(len(master.layout.annotations) - 1)
        if ann_dict.get("x") is not None:
            x_values.append(ann_dict["x"])
        if ann_dict.get("y") is not None:
            y_values.append(ann_dict["y"])

    return trace_indices, annotation_indices, x_values, y_values


def _zip_summary(panel, id_col=None):
    df = panel.data
    start = int(df[START_COL].min()) if START_COL in df.columns and len(df) else 0
    end = int(df[END_COL].max()) if END_COL in df.columns and len(df) else 1
    n_intervals = len(df)
    n_groups = None
    if id_col is not None:
        cols = [id_col] if isinstance(id_col, str) else list(id_col)
        if all(col in df.columns for col in cols):
            n_groups = len(df[cols].drop_duplicates())
    group_text = f"; {n_groups} groups" if n_groups is not None else ""
    text = f"{panel.label}: {n_intervals} intervals{group_text}<br>{start:,}–{end:,}"
    return start, end, text


def _add_zip_mode(master, panel, row, visible, id_col=None):
    start, end, text = _zip_summary(panel, id_col=id_col)
    if end <= start:
        end = start + 1
    y0, y1 = 0.25, 0.75
    trace = go.Scatter(
        x=[start, end, end, start, start],
        y=[y0, y0, y1, y1, y0],
        mode="lines",
        fill="toself",
        fillcolor="rgba(180, 190, 205, 0.25)",
        line=dict(color="rgba(80, 90, 110, 0.8)", width=1),
        hoverinfo="skip",
        showlegend=False,
        visible=visible,
    )
    master.add_trace(trace, row=row, col=1)
    trace_indices = [len(master.data) - 1]
    master.add_annotation(
        x=(start + end) / 2,
        y=0.5,
        text=text,
        showarrow=False,
        visible=visible,
        xanchor="center",
        yanchor="middle",
        row=row,
        col=1,
    )
    annotation_indices = [len(master.layout.annotations) - 1]
    return trace_indices, annotation_indices, [start, end], [0, 1]


def browser(
    data,
    adapter=None,
    *,
    modes=BROWSER_MODES,
    default_mode="packed",
    button_x=1.01,
    return_plot="fig",
    **kwargs,
):
    """Create a Plotly figure with per-panel view mode controls.

    Each panel gets a compact dropdown button on its right side. The control changes
    only that panel, so panels can independently use ``squish``, ``packed``,
    ``list`` or ``zip`` views.
    """
    if get_engine() not in {"ply", "plotly"}:
        raise ValueError("browser() is Plotly-only; call pre.set_engine('ply') first.")
    modes = tuple(modes)
    unknown_modes = set(modes) - set(BROWSER_MODES)
    if unknown_modes:
        raise ValueError(f"Unknown browser mode(s): {sorted(unknown_modes)}")
    if default_mode not in modes:
        raise ValueError("default_mode must be included in modes.")

    panels = _split_panels(data, adapter=adapter)
    if not panels:
        raise ValueError("No panels to render.")

    base_kwargs = deepcopy(kwargs)
    base_text = base_kwargs.pop("text", None)
    base_interval_height = float(base_kwargs.get("interval_height", 0.8))
    id_col = base_kwargs.get("id_col")

    master = make_subplots(
        rows=len(panels),
        cols=1,
        shared_xaxes=False,
        vertical_spacing=min(0.08, 0.35 / max(len(panels), 1)),
        subplot_titles=[panel.label for panel in panels],
    )

    panel_trace_indices = []
    panel_annotation_indices = []
    panel_mode_trace_indices = []
    panel_mode_annotation_indices = []
    panel_ranges = []

    for row_ix, panel in enumerate(panels, start=1):
        mode_trace = {}
        mode_ann = {}
        all_x = []
        all_y = []
        for mode in modes:
            visible = mode == default_mode
            if mode == "zip":
                trace_indices, ann_indices, x_vals, y_vals = _add_zip_mode(
                    master, panel, row_ix, visible, id_col=id_col
                )
            else:
                mode_specific_kwargs = _mode_kwargs(
                    mode, base_kwargs, base_text, base_interval_height
                )
                source = plot(
                    panel.data,
                    adapter=panel.adapter,
                    return_plot="fig",
                    warnings=False,
                    **mode_specific_kwargs,
                )
                trace_indices, ann_indices, x_vals, y_vals = _add_source_figure(
                    master, source, row_ix, visible
                )
            mode_trace[mode] = trace_indices
            mode_ann[mode] = ann_indices
            all_x.extend(x_vals)
            all_y.extend(y_vals)

        trace_indices = [ix for mode in modes for ix in mode_trace[mode]]
        ann_indices = [ix for mode in modes for ix in mode_ann[mode]]
        panel_trace_indices.append(trace_indices)
        panel_annotation_indices.append(ann_indices)
        panel_mode_trace_indices.append(mode_trace)
        panel_mode_annotation_indices.append(mode_ann)
        if all_x:
            xmin, xmax = min(all_x), max(all_x)
            xpad = (xmax - xmin) * 0.02 or 1
            master.update_xaxes(range=[xmin - xpad, xmax + xpad], row=row_ix, col=1)
        if all_y:
            ymin, ymax = min(all_y), max(all_y)
            ypad = (ymax - ymin) * 0.08 or 0.2
            master.update_yaxes(range=[ymin - ypad, ymax + ypad], row=row_ix, col=1)
        panel_ranges.append((min(all_y) if all_y else 0, max(all_y) if all_y else 1))

    updatemenus = []
    for panel_ix, panel in enumerate(panels):
        yaxis_name = "yaxis" if panel_ix == 0 else f"yaxis{panel_ix + 1}"
        domain = master.layout[yaxis_name].domain
        y = (domain[0] + domain[1]) / 2
        buttons = []
        panel_traces = panel_trace_indices[panel_ix]
        panel_annotations = panel_annotation_indices[panel_ix]
        for mode in modes:
            trace_visible = [
                ix in panel_mode_trace_indices[panel_ix][mode] for ix in panel_traces
            ]
            layout_update = {
                f"annotations[{ix}].visible": ix
                in panel_mode_annotation_indices[panel_ix][mode]
                for ix in panel_annotations
            }
            buttons.append(
                dict(
                    label=mode,
                    method="update",
                    args=[{"visible": trace_visible}, layout_update, panel_traces],
                )
            )
        updatemenus.append(
            dict(
                type="dropdown",
                buttons=buttons,
                direction="down",
                showactive=True,
                active=list(modes).index(default_mode),
                x=button_x,
                xanchor="left",
                y=y,
                yanchor="middle",
                pad={"r": 0, "t": 0},
            )
        )

    height = max(300, 260 * len(panels))
    master.update_layout(
        updatemenus=updatemenus,
        height=height,
        margin=dict(l=60, r=150, t=40, b=40),
        showlegend=False,
    )

    if return_plot == "fig":
        return master
    return master
