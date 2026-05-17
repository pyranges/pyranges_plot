"""Interactive Plotly browser controls for pyrangeyes panels."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import pandas as pd
import plotly.graph_objects as go
from pyranges1.core.names import CHROM_COL, START_COL, END_COL

from .core import get_engine
from .plot_main import plot

BROWSER_MODES = ("zip", "squish", "packed", "full")
_MODE_ALIASES = {"list": "full"}
_MODE_LABELS = {
    "zip": "Zip",
    "squish": "Squish",
    "packed": "Packed",
    "full": "Full",
}


@dataclass(frozen=True)
class BrowserPanel:
    """A single browser panel matching a Plotly y-axis/subplot."""

    index: int
    chrom: object
    data: object

    @property
    def label(self):
        return str(self.chrom) if self.chrom is not None else f"Panel {self.index + 1}"

    @property
    def xaxis_ref(self):
        return "x" if self.index == 0 else f"x{self.index + 1}"

    @property
    def yaxis_ref(self):
        return "y" if self.index == 0 else f"y{self.index + 1}"

    @property
    def yaxis_name(self):
        return "yaxis" if self.index == 0 else f"yaxis{self.index + 1}"


def _as_list(value):
    return value if isinstance(value, list) else [value]


def _normalize_modes(modes, default_mode):
    modes = tuple(_MODE_ALIASES.get(mode, mode) for mode in modes)
    default_mode = _MODE_ALIASES.get(default_mode, default_mode)
    unknown_modes = set(modes) - set(BROWSER_MODES)
    if unknown_modes:
        raise ValueError(f"Unknown browser mode(s): {sorted(unknown_modes)}")
    if default_mode not in modes:
        raise ValueError("default_mode must be included in modes.")
    return modes, default_mode


def _mode_kwargs(mode, base_kwargs, base_text, base_interval_height):
    kwargs = deepcopy(base_kwargs)
    if mode == "zip":
        return kwargs
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
    elif mode == "full":
        kwargs.update(packed=False, text=False)
    else:
        raise ValueError(
            f"Unknown browser mode {mode!r}; expected one of {BROWSER_MODES}."
        )
    return kwargs


def _data_for_chrom(data, chrom):
    if chrom is None:
        return data
    if isinstance(data, list):
        return [pr_obj[pr_obj[CHROM_COL] == chrom] for pr_obj in data]
    return data[data[CHROM_COL] == chrom]


def _chroms_from_data(data):
    data_l = _as_list(data)
    if not all(CHROM_COL in pr_obj.columns for pr_obj in data_l):
        return [None]
    chroms = []
    for pr_obj in data_l:
        for chrom in list(pr_obj[CHROM_COL].drop_duplicates()):
            if chrom not in chroms:
                chroms.append(chrom)
    return chroms


def _panels_for_figure(data, fig):
    """Return panels in the same row order as plot() uses for chromosome subplots."""
    chroms = _chroms_from_data(data)
    panels = [
        BrowserPanel(ix, chrom, _data_for_chrom(data, chrom))
        for ix, chrom in enumerate(chroms)
    ]
    yaxes = _layout_yaxis_names(fig)
    if len(panels) != len(yaxes):
        # Fallback for future plot() layouts that do not map 1:1 to chromosome rows.
        panels = [
            BrowserPanel(ix, chroms[ix] if ix < len(chroms) else None, data)
            for ix in range(len(yaxes))
        ]
    return panels


def _plot_mode_figure(
    data, adapter, mode, base_kwargs, base_text, base_interval_height
):
    kwargs = _mode_kwargs(mode, base_kwargs, base_text, base_interval_height)
    return plot(
        data,
        adapter=adapter,
        return_plot="fig",
        warnings=False,
        **kwargs,
    )


def _trace_axis_ref(trace, axis):
    value = getattr(trace, axis, None)
    if value is None:
        return axis[0]
    return value


def _annotation_axis_ref(annotation, axis):
    value = getattr(annotation, axis, None)
    if value is None:
        return axis[0]
    return value


def _is_panel_annotation(annotation):
    xref = getattr(annotation, "xref", None)
    yref = getattr(annotation, "yref", None)
    return xref not in {None, "paper"} and yref not in {None, "paper"}


def _is_panel_shape(shape):
    yref = getattr(shape, "yref", None)
    return yref not in {None, "paper"}


def _shape_yaxis_ref(shape):
    yref = getattr(shape, "yref", None) or "y"
    return yref.split()[0]


def _layout_yaxis_names(fig):
    names = []
    ix = 1
    while True:
        name = "yaxis" if ix == 1 else f"yaxis{ix}"
        if not hasattr(fig.layout, name):
            break
        axis = getattr(fig.layout, name)
        if axis is not None:
            names.append(name)
        ix += 1
    return names


def _axis_range(fig, yaxis_name):
    axis = getattr(fig.layout, yaxis_name)
    if axis.range is not None:
        return list(axis.range)
    return None


def _trace_values(values):
    if values is None:
        return []
    try:
        return [v for v in values if v is not None]
    except TypeError:
        return []


def _axis_range_from_trace_values(traces, yaxis_ref):
    values = []
    for trace in traces:
        if _trace_axis_ref(trace, "yaxis") == yaxis_ref:
            values.extend(_trace_values(getattr(trace, "y", None)))
    if not values:
        return [0, 1]
    ymin, ymax = min(values), max(values)
    ypad = (ymax - ymin) * 0.08 or 0.2
    return [ymin - ypad, ymax + ypad]


def _zip_summary(data, label, id_col=None):
    df = pd.concat(data, ignore_index=True) if isinstance(data, list) else data
    start = int(df[START_COL].min()) if START_COL in df.columns and len(df) else 0
    end = int(df[END_COL].max()) if END_COL in df.columns and len(df) else 1
    n_intervals = len(df)
    n_groups = None
    if id_col is not None:
        cols = [id_col] if isinstance(id_col, str) else list(id_col)
        if all(col in df.columns for col in cols):
            n_groups = len(df[cols].drop_duplicates())
    group_text = f"; {n_groups} groups" if n_groups is not None else ""
    text = f"{label}: {n_intervals} intervals{group_text}<br>{start:,}–{end:,}"
    return start, end, text


def _zip_trace_and_annotation(panel, id_col=None):
    start, end, text = _zip_summary(panel.data, panel.label, id_col=id_col)
    if end <= start:
        end = start + 1
    trace = go.Scatter(
        x=[start, end, end, start, start],
        y=[0.25, 0.25, 0.75, 0.75, 0.25],
        mode="lines",
        fill="toself",
        fillcolor="rgba(180, 190, 205, 0.25)",
        line=dict(color="rgba(80, 90, 110, 0.8)", width=1),
        hoverinfo="skip",
        showlegend=False,
        xaxis=panel.xaxis_ref,
        yaxis=panel.yaxis_ref,
    )
    annotation = dict(
        x=(start + end) / 2,
        y=0.5,
        text=text,
        showarrow=False,
        xanchor="center",
        yanchor="middle",
        xref=panel.xaxis_ref,
        yref=panel.yaxis_ref,
    )
    return trace, annotation


def _copy_layout_from_mode(base_fig, mode_fig, panel, mode):
    yaxis_name = panel.yaxis_name
    axis_range = _axis_range(mode_fig, yaxis_name)
    if axis_range is None:
        axis_range = _axis_range_from_trace_values(mode_fig.data, panel.yaxis_ref)
    return {mode: axis_range}


def browser(
    data,
    adapter=None,
    *,
    modes=BROWSER_MODES,
    default_mode="packed",
    button_x=0.98,
    return_plot="fig",
    **kwargs,
):
    """Create a Plotly figure with per-panel view controls.

    ``browser()`` deliberately reuses the regular Plotly ``plot()`` output for
    layout, axes, subplot structure, multi-object stacking, titles, and default
    interactivity. It layers alternate view traces and per-panel mode controls on
    top of that figure instead of rebuilding the plot layout from scratch.
    """
    if get_engine() not in {"ply", "plotly"}:
        raise ValueError("browser() is Plotly-only; call pre.set_engine('ply') first.")

    modes, default_mode = _normalize_modes(modes, default_mode)
    base_kwargs = deepcopy(kwargs)
    base_text = base_kwargs.pop("text", None)
    base_interval_height = float(base_kwargs.get("interval_height", 0.8))
    id_col = base_kwargs.get("id_col")

    mode_figures = {}
    for mode in modes:
        if mode != "zip":
            mode_figures[mode] = _plot_mode_figure(
                data, adapter, mode, base_kwargs, base_text, base_interval_height
            )

    if default_mode == "zip":
        layout_mode = next((mode for mode in modes if mode != "zip"), "packed")
        base_fig = deepcopy(
            mode_figures.get(layout_mode)
            or _plot_mode_figure(
                data, adapter, "packed", base_kwargs, base_text, base_interval_height
            )
        )
        for trace in base_fig.data:
            trace.visible = False
        for annotation in base_fig.layout.annotations or []:
            if _is_panel_annotation(annotation):
                annotation.visible = False
        for shape in base_fig.layout.shapes or []:
            if _is_panel_shape(shape):
                shape.visible = False
    else:
        base_fig = deepcopy(mode_figures[default_mode])

    panels = _panels_for_figure(data, base_fig)
    panel_trace_indices = {panel.yaxis_ref: [] for panel in panels}
    panel_annotation_indices = {panel.yaxis_ref: [] for panel in panels}
    panel_shape_indices = {panel.yaxis_ref: [] for panel in panels}
    panel_mode_traces = {
        panel.yaxis_ref: {mode: [] for mode in modes} for panel in panels
    }
    panel_mode_annotations = {
        panel.yaxis_ref: {mode: [] for mode in modes} for panel in panels
    }
    panel_mode_shapes = {
        panel.yaxis_ref: {mode: [] for mode in modes} for panel in panels
    }
    panel_mode_ranges = {panel.yaxis_ref: {} for panel in panels}

    for ix, trace in enumerate(base_fig.data):
        yref = _trace_axis_ref(trace, "yaxis")
        if yref in panel_trace_indices:
            panel_trace_indices[yref].append(ix)
            if default_mode != "zip":
                panel_mode_traces[yref][default_mode].append(ix)

    for ix, annotation in enumerate(base_fig.layout.annotations or []):
        if not _is_panel_annotation(annotation):
            continue
        yref = _annotation_axis_ref(annotation, "yref")
        if yref in panel_annotation_indices:
            panel_annotation_indices[yref].append(ix)
            if default_mode != "zip":
                panel_mode_annotations[yref][default_mode].append(ix)

    for ix, shape in enumerate(base_fig.layout.shapes or []):
        if not _is_panel_shape(shape):
            continue
        yref = _shape_yaxis_ref(shape)
        if yref in panel_shape_indices:
            panel_shape_indices[yref].append(ix)
            if default_mode != "zip":
                panel_mode_shapes[yref][default_mode].append(ix)

    for panel in panels:
        if default_mode == "zip":
            panel_mode_ranges[panel.yaxis_ref][default_mode] = [0, 1]
        else:
            panel_mode_ranges[panel.yaxis_ref].update(
                _copy_layout_from_mode(
                    base_fig, mode_figures[default_mode], panel, default_mode
                )
            )

    for mode in modes:
        if mode == default_mode:
            continue
        for panel in panels:
            if mode == "zip":
                trace, annotation = _zip_trace_and_annotation(panel, id_col=id_col)
                trace.visible = default_mode == "zip"
                base_fig.add_trace(trace)
                trace_ix = len(base_fig.data) - 1
                base_fig.add_annotation(
                    {**annotation, "visible": default_mode == "zip"}
                )
                ann_ix = len(base_fig.layout.annotations) - 1
                panel_trace_indices[panel.yaxis_ref].append(trace_ix)
                panel_annotation_indices[panel.yaxis_ref].append(ann_ix)
                panel_mode_traces[panel.yaxis_ref][mode].append(trace_ix)
                panel_mode_annotations[panel.yaxis_ref][mode].append(ann_ix)
                panel_mode_ranges[panel.yaxis_ref][mode] = [0, 1]
                continue

            mode_fig = mode_figures[mode]
            panel_mode_ranges[panel.yaxis_ref].update(
                _copy_layout_from_mode(base_fig, mode_fig, panel, mode)
            )
            for trace in mode_fig.data:
                if _trace_axis_ref(trace, "yaxis") != panel.yaxis_ref:
                    continue
                new_trace = deepcopy(trace)
                new_trace.visible = False
                base_fig.add_trace(new_trace)
                trace_ix = len(base_fig.data) - 1
                panel_trace_indices[panel.yaxis_ref].append(trace_ix)
                panel_mode_traces[panel.yaxis_ref][mode].append(trace_ix)
            for annotation in mode_fig.layout.annotations or []:
                if not _is_panel_annotation(annotation):
                    continue
                if _annotation_axis_ref(annotation, "yref") != panel.yaxis_ref:
                    continue
                ann_dict = annotation.to_plotly_json()
                ann_dict["visible"] = False
                base_fig.add_annotation(ann_dict)
                ann_ix = len(base_fig.layout.annotations) - 1
                panel_annotation_indices[panel.yaxis_ref].append(ann_ix)
                panel_mode_annotations[panel.yaxis_ref][mode].append(ann_ix)
            for shape in mode_fig.layout.shapes or []:
                if not _is_panel_shape(shape):
                    continue
                if _shape_yaxis_ref(shape) != panel.yaxis_ref:
                    continue
                shape_dict = shape.to_plotly_json()
                shape_dict["visible"] = False
                shapes = list(base_fig.layout.shapes or [])
                shapes.append(shape_dict)
                base_fig.update_layout(shapes=shapes)
                shape_ix = len(base_fig.layout.shapes) - 1
                panel_shape_indices[panel.yaxis_ref].append(shape_ix)
                panel_mode_shapes[panel.yaxis_ref][mode].append(shape_ix)

    if default_mode == "zip":
        for panel in panels:
            trace, annotation = _zip_trace_and_annotation(panel, id_col=id_col)
            trace.visible = True
            base_fig.add_trace(trace)
            trace_ix = len(base_fig.data) - 1
            base_fig.add_annotation({**annotation, "visible": True})
            ann_ix = len(base_fig.layout.annotations) - 1
            panel_trace_indices[panel.yaxis_ref].append(trace_ix)
            panel_annotation_indices[panel.yaxis_ref].append(ann_ix)
            panel_mode_traces[panel.yaxis_ref]["zip"].append(trace_ix)
            panel_mode_annotations[panel.yaxis_ref]["zip"].append(ann_ix)

    updatemenus = list(base_fig.layout.updatemenus or [])
    for panel in panels:
        axis = getattr(base_fig.layout, panel.yaxis_name)
        domain = axis.domain or [0, 1]
        buttons = []
        panel_traces = panel_trace_indices[panel.yaxis_ref]
        panel_annotations = panel_annotation_indices[panel.yaxis_ref]
        panel_shapes = panel_shape_indices[panel.yaxis_ref]
        for mode in modes:
            trace_visible = [
                ix in panel_mode_traces[panel.yaxis_ref][mode] for ix in panel_traces
            ]
            layout_update = {
                f"annotations[{ix}].visible": ix
                in panel_mode_annotations[panel.yaxis_ref][mode]
                for ix in panel_annotations
            }
            layout_update.update(
                {
                    f"shapes[{ix}].visible": ix
                    in panel_mode_shapes[panel.yaxis_ref][mode]
                    for ix in panel_shapes
                }
            )
            layout_update[f"{panel.yaxis_name}.range"] = panel_mode_ranges[
                panel.yaxis_ref
            ][mode]
            buttons.append(
                dict(
                    label=_MODE_LABELS[mode],
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
                xanchor="right",
                y=min(domain[1] + 0.02, 1.0),
                yanchor="top",
                pad={"r": 0, "t": 0},
            )
        )

    base_fig.update_layout(
        updatemenus=updatemenus,
        showlegend=False,
        dragmode="select",
        selectdirection="h",
    )

    if return_plot == "fig":
        return base_fig
    return base_fig
