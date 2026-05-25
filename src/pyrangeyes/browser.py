"""Interactive Plotly browser controls for pyrangeyes panels."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import pandas as pd
import plotly.graph_objects as go
from pyranges1.core.names import CHROM_COL, START_COL, END_COL

from .core import get_engine
from .plot_main import plot
from .track import Track

BROWSER_MODES = ("zip", "squish", "pack", "full")
_MODE_ALIASES = {"list": "full", "packed": "pack"}
_MODE_LABELS = {
    "zip": "Zip",
    "squish": "Squish",
    "pack": "Packed",
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
    def xaxis_name(self):
        return "xaxis" if self.index == 0 else f"xaxis{self.index + 1}"

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


def _mode_tracks(data, mode, base_label):
    tracks = data if isinstance(data, list) else [data]
    out = []
    for item in tracks:
        track = item if isinstance(item, Track) else Track(item)
        options = dict(track.options)
        if mode == "squish":
            options.update(pack=True, label=False, squish=True)
        elif mode == "pack":
            options.update(pack=True, label=base_label)
            options.pop("squish", None)
        elif mode == "full":
            options.update(pack=False, label=False)
            options.pop("squish", None)
        out.append(Track(track.data, track.adapter, **options))
    return out if isinstance(data, list) else out[0]


def _mode_kwargs(mode, base_kwargs, base_interval_height):
    kwargs = deepcopy(base_kwargs)
    if mode == "zip":
        return kwargs
    if mode == "squish":
        kwargs.update(
            interval_height=min(base_interval_height, 0.25),
            v_spacer=min(kwargs.get("v_spacer", 0.2), 0.08),
        )
    else:
        if mode not in {"pack", "full"}:
            raise ValueError(
                f"Unknown browser mode {mode!r}; expected one of {BROWSER_MODES}."
            )
    return kwargs


def _data_for_chrom(data, chrom):
    if chrom is None:
        return data
    if isinstance(data, list):
        result = []
        for item in data:
            if isinstance(item, Track):
                subset = item.data[item.data[CHROM_COL] == chrom]
                result.append(Track(subset, item.adapter, **item.options))
            else:
                result.append(item[item[CHROM_COL] == chrom])
        return result
    if isinstance(data, Track):
        return Track(
            data.data[data.data[CHROM_COL] == chrom], data.adapter, **data.options
        )
    return data[data[CHROM_COL] == chrom]


def _chroms_from_data(data):
    data_l = _as_list(data)
    data_l = [item.data if isinstance(item, Track) else item for item in data_l]
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


def _plot_mode_figure(data, mode, base_kwargs, base_label, base_interval_height):
    kwargs = _mode_kwargs(mode, base_kwargs, base_interval_height)
    mode_data = _mode_tracks(data, mode, base_label)
    return plot(
        mode_data,
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


def _jsonish(value):
    if value is None:
        return None
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, list):
        return [_jsonish(item) for item in value]
    return value


def _axis_range(fig, yaxis_name):
    axis = getattr(fig.layout, yaxis_name)
    if axis.range is not None:
        return list(axis.range)
    return None


def _axis_layout(fig, yaxis_name):
    axis = getattr(fig.layout, yaxis_name)
    axis_json = axis.to_plotly_json()
    layout = {}
    for key in (
        "range",
        "domain",
        "tickmode",
        "tickvals",
        "ticktext",
        "showticklabels",
        "showgrid",
        "zeroline",
        "showline",
        "linewidth",
        "linecolor",
        "mirror",
        "color",
        "fixedrange",
        "visible",
    ):
        if key in axis_json:
            layout[key] = _jsonish(axis_json[key])
    layout.setdefault("visible", True)

    # Plotly does not remove old tick labels unless they are explicitly
    # overwritten, so pack/squish/zip switches must clear labels left by full.
    layout.setdefault("tickvals", [])
    layout.setdefault("ticktext", [])
    return layout


def _mode_weight(axis_layout, mode):
    if mode == "zip":
        return 0.12
    axis_range = axis_layout.get("range") or [0, 1]
    if len(axis_range) < 2:
        return 1.0
    span = abs(float(axis_range[1]) - float(axis_range[0]))
    return span or 1.0


def _panel_pixel_height(axis_layout, mode):
    if mode == "zip":
        return 58
    span = _mode_weight(axis_layout, mode)
    if mode == "squish":
        return max(20, round(span / 0.25 * 8))
    if mode in {"pack", "full"}:
        return max(48, round(span / 0.6 * 18))
    return max(48, round(span * 30))


def _figure_height(panel_heights, gap_px=60, margin_px=120):
    if not panel_heights:
        return margin_px
    return margin_px + sum(panel_heights) + gap_px * (len(panel_heights) - 1)


def _mixed_domains(
    panels,
    base_fig,
    panel_mode_axis_layouts,
    active_panel,
    mode,
    default_mode,
    panel_title_indices=(),
    menu_indices=(),
):
    panel_modes = [
        mode if panel.yaxis_ref == active_panel.yaxis_ref else default_mode
        for panel in panels
    ]
    panel_heights = [
        _panel_pixel_height(
            panel_mode_axis_layouts[panel.yaxis_ref][panel_mode], panel_mode
        )
        for panel, panel_mode in zip(panels, panel_modes)
    ]
    gap_px = 60 if len(panels) > 1 else 0
    plot_height = sum(panel_heights) + gap_px * max(0, len(panel_heights) - 1)
    if plot_height <= 0:
        return {}

    domains = {"height": _figure_height(panel_heights, gap_px=gap_px)}
    top_px = plot_height
    for ix, (panel, panel_height) in enumerate(zip(panels, panel_heights)):
        top = top_px / plot_height
        bottom = (top_px - panel_height) / plot_height
        domains[f"{panel.yaxis_name}.domain"] = [bottom, top]
        if ix < len(panel_title_indices):
            domains[f"annotations[{panel_title_indices[ix]}].y"] = (
                top + 20 / plot_height
            )
        if ix < len(menu_indices):
            domains[f"updatemenus[{menu_indices[ix]}].y"] = top + 20 / plot_height
        top_px -= panel_height + gap_px
    return domains


def _panel_title_indices(fig, panels):
    indices = []
    for ix, annotation in enumerate(fig.layout.annotations or []):
        if len(indices) >= len(panels):
            break
        if _is_panel_annotation(annotation):
            continue
        xref = getattr(annotation, "xref", None)
        yref = getattr(annotation, "yref", None)
        if xref == "paper" and yref == "paper":
            indices.append(ix)
    return indices


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
    if isinstance(data, list):
        dfs = [item.data if isinstance(item, Track) else item for item in data]
        df = pd.concat(dfs, ignore_index=True)
    else:
        df = data.data if isinstance(data, Track) else data
    start = int(df[START_COL].min()) if START_COL in df.columns and len(df) else 0
    end = int(df[END_COL].max()) if END_COL in df.columns and len(df) else 1
    n_intervals = len(df)
    n_groups = None
    if id_col is not None:
        cols = [id_col] if isinstance(id_col, str) else list(id_col)
        if all(col in df.columns for col in cols):
            n_groups = len(df[cols].drop_duplicates())
    group_label = f"; {n_groups} groups" if n_groups is not None else ""
    text = f"{label}: {n_intervals} intervals{group_label}<br>{start:,}–{end:,}"
    return start, end, text


def _zip_trace_and_annotation(panel, id_col=None):
    start, end, label = _zip_summary(panel.data, panel.label, id_col=id_col)
    if end <= start:
        end = start + 1
    trace = go.Scatter(
        x=[],
        y=[],
        mode="markers",
        hoverinfo="skip",
        showlegend=False,
        xaxis=panel.xaxis_ref,
        yaxis=panel.yaxis_ref,
    )
    annotation = dict(
        x=(start + end) / 2,
        y=0.5,
        text=label,
        showarrow=False,
        xanchor="center",
        yanchor="middle",
        xref=panel.xaxis_ref,
        yref=panel.yaxis_ref,
    )
    return trace, annotation


def _copy_layout_from_mode(base_fig, mode_fig, panel, mode):
    yaxis_name = panel.yaxis_name
    axis_layout = _axis_layout(mode_fig, yaxis_name)
    if "range" not in axis_layout or axis_layout["range"] is None:
        axis_layout["range"] = _axis_range_from_trace_values(
            mode_fig.data, panel.yaxis_ref
        )
    return {mode: axis_layout}


def browse(
    data,
    *,
    modes=BROWSER_MODES,
    default_mode="pack",
    button_x=0.98,
    return_plot="fig",
    **kwargs,
):
    """Build an interactive Plotly interval browser.

    ``browse()`` starts from the regular ``plot()`` Plotly figure, preserving its
    layout, axes, subplot structure, titles, multi-object stacking, and standard
    interactivity. It then adds alternate per-panel interval views and compact
    mode selectors so each chromosome or input panel can switch independently
    between ``"squish"``, ``"pack"``, ``"list"``, and ``"zip"`` views.

    Parameters are the same as ``plot()`` unless listed here. ``modes`` selects
    the available views, ``default_mode`` chooses the initial view, and
    ``button_x`` positions the panel selectors. The function is Plotly-only and
    returns a Plotly figure by default; pass ``return_plot`` as for ``plot()`` to
    request another supported Plotly return type.
    """
    if get_engine() not in {"ply", "plotly"}:
        raise ValueError("browse() is Plotly-only; call pe.set_engine('ply') first.")

    modes, default_mode = _normalize_modes(modes, default_mode)
    base_kwargs = deepcopy(kwargs)
    base_label = base_kwargs.pop("label", None)
    base_interval_height = float(base_kwargs.get("interval_height", 0.8))
    id_col = base_kwargs.get("id_col")

    mode_figures = {}
    for mode in modes:
        if mode != "zip":
            mode_figures[mode] = _plot_mode_figure(
                data, mode, base_kwargs, base_label, base_interval_height
            )

    if default_mode == "zip":
        layout_mode = next((mode for mode in modes if mode != "zip"), "pack")
        base_fig = deepcopy(
            mode_figures.get(layout_mode)
            or _plot_mode_figure(
                data, "pack", base_kwargs, base_label, base_interval_height
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
    panel_mode_axis_layouts = {panel.yaxis_ref: {} for panel in panels}
    panel_mode_xaxis_layouts = {panel.yaxis_ref: {} for panel in panels}

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
            panel_mode_axis_layouts[panel.yaxis_ref][default_mode] = {
                "range": [0, 1],
                "visible": False,
                "tickvals": [],
                "ticktext": [],
            }
            panel_mode_xaxis_layouts[panel.yaxis_ref][default_mode] = {"visible": False}
        else:
            panel_mode_axis_layouts[panel.yaxis_ref].update(
                _copy_layout_from_mode(
                    base_fig, mode_figures[default_mode], panel, default_mode
                )
            )
            panel_mode_xaxis_layouts[panel.yaxis_ref][default_mode] = {"visible": True}

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
                panel_mode_axis_layouts[panel.yaxis_ref][mode] = {
                    "range": [0, 1],
                    "visible": False,
                    "tickvals": [],
                    "ticktext": [],
                }
                panel_mode_xaxis_layouts[panel.yaxis_ref][mode] = {"visible": False}
                continue

            mode_fig = mode_figures[mode]
            panel_mode_axis_layouts[panel.yaxis_ref].update(
                _copy_layout_from_mode(base_fig, mode_fig, panel, mode)
            )
            panel_mode_xaxis_layouts[panel.yaxis_ref][mode] = {"visible": True}
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

    panel_title_indices = _panel_title_indices(base_fig, panels)
    updatemenus = list(base_fig.layout.updatemenus or [])
    base_menu_count = len(updatemenus)
    menu_indices = list(range(base_menu_count, base_menu_count + len(panels)))
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
            for axis_prop, axis_value in panel_mode_axis_layouts[panel.yaxis_ref][
                mode
            ].items():
                layout_update[f"{panel.yaxis_name}.{axis_prop}"] = axis_value
            for axis_prop, axis_value in panel_mode_xaxis_layouts[panel.yaxis_ref][
                mode
            ].items():
                layout_update[f"{panel.xaxis_name}.{axis_prop}"] = axis_value
            layout_update.update(
                _mixed_domains(
                    panels,
                    base_fig,
                    panel_mode_axis_layouts,
                    panel,
                    mode,
                    default_mode,
                    panel_title_indices,
                    menu_indices,
                )
            )
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
                showactive=False,
                active=-1,
                font={"size": 10},
                bgcolor="rgba(255,255,255,0.75)",
                x=button_x,
                xanchor="right",
                y=min(domain[1] + 0.02, 1.0),
                yanchor="top",
                pad={"r": 0, "t": 0},
            )
        )

    default_panel_heights = [
        _panel_pixel_height(
            panel_mode_axis_layouts[panel.yaxis_ref][default_mode], default_mode
        )
        for panel in panels
    ]
    base_fig.update_layout(
        height=_figure_height(
            default_panel_heights, gap_px=60 if len(panels) > 1 else 0
        ),
        updatemenus=updatemenus,
        showlegend=False,
        dragmode="select",
        selectdirection="h",
    )
    if default_mode == "zip":
        for panel in panels:
            base_fig.update_xaxes(visible=False, row=panel.index + 1, col=1)
            base_fig.update_yaxes(visible=False, row=panel.index + 1, col=1)

    if return_plot == "fig":
        return base_fig
    return base_fig
