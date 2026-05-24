"""Legend helpers shared by Matplotlib and Plotly renderers."""

import pandas as pd
import matplotlib.colors as mcolors

from .names import (
    BORDER_COLOR_COL,
    COLOR_INFO,
    COLOR_LEGEND_KIND_COL,
    COLOR_LEGEND_TITLE_COL,
    COLOR_TAG_COL,
    OUTLINE_LEGEND_KIND_COL,
    OUTLINE_LEGEND_TITLE_COL,
    OUTLINE_TAG_COL,
)


def _first_present(df, column, default=None):
    if column not in df.columns or df.empty:
        return default
    return df[column].iloc[0]


def _as_label(value):
    return str(value)


def categorical_fill_items(df):
    """Return ``[(label, fill_color, outline_color), ...]`` for fill legend."""
    if _first_present(df, COLOR_LEGEND_KIND_COL) == "quantitative":
        return []
    items = []
    seen = set()
    for _, row in (
        df[[COLOR_TAG_COL, COLOR_INFO, BORDER_COLOR_COL]].drop_duplicates().iterrows()
    ):
        label = _as_label(row[COLOR_TAG_COL])
        if label in seen:
            continue
        seen.add(label)
        items.append((label, row[COLOR_INFO], row[BORDER_COLOR_COL]))
    return items


def categorical_outline_items(df):
    """Return ``[(label, outline_color), ...]`` for mapped outline legend."""
    if OUTLINE_TAG_COL not in df.columns:
        return []
    if _first_present(df, OUTLINE_LEGEND_KIND_COL) == "quantitative":
        return []
    items = []
    seen = set()
    for _, row in df[[OUTLINE_TAG_COL, BORDER_COLOR_COL]].drop_duplicates().iterrows():
        label = _as_label(row[OUTLINE_TAG_COL])
        if label in seen:
            continue
        seen.add(label)
        items.append((label, row[BORDER_COLOR_COL]))
    return items


def _quantitative_info(df, tag_col, fill_col, title_col):
    if tag_col not in df.columns:
        return None
    values = pd.to_numeric(df[tag_col], errors="coerce")
    mask = values.notna()
    if not mask.any():
        return None

    qdf = pd.DataFrame({"value": values[mask], "fill": df.loc[mask, fill_col]})
    qdf = qdf.drop_duplicates("value").sort_values("value")
    vmin = float(qdf["value"].min())
    vmax = float(qdf["value"].max())
    if vmin == vmax:
        vmax = vmin + 1

    colors = [mcolors.to_hex(color) for color in qdf["fill"]]
    if len(colors) == 1:
        colors = colors * 2
    title = _first_present(df, title_col, "")
    return {"title": title, "vmin": vmin, "vmax": vmax, "colors": colors}


def quantitative_fill_info(df):
    if _first_present(df, COLOR_LEGEND_KIND_COL) != "quantitative":
        return None
    return _quantitative_info(df, COLOR_TAG_COL, COLOR_INFO, COLOR_LEGEND_TITLE_COL)


def quantitative_outline_info(df):
    if _first_present(df, OUTLINE_LEGEND_KIND_COL) != "quantitative":
        return None
    return _quantitative_info(
        df, OUTLINE_TAG_COL, BORDER_COLOR_COL, OUTLINE_LEGEND_TITLE_COL
    )
