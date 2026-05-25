import pandas as pd
import numpy as np


# import pyranges as pr
from .core import (
    get_id_col,
    get_engine,
    print_options,
    get_options,
    get_warnings,
    set_theme,
    get_theme,
    set_options,
)
from .plot_features import prp_cmap
from .track import Track
from . import adapters
from .data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
    compute_thresh,
    subdf_assigncolor,
    assign_label_rows,
    _assign_color_channel,
    _channel_colormap,
    _is_channel_colormap,
    _is_quantitative_colormap,
    _normalize_limits_to_panels,
    _normalize_regions_to_panels,
)
from .introns_off import introns_resize, recalc_axis
from pyranges1.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL
from .names import (
    PR_INDEX_COL,
    ORISTART_COL,
    ORIEND_COL,
    SHRTHRES_COL,
    ADJSTART_COL,
    ADJEND_COL,
    CUM_DELTA_COL,
    EXON_IX_COL,
    TEXT_PAD_COL,
    TEXT_LABEL_COL,
    TEXT_START_COL,
    TEXT_END_COL,
    TEXT_MID_COL,
    TEXT_PAD_FRAC_COL,
    TEXT_PAD_Y_COL,
    TEXT_HEIGHT_COL,
    TEXT_EXTRA_Y_COL,
    TEXT_POSITION_COL,
    THICK_COL,
    SHAPE_COL,
    MARKER_SIZE_COL,
    REVERSE_COL,
    SQUISH_COL,
    SQUISH_FACTOR_COL,
    COLOR_TAG_COL,
    COLOR_INFO,
)


def _normalize_track_flags(value, n_tracks, *, name):
    """Normalize a bool or per-track bool list into ``list[bool]``."""
    if isinstance(value, bool):
        return [value] * n_tracks
    if isinstance(value, (list, tuple)):
        if len(value) != n_tracks:
            raise ValueError(
                f"{name} must be a bool or a list with one bool per track "
                f"({n_tracks} expected, got {len(value)})."
            )
        if not all(isinstance(item, bool) for item in value):
            raise TypeError(f"{name} list entries must all be bool values.")
        return list(value)
    raise TypeError(f"{name} must be a bool or a list of bool values.")


def _normalize_label_spec(label, pack, *, label_position, label_fit, label_angle):
    """Normalize the public ``label=`` argument into a rendering spec."""
    defaulted = label is None
    if label is None:
        label = bool(pack)
    if not isinstance(label, (bool, str)):
        raise TypeError("label must be None, bool, or a format string.")

    position_aliases = {"top": "above", "bottom": "below"}
    label_position = position_aliases.get(label_position, label_position)
    allowed_positions = {"left", "right", "center", "above", "below"}
    if label_position not in allowed_positions:
        public_positions = sorted(allowed_positions | set(position_aliases))
        raise ValueError(
            f"label_position must be one of {public_positions}; got {label_position!r}."
        )

    enabled = bool(label)
    return {
        "enabled": enabled,
        "defaulted": defaulted,
        "label": label if isinstance(label, str) else None,
        "position": label_position,
        "angle": label_angle,
        "fit": bool(label_fit),
        "use_label_for_fit": isinstance(label, str),
    }


def _format_text_label(row, text_spec, genename):
    if not text_spec["enabled"]:
        return ""
    label = text_spec.get("label")
    if label is None:
        return str(genename)
    return str(label).format_map(row.to_dict())


def _assign_text_group_spans(subdf, id_col, interval_height):
    """Attach post-transform group spans used to position text labels."""
    group_cols = [CHROM_COL, PR_INDEX_COL] + id_col
    grouped = subdf.groupby(group_cols, observed=True, sort=False)
    subdf[TEXT_START_COL] = grouped[START_COL].transform("min")
    subdf[TEXT_END_COL] = grouped[END_COL].transform("max")
    subdf[TEXT_MID_COL] = (subdf[TEXT_START_COL] + subdf[TEXT_END_COL]) / 2
    # Vertical text at pad=0 should sit outside the full allocated interval row,
    # not the possibly shorter rendered rectangle from height_col/adapters. This
    # matches the intron-arrow envelope and keeps UTR-only labels from appearing
    # inside the row's visual space. Squished tracks use the squished row height.
    if SQUISH_FACTOR_COL in subdf.columns:
        subdf[TEXT_HEIGHT_COL] = interval_height * subdf[SQUISH_FACTOR_COL]
    else:
        subdf[TEXT_HEIGHT_COL] = interval_height
    return subdf


def _attach_panel_y_height(chrmd_df, genesmd_df, interval_height):
    """Attach per-PyRanges-object panel height for percentage text padding."""
    grouped = genesmd_df.groupby([CHROM_COL, PR_INDEX_COL], observed=True)["ycoord"]
    panel_y_height = grouped.max() - grouped.min() + 0.5 + interval_height
    return chrmd_df.join(panel_y_height.rename("y_height"))


TRACK_ID_COL = "__pe_track_id__"
TRACK_FILL_COL = "__pe_track_fill__"
TRACK_OUTLINE_COL = "__pe_track_outline__"
TRACK_LABEL_COLOR_COL = "__pe_track_label_color__"
TRACK_DEPTH_COL = "__pe_track_depth__"

TRACK_OPTION_KEYS = {
    "name",
    "squish",
    "id_col",
    "pack",
    "fill_col",
    "outline_col",
    "label_color_col",
    "height_col",
    "depth_col",
    "shape_col",
    "colormap",
    "outline_color",
    "label",
    "label_color",
    "label_position",
    "label_fit",
    "label_angle",
}


def _as_list(value):
    if value is None or isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def _make_tag_column(df, columns, output_col):
    columns = _as_list(columns)
    if columns is None:
        return None
    for column in columns:
        if column not in df.columns:
            raise ValueError(
                f"The provided column {column!r} is not present in the given data."
            )
    if len(columns) == 1:
        df[output_col] = df[columns[0]]
    else:
        df[output_col] = list(zip(*[df[column] for column in columns]))
    return [output_col]


def _make_id_column(df, id_col):
    id_col = _as_list(id_col)
    if id_col is None:
        id_col = get_id_col() or "__interval_index__"
    id_col = _as_list(id_col)
    for id_str in id_col:
        if id_str not in df.columns and id_str != "__interval_index__":
            raise ValueError(
                "Please define a correct name of the ID column using id_col=... "
                "or Track(..., id_col=...)."
            )
    if id_col == ["__interval_index__"]:
        df[TRACK_ID_COL] = [str(i) for i in range(len(df))]
    elif len(id_col) == 1:
        df[TRACK_ID_COL] = df[id_col[0]]
    else:
        df[TRACK_ID_COL] = list(zip(*[df[column] for column in id_col]))
    return [TRACK_ID_COL]


def _legend_title_from_public_cols(columns, fallback_columns):
    columns = _as_list(columns) or _as_list(fallback_columns)
    if columns is None:
        return None
    if columns == ["__interval_index__"]:
        return "interval"
    return ", ".join(columns)


def _shared_color_key(legend_title, value):
    if legend_title is None:
        return str(value)
    return f"{legend_title}: {value}"


def _can_share_global_colormap(colormap):
    fill_cmap = _channel_colormap(colormap, "fill", fallback=prp_cmap)
    return fill_cmap not in ["direct", None, False] and not isinstance(
        fill_cmap, dict
    ) and not _is_quantitative_colormap(fill_cmap)


def _build_shared_fill_colormap(prepared_tracks, global_colormap, warnings):
    if not _can_share_global_colormap(global_colormap):
        return None
    tags = []
    seen = set()
    for meta in prepared_tracks.values():
        if not meta["uses_global_colormap"]:
            continue
        legend_title = meta["fill_legend_title"]
        for value in meta["df"][TRACK_FILL_COL].drop_duplicates():
            key = _shared_color_key(legend_title, value)
            if key not in seen:
                tags.append(key)
                seen.add(key)
    if not tags:
        return None

    fill_cmap = _channel_colormap(global_colormap, "fill", fallback=prp_cmap)
    color_df = pd.DataFrame({"tag": tags})
    color_df = _assign_color_channel(color_df, "tag", COLOR_INFO, fill_cmap, warnings)
    return dict(zip(color_df["tag"], color_df[COLOR_INFO], strict=False))


def _replace_fill_channel(colormap, fill_cmap):
    if _is_channel_colormap(colormap):
        shared = dict(colormap)
        shared["fill"] = fill_cmap
        return shared
    return fill_cmap


def _normalize_tracks(data):
    if not isinstance(data, list):
        data = [data]
    return [item if isinstance(item, Track) else Track(item) for item in data]


def _validate_track_options(track, feature_keys, adapter_keys=()):
    allowed = TRACK_OPTION_KEYS | set(feature_keys) | set(adapter_keys)
    unknown = sorted(set(track.options) - allowed)
    if unknown:
        raise Exception(
            f"Unknown Track option(s) {unknown!r}. Use plot(...) for global options, "
            "Track(..., name=..., squish=..., adapter-specific options), or print_options()."
        )


def _track_value(track, key, global_value):
    return track.options[key] if key in track.options else global_value


def _assign_label_pad_fraction(df, chrmd_df):
    """Attach text pad as a fraction of the visible span on each axis."""
    chrom = df[CHROM_COL].iloc[0]
    pr_ix = df[PR_INDEX_COL].iloc[0]
    chrmd = chrmd_df.loc[(chrom, pr_ix)]
    x_range = df[END_COL].max() - df[START_COL].min()
    if TEXT_PAD_FRAC_COL in df.columns:
        frac = df[TEXT_PAD_FRAC_COL].iloc[0]
    else:
        frac = 0 if x_range == 0 else df[TEXT_PAD_COL].iloc[0] / x_range
    df[TEXT_PAD_FRAC_COL] = [frac] * len(df)
    df[TEXT_PAD_COL] = [frac * x_range] * len(df)
    df[TEXT_PAD_Y_COL] = [frac * chrmd["y_height"]] * len(df)
    return df


# Check for matplotlib
try:
    from .matplotlib_base.plot_exons_plt import plot_exons_plt

    missing_plt_flag = 0
except ImportError:
    missing_plt_flag = 1

# Check for plotly
try:
    from .plotly_base.plot_exons_ply import plot_exons_ply

    missing_ply_flag = 0
except ImportError:
    missing_ply_flag = 1


def _attach_panel_display(chrmd_df_grouped, panel_display):
    """Attach per-panel display info (real chromosome + window) to
    chrmd_df_grouped as columns the rendering layer can read.

    When ``panel_display`` is None (no exploding happened) we still populate
    sensible defaults so the title formatter has uniform data: display_chrom
    is the index value (the real chromosome), and display_start/display_end
    are taken from the dynamic ``min_max`` tuple.
    """
    if panel_display:
        chrmd_df_grouped["display_chrom"] = [
            panel_display.get(ix, {}).get("chrom", ix) for ix in chrmd_df_grouped.index
        ]
        # Prefer the explicit window edges from panel_display when given;
        # fall back to the dynamic min_max when an edge was None.
        ds, de = [], []
        for ix in chrmd_df_grouped.index:
            info = panel_display.get(ix, {})
            mn, mx = chrmd_df_grouped.loc[ix]["min_max"]
            ds.append(info.get("start") if info.get("start") is not None else mn)
            de.append(info.get("end") if info.get("end") is not None else mx)
        chrmd_df_grouped["display_start"] = ds
        chrmd_df_grouped["display_end"] = de
    else:
        chrmd_df_grouped["display_chrom"] = list(chrmd_df_grouped.index)
        chrmd_df_grouped["display_start"] = [
            mm[0] for mm in chrmd_df_grouped["min_max"]
        ]
        chrmd_df_grouped["display_end"] = [mm[1] for mm in chrmd_df_grouped["min_max"]]


def _panel_matches_selector(panel_id, row, selector):
    """Return whether a reverse selector refers to a panel."""
    display_chrom = row.get("display_chrom", panel_id)
    display_start = row.get("display_start")
    display_end = row.get("display_end")

    if selector == panel_id or selector == display_chrom:
        return True
    if isinstance(selector, tuple) and len(selector) == 3:
        return (
            selector[0] == display_chrom
            and selector[1] == display_start
            and selector[2] == display_end
        )
    return False


def _auto_reverse_flags(subdf, panels):
    """Reverse panels where all intervals with known strand are negative."""
    flags = {}
    for panel_id in panels:
        if STRAND_COL not in subdf.columns:
            flags[panel_id] = False
            continue
        strands = subdf.loc[subdf[CHROM_COL] == panel_id, STRAND_COL].dropna()
        strands = [s for s in strands if s in {"+", "-"}]
        flags[panel_id] = bool(strands) and all(s == "-" for s in strands)
    return flags


def _normalize_reverse(reverse, subdf, chrmd_df_grouped):
    """Normalize public reverse= input into one bool per panel."""
    panels = list(chrmd_df_grouped.index)
    if isinstance(reverse, bool):
        return {panel_id: reverse for panel_id in panels}
    if reverse is None:
        return {panel_id: False for panel_id in panels}
    if isinstance(reverse, str) and reverse == "auto":
        return _auto_reverse_flags(subdf, panels)

    if isinstance(reverse, dict):
        flags = {panel_id: False for panel_id in panels}
        for selector, value in reverse.items():
            matched = False
            for panel_id, row in chrmd_df_grouped.iterrows():
                if _panel_matches_selector(panel_id, row, selector):
                    flags[panel_id] = bool(value)
                    matched = True
            if not matched:
                raise ValueError(
                    f"reverse selector {selector!r} did not match any panel."
                )
        return flags

    if isinstance(reverse, tuple) and len(reverse) == 3:
        selectors = [reverse]
    elif isinstance(reverse, (list, tuple, set)):
        selectors = list(reverse)
    else:
        selectors = [reverse]

    flags = {panel_id: False for panel_id in panels}
    for selector in selectors:
        matched = False
        for panel_id, row in chrmd_df_grouped.iterrows():
            if _panel_matches_selector(panel_id, row, selector):
                flags[panel_id] = True
                matched = True
        if not matched:
            raise ValueError(f"reverse selector {selector!r} did not match any panel.")
    return flags


def _reverse_limits(limits, reverse_flags):
    """Mirror per-panel limits for reversed panels."""
    if limits is None:
        return None
    if isinstance(limits, dict):
        out = {}
        for panel, value in limits.items():
            if value is None or not reverse_flags.get(panel, False):
                out[panel] = value
            else:
                start, end = value
                out[panel] = (
                    None if end is None else -end,
                    None if start is None else -start,
                )
        return out
    if isinstance(limits, tuple):
        # Tuple limits apply globally. Mirror only when every panel is reversed;
        # mixed manual reversal with tuple limits would be ambiguous.
        if all(reverse_flags.values()):
            start, end = limits
            return (None if end is None else -end, None if start is None else -start)
    return limits


def _apply_reverse_transform(subdf, reverse_flags):
    """Mirror coordinates and strand for panels selected by reverse=."""
    if not any(reverse_flags.values()):
        return subdf

    subdf = subdf.copy()
    mask = subdf[CHROM_COL].map(reverse_flags).fillna(False).astype(bool)
    if not mask.any():
        return subdf

    old_start = subdf.loc[mask, START_COL].copy()
    old_end = subdf.loc[mask, END_COL].copy()
    subdf.loc[mask, START_COL] = -old_end
    subdf.loc[mask, END_COL] = -old_start

    if STRAND_COL in subdf.columns:
        subdf.loc[mask, "__oriStrand__"] = subdf.loc[mask, STRAND_COL]
        subdf.loc[mask, STRAND_COL] = subdf.loc[mask, STRAND_COL].replace(
            {"+": "__pyrangeyes_plus__", "-": "+"}
        )
        subdf.loc[mask, STRAND_COL] = subdf.loc[mask, STRAND_COL].replace(
            {"__pyrangeyes_plus__": "-"}
        )
    return subdf


def _reverse_ts_data_for_plot(ts_data, reverse_flags):
    """Mirror shrink metadata for selected panels after shrink coordinates exist."""
    if not any(reverse_flags.values()):
        return ts_data
    out = {}
    for chrom, df in ts_data.items():
        if df.empty or not reverse_flags.get(chrom, False):
            out[chrom] = df
            continue
        r = df.copy()
        r[ORISTART_COL] = r[START_COL]
        r[ORIEND_COL] = r[END_COL]
        cum_start = r[CUM_DELTA_COL].shift(periods=1, fill_value=0)
        adj_start = r[START_COL] - cum_start
        adj_end = r[END_COL] - r[CUM_DELTA_COL]
        r[START_COL] = -adj_end
        r[END_COL] = -adj_start
        r[ADJSTART_COL] = r[START_COL]
        r[ADJEND_COL] = r[END_COL]
        r[CUM_DELTA_COL] = 0
        out[chrom] = r
    return out


def _assign_custom_tooltip_column(subdf, tooltip):
    """Render user tooltip templates before coordinate transforms are applied."""
    if tooltip is None:
        return subdf, tooltip

    rendered = []
    for _, row in subdf.iterrows():
        template = tooltip
        if isinstance(template, str) and template.startswith("$"):
            col = template[1:]
            if col in row.index:
                template = row[col]
        rendered.append(str(template).format_map(row.to_dict()))
    subdf = subdf.copy()
    subdf["__tooltip__"] = rendered
    return subdf, "{__tooltip__}"


def plot(
    data,
    *,
    id_col=None,
    warnings=None,
    max_shown=25,
    pack=True,
    return_plot=None,
    add_aligned_plots=None,
    fill_col=None,
    outline_col=None,
    label_color_col=None,
    shrink=False,
    limits=None,
    regions=None,
    reverse=False,
    label=None,
    legend=False,
    panel_title=None,
    tooltip=None,
    to_file=None,
    theme=None,
    sort_ranges=False,
    height_col=None,
    depth_col=None,
    shape_col=None,
    **kargs,
):
    """
    Create genes plot from 1/+ PyRanges objects.

    Parameters
    ----------
    data: {pyranges.PyRanges, Track, or list}
        One PyRanges/Track object, or a list of PyRanges/Track objects displayed as separate tracks.
        Use ``Track(data, adapter=None, name=None, **options)`` for per-track configuration.

    id_col: str, default None
        Name of the column containing gene ID.

    warnings: bool, default True
        Whether the warnings should be shown or not.

    max_shown: int, default 20
        Maximum number of genes plotted in the dataframe order.

    pack: bool, default True
        Disposition of the genes in the plot. Use True for a pack disposition (genes in the same line if
        they do not overlap) and False for unpack (one row per gene). Per-track ``Track(..., pack=...)``
        overrides this default.

    return_plot: {None, "fig", "app"}, default None
        Return the backend figure/app instead of only displaying or saving it.

    add_aligned_plots: list, default None
        Extra backend traces/axes aligned below the genomic x-axis. Currently accepts one panel/chromosome.

    fill_col: str, default None
        Name of the column used to color the interval fill. If not specified, id_col will be used.
        Values are mapped through ``colormap`` unless ``colormap="direct"`` is used.

    outline_col: str, default None
        Name of the column used to color interval outlines. If not specified, interval outlines use the
        resolved fill colors. For one fixed outline color, use ``outline_color="black"``.

    label_color_col: str, default None
        Name of the column used to color labels. If provided, values are mapped through
        ``colormap["label"]`` when present, otherwise through the fill colormap. This overrides
        the fixed ``label_color`` option.

    colormap: str, list, dict, or "direct", default "popart"
        Colors used for interval fills and, optionally, mapped outlines and labels.

        If ``"direct"``, values in ``fill_col`` and ``outline_col`` are interpreted as literal colors.
        If a string, use the named Matplotlib/Plotly colormap or color sequence.
        If a list, assign colors from the list to distinct values.
        If a dict, use a channel mapping with required ``"fill"`` and optional ``"outline"``
        and ``"label"`` entries. ``"outline": "fill"`` reuses the fill mapping.
        ``"label": None`` or an omitted ``"label"`` entry uses fixed ``label_color`` unless
        ``label_color_col`` is provided; ``"label": "fill"`` and ``"label": "outline"`` reuse
        those channels. Other label colormap specs require ``label_color_col``::

            colormap={
                "fill": {"exon": "skyblue", "CDS": "orange"},
                "outline": "fill",
                "label": {"low": "black", "high": "white"},
            }

        For quantitative coloring, use ``type="quantitative"``. Values are normalized to the observed
        min/max by default; set ``range=(min, max)`` to choose the normalization range manually::

            colormap={"type": "quantitative", "colors": "viridis"}
            colormap={"type": "quantitative", "colors": ["blue", "white", "red"], "range": (-1, 1)}

        Quantitative ``colors`` may be a named continuous colormap, a list of gradient colors, or
        normalized stops such as ``[(0, "blue"), (0.5, "white"), (1, "red")]``.

    shrink: bool, default False
        Whether to compress the intron ranges to facilitate visualization or not.

    limits: {None, dict, tuple, pyranges.PyRanges}, default None
        Customization of coordinates for the chromosome plots.

        - None: minimum and maximum exon coordinate plotted plus a 5% of the range on each side.
        - dict: {chr_name1: (min_coord, max_coord), chr_name2: (min_coord, max_coord), ...}.
          Not all the plotted chromosomes need to be specified in the dictionary and some coordinates
          can be indicated as None, both cases lead to the use of the default value.
        - tuple: the coordinate limits of all chromosomes will be defined as indicated.
        - pyranges.PyRanges: for each matching chromosome between the plotted data
          and the limits data, the limits will be defined by the minimum and maximum coordinates
          in the pyranges object defined as limits. If some plotted chromosomes are not present they
          will be left as default.

    regions: {None, list, str, pyranges.PyRanges}, default None
        Optional panel layout replacing the default one-panel-per-chromosome layout.
        If provided, panels are exactly these regions in order and ``limits`` is ignored.
        Use a list of ``(chromosome, start, end)`` tuples and/or PyRanges objects,
        a PyRanges object (one row per panel), or a column name whose values define panels.

    reverse: {bool, "auto", str, tuple, list, dict}, default False
        Mirror selected panels for transcript-direction views while keeping tick labels,
        titles, and tooltips in original genomic coordinates. Use ``True`` to reverse all
        panels; ``"auto"`` to reverse panels whose known strands are all negative; a panel
        name, ``(chromosome, start, end)`` region tuple, or list of these to reverse selected
        panels; or a dict mapping selectors to booleans.

    label: {None, bool, str}, default None
        Controls interval labels. If None, labels are enabled for pack
        plots and disabled for unpack plots to avoid duplicated row labels. If
        a track sets ``pack=False``, the default is disabled for that track.
        If True, the id/index is used; if False, labels are disabled. A string
        is interpreted as a row-value format template such as ``"{Feature}: {id}"``.
        Use ``label_pad``, ``label_size``, ``label_color``, ``label_angle``,
        ``label_position``, and ``label_fit`` to control label appearance and layout.

    legend: bool, default False
        Whether the legend should appear in the plot.

    panel_title: {None, str}, default None
        Subplot title template. Available placeholders: ``{chrom}``, ``{start}``, ``{end}``,
        ``{orientation}`` (``"fwd"``/``"rev"``), and ``{rev_flag}`` (``""``/``" (rev)"``).
        If None, pyrangeyes chooses ``"Chromosome {chrom}{rev_flag}"`` normally
        (identical to the old default unless reversed), ``"{chrom}:{start}-{end}"``
        for explicit ``regions``, and ``"{chrom}"`` when ``regions`` is a column name.

    tooltip: str, default None
        Dataframe information to show in a tooltip when placing the mouse over a gene, the given
        information will be added to the default: strand, start-end coordinates and id. This must be
        provided as a string containing the column names of the values to be shown within curly brackets.
        For example if you want to show the value of the pointed gene for the column "col1" a valid tooltip
        string could be: "Value of col1: {col1}". Note that the values in the curly brackets are not
        strings. If you want to introduce a newline you can use a newline character "\" + "n".

    to_file: {str, tuple}, default None
        Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
        Optionally, a tuple can be privided where the file name is specified as a str in the first position and in the
        second position there is a tuple specifying the height and width of the figure in px.

    theme: str, default "light"
        General color appearance of the plot. Available modes: "light", "dark", "pastel", "swimming_pool".

    sort_ranges: bool, default False
        Whether to sort interval groups by genomic coordinates before plotting.
        If False, the default, unpack plots preserve the first-seen order of rows/groups in the input.
        If True, interval groups are ordered by the internal genomic sorting behavior.

    height_col: str, default None
        Numeric column defining interval heights. Values must range from 0 to 1,
        where 1 uses the full ``interval_height`` and smaller values are rendered
        proportionally shorter. If provided, this parameter overrides the default
        uniform interval height. Usually this is set by adapters rather than
        provided directly.

    depth_col: str, default None
        Numeric column defining interval draw order for overlapping intervals. Lower values are drawn first;
        higher values are drawn later, on top of lower-depth intervals. No range constraint is applied.
        Usually this is set by adapters rather than provided directly.

    shape_col: str, default None
        Column defining interval shapes. Supported values are ``"rectangle"``,
        ``"diamond"``, ``"triangle-up"``, ``"triangle-down"``, and ``"circle"``.
        Usually this is set by adapters rather than provided directly.

    kwargs
        Customizable plot features can be defined using keyword arguments. Use print_options() function to check the variables'
        nomenclature, description and default values. Adapter-specific options are passed via
        ``Track(data, adapter, **options)``; inspect them with, for example, ``print_options(adapter="mRNA")``.



    Examples
    --------

    >>> import pyranges1 as pr, pyrangeyes as pe

    >>> pe.set_engine('plotly')

    >>> p = pr.PyRanges({
    ...     "Chromosome": ["1"] * 5,
    ...     "Strand": ["+"] * 3 + ["-"] * 2,
    ...     "Start": [10, 20, 30, 25, 40],
    ...     "End": [15, 25, 35, 30, 50],
    ...     "transcript_id": ["t1"] * 3 + ["t2"] * 2,
    ...     "Feature": ["exon", "CDS", "exon", "CDS", "exon"],
    ...     "feature1": ["A", "B", "C", "A", "B"],
    ...     "fill_hex": ["#1f77b4"] * 5,
    ...     "outline_hex": ["#333333"] * 5,
    ... })

    >>> pe.plot(p, id_col="transcript_id", max_shown=25, colormap='Set3', label=False)

    >>> pe.plot(p, id_col="transcript_id", fill_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> pe.plot(p, id_col="transcript_id", fill_col='fill_hex', outline_col='outline_hex', colormap='direct')

    >>> pe.plot(
    ...     p,
    ...     id_col="transcript_id",
    ...     fill_col='Strand',
    ...     outline_col='Feature',
    ...     colormap={'fill': {'+': 'green', '-': 'red'}, 'outline': {'exon': 'black', 'CDS': 'gold'}},
    ... )

    >>> pe.plot(p, limits={'1': (1000, 50000)}, panel_title="Chrom: {chrom}")

    >>> # Two windows on chromosome 1 shown as separate panels:
    >>> pe.plot(p, regions=[('1', 10, 30), ('1', 30, 60)])

    >>> # Or use a column to define region panels:
    >>> pe.plot(p, regions='transcript_id')

    >>> pe.plot([p, p], id_col="transcript_id", shrink=True, tooltip="Feature1: {feature1}")

    >>> pe.plot([pe.Track(p, name="first_p"), pe.Track(p, name="second_p")], id_col="transcript_id", pack=False, to_file='my_plot.pdf')
    """

    tracks = _normalize_tracks(data)
    data = [track.data for track in tracks]

    track_names = [track.options.get("name") for track in tracks]
    track_names = (
        track_names if any(name is not None for name in track_names) else False
    )

    # Deal with export
    if to_file:
        # given str file name
        if isinstance(to_file, str):
            ext = to_file[-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = (1600, 800)
        # given tuple (name, size)
        else:
            ext = to_file[0][-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = to_file[1]
            to_file = to_file[0]
    # not given to_file, store default size
    else:
        file_size = (1600, 800)

    ID_COL = [TRACK_ID_COL]

    # Deal with warnings
    if warnings is None:
        warnings = get_warnings()

    # Deal with engine
    engine = get_engine()

    # PREPARE DATA for plot
    # Deal with plot features as kargs
    feature_keys = print_options(return_keys=True)
    wrong_keys = [k for k in kargs if k not in feature_keys]
    if wrong_keys:
        raise Exception(
            f"The following keys do not match any customizable features: {wrong_keys}.\nCheck the customizable variable names using the print_options function."
        )

    def getvalue(key):
        if key in kargs:
            value = kargs[key]
            return value  ## add invalid data type??
        else:
            return get_options(key)

    # Get default plot features
    # store old options to reset them after the plot
    oldtheme = get_theme()
    oldfeat_dict = get_options("values")

    # check option modifications in params
    if theme is None:  # not specified in params, check if it was set
        theme = get_theme()
    set_theme(theme)

    feat_dict = {
        "colormap": getvalue("colormap"),
        "intron_color": getvalue("intron_color"),
        "tag_bkg": getvalue("tag_bkg"),
        "figure_bg": getvalue("figure_bg"),
        "track_bg": getvalue("track_bg"),
        "plot_border": getvalue("plot_border"),
        "title_dict_plt": {
            "family": "sans-serif",
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")) - 5,
        },
        "title_dict_ply": {
            "family": getvalue("title_font"),
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")),
        },
        "grid_color": getvalue("grid_color"),
        "outline_color": getvalue("outline_color"),
        "interval_height": float(getvalue("interval_height")),
        "transcript_utr_width": 0.3 * float(getvalue("interval_height")),
        "v_spacer": getvalue("v_spacer"),
        "label_size": float(getvalue("label_size")),
        "label_color": getvalue("label_color"),
        "label_angle": float(getvalue("label_angle")),
        "label_position": getvalue("label_position"),
        "label_fit": getvalue("label_fit"),
        "label_pad": float(getvalue("label_pad")) / 100,
        "plotly_port": getvalue("plotly_port"),
        "arrow_line_width": float(getvalue("arrow_line_width")),
        "arrow_color": getvalue("arrow_color"),
        "arrow_size": getvalue("arrow_size"),
        "shrink_threshold": getvalue("shrink_threshold"),
        "shrunk_bg": getvalue("shrunk_bg"),
        "squish_factor": float(getvalue("squish_factor")),
        "x_ticks": getvalue("x_ticks"),
    }
    if not 0 < feat_dict["squish_factor"] <= 1:
        raise ValueError("squish_factor must be > 0 and <= 1.")
    feat_dict["track_bg_by_pr"] = (
        {
            pr_ix: _track_value(track, "track_bg", feat_dict["track_bg"])
            for pr_ix, track in enumerate(tracks)
        }
        if any("track_bg" in track.options for track in tracks)
        else {}
    )

    squish_flags = [bool(track.options.get("squish", False)) for track in tracks]
    pack_flags = [bool(_track_value(track, "pack", pack)) for track in tracks]
    track_scales = {
        pr_ix: feat_dict["squish_factor"] if flag else 1.0
        for pr_ix, flag in enumerate(squish_flags)
    }
    label_specs = {
        pr_ix: _normalize_label_spec(
            _track_value(track, "label", label),
            pack_flags[pr_ix],
            label_position=_track_value(
                track, "label_position", feat_dict["label_position"]
            ),
            label_fit=_track_value(track, "label_fit", feat_dict["label_fit"]),
            label_angle=_track_value(track, "label_angle", feat_dict["label_angle"]),
        )
        for pr_ix, track in enumerate(tracks)
    }
    label = {
        "enabled": any(spec["enabled"] for spec in label_specs.values()),
        "fit": any(spec["enabled"] and spec["fit"] for spec in label_specs.values()),
        "use_label_for_fit": any(
            spec.get("use_label_for_fit") for spec in label_specs.values()
        ),
        "position": feat_dict["label_position"],
        "angle": feat_dict["label_angle"],
    }
    pack_by_track = dict(enumerate(pack_flags))
    pack_for_axes = all(pack_flags)
    shrink_threshold = feat_dict["shrink_threshold"]
    colormap = feat_dict["colormap"]
    if colormap == "popart":
        colormap = prp_cmap

    # restore options set before plot is called
    set_theme(oldtheme)
    set_options(oldfeat_dict)

    # Make DataFrame subset if needed
    df_d = {}
    prepared_tracks = {}
    tot_ngenes_l = []
    depth_cols = {}
    height_cols = {}
    shape_cols = {}
    for pr_ix, track in enumerate(tracks):
        df_item = track.data
        track_id_col = _track_value(track, "id_col", id_col)
        track_fill_col = _track_value(track, "fill_col", fill_col)
        track_outline_col = _track_value(track, "outline_col", outline_col)
        track_label_color_col = _track_value(track, "label_color_col", label_color_col)
        legend_id_col = track_id_col or get_id_col() or "__interval_index__"
        fill_legend_title = _legend_title_from_public_cols(
            track_fill_col, legend_id_col
        )
        outline_legend_title = _legend_title_from_public_cols(track_outline_col, None)
        track_height_col = _track_value(track, "height_col", height_col)
        track_depth_col = _track_value(track, "depth_col", depth_col)
        track_shape_col = _track_value(track, "shape_col", shape_col)
        adapter_name = track.adapter

        if adapter_name is not None:
            adapter_func = adapters.get(adapter_name)
            adapter_kwargs = adapters.get_options(adapter_name, "values")
            accepted_adapter_kwargs = adapters.accepted_kwargs(adapter_name)
            _validate_track_options(track, feature_keys, accepted_adapter_kwargs)
            for arg_name, arg_value in {
                "id_col": track_id_col,
                "height_col": track_height_col,
                "depth_col": track_depth_col,
                "shape_col": track_shape_col,
            }.items():
                if arg_name in accepted_adapter_kwargs and arg_value is not None:
                    adapter_kwargs[arg_name] = arg_value
            for arg_name, arg_value in track.options.items():
                if arg_name in accepted_adapter_kwargs:
                    adapter_kwargs[arg_name] = arg_value
            df_item = adapter_func(df_item, **adapter_kwargs)
            defaults = adapters.default_plot_args(adapter_name, adapter_kwargs)
            if track_id_col is None:
                track_id_col = defaults.get("id_col")
            if track_height_col is None:
                track_height_col = defaults.get("height_col")
            if track_depth_col is None:
                track_depth_col = defaults.get("depth_col")
            if track_shape_col is None:
                track_shape_col = defaults.get("shape_col")
        else:
            _validate_track_options(track, feature_keys)

        # deal with empty PyRanges
        if df_item.empty:
            continue
        df_item = df_item.copy()

        _make_id_column(df_item, track_id_col)
        if track_fill_col is None:
            df_item[TRACK_FILL_COL] = df_item[TRACK_ID_COL]
        else:
            _make_tag_column(df_item, track_fill_col, TRACK_FILL_COL)
        outline_cols = _make_tag_column(df_item, track_outline_col, TRACK_OUTLINE_COL)
        label_color_cols = _make_tag_column(
            df_item, track_label_color_col, TRACK_LABEL_COLOR_COL
        )

        df_d[pr_ix], tot_ngenes = make_subset(df_item, ID_COL, max_shown)
        tot_ngenes_l.append(tot_ngenes)

        track_colormap = _track_value(track, "colormap", colormap)
        if track_colormap == "popart":
            track_colormap = prp_cmap
        prepared_tracks[pr_ix] = {
            "df": df_d[pr_ix],
            "colormap": track_colormap,
            "uses_global_colormap": "colormap" not in track.options,
            "outline_cols": outline_cols,
            "label_color_cols": label_color_cols,
            "outline_color": _track_value(
                track, "outline_color", feat_dict["outline_color"]
            ),
            "label_color": _track_value(track, "label_color", feat_dict["label_color"]),
            "fill_legend_title": fill_legend_title,
            "outline_legend_title": outline_legend_title,
        }
        depth_cols[pr_ix] = track_depth_col
        height_cols[pr_ix] = track_height_col
        shape_cols[pr_ix] = track_shape_col

    shared_fill_colormap = _build_shared_fill_colormap(
        prepared_tracks, colormap, warnings
    )
    for pr_ix, meta in prepared_tracks.items():
        fill_cols = [TRACK_FILL_COL]
        track_colormap = meta["colormap"]
        df_to_color = meta["df"]
        if meta["uses_global_colormap"] and shared_fill_colormap is not None:
            df_to_color = df_to_color.copy()
            shared_key_col = "__pe_shared_color_key__"
            df_to_color[shared_key_col] = [
                _shared_color_key(meta["fill_legend_title"], value)
                for value in df_to_color[TRACK_FILL_COL]
            ]
            fill_cols = [shared_key_col]
            track_colormap = _replace_fill_channel(track_colormap, shared_fill_colormap)

        df_d[pr_ix] = subdf_assigncolor(
            df_to_color,
            track_colormap,
            fill_cols,
            meta["outline_cols"],
            meta["outline_color"],
            meta["label_color"],
            meta["label_color_cols"],
            warnings,
            fill_legend_title=meta["fill_legend_title"],
            outline_legend_title=meta["outline_legend_title"],
        )
        if fill_cols != [TRACK_FILL_COL]:
            df_d[pr_ix][COLOR_TAG_COL] = df_d[pr_ix][TRACK_FILL_COL]

    for tot_ngenes in tot_ngenes_l:
        if tot_ngenes > max_shown:
            subset_warn = 1
            break
        else:
            subset_warn = 0

    # concat subset dataframes and create new column with input list index
    if not df_d:
        raise Exception("The provided PyRanges object/s are empty.")
    subdf = pd.concat(df_d, names=[PR_INDEX_COL]).reset_index(
        level=PR_INDEX_COL
    )  ### change to pr but doesn't work yet!!
    subdf[SQUISH_COL] = subdf[PR_INDEX_COL].map(dict(enumerate(squish_flags)))
    subdf[SQUISH_FACTOR_COL] = subdf[PR_INDEX_COL].map(track_scales)

    # If `regions` is provided, replace the default chromosome layout with the
    # requested region panels and ignore `limits` for this call. Otherwise,
    # keep legacy `limits` behavior as coordinate customization for chromosome
    # panels.
    if panel_title is None:
        if isinstance(regions, str):
            panel_title = "{chrom}"
        elif regions is not None:
            panel_title = "{chrom}:{start}-{end}"
        else:
            panel_title = "Chromosome {chrom}{rev_flag}"

    if regions is not None:
        subdf, limits, panel_display = _normalize_regions_to_panels(subdf, regions)
    else:
        subdf, limits, panel_display = _normalize_limits_to_panels(subdf, limits)

    # group id_cols in one column to count genes in chrmd
    # if len(ID_COL) > 1:
    #   subdf["__id_col_2count__"] = list(zip(*[subdf[c] for c in ID_COL+[PR_INDEX_COL]+[CHROM_COL]]))
    # else:
    subdf["__id_col_2count__"] = list(
        zip(*[subdf[c] for c in [CHROM_COL] + [PR_INDEX_COL] + ID_COL])
    )
    if label["enabled"]:

        def _row_genename(row):
            vals = [row[c] for c in ID_COL]
            return vals[0] if len(vals) == 1 else tuple(vals)

        subdf[TEXT_LABEL_COL] = [
            ""
            if bool(row.get(SQUISH_COL, False))
            else _format_text_label(
                row,
                label_specs[int(row[PR_INDEX_COL])],
                _row_genename(row),
            )
            for _, row in subdf.iterrows()
        ]

    # Validate per-track depth columns before rendering. Higher depth values are
    # drawn later, so they appear on top of lower-depth intervals.
    depth_col_for_render = None
    if any(col is not None for col in depth_cols.values()):
        subdf[TRACK_DEPTH_COL] = 0
        for pr_ix, col in depth_cols.items():
            if col is None:
                continue
            mask = subdf[PR_INDEX_COL] == pr_ix
            if col not in subdf.columns:
                raise ValueError(
                    f"The provided depth_col {col!r} is not present in the given data."
                )
            depth_values = pd.to_numeric(subdf.loc[mask, col], errors="coerce")
            if depth_values.isna().any():
                raise ValueError(
                    f"depth_col {col!r} must contain only numeric, non-missing values."
                )
            subdf.loc[mask, TRACK_DEPTH_COL] = depth_values
        depth_col_for_render = TRACK_DEPTH_COL

    # Deal with height_col
    # set proper height values
    subdf[THICK_COL] = [feat_dict["interval_height"]] * len(subdf)
    for pr_ix, col in height_cols.items():
        if col is None:
            continue
        mask = subdf[PR_INDEX_COL] == pr_ix
        if col not in subdf.columns:
            raise ValueError(
                f"The provided height_col {col!r} is not present in the given data."
            )
        height_values = pd.to_numeric(subdf.loc[mask, col], errors="coerce")
        if height_values.isna().any():
            raise ValueError(
                f"height_col {col!r} must contain only numeric, non-missing values."
            )
        if ((height_values < 0) | (height_values > 1)).any():
            raise ValueError(
                f"height_col {col!r} values must range from 0 to 1; 1 is rendered at the full interval_height."
            )
        subdf.loc[mask, THICK_COL] = height_values * feat_dict["interval_height"]

    if any(squish_flags):
        subdf[THICK_COL] = subdf[THICK_COL] * subdf[SQUISH_FACTOR_COL]
    feat_dict["layout_interval_height"] = (
        float(subdf[THICK_COL].max())
        if any(squish_flags)
        else feat_dict["interval_height"]
    )

    subdf[SHAPE_COL] = "rectangle"
    accepted_shapes = {"rectangle", "diamond", "triangle-up", "triangle-down", "circle"}
    for pr_ix, col in shape_cols.items():
        if col is None:
            continue
        mask = subdf[PR_INDEX_COL] == pr_ix
        if col not in subdf.columns:
            raise ValueError(
                f"The provided shape_col {col!r} is not present in the given data."
            )
        subdf.loc[mask, SHAPE_COL] = subdf.loc[mask, col].fillna("rectangle")
        unknown_shapes = set(subdf.loc[mask, SHAPE_COL]) - accepted_shapes
        if unknown_shapes:
            raise ValueError(
                f"shape_col {col!r} contains unsupported shapes: {sorted(unknown_shapes)}. "
                f"Supported shapes are: {sorted(accepted_shapes)}."
            )

    if adapters.ADAPTER_MARKER_SIZE_COL in subdf.columns:
        subdf[MARKER_SIZE_COL] = subdf[adapters.ADAPTER_MARKER_SIZE_COL]
        if any(squish_flags):
            subdf[MARKER_SIZE_COL] = subdf[MARKER_SIZE_COL] * subdf[SQUISH_FACTOR_COL]

    fill_col = [TRACK_FILL_COL]

    # This is needed to maintain the order of the rows when adding multiple pr
    if len(ID_COL) == 1:
        # Only one column
        order = subdf[ID_COL[0]].drop_duplicates().tolist()
    else:
        # Multi-column, we use tuples
        order = (
            subdf[ID_COL]
            .drop_duplicates()
            .apply(lambda row: tuple(row), axis=1)
            .tolist()
        )

    # Create genes metadata DataFrame
    genesmd_df = get_genes_metadata(
        subdf,
        ID_COL,
        fill_col,
        pack_for_axes,
        feat_dict["interval_height"],
        feat_dict["v_spacer"],
        order,
        sort_ranges,
    )

    if label["enabled"]:
        labels = subdf.groupby(
            [CHROM_COL, PR_INDEX_COL] + ID_COL, observed=True, sort=False
        )[TEXT_LABEL_COL].first()
        genesmd_df = genesmd_df.join(labels.rename(TEXT_LABEL_COL))

    genesmd_df = assign_label_rows(
        genesmd_df,
        ID_COL,
        PR_INDEX_COL,
        label_pad=feat_dict["label_pad"],
        pack=pack_for_axes,
        sort_ranges=sort_ranges,
        interval_height=feat_dict["interval_height"],
        v_spacer=feat_dict["v_spacer"],
        plot_limits=None,  # You can pass limits if needed
        text_label_col=TEXT_LABEL_COL if label.get("use_label_for_fit") else None,
        text_avoid=label["enabled"] and label["fit"],
        track_scales=track_scales,
        pack_by_track=pack_by_track,
        text_position_by_track={
            pr_ix: spec["position"] for pr_ix, spec in label_specs.items()
        },
        label_size=feat_dict["label_size"],
    )

    # Create chromosome metadata DataFrame
    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf,
        limits,
        genesmd_df,
        pack_for_axes,
        feat_dict["v_spacer"],
        feat_dict["layout_interval_height"],
    )
    chrmd_df = _attach_panel_y_height(
        chrmd_df, genesmd_df, feat_dict["layout_interval_height"]
    )
    _attach_panel_display(chrmd_df_grouped, panel_display)
    reverse_flags = _normalize_reverse(reverse, subdf, chrmd_df_grouped)
    plot_limits = limits

    # Deal with introns off
    # adapt coordinates to shrunk
    ts_data = {}
    subdf[ORISTART_COL] = subdf[START_COL]
    subdf[ORIEND_COL] = subdf[END_COL]
    subdf, tooltip = _assign_custom_tooltip_column(subdf, tooltip)
    tick_pos_d = {}
    ori_tick_pos_d = {}

    if shrink:
        # compute threshold
        if isinstance(shrink_threshold, int):
            subdf[SHRTHRES_COL] = [shrink_threshold] * len(subdf)
        elif isinstance(shrink_threshold, float):
            subdf[SHRTHRES_COL] = [shrink_threshold] * len(subdf)
            subdf = subdf.groupby(CHROM_COL, group_keys=False, observed=True)[
                subdf.columns
            ].apply(
                lambda x: compute_thresh(x, chrmd_df_grouped) if not x.empty else None,
            )

        subdf = subdf.groupby(CHROM_COL, group_keys=False, observed=True)[
            subdf.columns
        ].apply(
            lambda x: introns_resize(x, ts_data, ID_COL),
        )  # empty rows when subset
        subdf[START_COL] = subdf[ADJSTART_COL]
        subdf[END_COL] = subdf[ADJEND_COL]

    else:
        subdf[CUM_DELTA_COL] = [0] * len(subdf)

    if any(reverse_flags.values()):
        subdf = _apply_reverse_transform(subdf, reverse_flags)
        ts_data = _reverse_ts_data_for_plot(ts_data, reverse_flags)
        plot_limits = _reverse_limits(limits, reverse_flags)

    genesmd_df = get_genes_metadata(
        subdf,
        ID_COL,
        fill_col,
        pack_for_axes,
        feat_dict["interval_height"],
        feat_dict["v_spacer"],
        order,
        sort_ranges,
    )

    if label["enabled"]:
        labels = subdf.groupby(
            [CHROM_COL, PR_INDEX_COL] + ID_COL, observed=True, sort=False
        )[TEXT_LABEL_COL].first()
        genesmd_df = genesmd_df.join(labels.rename(TEXT_LABEL_COL))

    genesmd_df = assign_label_rows(
        genesmd_df,
        ID_COL,
        PR_INDEX_COL,
        label_pad=feat_dict["label_pad"],
        pack=pack_for_axes,
        sort_ranges=sort_ranges,
        interval_height=feat_dict["interval_height"],
        v_spacer=feat_dict["v_spacer"],
        plot_limits=None,  # You can pass limits if needed
        text_label_col=TEXT_LABEL_COL if label.get("use_label_for_fit") else None,
        text_avoid=label["enabled"] and label["fit"],
        track_scales=track_scales,
        pack_by_track=pack_by_track,
        text_position_by_track={
            pr_ix: spec["position"] for pr_ix, spec in label_specs.items()
        },
        label_size=feat_dict["label_size"],
    )

    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf,
        plot_limits,
        genesmd_df,
        pack_for_axes,
        feat_dict["v_spacer"],
        feat_dict["layout_interval_height"],
        ts_data=ts_data if shrink else None,
    )
    chrmd_df = _attach_panel_y_height(
        chrmd_df, genesmd_df, feat_dict["layout_interval_height"]
    )
    _attach_panel_display(chrmd_df_grouped, panel_display)
    chrmd_df_grouped[REVERSE_COL] = [
        reverse_flags[chrom] for chrom in chrmd_df_grouped.index
    ]
    chrmd_df[REVERSE_COL] = [reverse_flags[chrom] for chrom, _ in chrmd_df.index]

    if shrink:
        # compute new axis values and positions if needed
        if ts_data:
            tick_pos_d, ori_tick_pos_d = recalc_axis(
                ts_data, tick_pos_d, ori_tick_pos_d
            )

    subdf.sort_values([CHROM_COL, PR_INDEX_COL] + ID_COL + [START_COL], inplace=True)
    chrmd_df.sort_values([CHROM_COL, PR_INDEX_COL], inplace=True)
    subdf[EXON_IX_COL] = subdf.groupby(
        [CHROM_COL, PR_INDEX_COL] + ID_COL, group_keys=False, observed=True
    ).cumcount()
    if label["enabled"]:
        subdf = _assign_text_group_spans(subdf, ID_COL, feat_dict["interval_height"])
    genesmd_df.sort_values([CHROM_COL, PR_INDEX_COL] + [START_COL], inplace=True)

    # Deal with label_pad
    label_pad = feat_dict["label_pad"]
    subdf[TEXT_PAD_FRAC_COL] = [label_pad] * len(subdf)
    subdf[TEXT_PAD_COL] = [label_pad] * len(subdf)
    subdf = subdf.groupby([CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True)[
        subdf.columns
    ].apply(lambda x: _assign_label_pad_fraction(x, chrmd_df))

    # Deal with added plots
    if (len(chrmd_df_grouped) > 1) and add_aligned_plots:
        raise Exception(
            f"The parameter add_aligned_plots accepts only one chromosome in the input data. The provided data contains {len(chrmd_df_grouped)}"
        )

    if "REF" in subdf.columns:
        subdf["REF"] = subdf["REF"].astype(str)
        subdf["REF"] = subdf["REF"].replace(["nan", "NaN", "None"], np.nan)

    if tooltip is None:
        # Create a list to store the updated tooltips
        updated_tooltips = []
        subdf["__tooltip__"] = ""
        for index, row in subdf.iterrows():
            if STRAND_COL in subdf.columns:
                strand = row.get("__oriStrand__", row.get(STRAND_COL))
            else:
                strand = ""
            if "REF" in subdf.columns:
                if pd.notna(row.get("REF", None)):
                    tool_str = row["REF"] + ">" + row["ALT"]
                    geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}<br>{tool_str}"
                else:
                    if strand:
                        geneinfo = f"[{strand}] ({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default with strand
                    else:
                        geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default without strand
            else:
                if strand:
                    geneinfo = f"[{strand}] ({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default with strand
                else:
                    geneinfo = f"({(row.__oriStart__)}, {(row.__oriEnd__)})<br>ID: {row['__id_col_2count__'][2]}"  # default without strand

            updated_tooltips.append(geneinfo)
        # Assign the updated tooltips back to the DataFrame
        subdf["__tooltip__"] = updated_tooltips

    if tooltip is None:
        tooltip = "{__tooltip__}"

    if return_plot is not None:
        # deal with engine and call proper plot
        if engine in ["plt", "matplotlib"]:
            if not missing_plt_flag:
                return plot_exons_plt(
                    subdf=subdf,
                    depth_col=depth_col_for_render,
                    tot_ngenes_l=tot_ngenes_l,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=False,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    track_names=track_names,
                    label=label,
                    panel_title=panel_title,
                    pack=pack_for_axes,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                )

            else:
                raise Exception(
                    "Make sure to install matplotlib dependecies by running `pip install pyranges-plot[plt]`"
                )

        elif engine in ["ply", "plotly"]:
            if not missing_ply_flag:
                return plot_exons_ply(
                    subdf=subdf,
                    depth_col=depth_col_for_render,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=False,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    track_names=track_names,
                    label=label,
                    panel_title=panel_title,
                    pack=pack_for_axes,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                    subset_warn=subset_warn,
                )
            else:
                raise Exception(
                    "Make sure to install plotly dependecies by running `pip install pyranges-plot[plotly]`"
                )

        else:
            raise Exception("Please define engine with set_engine().")
    else:
        if engine in ["plt", "matplotlib"]:
            if not missing_plt_flag:
                plot_exons_plt(
                    subdf=subdf,
                    depth_col=depth_col_for_render,
                    tot_ngenes_l=tot_ngenes_l,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=False,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    track_names=track_names,
                    label=label,
                    panel_title=panel_title,
                    pack=pack_for_axes,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                )
            else:
                raise Exception(
                    "Make sure to install matplotlib dependecies by running `pip install pyrangeyes[plt]`"
                )
        elif engine in ["ply", "plotly"]:
            if not missing_ply_flag:
                plot_exons_ply(
                    subdf=subdf,
                    depth_col=depth_col_for_render,
                    feat_dict=feat_dict,
                    genesmd_df=genesmd_df,
                    chrmd_df=chrmd_df,
                    chrmd_df_grouped=chrmd_df_grouped,
                    ts_data=ts_data,
                    max_shown=max_shown,
                    id_col=ID_COL,
                    transcript_str=False,
                    tooltip=tooltip,
                    legend=legend,
                    return_plot=return_plot,
                    add_aligned_plots=add_aligned_plots,
                    track_names=track_names,
                    label=label,
                    panel_title=panel_title,
                    pack=pack_for_axes,
                    to_file=to_file,
                    file_size=file_size,
                    warnings=warnings,
                    tick_pos_d=tick_pos_d,
                    ori_tick_pos_d=ori_tick_pos_d,
                    subset_warn=subset_warn,
                )
            else:
                raise Exception(
                    "Make sure to install plotly dependecies by running `pip install pyrangeyes[plotly]`"
                )
        else:
            raise Exception("Please define engine with set_engine().")
