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
from . import adapters
from .data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
    compute_thresh,
    subdf_assigncolor,
    assign_label_rows,
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
    THICK_COL,
    SHAPE_COL,
    MARKER_SIZE_COL,
    REVERSE_COL,
)


def _normalize_text_spec(text, packed, *, text_position, text_fit, text_angle):
    """Normalize the public ``text=`` argument into a rendering spec."""
    if text is None:
        text = bool(packed)
    if not isinstance(text, (bool, str)):
        raise TypeError("text must be None, bool, or a format string.")

    allowed_positions = {"left", "right", "center", "above", "below"}
    if text_position not in allowed_positions:
        raise ValueError(
            f"text_position must be one of {sorted(allowed_positions)}; "
            f"got {text_position!r}."
        )

    enabled = bool(text)
    return {
        "enabled": enabled,
        "label": text if isinstance(text, str) else None,
        "position": text_position,
        "angle": text_angle,
        "fit": bool(text_fit),
        "use_label_for_fit": isinstance(text, str),
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
    # inside the row's visual space.
    subdf[TEXT_HEIGHT_COL] = interval_height
    return subdf


def _attach_panel_y_height(chrmd_df, genesmd_df, interval_height):
    """Attach per-PyRanges-object panel height for percentage text padding."""
    grouped = genesmd_df.groupby([CHROM_COL, PR_INDEX_COL], observed=True)["ycoord"]
    panel_y_height = grouped.max() - grouped.min() + 0.5 + interval_height
    return chrmd_df.join(panel_y_height.rename("y_height"))


def _assign_text_pad_fraction(df, chrmd_df):
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
    adapter=None,
    *,
    id_col=None,
    warnings=None,
    max_shown=25,
    packed=True,
    return_plot=None,
    add_aligned_plots=None,
    color_col=None,
    outline_col=None,
    text_color_col=None,
    shrink=False,
    limits=None,
    regions=None,
    reverse=False,
    text=None,
    legend=False,
    title_chr=None,
    track_labels=None,
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
    data: {pyranges.PyRanges or list of pyranges.PyRanges}
        One PyRanges object, or a list of PyRanges objects displayed as separate tracks.

    adapter: {None, str, list}, default None
        Optional shortcut for a pre-configured visualization. For example,
        ``plot(annotation, "mRNA")`` renders GTF/GFF-like mRNA structure with
        thin exon/UTR regions and thick CDS regions. For a list of PyRanges,
        pass one adapter per track, such as ``plot([transcripts, variants],
        adapter=["mRNA", "SNP"])``. Use ``pe.adapters.describe()`` to list
        available adapters.

    id_col: str, default None
        Name of the column containing gene ID.

    warnings: bool, default True
        Whether the warnings should be shown or not.

    max_shown: int, default 20
        Maximum number of genes plotted in the dataframe order.

    packed: bool, default True
        Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).

    return_plot: {None, "fig", "app"}, default None
        Return the backend figure/app instead of only displaying or saving it.

    add_aligned_plots: list, default None
        Extra backend traces/axes aligned below the genomic x-axis. Currently accepts one panel/chromosome.

    color_col: str, default None
        Name of the column used to color the interval fill. If not specified, id_col will be used.
        Values are mapped through ``colormap`` unless ``colormap="direct"`` is used.

    outline_col: str, default None
        Name of the column used to color interval outlines. If not specified, interval outlines use the
        resolved fill colors. For one fixed outline color, use ``outline_color="black"``.

    text_color_col: str, default None
        Name of the column used to color text labels. If provided, values are mapped through
        ``colormap["text"]`` when present, otherwise through the fill colormap. This overrides
        the fixed ``text_color`` option.

    colormap: str, list, dict, or "direct", default "popart"
        Colors used for interval fills and, optionally, mapped outlines and text.

        If ``"direct"``, values in ``color_col`` and ``outline_col`` are interpreted as literal colors.
        If a string, use the named Matplotlib/Plotly colormap or color sequence.
        If a list, assign colors from the list to distinct values.
        If a dict, use a channel mapping with required ``"color"`` and optional ``"outline"``
        and ``"text"`` entries. ``"outline": "color"`` reuses the fill mapping.
        ``"text": None`` or an omitted ``"text"`` entry uses fixed ``text_color`` unless
        ``text_color_col`` is provided; ``"text": "color"`` and ``"text": "outline"`` reuse
        those channels. Other text colormap specs require ``text_color_col``::

            colormap={
                "color": {"exon": "skyblue", "CDS": "orange"},
                "outline": "color",
                "text": {"low": "black", "high": "white"},
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

    text: {None, bool, str}, default None
        Controls interval text annotations. If None, text is enabled for packed
        plots and disabled for unpacked plots to avoid duplicated row labels.
        If True, the id/index is used; if False, labels are disabled. A string
        is interpreted as a row-value format template such as ``"{Feature}: {id}"``.
        Use ``text_pad``, ``text_size``, ``text_color``, ``text_angle``,
        ``text_position``, and ``text_fit`` to control label appearance and layout.

    legend: bool, default False
        Whether the legend should appear in the plot.

    title_chr: {None, str}, default None
        Subplot title template. Available placeholders: ``{chrom}``, ``{start}``, ``{end}``,
        ``{orientation}`` (``"fwd"``/``"rev"``), and ``{rev_flag}`` (``""``/``" (rev)"``).
        If None, pyrangeyes chooses ``"Chromosome {chrom}{rev_flag}"`` normally
        (identical to the old default unless reversed), ``"{chrom}:{start}-{end}"``
        for explicit ``regions``, and ``"{chrom}"`` when ``regions`` is a column name.

    track_labels: list, default None
        Track labels shown when plotting multiple PyRanges objects.


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
        If False, the default, unpacked plots preserve the first-seen order of rows/groups in the input.
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
        nomenclature, description and default values. Adapter options can also be passed here when ``adapter``
        is set; inspect them with, for example, ``print_options(adapter="mRNA")``.



    Examples
    --------

    >>> import pyranges as pr, pyrangeyes as pe

    >>> pe.set_engine('plotly')

    >>> p = pr.PyRanges({"Chromosome": [1]*5, "Strand": ["+"]*3 + ["-"]*2, "Start": [10,20,30,25,40], "End": [15,25,35,30,50], "transcript_id": ["t1"]*3 + ["t2"]*2}, "feature1": ["A", "B", "C", "A", "B"])

    >>> plot(p, id_col="transcript_id",  max_shown=25, colormap='Set3', text=False)

    >>> plot(p, id_col="transcript_id", color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot(p, id_col="transcript_id", color_col='fill_hex', outline_col='outline_hex', colormap='direct')

    >>> plot(
    ...     p,
    ...     id_col="transcript_id",
    ...     color_col='Strand',
    ...     outline_col='Feature',
    ...     colormap={'color': {'+': 'green', '-': 'red'}, 'outline': {'exon': 'black', 'CDS': 'gold'}},
    ... )

    >>> plot(p, limits = {'1': (1000, 50000), '2': None, '3': (10000, None)}, title_chr="Chrom: {chrom}")

    >>> # Two windows on chromosome 1 shown as separate panels:
    >>> plot(p, regions = [('1', 1_000, 5_000), ('1', 50_000, 60_000)])

    >>> # Or use a column to define region panels:
    >>> plot(p, regions = 'transcript_id')

    >>> plot([p, p], id_col="transcript_id", shrink=True, tooltip = "Feature1: {feature1}")

    >>> plot([p, p], id_col="transcript_id", track_labels=["first_p", "second_p"], packed=False, to_file='my_plot.pdf')
    """

    # Treat input data as list
    if not isinstance(data, list):
        data = [data]

    if adapter is not None:
        if isinstance(adapter, str):
            adapter_names = [adapter] * len(data)
        else:
            adapter_names = list(adapter)
            if len(adapter_names) != len(data):
                raise ValueError(
                    "When adapter is a list, provide exactly one adapter per track."
                )

        plot_arg_values = {
            "id_col": id_col,
            "height_col": height_col,
            "depth_col": depth_col,
            "shape_col": shape_col,
        }
        adapted_data = []
        adapter_defaults = []
        consumed_kargs = set()
        for adapter_name, df_item in zip(adapter_names, data):
            adapter_func = adapters.get(adapter_name)
            adapter_kwargs = adapters.get_options(adapter_name, "values")
            accepted_adapter_kwargs = adapters.accepted_kwargs(adapter_name)
            for arg_name, arg_value in plot_arg_values.items():
                if arg_name in accepted_adapter_kwargs and arg_value is not None:
                    adapter_kwargs[arg_name] = arg_value
            for arg_name, arg_value in kargs.items():
                if arg_name in accepted_adapter_kwargs:
                    adapter_kwargs[arg_name] = arg_value
                    consumed_kargs.add(arg_name)
            adapted_data.append(adapter_func(df_item, **adapter_kwargs))
            adapter_defaults.append(
                adapters.default_plot_args(adapter_name, adapter_kwargs)
            )
        for arg_name in consumed_kargs:
            kargs.pop(arg_name)
        data = adapted_data
        default_plot_args = {}
        for defaults in adapter_defaults:
            for arg_name, default_value in defaults.items():
                if arg_name not in default_plot_args:
                    default_plot_args[arg_name] = default_value
                elif default_plot_args[arg_name] != default_value:
                    default_plot_args.pop(arg_name, None)
        if id_col is None and "id_col" in default_plot_args:
            id_col = default_plot_args["id_col"]
        if height_col is None and "height_col" in default_plot_args:
            height_col = default_plot_args["height_col"]
        if depth_col is None and "depth_col" in default_plot_args:
            depth_col = default_plot_args["depth_col"]
        if shape_col is None and "shape_col" in default_plot_args:
            shape_col = default_plot_args["shape_col"]

    # Ensure correct track_labels
    if track_labels:
        if len(track_labels) != len(data):
            raise Exception(
                f"The number of provided track_labels {track_labels} does not match the number of tracks ({len(data)})."
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

    # Deal with id column
    if id_col is None:
        ID_COL = get_id_col()
        if not ID_COL:
            ID_COL = ["__interval_index__"]
    else:
        ID_COL = id_col
    # treat as list
    if isinstance(ID_COL, str):
        ID_COL = [ID_COL]

    for df_item in data:
        for id_str in ID_COL:
            # Ensure correct names
            if (
                id_str is not None
                and id_str not in df_item.columns
                and id_str != "__interval_index__"
            ):
                raise Exception(
                    "Please define a correct name of the ID column using either set_id_col() function or plot_generic parameter as plot_generic(..., id_col = 'your_id_col')"
                )
            # Avoid Nan in id column

    # Deal with warnings
    if warnings is None:
        warnings = get_warnings()

    # Deal with engine
    engine = get_engine()

    # PREPARE DATA for plot
    # Deal with plot features as kargs
    wrong_keys = [k for k in kargs if k not in print_options(return_keys=True)]
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
        "fig_bkg": getvalue("fig_bkg"),
        "plot_bkg": getvalue("plot_bkg"),
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
        "text_size": float(getvalue("text_size")),
        "text_color": getvalue("text_color"),
        "text_angle": float(getvalue("text_angle")),
        "text_position": getvalue("text_position"),
        "text_fit": getvalue("text_fit"),
        "text_pad": float(getvalue("text_pad")) / 100,
        "plotly_port": getvalue("plotly_port"),
        "arrow_line_width": float(getvalue("arrow_line_width")),
        "arrow_color": getvalue("arrow_color"),
        "arrow_size": getvalue("arrow_size"),
        "shrink_threshold": getvalue("shrink_threshold"),
        "shrunk_bkg": getvalue("shrunk_bkg"),
        "x_ticks": getvalue("x_ticks"),
    }
    text = _normalize_text_spec(
        text,
        packed,
        text_position=feat_dict["text_position"],
        text_fit=feat_dict["text_fit"],
        text_angle=feat_dict["text_angle"],
    )
    shrink_threshold = feat_dict["shrink_threshold"]
    colormap = feat_dict["colormap"]
    if colormap == "popart":
        colormap = prp_cmap

    # restore options set before plot is called
    set_theme(oldtheme)
    set_options(oldfeat_dict)

    # Make DataFrame subset if needed
    df_d = {}
    tot_ngenes_l = []
    for pr_ix, df_item in enumerate(data):
        # deal with empty PyRanges
        if df_item.empty:
            continue
        df_item = df_item.copy()

        # consider not known id_col, plot each interval individually
        if ID_COL == ["__interval_index__"]:
            df_item["__interval_index__"] = [str(i) for i in range(len(df_item))]
            df_d[pr_ix], tot_ngenes = make_subset(
                df_item, "__interval_index__", max_shown
            )
            tot_ngenes_l.append(tot_ngenes)

        # known id_col
        else:
            df_d[pr_ix], tot_ngenes = make_subset(df_item, ID_COL, max_shown)
            tot_ngenes_l.append(tot_ngenes)

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

    # If `regions` is provided, replace the default chromosome layout with the
    # requested region panels and ignore `limits` for this call. Otherwise,
    # keep legacy `limits` behavior as coordinate customization for chromosome
    # panels.
    if title_chr is None:
        if isinstance(regions, str):
            title_chr = "{chrom}"
        elif regions is not None:
            title_chr = "{chrom}:{start}-{end}"
        else:
            title_chr = "Chromosome {chrom}{rev_flag}"

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
    if text["enabled"]:

        def _row_genename(row):
            vals = [row[c] for c in ID_COL]
            return vals[0] if len(vals) == 1 else tuple(vals)

        subdf[TEXT_LABEL_COL] = [
            _format_text_label(row, text, _row_genename(row))
            for _, row in subdf.iterrows()
        ]

    # Validate depth_col before rendering. Higher depth values are drawn later,
    # so they appear on top of lower-depth intervals when intervals overlap.
    if depth_col is not None:
        if depth_col not in subdf.columns:
            raise ValueError(
                f"The provided depth_col {depth_col!r} is not present in the given data."
            )
        depth_values = pd.to_numeric(subdf[depth_col], errors="coerce")
        if depth_values.isna().any():
            raise ValueError(
                f"depth_col {depth_col!r} must contain only numeric, non-missing values."
            )
        subdf[depth_col] = depth_values

    # Deal with height_col
    # set proper height values
    if height_col:
        # Is it present in data?
        if height_col not in subdf.columns:
            raise ValueError(
                f"The provided height_col {height_col!r} is not present in the given data."
            )

        height_values = pd.to_numeric(subdf[height_col], errors="coerce")
        if height_values.isna().any():
            raise ValueError(
                f"height_col {height_col!r} must contain only numeric, non-missing values."
            )
        if ((height_values < 0) | (height_values > 1)).any():
            raise ValueError(
                f"height_col {height_col!r} values must range from 0 to 1; "
                "1 is rendered at the full interval_height."
            )
        subdf[THICK_COL] = height_values * feat_dict["interval_height"]

    else:
        subdf[THICK_COL] = [feat_dict["interval_height"]] * len(subdf)

    if shape_col is not None:
        if shape_col not in subdf.columns:
            raise ValueError(
                f"The provided shape_col {shape_col!r} is not present in the given data."
            )
        subdf[SHAPE_COL] = subdf[shape_col].fillna("rectangle")
        accepted_shapes = {
            "rectangle",
            "diamond",
            "triangle-up",
            "triangle-down",
            "circle",
        }
        unknown_shapes = set(subdf[SHAPE_COL]) - accepted_shapes
        if unknown_shapes:
            raise ValueError(
                f"shape_col {shape_col!r} contains unsupported shapes: {sorted(unknown_shapes)}. "
                f"Supported shapes are: {sorted(accepted_shapes)}."
            )
    else:
        subdf[SHAPE_COL] = "rectangle"

    if adapters.ADAPTER_MARKER_SIZE_COL in subdf.columns:
        subdf[MARKER_SIZE_COL] = subdf[adapters.ADAPTER_MARKER_SIZE_COL]

    # Store color information in data
    # color_col as list
    if color_col is None:
        color_col = ID_COL
    elif isinstance(color_col, str):
        color_col = [color_col]

    if outline_col is not None:
        if isinstance(outline_col, str):
            outline_col = [outline_col]
        for outline_str in outline_col:
            if outline_str not in subdf.columns:
                raise Exception(
                    f"The provided outline_col {outline_str} is not present in the given data."
                )

    if text_color_col is not None:
        if isinstance(text_color_col, str):
            text_color_col = [text_color_col]
        for text_color_str in text_color_col:
            if text_color_str not in subdf.columns:
                raise Exception(
                    f"The provided text_color_col {text_color_str} is not present in the given data."
                )

    subdf = subdf_assigncolor(
        subdf,
        colormap,
        color_col,
        outline_col,
        feat_dict["outline_color"],
        feat_dict["text_color"],
        text_color_col,
        warnings,
    )

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
        color_col,
        packed,
        feat_dict["interval_height"],
        feat_dict["v_spacer"],
        order,
        sort_ranges,
    )

    if text["enabled"]:
        labels = subdf.groupby(
            [CHROM_COL, PR_INDEX_COL] + ID_COL, observed=True, sort=False
        )[TEXT_LABEL_COL].first()
        genesmd_df = genesmd_df.join(labels.rename(TEXT_LABEL_COL))

    genesmd_df = assign_label_rows(
        genesmd_df,
        ID_COL,
        PR_INDEX_COL,
        text_pad=feat_dict["text_pad"],
        packed=packed,
        sort_ranges=sort_ranges,
        plot_limits=None,  # You can pass limits if needed
        text_label_col=TEXT_LABEL_COL if text.get("use_label_for_fit") else None,
        text_avoid=text["enabled"] and text["fit"],
    )

    # Create chromosome metadata DataFrame
    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf,
        limits,
        genesmd_df,
        packed,
        feat_dict["v_spacer"],
        feat_dict["interval_height"],
    )
    chrmd_df = _attach_panel_y_height(
        chrmd_df, genesmd_df, feat_dict["interval_height"]
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
        color_col,
        packed,
        feat_dict["interval_height"],
        feat_dict["v_spacer"],
        order,
        sort_ranges,
    )

    if text["enabled"]:
        labels = subdf.groupby(
            [CHROM_COL, PR_INDEX_COL] + ID_COL, observed=True, sort=False
        )[TEXT_LABEL_COL].first()
        genesmd_df = genesmd_df.join(labels.rename(TEXT_LABEL_COL))

    genesmd_df = assign_label_rows(
        genesmd_df,
        ID_COL,
        PR_INDEX_COL,
        text_pad=feat_dict["text_pad"],
        packed=packed,
        sort_ranges=sort_ranges,
        plot_limits=None,  # You can pass limits if needed
        text_label_col=TEXT_LABEL_COL if text.get("use_label_for_fit") else None,
        text_avoid=text["enabled"] and text["fit"],
    )

    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf,
        plot_limits,
        genesmd_df,
        packed,
        feat_dict["v_spacer"],
        feat_dict["interval_height"],
        ts_data=ts_data if shrink else None,
    )
    chrmd_df = _attach_panel_y_height(
        chrmd_df, genesmd_df, feat_dict["interval_height"]
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
    if text["enabled"]:
        subdf = _assign_text_group_spans(subdf, ID_COL, feat_dict["interval_height"])
    genesmd_df.sort_values([CHROM_COL, PR_INDEX_COL] + [START_COL], inplace=True)

    # Deal with text_pad
    text_pad = feat_dict["text_pad"]
    subdf[TEXT_PAD_FRAC_COL] = [text_pad] * len(subdf)
    subdf[TEXT_PAD_COL] = [text_pad] * len(subdf)
    subdf = subdf.groupby([CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True)[
        subdf.columns
    ].apply(lambda x: _assign_text_pad_fraction(x, chrmd_df))

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
                    depth_col=depth_col,
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
                    track_labels=track_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
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
                    depth_col=depth_col,
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
                    track_labels=track_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
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
                    depth_col=depth_col,
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
                    track_labels=track_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
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
                    depth_col=depth_col,
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
                    track_labels=track_labels,
                    text=text,
                    title_chr=title_chr,
                    packed=packed,
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
