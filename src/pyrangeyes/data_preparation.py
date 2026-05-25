import warnings

import numpy as np
from intervaltree import IntervalTree
import pyranges1 as pr
from pyranges1.core.names import CHROM_COL, START_COL, END_COL
import pandas as pd

# Check for matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    missing_plt_flag = 0
except ImportError:
    missing_plt_flag = 1
# Check for plotly
try:
    import plotly.colors as pc

    missing_ply_flag = 0
except ImportError:
    missing_ply_flag = 1


from .names import (
    PR_INDEX_COL,
    SHRTHRES_COL,
    TEXT_PAD_COL,
    COLOR_INFO,
    COLOR_TAG_COL,
    BORDER_COLOR_COL,
    OUTLINE_TAG_COL,
    COLOR_LEGEND_KIND_COL,
    OUTLINE_LEGEND_KIND_COL,
    COLOR_LEGEND_TITLE_COL,
    OUTLINE_LEGEND_TITLE_COL,
    TEXT_COLOR_COL,
    TEXT_COLOR_TAG_COL,
    TEXT_EXTRA_Y_COL,
    TEXT_POSITION_COL,
    SQUISH_FACTOR_COL,
    THICK_COL,
    PANEL_SEP,
)
from .core import cumdelting, get_engine, get_warnings, check4dependency
from .matplotlib_base.core import plt_popup_warning
from .plot_features import prp_cmap


############ COMPUTE INTRONS OFF THRESHOLD
def compute_thresh(df, chrmd_df_grouped):
    """Get shrink threshold from limits"""

    chrom = df[CHROM_COL].iloc[0]
    chrmd = chrmd_df_grouped.loc[chrom]
    limit_range = chrmd["max"] - chrmd["min"]
    df[SHRTHRES_COL] = [int(df[SHRTHRES_COL].iloc[0] * limit_range)] * len(df)

    return df


############ COMPUTE TEXT PAD SIZE
def compute_tpad(df, chrmd_df_grouped):
    """Get text pad size from limits"""

    chrom = df[CHROM_COL].iloc[0]
    chrmd = chrmd_df_grouped.loc[chrom]
    limit_range = chrmd["max"] - chrmd["min"]
    df[TEXT_PAD_COL] = [df[TEXT_PAD_COL].iloc[0] * limit_range] * len(df)

    return df


############ SUBSET
def make_subset(df, id_col, max_shown):
    """Reduce the number of genes to work with."""

    # create a column indexing all the genes in the df
    df["gene_index"] = df.groupby(id_col, group_keys=False).ngroup()
    tot_ngenes = max(df["gene_index"])

    # select maximum number of genes
    if max(df.gene_index) + 1 <= max_shown:
        subdf = df
    else:
        subdf = df[df.gene_index < max_shown]

    # remove the gene_index column from the original df
    df.drop("gene_index", axis=1, inplace=True)

    return subdf, tot_ngenes


############ GENESMD_DF


###pack
def genesmd_pack(genesmd_df):
    """xxx"""

    # Initialize IntervalTree and used y-coordinates list
    trees = [IntervalTree()]

    def find_tree(row):
        for tree in trees:
            if not tree.overlaps(row[START_COL], row[END_COL]):
                return tree
        trees.append(IntervalTree())
        return trees[-1]

    # Assign y-coordinates
    for idx, row in genesmd_df.iterrows():
        tree = find_tree(row)
        tree.addi(row[START_COL], row[END_COL], idx)
        genesmd_df.at[idx, "ycoord"] = trees.index(tree)

    return genesmd_df


def update_y(genesmd_df, exon_height, v_spacer):
    """Update y coords according to previous prs"""

    # Consider pr dividing lines spot and the height of the previous pr to update y coords
    y_prev_df = (
        genesmd_df.groupby(PR_INDEX_COL)["ycoord"]
        .max()
        .shift(-1, fill_value=-(exon_height + v_spacer * 2))
        .apply(lambda x: x + (exon_height + v_spacer * 2))
        .loc[::-1]
        .cumsum()[::-1]
    )
    y_prev_df.name = "update_y_prev"
    genesmd_df = genesmd_df.join(y_prev_df, on=PR_INDEX_COL)
    genesmd_df["ycoord"] += genesmd_df["update_y_prev"]

    return genesmd_df


###colors for genes
def is_pltcolormap(colormap_string):
    """Checks whether the string given is a valid plt colormap name."""

    if check4dependency("matplotlib"):
        try:
            colormap = plt.colormaps[colormap_string]
            if colormap is not None and isinstance(colormap, mcolors.Colormap):
                return True
            else:
                return False

        except KeyError:
            return False

    else:
        return False


def is_plycolormap(colormap_string):
    """Checks whether the string given is a valid plotly color object name."""

    if check4dependency("plotly"):
        if hasattr(pc.sequential, colormap_string):
            return True
        elif hasattr(pc.diverging, colormap_string):
            return True
        elif hasattr(pc.cyclical, colormap_string):
            return True
        elif hasattr(pc.qualitative, colormap_string):
            return True


def get_plycolormap(colormap_string):
    """Provides the plotly color object corresponding to the string given."""

    if hasattr(pc.sequential, colormap_string):
        return getattr(pc.sequential, colormap_string)
    elif hasattr(pc.diverging, colormap_string):
        return getattr(pc.diverging, colormap_string)
    elif hasattr(pc.cyclical, colormap_string):
        return getattr(pc.cyclical, colormap_string)
    elif hasattr(pc.qualitative, colormap_string):
        return getattr(pc.qualitative, colormap_string)


def _is_channel_colormap(colormap):
    return isinstance(colormap, dict) and any(
        channel in colormap for channel in ("fill", "outline", "label")
    )


def _is_quantitative_colormap(colormap):
    return isinstance(colormap, dict) and colormap.get("type") == "quantitative"


def _validate_color_sequence(colors, context):
    if not isinstance(colors, list) or len(colors) == 0:
        raise ValueError(f"{context} must be a non-empty list of colors.")


def _validate_colormap_spec(colormap, context="colormap"):
    """Validate one colormap specification and fail with actionable errors."""
    if colormap in ["direct", None, False]:
        return

    if isinstance(colormap, str):
        if colormap == "popart" or is_pltcolormap(colormap) or is_plycolormap(colormap):
            return
        raise ValueError(
            f"{context}={colormap!r} is not a known Matplotlib/Plotly colormap. "
            "For a fixed outline color, use outline_color='black' instead of "
            "colormap={'outline': 'black'}."
        )

    if isinstance(colormap, list):
        _validate_color_sequence(colormap, context)
        return

    if not missing_plt_flag and isinstance(colormap, mcolors.Colormap):
        return

    if isinstance(colormap, dict):
        if _is_channel_colormap(colormap):
            allowed = {"fill", "outline", "label"}
            extra = set(colormap) - allowed
            if extra:
                raise ValueError(
                    f"{context} with channel-specific colors only accepts keys "
                    f"'fill', 'outline', and 'label'; found {sorted(extra)!r}."
                )
            if "fill" not in colormap:
                raise ValueError(f"{context} channel mapping requires a 'fill' entry.")
            _validate_colormap_spec(colormap["fill"], f"{context}['fill']")
            if "outline" in colormap:
                outline = colormap["outline"]
                if outline != "fill":
                    _validate_colormap_spec(outline, f"{context}['outline']")
            if "label" in colormap:
                label = colormap["label"]
                if label not in (None, "fill", "outline"):
                    _validate_colormap_spec(label, f"{context}['label']")
            return

        if "type" in colormap:
            allowed = {"type", "colors", "range", "na_color"}
            extra = set(colormap) - allowed
            if extra:
                raise ValueError(
                    f"{context} with type='quantitative' only accepts keys "
                    f"'type', 'colors', 'range', and 'na_color'; found {sorted(extra)!r}."
                )
            if colormap["type"] != "quantitative":
                raise ValueError(
                    f"{context} type must be 'quantitative'; found "
                    f"{colormap['type']!r}."
                )
            if "colors" not in colormap:
                raise ValueError(
                    f"{context} with type='quantitative' requires a 'colors' entry."
                )
            colors = colormap["colors"]
            if isinstance(colors, dict):
                raise ValueError(
                    f"{context} quantitative 'colors' must be a named colormap, "
                    "a list of colors, or normalized color stops; dict mappings "
                    "are only valid for categorical colors."
                )
            if isinstance(colors, str):
                if not (is_pltcolormap(colors) or is_plycolormap(colors)):
                    raise ValueError(
                        f"{context} quantitative colors={colors!r} is not a "
                        "known Matplotlib/Plotly colormap."
                    )
            elif isinstance(colors, list):
                _validate_color_sequence(colors, f"{context}['colors']")
                if all(isinstance(item, tuple) and len(item) == 2 for item in colors):
                    stops = [item[0] for item in colors]
                    if stops != sorted(stops) or stops[0] < 0 or stops[-1] > 1:
                        raise ValueError(
                            f"{context} quantitative color stops must be sorted "
                            "and normalized between 0 and 1."
                        )
                elif any(isinstance(item, tuple) for item in colors):
                    raise ValueError(
                        f"{context} quantitative color stops must all be "
                        "(position, color) tuples."
                    )
            else:
                raise ValueError(
                    f"{context} quantitative 'colors' must be a string or list."
                )

            value_range = colormap.get("range")
            if value_range is not None:
                if not isinstance(value_range, tuple) or len(value_range) != 2:
                    raise ValueError(
                        f"{context} quantitative range must be a (min, max) tuple."
                    )
                low, high = value_range
                if low is not None and high is not None and low >= high:
                    raise ValueError(
                        f"{context} quantitative range min must be smaller than max."
                    )
            return

        # Plain dictionaries are categorical value-to-color mappings.
        return

    raise ValueError(
        f"{context} must be 'direct', a colormap name, a color list, a mapping, "
        "or a quantitative colormap spec."
    )


def _channel_colormap(colormap, channel, fallback=None):
    """Return the colormap configuration for a style channel."""
    if _is_channel_colormap(colormap):
        return colormap.get(channel, fallback)
    return colormap


def _colors_to_quantitative_cmap(colors):
    if missing_plt_flag:
        raise ImportError("Quantitative colormaps require matplotlib.colors.")

    if isinstance(colors, str):
        if is_pltcolormap(colors):
            return plt.get_cmap(colors)
        return mcolors.LinearSegmentedColormap.from_list(
            colors, get_plycolormap(colors)
        )

    if all(isinstance(item, tuple) and len(item) == 2 for item in colors):
        return mcolors.LinearSegmentedColormap.from_list("pyrangeyes_quant", colors)

    return mcolors.LinearSegmentedColormap.from_list("pyrangeyes_quant", colors)


def _assign_quantitative_color_channel(subdf, tag_col, output_col, colormap):
    values = pd.to_numeric(subdf[tag_col], errors="coerce")
    invalid = subdf[tag_col].notna() & values.isna()
    if invalid.any():
        examples = subdf.loc[invalid, tag_col].drop_duplicates().head(3).tolist()
        raise ValueError(
            f"Quantitative colormap requires numeric values in {tag_col!r}; "
            f"found non-numeric value(s): {examples!r}."
        )

    na_color = colormap.get("na_color", "black")
    low, high = colormap.get("range", (None, None))
    finite_values = values[np.isfinite(values)]
    if low is None:
        low = finite_values.min() if len(finite_values) else 0
    if high is None:
        high = finite_values.max() if len(finite_values) else 1
    if low >= high:
        raise ValueError("Quantitative colormap range min must be smaller than max.")

    cmap = _colors_to_quantitative_cmap(colormap["colors"])
    norm = mcolors.Normalize(vmin=low, vmax=high, clip=True)
    subdf[output_col] = [
        na_color if pd.isna(value) else mcolors.to_hex(cmap(norm(value)))
        for value in values
    ]
    return subdf


def _assign_color_channel(subdf, tag_col, output_col, colormap, warnings):
    """Resolve a tag column into concrete colors."""
    _validate_colormap_spec(colormap)
    if _is_quantitative_colormap(colormap):
        return _assign_quantitative_color_channel(subdf, tag_col, output_col, colormap)

    color_tags = subdf[tag_col].drop_duplicates()
    n_color_tags = len(color_tags)

    if colormap in ["direct", None, False]:
        subdf[output_col] = subdf[tag_col]
        return subdf
    if colormap == "popart":
        colormap = prp_cmap

    # 0-string to colormap object if possible
    if isinstance(colormap, str):
        if is_pltcolormap(colormap):
            colormap = plt.get_cmap(colormap)
        elif is_plycolormap(colormap):
            colormap = get_plycolormap(colormap)
        else:
            raise Exception(
                "The provided string does not match any installed dependency colormap."
            )

    # 1-plt colormap to list
    if not missing_plt_flag and isinstance(colormap, mcolors.ListedColormap):
        # Many continuous Matplotlib colormaps (e.g. viridis) are implemented
        # as 256-color ListedColormaps.  Taking the first N entries makes all
        # early categories nearly identical, so sample large palettes evenly.
        if getattr(colormap, "N", len(colormap.colors)) > 20:
            positions = np.linspace(0, 1, max(n_color_tags, 1))
            colormap = [mcolors.to_hex(colormap(position)) for position in positions]
        else:
            colormap = list(colormap.colors)  # colors of plt object
            colormap = [mcolors.to_hex(color) for color in colormap]
    elif not missing_plt_flag and isinstance(colormap, mcolors.Colormap):
        positions = np.linspace(0, 1, max(n_color_tags, 1))
        colormap = [mcolors.to_hex(colormap(position)) for position in positions]

    # 2-list to dict
    if isinstance(colormap, list):
        # adjust number of colors
        if n_color_tags > len(colormap):
            engine = get_engine()
            if warnings is None:
                warnings = get_warnings()
            # A one-color list is an intentional uniform-color shorthand, not
            # an accidental palette wrap.
            if len(colormap) > 1 and engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "The genes are colored by iterating over the given color list."
                )
            elif len(colormap) > 1 and engine in ["ply", "plotly"] and warnings:
                subdf["_iterwarning!"] = [1] * len(subdf)
        else:
            colormap = colormap[:n_color_tags]
        colormap = {
            str(color_tags.iloc[i]): colormap[i % len(colormap)]
            for i in range(n_color_tags)
        }

    # 3- Use dict to assign color to gene
    if isinstance(colormap, dict):
        subdf[tag_col] = subdf[tag_col].astype(str)
        subdf[output_col] = subdf[tag_col].map(colormap)

        # add black genes warning if needed
        if subdf[output_col].isna().any():
            engine = get_engine()
            if warnings is None:
                warnings = get_warnings()
            if engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "Some genes do not have a color assigned so they are colored in black."
                )
            elif engine in ["ply", "plotly"] and warnings:
                subdf["_blackwarning!"] = [1] * len(subdf)
            subdf.fillna({output_col: "black"}, inplace=True)

    return subdf


def _legend_kind(colormap):
    if _is_quantitative_colormap(colormap):
        return "quantitative"
    return "categorical"


def _legend_title(cols, fallback):
    if cols is None:
        return fallback
    return ", ".join(cols)


def subdf_assigncolor(
    subdf,
    colormap,
    fill_col,
    outline_col,
    outline_color,
    label_color,
    label_color_col,
    warnings,
    fill_legend_title=None,
    outline_legend_title=None,
):
    """Add fill, outline, and text color information to data."""
    _validate_colormap_spec(colormap)

    # Create COLOR_COL column
    if len(fill_col) > 1:
        subdf[COLOR_TAG_COL] = list(zip(*[subdf[c] for c in fill_col]))
    else:
        subdf[COLOR_TAG_COL] = subdf[fill_col[0]]

    color_cmap = _channel_colormap(colormap, "fill", fallback=prp_cmap)
    subdf = _assign_color_channel(
        subdf, COLOR_TAG_COL, COLOR_INFO, color_cmap, warnings
    )
    subdf[COLOR_LEGEND_KIND_COL] = _legend_kind(color_cmap)
    subdf[COLOR_LEGEND_TITLE_COL] = fill_legend_title or _legend_title(fill_col, "fill")

    outline_spec = (
        _channel_colormap(colormap, "outline", fallback="fill")
        if _is_channel_colormap(colormap)
        else "fill"
    )
    if (
        outline_color is not None
        and _is_channel_colormap(colormap)
        and "outline" in colormap
    ):
        raise ValueError(
            "Do not provide both colormap['outline'] and outline_color. "
            "Use outline_color for one fixed outline color, or outline_col "
            "with colormap['outline'] for mapped outlines."
        )

    if outline_color is not None:
        subdf[BORDER_COLOR_COL] = [outline_color] * len(subdf)
        subdf[OUTLINE_LEGEND_KIND_COL] = "fixed"
        subdf[OUTLINE_LEGEND_TITLE_COL] = "outline"
        resolved_outline_cmap = None
    elif outline_col is None:
        if outline_spec != "fill":
            raise ValueError(
                "Mapped colormap['outline'] requires outline_col. For one fixed "
                "outline color, use outline_color='black'."
            )
        subdf[BORDER_COLOR_COL] = subdf[COLOR_INFO]
        subdf[OUTLINE_LEGEND_KIND_COL] = "same"
        subdf[OUTLINE_LEGEND_TITLE_COL] = subdf[COLOR_LEGEND_TITLE_COL]
        resolved_outline_cmap = color_cmap
    else:
        if len(outline_col) > 1:
            subdf[OUTLINE_TAG_COL] = list(zip(*[subdf[c] for c in outline_col]))
        else:
            subdf[OUTLINE_TAG_COL] = subdf[outline_col[0]]
        resolved_outline_cmap = color_cmap if outline_spec == "fill" else outline_spec
        subdf = _assign_color_channel(
            subdf, OUTLINE_TAG_COL, BORDER_COLOR_COL, resolved_outline_cmap, warnings
        )
        subdf[OUTLINE_LEGEND_KIND_COL] = _legend_kind(resolved_outline_cmap)
        subdf[OUTLINE_LEGEND_TITLE_COL] = outline_legend_title or _legend_title(
            outline_col, "outline"
        )

    label_spec = (
        _channel_colormap(colormap, "label", fallback=None)
        if _is_channel_colormap(colormap)
        else None
    )
    if label_color_col is None:
        if label_spec is None:
            subdf[TEXT_COLOR_COL] = [label_color] * len(subdf)
        elif label_spec == "fill":
            subdf[TEXT_COLOR_COL] = subdf[COLOR_INFO]
        elif label_spec == "outline":
            subdf[TEXT_COLOR_COL] = subdf[BORDER_COLOR_COL]
        else:
            raise ValueError(
                "colormap['label'] requires label_color_col unless it is None, "
                "'fill', or 'outline'."
            )
    else:
        if len(label_color_col) > 1:
            subdf[TEXT_COLOR_TAG_COL] = list(zip(*[subdf[c] for c in label_color_col]))
        else:
            subdf[TEXT_COLOR_TAG_COL] = subdf[label_color_col[0]]
        if label_spec in (None, "fill"):
            label_cmap = color_cmap
        elif label_spec == "outline":
            label_cmap = resolved_outline_cmap or color_cmap
        else:
            label_cmap = label_spec
        subdf = _assign_color_channel(
            subdf, TEXT_COLOR_TAG_COL, TEXT_COLOR_COL, label_cmap, warnings
        )

    return subdf


def codes(vals, desc=False):
    """Function for ordering multiindex df"""
    c, _ = pd.factorize(vals)
    return (c.max() - c) if desc else c


def get_genes_metadata(
    df, id_col, fill_col, pack, exon_height, v_spacer, order, sort_ranges
):
    """Create genes metadata df."""

    # Check if Chromosome column has mixed types
    chrom_dtype = pd.api.types.infer_dtype(df[CHROM_COL], skipna=True)
    if "mixed" in chrom_dtype:
        warnings.warn(
            "The Chromosome column contains mixed data types. Please ensure all values are of the same type.",
            UserWarning,
        )

    # Start df with chromosome and the column defining color
    # Define the aggregation functions for each column
    agg_funcs = {
        col: "first"
        for col in id_col + fill_col
        # if col not in [START_COL, END_COL, PR_INDEX_COL]
    }
    agg_funcs[START_COL] = "min"
    agg_funcs[END_COL] = "max"
    # workaround for Chromosome in fill_col list
    if CHROM_COL in fill_col:
        genesmd_df = (
            df.groupby(
                [CHROM_COL, PR_INDEX_COL] + id_col,
                group_keys=False,
                observed=True,
                sort=sort_ranges,
            ).agg(agg_funcs)
            # .reset_index(level=[PR_INDEX_COL, CHROM_COL])
        )
        genesmd_df["chromosome"] = genesmd_df[CHROM_COL]
        for i in range(len(fill_col)):
            if fill_col[i] == CHROM_COL:
                fill_col[i] = "chromosome"

    else:
        genesmd_df = (
            df.groupby(
                [CHROM_COL, PR_INDEX_COL] + id_col,
                group_keys=False,
                observed=True,
                sort=sort_ranges,
            ).agg(agg_funcs)
            # .reset_index(level=[PR_INDEX_COL, CHROM_COL])
        )

    genesmd_df["chrix"] = genesmd_df.groupby(
        CHROM_COL, group_keys=False, observed=True
    ).ngroup()

    # Sort by pr_ix and chromosome / If user wants to sort the df
    if sort_ranges:
        genesmd_df.sort_values(by=[PR_INDEX_COL, "chrix"], inplace=True)

    else:
        # With y increasing upwards, reverse first-seen input order in the
        # metadata so the plotted rows read top-to-bottom like the input rows.
        order_map = {v: i for i, v in enumerate(order)}
        if len(id_col) == 1:
            rank = [
                order_map.get(v, len(order))
                for v in genesmd_df.index.get_level_values(id_col[0])
            ]
        else:
            id_values = zip(*[genesmd_df.index.get_level_values(col) for col in id_col])
            rank = [order_map.get(tuple(v), len(order)) for v in id_values]

        genesmd_df = genesmd_df.assign(__input_order_rank__=rank).sort_values(
            "__input_order_rank__", ascending=False, kind="stable"
        )
        genesmd_df.drop(columns="__input_order_rank__", inplace=True)

    genesmd_df["gene_ix_xchrom"] = genesmd_df.groupby(
        ["chrix", PR_INDEX_COL], group_keys=False, observed=True, sort=False
    ).cumcount()

    return genesmd_df


############ CHRMD_DF


##limits
def _region_specs_from_pyranges(regions):
    """Return ``(chrom, start, end)`` tuples from a PyRanges object."""
    specs = []
    for _, row in regions.iterrows():
        specs.append((row[CHROM_COL], int(row[START_COL]), int(row[END_COL])))
    return specs


def _normalize_region_specs(regions):
    """Normalize explicit ``regions=`` input into ``(chrom, start, end)`` tuples."""
    if isinstance(regions, pr.PyRanges):
        return _region_specs_from_pyranges(regions)

    if not isinstance(regions, list):
        raise Exception(
            "regions must be either a column name, a PyRanges object, or a list of (chrom, start, end) tuples/PyRanges objects."
        )

    specs = []
    for region in regions:
        if isinstance(region, pr.PyRanges):
            specs.extend(_region_specs_from_pyranges(region))
        elif isinstance(region, tuple) and len(region) == 3:
            specs.append(region)
        else:
            raise Exception(
                f"regions entries must be (chrom, start, end) tuples or PyRanges objects; got {region!r}"
            )
    return specs


def _normalize_regions_to_panels(subdf, regions):
    """Apply the dedicated ``regions=`` layout.

    When ``regions`` is provided, chromosome-grouped layout is replaced by the
    exact listed/grouped regions, in order. Synthetic chromosome identifiers are
    used internally so the existing plotting pipeline can render each region as
    a separate panel.

    Returns
    -------
    new_subdf : pandas.DataFrame
        Either the original ``subdf`` (when there is at most one window per
        chromosome) or a row-exploded copy with ``Chromosome`` relabeled to
        synthetic panel ids of the form ``f"{chrom}{PANEL_SEP}{ix}"``.
    panel_limits : dict | None
        Object to pass to :func:`chrmd_limits`, keyed by synthetic panel id.
    panel_display : dict[str, dict] | None
        Map ``synthetic_chrom -> {"chrom": real_chrom, "start": s|None, "end": e|None}``
        used by the rendering layer to build subplot titles. Returns ``None``
        when no exploding was needed (the rendering code then uses the real
        chromosome value and the dynamic min/max).
    """
    panel_limits = {}
    panel_display = {}
    new_frames = []

    if isinstance(regions, str):
        if regions not in subdf.columns:
            raise Exception(f"regions column {regions!r} is not present in the data.")
        for ix, region_name in enumerate(pd.unique(subdf[regions])):
            panel_id = f"{region_name}{PANEL_SEP}{ix}"
            block = subdf[subdf[regions] == region_name].copy()
            block[CHROM_COL] = panel_id
            new_frames.append(block)
            panel_limits[panel_id] = (None, None)
            panel_display[panel_id] = {
                "chrom": region_name,
                "start": None,
                "end": None,
            }
    else:
        for ix, (chrom, start, end) in enumerate(_normalize_region_specs(regions)):
            panel_id = f"{chrom}{PANEL_SEP}{ix}"
            block = subdf[subdf[CHROM_COL] == chrom]
            if start is not None:
                block = block[block[END_COL] > start]
            if end is not None:
                block = block[block[START_COL] < end]
            block = block.copy()
            block[CHROM_COL] = panel_id
            new_frames.append(block)
            panel_limits[panel_id] = (start, end)
            panel_display[panel_id] = {
                "chrom": chrom,
                "start": start,
                "end": end,
            }

    new_subdf = pd.concat(new_frames) if new_frames else subdf.iloc[0:0].copy()
    panel_order = list(panel_display)
    present_panel_order = [p for p in panel_order if p in set(new_subdf[CHROM_COL])]
    new_subdf[CHROM_COL] = new_subdf[CHROM_COL].astype(
        pd.CategoricalDtype(categories=present_panel_order, ordered=True)
    )
    panel_limits = {p: panel_limits[p] for p in present_panel_order}
    panel_display = {p: panel_display[p] for p in present_panel_order}
    if not panel_display:
        raise Exception("regions did not select any intervals to plot.")
    return new_subdf, panel_limits, panel_display


def _normalize_limits_to_panels(subdf, limits):
    """Keep legacy ``limits`` behavior: it customizes existing chromosome panels."""
    if isinstance(limits, dict):
        for chrom, val in limits.items():
            if val is not None and not isinstance(val, tuple):
                raise Exception(
                    f"limits[{chrom!r}] must be a (start, end) tuple or None; use regions= for multiple panels/windows."
                )
    return subdf, limits, None


def chrmd_limits(chrmd_df, limits):
    """Compute 'min_max' column for chromosome metadata"""

    # 1- create min_max column containing (plot min, plot max)

    # no limits no info
    if limits is None:
        chrmd_df["min_max"] = [(np.nan, np.nan)] * len(chrmd_df)

    # one tuple for all chromosomes
    elif type(limits) is tuple:
        chrmd_df["min_max"] = [limits] * len(chrmd_df)

    # pyranges object
    elif type(limits) is pr.PyRanges:
        # create dict to map limits
        limits_chrmd_df = limits.groupby(
            CHROM_COL, group_keys=False, observed=True
        ).agg({START_COL: "min", END_COL: "max"})
        # limits_chrmd_dict = limits_chrmd_df.to_dict(orient="index")

        # function to get matching values from limits_chrmd_df
        def make_min_max(row):
            chromosome = row.name[0]
            if chromosome in limits_chrmd_df.index:
                limits = limits_chrmd_df.loc[chromosome]

                return (
                    limits[START_COL],
                    limits[END_COL],
                )  # chromosome in both sets of data
            else:
                return (np.nan, np.nan)  # chromosome does not match

        # create limits column in plotting data
        chrmd_df["min_max"] = chrmd_df.apply(make_min_max, axis=1)

    # dictionary as limits
    else:
        chrmd_df["min_max"] = [
            limits.get(index)
            for index in list(chrmd_df.index.get_level_values(CHROM_COL))
        ]  # fills with None the chromosomes not specified


def fill_min_max(row, ts_data):
    """Complete min_max column for chromosome metadata if needed."""

    minmax_t = row["min_max"]
    # deal with empty rows
    if minmax_t is None:
        minmax_t = (np.nan, np.nan)

    # check both items and put default if necessary
    minmax_l = list(minmax_t)

    # add default to lower limit
    if minmax_l[0] is None or np.isnan(minmax_l[0]):
        minmax_l[0] = row["min"]
    # add default to higher limit
    if minmax_l[1] is None or np.isnan(minmax_l[1]):
        minmax_l[1] = row["max"]
    # consider introns off for higher limit
    else:
        if len(row) == 5:
            new_upper_lim = cumdelting([minmax_l[1]], ts_data, row.name[0], row.name[1])
            minmax_l[1] = new_upper_lim[0]

    # put plot coordinates in min_max
    row["min_max"] = minmax_l
    return row


def get_chromosome_metadata(
    df, limits, genesmd_df, pack, v_spacer, exon_height, ts_data=None
):
    """Create chromosome metadata df."""

    # Start df
    agg_funcs = {
        START_COL: "min",
        END_COL: "max",
        "__id_col_2count__": "nunique",
    }

    chrmd_df = df.groupby([CHROM_COL, PR_INDEX_COL], observed=True).agg(agg_funcs)
    chrmd_df.rename(
        columns={START_COL: "min", END_COL: "max", "__id_col_2count__": "n_genes"},
        inplace=True,
    )

    # Adjust limits in case +1 pr
    if len(df[PR_INDEX_COL].drop_duplicates()) > 1:
        chrmd_df["min"] = chrmd_df.groupby(CHROM_COL, group_keys=False, observed=True)[
            "min"
        ].transform("min")
        chrmd_df["max"] = chrmd_df.groupby(CHROM_COL, group_keys=False, observed=True)[
            "max"
        ].transform("max")

    # Add limits
    chrmd_limits(chrmd_df, limits)  # unknown limits are nan
    chrmd_df = chrmd_df.apply(lambda x: fill_min_max(x, ts_data), axis=1)

    # Store per-pr top y and order prs by visual position (top to bottom)
    pr_top_y = genesmd_df.groupby(
        [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
    )["ycoord"].max()
    chrmd_df = chrmd_df.join(pr_top_y.rename("pr_top_y"))
    if SQUISH_FACTOR_COL in df.columns:
        pr_scale = df.groupby(
            [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
        )[SQUISH_FACTOR_COL].first()
    else:
        pr_scale = pr_top_y.copy() * 0 + 1.0
    chrmd_df = chrmd_df.join(pr_scale.rename("pr_scale"))
    if THICK_COL in df.columns:
        pr_render_height = df.groupby(
            [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
        )[THICK_COL].max()
        pr_render_height = pr_render_height.where(pr_scale < 1, exon_height)
    else:
        pr_render_height = pr_top_y.copy() * 0 + exon_height
    chrmd_df = chrmd_df.join(pr_render_height.rename("pr_render_height"))
    if TEXT_EXTRA_Y_COL in genesmd_df.columns:
        text_extra_y = genesmd_df.groupby(
            [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
        )[TEXT_EXTRA_Y_COL].max()
        text_position = genesmd_df.groupby(
            [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
        )[TEXT_POSITION_COL].first()
    else:
        text_extra_y = pr_top_y.copy() * 0
        text_position = pr_top_y.copy().astype(object)
        text_position[:] = "left"
    chrmd_df = chrmd_df.join(text_extra_y.rename(TEXT_EXTRA_Y_COL))
    chrmd_df = chrmd_df.join(text_position.rename(TEXT_POSITION_COL))

    pr_bottom_y = genesmd_df.groupby(
        [CHROM_COL, PR_INDEX_COL], group_keys=False, observed=True
    )["ycoord"].min()
    top_label_extra = chrmd_df[TEXT_EXTRA_Y_COL].where(
        chrmd_df[TEXT_POSITION_COL] == "above", 0
    )
    bottom_label_extra = chrmd_df[TEXT_EXTRA_Y_COL].where(
        chrmd_df[TEXT_POSITION_COL] == "below", 0
    )
    chrmd_df["__render_top__"] = (
        chrmd_df["pr_top_y"] + 0.5 + pr_render_height / 2 + top_label_extra
    )
    chrmd_df["__render_bottom__"] = (
        pr_bottom_y + 0.5 - pr_render_height / 2 - bottom_label_extra
    )
    chrmd_df["track_gap"] = v_spacer * chrmd_df["pr_scale"]
    chrmd_df["__track_top__"] = chrmd_df["__render_top__"] + chrmd_df["track_gap"]
    chrmd_df["__track_bottom__"] = chrmd_df["__render_bottom__"] - chrmd_df["track_gap"]
    chrmd_df = (
        chrmd_df.reset_index()
        .sort_values(
            [CHROM_COL, "pr_top_y", PR_INDEX_COL], ascending=[True, False, True]
        )
        .set_index([CHROM_COL, PR_INDEX_COL])
    )

    chrmd_df_grouped = (
        chrmd_df.reset_index(level=PR_INDEX_COL)
        .groupby(CHROM_COL, group_keys=False, observed=True)
        .agg(
            {
                "min": "first",
                "max": "first",
                "min_max": "first",
                PR_INDEX_COL: ["size", list],
            }
        )
    )
    chrmd_df_grouped.columns = ["min", "max", "min_max", "n_pr_ix", "present_pr"]

    # Store plot y height
    chrmd_df_grouped = chrmd_df_grouped.join(
        genesmd_df.groupby([CHROM_COL], group_keys=False, observed=True)["ycoord"].max()
    )
    chrmd_df_grouped.rename(columns={"ycoord": "y_height"}, inplace=True)
    chrmd_df_grouped["y_height"] += (
        0.5 + exon_height / 2
    )  # the middle of the rectangle is +.5 of ycoord

    panel_edges = chrmd_df.groupby(CHROM_COL, group_keys=False, observed=True).agg(
        y_min=("__track_bottom__", "min"),
        y_max=("__track_top__", "max"),
    )
    chrmd_df_grouped = chrmd_df_grouped.join(panel_edges)
    chrmd_df_grouped["use_render_y_limits"] = True
    panel_scale = chrmd_df.groupby(CHROM_COL, observed=True)["pr_scale"].min()
    chrmd_df_grouped = chrmd_df_grouped.join(panel_scale.rename("y_pad_scale"))
    chrmd_df_grouped["y_pad"] = 0
    chrmd_df_grouped["y_height"] = chrmd_df_grouped["y_max"] - chrmd_df_grouped["y_min"]

    # Obtain the positions of lines separating pr objects. For squished tracks,
    # use the actual rendered edges of adjacent tracks so dividers stay between
    # condensed glyphs instead of crossing them.
    next_track_top = chrmd_df.groupby(CHROM_COL, observed=True)["__track_top__"].shift(
        -1
    )
    chrmd_df["pr_line"] = next_track_top
    chrmd_df["pr_line"] = chrmd_df["pr_line"].fillna(0).round(10)
    chrmd_df.drop(
        columns=[
            "pr_top_y",
            "pr_render_height",
            "pr_scale",
            TEXT_EXTRA_Y_COL,
            TEXT_POSITION_COL,
            "track_gap",
            "__render_top__",
            "__render_bottom__",
            "__track_top__",
            "__track_bottom__",
        ],
        inplace=True,
    )

    # Set chrom_ix to get the right association to the plot index
    chrmd_df_grouped["chrom_ix"] = chrmd_df_grouped.groupby(
        CHROM_COL, group_keys=False, observed=True
    ).ngroup()

    return chrmd_df, chrmd_df_grouped


def no_overlap(a, b, pad=2, pw=None):
    """Check if two intervals a and b overlap, considering a padding."""
    if pw is not None:
        if pw > 10000:
            pad = 50
        elif pw < 10000 and pw >= 1000:
            pad = 20
        elif pw < 1000 and pw >= 200:
            pad = 10
        elif pw < 200 and pw >= 50:
            pad = 5
        elif pw <= 50:
            pad = 0
    return a[1] + pad <= b[0] or a[0] >= b[1] + pad


def assign_label_rows(
    subdf,
    id_col,
    PR_INDEX_COL,
    label_pad,
    pack,
    sort_ranges,
    interval_height=None,
    v_spacer=0.3,
    plot_limits=None,
    text_label_col=None,
    text_avoid=True,
    track_scales=None,
    pack_by_track=None,
    text_position_by_track=None,
    label_size=None,
):
    """
    Assign non-overlapping ycoord rows to groups defined by (PR_INDEX_COL, id_col).

    - Does NOT reset the index of `subdf` (we preserve its index).
    - If PR_INDEX_COL or id_col are not regular columns, they are extracted
      from the index levels and used internally (no permanent index reset).
    - Returns the original `df` with updated ycoord for the rows present in subdf.

    Parameters
    ----------
    subdf : pd.DataFrame
        Subset of df to compute label rows for (same index as corresponding rows in df).
    id_col : str
        Column name with the label text (e.g. "ID").
    PR_INDEX_COL : str
        Name of index level or column that contains the pyranges index (default "__pr_ix__").
    label_pad : float
        Fractional padding (relative to plot width). Default 0.005.
    plot_limits : tuple(xmin, xmax) or None
        If provided, used to compute padding scale; otherwise taken from subdf Start/End.
    """
    s = subdf.copy()
    s[PR_INDEX_COL] = s.index.get_level_values(PR_INDEX_COL)
    s = s.reset_index(level=id_col, drop=True)
    s = s.reset_index(level=PR_INDEX_COL, drop=True)

    ycoord_map = {}
    interval_height = 0.6 if interval_height is None else interval_height
    s[TEXT_EXTRA_Y_COL] = 0.0
    s[TEXT_POSITION_COL] = "left"

    # Iterating in sorted cromosomes if sort_ranges == True else in original order
    chrom_iter = sorted(s["chrix"].unique()) if sort_ranges else pd.unique(s["chrix"])

    for chrom in chrom_iter:
        next_track_bottom = 0
        chrom_df = s[s["chrix"] == chrom]

        # Compute plot limits
        if plot_limits is None:
            xmin = float(chrom_df["Start"].min())
            xmax = float(chrom_df["End"].max())
        else:
            xmin, xmax = plot_limits

        plot_width = float(xmax - xmin) if xmax - xmin != 0 else 1.0
        pad_unit = label_pad * plot_width

        # Visual interval, all ranges occupies by groups connected by id_col
        visual_spans = (
            chrom_df.groupby([PR_INDEX_COL] + id_col, observed=True)
            .agg(
                VStart=("Start", "min"),
                VEnd=("End", "max"),
            )
            .reset_index()
        )

        pr_values = list(pd.unique(chrom_df[PR_INDEX_COL]))
        pr_iter = sorted(pr_values, reverse=True)
        for pr_val in pr_iter:
            track_scale = (
                1.0 if track_scales is None else float(track_scales.get(pr_val, 1.0))
            )
            track_pack = (
                pack if pack_by_track is None else bool(pack_by_track.get(pr_val, pack))
            )
            text_position = (
                "left"
                if text_position_by_track is None
                else text_position_by_track.get(pr_val, "left")
            )
            effective_text_avoid = text_avoid and track_scale >= 1
            # Font size is specified in points by the rendering engines, while
            # layout is expressed in pyrangeyes' row units.  The default 12 pt
            # label visually occupies about one default interval-height row
            # (0.6), so use that conversion to reserve vertical label room for
            # above/below labels without making horizontal pack artificially
            # sparse.
            text_extra_y = (
                0.0
                if not (effective_text_avoid and text_position in {"above", "below"})
                else float(label_size or 12) / 20.0
            )
            track_height = interval_height * track_scale
            track_gap = v_spacer * track_scale
            sub = chrom_df[chrom_df[PR_INDEX_COL] == pr_val]
            # In case sort_ranges is true we reorder the df by start
            if sort_ranges:
                agg = {"Start_min": ("Start", "min"), "End_max": ("End", "max")}
                if text_label_col and text_label_col in sub.columns:
                    agg[text_label_col] = (text_label_col, "first")
                gdf = (
                    sub.groupby(
                        [PR_INDEX_COL] + id_col, observed=True, sort=sort_ranges
                    )
                    .agg(**agg)
                    .reset_index()
                )
            else:
                # maintaining original order
                seen = []
                records = []

                for _, r in sub.iterrows():
                    key = (r[PR_INDEX_COL], tuple(r[id_col]))
                    if key not in seen:
                        seen.append(key)

                for pr_ix, id_vals in seen:
                    g = sub[
                        (sub[PR_INDEX_COL] == pr_ix)
                        & (sub[id_col].apply(tuple, axis=1) == id_vals)
                    ]
                    record = {
                        PR_INDEX_COL: pr_ix,
                        **{c: v for c, v in zip(id_col, id_vals)},
                        "Start_min": g["Start"].min(),
                        "End_max": g["End"].max(),
                    }
                    if text_label_col and text_label_col in g.columns:
                        record[text_label_col] = g[text_label_col].iloc[0]
                    records.append(record)

                gdf = pd.DataFrame(records)

            rows = []  # each element = row interval
            for _, g in gdf.iterrows():
                if track_pack:
                    label_len = 0
                    if effective_text_avoid and text_position in {"left", "right"}:
                        if text_label_col and text_label_col in g:
                            label_len = len(str(g[text_label_col]))
                        else:
                            label_len = len(str(tuple(g[id_col])))
                    # agafem l'interval VISUAL del grup connectat
                    vsp = visual_spans[
                        (visual_spans[PR_INDEX_COL] == g[PR_INDEX_COL])
                        & (
                            visual_spans[id_col].apply(tuple, axis=1)
                            == tuple(g[id_col])
                        )
                    ].iloc[0]

                    interval = (
                        vsp["VStart"] - label_len * pad_unit,
                        vsp["VEnd"],
                    )

                    assigned_row = None
                    for rid, row_intervals in enumerate(rows):
                        if all(
                            no_overlap(
                                interval,
                                (s0, e0),
                                pad=0 if not effective_text_avoid else 2,
                                pw=None if not effective_text_avoid else plot_width,
                            )
                            for s0, e0 in row_intervals
                        ):
                            assigned_row = rid
                            row_intervals.append(interval)
                            break

                    if assigned_row is None:
                        assigned_row = len(rows)
                        rows.append([interval])

                    key = (chrom, g[PR_INDEX_COL], tuple(g[id_col]))
                    if key not in ycoord_map:
                        ycoord = (
                            next_track_bottom
                            + track_gap
                            + track_height / 2
                            + assigned_row * (track_height + track_gap)
                            - 0.5
                        )
                        ycoord_map[key] = round(ycoord, 10)
                        s.loc[
                            (s["chrix"] == chrom)
                            & (s[PR_INDEX_COL] == g[PR_INDEX_COL])
                            & (s[id_col].apply(tuple, axis=1) == tuple(g[id_col])),
                            [TEXT_EXTRA_Y_COL, TEXT_POSITION_COL],
                        ] = [text_extra_y, text_position]

                else:
                    interval = (
                        g["Start_min"],
                        g["End_max"],
                    )
                    assigned_row = len(rows)
                    rows.append([interval])

                    key = (chrom, g[PR_INDEX_COL], tuple(g[id_col]))
                    ycoord = (
                        next_track_bottom
                        + track_gap
                        + track_height / 2
                        + assigned_row * (track_height + track_gap)
                        - 0.5
                    )
                    ycoord_map[key] = round(ycoord, 10)
                    s.loc[
                        (s["chrix"] == chrom)
                        & (s[PR_INDEX_COL] == g[PR_INDEX_COL])
                        & (s[id_col].apply(tuple, axis=1) == tuple(g[id_col])),
                        [TEXT_EXTRA_Y_COL, TEXT_POSITION_COL],
                    ] = [text_extra_y, text_position]

            next_track_bottom += (
                max(len(rows), 1) * (track_height + text_extra_y)
                + (max(len(rows), 1) + 1) * track_gap
            )

    # Assign ycoord back to all rows
    def _assign_y(r):
        key = (r["chrix"], r[PR_INDEX_COL], tuple(r[id_col]))
        return ycoord_map[key]

    s["ycoord"] = s.apply(_assign_y, axis=1)

    # restore multi-index
    s.set_index([PR_INDEX_COL], append=True, inplace=True)
    s.set_index(id_col, append=True, inplace=True, drop=False)

    return s
