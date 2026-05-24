"""Adapters that prepare domain-specific inputs for :func:`pyrangeyes.plot`."""

from __future__ import annotations

import sys
import types

import pandas as pd
import pyranges1 as pr


class _DefaultOption:
    """Sentinel meaning "read this adapter argument from get_options()"."""

    def __repr__(self):
        return "DEFAULT"


DEFAULT = _DefaultOption()

_CDS_FEATURES = {"CDS"}
_EXON_FEATURES = {"exon"}
ADAPTER_ID_COL = "__adapter_id__"
ADAPTER_HEIGHT_COL = "__adapter_height__"
ADAPTER_DEPTH_COL = "__adapter_depth__"
ADAPTER_SHAPE_COL = "__adapter_shape__"
ADAPTER_MARKER_SIZE_COL = "__adapter_marker_size__"


def _split_parent(value):
    """Return the first parent ID from a GFF3 Parent-like value."""

    if pd.isna(value):
        return value
    value = str(value)
    if not value:
        return value
    return value.split(",", maxsplit=1)[0]


def _first_existing(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def _as_dataframe(annotation):
    if isinstance(annotation, pd.DataFrame):
        return annotation.copy()
    if isinstance(annotation, pr.PyRanges):
        return annotation.copy()
    raise TypeError("annotation must be a pyranges1.PyRanges or pandas.DataFrame.")


def mRNA(
    annotation,
    *,
    id_col=DEFAULT,
    feature_col=DEFAULT,
    height_col=DEFAULT,
    depth_col=DEFAULT,
    utr_height=DEFAULT,
    **_,
):
    """Prepare GTF/GFF mRNA annotations for plotting.

    The adapter returns a PyRanges object ready for ``plot(..., "mRNA")``:
    exon/UTR intervals receive a thinner relative height, CDS intervals receive
    full height, and CDS intervals draw on top.

    Arguments left as ``DEFAULT`` are read from the current adapter options at
    call time. Inspect them with ``pe.get_options(adapter="mRNA", varname="values")``
    or ``pe.print_options(adapter="mRNA")``.

    Parameters
    ----------
    annotation : pyranges.PyRanges or pandas.DataFrame
        Freshly loaded GTF/GFF-style annotation. GTF inputs usually provide a
        ``transcript_id`` column. GFF3 inputs usually provide ``ID`` and
        ``Parent`` columns; for exon/CDS rows, the first ``Parent`` value is
        used as the transcript identifier when needed.
    id_col : str, None, or DEFAULT, default DEFAULT
        Transcript/mRNA identifier column. ``None`` means ``transcript_id`` (or
        similar) is used when present, otherwise it is derived from GFF3
        ``Parent``/``ID`` columns and stored as ``transcript_id``.
    feature_col : str or DEFAULT, default DEFAULT
        Column containing feature types such as ``exon`` and ``CDS``.
    height_col : str or DEFAULT, default DEFAULT
        Output relative-height column.
    depth_col : str or DEFAULT, default DEFAULT
        Output draw-order column. Exons are ``0``; CDS intervals are ``1`` and
        draw later/on top.
    utr_height : float or DEFAULT, default DEFAULT
        Relative height assigned to exon/UTR intervals. Must be in ``[0, 1]``.

    Returns
    -------
    pyranges.PyRanges
        Exon/CDS intervals enriched with ID, height, and depth columns.
    """

    id_col, feature_col, height_col, depth_col, utr_height = _resolve_defaults(
        "mRNA",
        {
            "id_col": id_col,
            "feature_col": feature_col,
            "height_col": height_col,
            "depth_col": depth_col,
            "utr_height": utr_height,
        },
    ).values()

    if not 0 <= utr_height <= 1:
        raise ValueError("utr_height must be in the [0, 1] range.")

    df = _as_dataframe(annotation)

    if feature_col not in df.columns:
        raise ValueError(f"feature_col {feature_col!r} is not present in annotation.")

    df = df[df[feature_col].isin(_EXON_FEATURES | _CDS_FEATURES)].copy()
    if df.empty:
        raise ValueError("annotation does not contain exon or CDS features to plot.")

    resolved_id_col = id_col or _first_existing(
        df.columns, ["transcript_id", "transcript", "transcriptId"]
    )

    if resolved_id_col is None:
        if "Parent" in df.columns:
            resolved_id_col = "transcript_id"
            df[resolved_id_col] = df["Parent"].map(_split_parent)
        elif "ID" in df.columns:
            resolved_id_col = "transcript_id"
            df[resolved_id_col] = df["ID"].map(_split_parent)
        else:
            raise ValueError(
                "Could not infer a transcript identifier column. Provide "
                "id_col, or load a GTF/GFF with transcript_id or Parent."
            )
    elif resolved_id_col not in df.columns:
        raise ValueError(f"id_col {resolved_id_col!r} is not present in annotation.")

    missing_ids = df[resolved_id_col].isna() | (df[resolved_id_col].astype(str) == "")
    if missing_ids.any() and "Parent" in df.columns:
        df.loc[missing_ids, resolved_id_col] = df.loc[missing_ids, "Parent"].map(
            _split_parent
        )
        missing_ids = df[resolved_id_col].isna() | (
            df[resolved_id_col].astype(str) == ""
        )

    if missing_ids.any():
        raise ValueError(f"id_col {resolved_id_col!r} contains missing values.")

    df[height_col] = df[feature_col].map({"exon": utr_height, "CDS": 1.0}).astype(float)
    df[depth_col] = df[feature_col].map({"exon": 0, "CDS": 1}).astype(int)
    df[ADAPTER_ID_COL] = df[resolved_id_col]
    df[ADAPTER_HEIGHT_COL] = df[height_col]
    df[ADAPTER_DEPTH_COL] = df[depth_col]
    df[ADAPTER_SHAPE_COL] = "rectangle"

    return pr.PyRanges(df)


def SNP(
    variants,
    *,
    id_col=DEFAULT,
    ref_col=DEFAULT,
    alt_col=DEFAULT,
    height=DEFAULT,
    shape=DEFAULT,
    size=DEFAULT,
    shape_col=DEFAULT,
    **_,
):
    """Prepare SNP/VCF-like variants for plotting as fixed-size markers.

    Arguments left as ``DEFAULT`` are read from the current adapter options at
    call time. Inspect them with ``pe.get_options(adapter="SNP", varname="values")``
    or ``pe.print_options(adapter="SNP")``.

    Parameters
    ----------
    variants : pyranges.PyRanges or pandas.DataFrame
        Variant intervals, commonly from VCF-like data with ``REF`` and ``ALT``
        columns.
    id_col : str, None, or DEFAULT, default DEFAULT
        Variant identifier column. If no usable ID is available, a readable
        ``chrom:start REF>ALT`` identifier is created.
    ref_col, alt_col : str or DEFAULT, default DEFAULT
        Reference and alternate allele columns used for generated IDs and
        tooltips when present.
    height : float or DEFAULT, default DEFAULT
        Relative marker row height; must range from 0 to 1.
    shape : {"diamond", "triangle-up", "triangle-down", "circle"} or DEFAULT
        Marker shape.
    size : int or float or DEFAULT, default DEFAULT
        Marker size in screen points. Width and height are kept equal in the
        rendered visualization.
    shape_col : str or DEFAULT, default DEFAULT
        Output shape column used internally by ``plot``.

    Returns
    -------
    pyranges.PyRanges
        Variant intervals enriched with ID, height, and diamond-shape columns.
    """

    id_col, ref_col, alt_col, height, shape, size, shape_col = _resolve_defaults(
        "SNP",
        {
            "id_col": id_col,
            "ref_col": ref_col,
            "alt_col": alt_col,
            "height": height,
            "shape": shape,
            "size": size,
            "shape_col": shape_col,
        },
    ).values()

    if not 0 <= height <= 1:
        raise ValueError("height must be in the [0, 1] range.")
    if size <= 0:
        raise ValueError("size must be positive.")
    accepted_shapes = {"diamond", "triangle-up", "triangle-down", "circle"}
    if shape not in accepted_shapes:
        raise ValueError(
            f"shape must be one of {sorted(accepted_shapes)}; got {shape!r}."
        )

    df = _as_dataframe(variants)
    lengths = df["End"].astype(int) - df["Start"].astype(int)
    if not (lengths == 1).all():
        raise ValueError("SNP adapter expects single-position intervals of length 1.")

    resolved_id_col = id_col or _first_existing(
        df.columns, ["ID", "Name", "variant_id"]
    )
    has_usable_id = (
        resolved_id_col in df.columns
        and not df[resolved_id_col].isna().all()
        and not (df[resolved_id_col].astype(str) == ".").all()
    )
    if not has_usable_id:
        resolved_id_col = "variant_id"
        ref_values = df[ref_col].astype(str) if ref_col in df.columns else "?"
        alt_values = df[alt_col].astype(str) if alt_col in df.columns else "?"
        df[resolved_id_col] = (
            df["Chromosome"].astype(str)
            + ":"
            + df["Start"].astype(str)
            + " "
            + ref_values
            + ">"
            + alt_values
        )

    df[ADAPTER_ID_COL] = df[resolved_id_col]
    df[ADAPTER_HEIGHT_COL] = height
    df[ADAPTER_DEPTH_COL] = 0
    df[shape_col] = shape
    df[ADAPTER_SHAPE_COL] = df[shape_col]
    df[ADAPTER_MARKER_SIZE_COL] = size

    return pr.PyRanges(df)


_ADAPTERS = {"mRNA": mRNA, "SNP": SNP}
_ADAPTER_DESCRIPTIONS = {
    "mRNA": "Pre-configured mRNA/GTF/GFF visualization: exon/UTR intervals are thin and CDS intervals are thick.",
    "SNP": "Pre-configured SNP/VCF-like visualization: single-position variants are drawn as fixed-size markers.",
}
_ADAPTER_OPTIONS = {
    "mRNA": {
        "id_col": (
            None,
            "Transcript/mRNA identifier column. None means infer from transcript_id or GFF3 Parent/ID.",
            " ",
        ),
        "feature_col": (
            "Feature",
            "Column containing feature types such as exon and CDS.",
            " ",
        ),
        "height_col": (
            "__mrna_height__",
            "Output relative-height column used by the adapter.",
            " ",
        ),
        "depth_col": (
            "__mrna_depth__",
            "Output draw-order column used by the adapter. CDS draws on top of exon.",
            " ",
        ),
        "utr_height": (
            0.3,
            "Relative height assigned to exon/UTR intervals; must range from 0 to 1.",
            " ",
        ),
    },
    "SNP": {
        "id_col": (
            None,
            "Variant identifier column. None means infer from ID/Name/variant_id, or create one from position and alleles.",
            " ",
        ),
        "ref_col": (
            "REF",
            "Reference allele column used for generated IDs/tooltips.",
            " ",
        ),
        "alt_col": (
            "ALT",
            "Alternate allele column used for generated IDs/tooltips.",
            " ",
        ),
        "height": (0.8, "Relative marker row height; must range from 0 to 1.", " "),
        "shape": (
            "diamond",
            "Marker shape: diamond, triangle-up, triangle-down, or circle.",
            " ",
        ),
        "size": (18, "Marker size in screen points; width and height are equal.", " "),
        "shape_col": (ADAPTER_SHAPE_COL, "Output shape column used by plot().", " "),
    },
}
_ADAPTER_OPTIONS_IN_USE = {
    adapter_name: dict(options) for adapter_name, options in _ADAPTER_OPTIONS.items()
}
_DEFAULT_PLOT_ARGS_FROM_OPTIONS = {
    "mRNA": {
        "id_col": ADAPTER_ID_COL,
        "height_col": ADAPTER_HEIGHT_COL,
        "depth_col": ADAPTER_DEPTH_COL,
        "shape_col": ADAPTER_SHAPE_COL,
    },
    "SNP": {
        "id_col": ADAPTER_ID_COL,
        "height_col": ADAPTER_HEIGHT_COL,
        "depth_col": ADAPTER_DEPTH_COL,
        "shape_col": ADAPTER_SHAPE_COL,
    },
}


def _resolve_defaults(name, values):
    """Replace DEFAULT sentinel values with current adapter option values."""

    option_values = get_options(name, "values")
    return {
        key: option_values[key] if value is DEFAULT else value
        for key, value in values.items()
    }


def get(name):
    """Return an adapter function by name."""

    try:
        return _ADAPTERS[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown adapter {name!r}. Available adapters are: {sorted(_ADAPTERS)}"
        ) from exc


def names():
    """Return available adapter names."""

    return sorted(_ADAPTERS)


def describe(name=None):
    """Return a short description of available adapters.

    Parameters
    ----------
    name : str, optional
        Adapter name. If omitted, descriptions for all adapters are returned.
    """

    if name is not None:
        get(name)  # validates name
        return _ADAPTER_DESCRIPTIONS[name]
    return {
        adapter_name: _ADAPTER_DESCRIPTIONS[adapter_name] for adapter_name in names()
    }


def _format_description(name=None):
    descriptions = describe(name)
    if isinstance(descriptions, str):
        descriptions = {name: descriptions}
    lines = ["Available pyrangeyes adapters:"]
    for adapter_name, description in descriptions.items():
        lines.append(f"- {adapter_name}: {description}")
    lines.append("Use pe.print_options('mRNA') to inspect adapter options.")
    return "\n".join(lines)


class _AdaptersModule(types.ModuleType):
    def __repr__(self):
        return _format_description()

    __str__ = __repr__


def get_options(name, varname="all"):
    """Return adapter options currently in use.

    ``varname="all"`` returns the full ``{name: (value, description,
    modified_tag)}`` mapping. ``varname="values"`` returns only option values
    as a dict. A single option name returns that option value, and a list of
    names returns the corresponding values in order.

    Adapter functions use the public ``DEFAULT`` sentinel in their signatures;
    when an argument is left as ``DEFAULT`` at runtime, the corresponding value
    is fetched through this function. This means ``set_options`` affects both
    ``Track(data, adapter=...)`` and direct adapter calls such as ``mRNA(...)``.
    """

    get(name)  # validates name
    options = _ADAPTER_OPTIONS_IN_USE[name]
    if varname == "all":
        return options
    if varname == "values":
        return {key: val[0] for key, val in options.items()}
    if isinstance(varname, list):
        return [options[var][0] for var in varname]
    if varname in options:
        return options[varname][0]
    raise ValueError(
        f"The variable {varname!r} is not customizable for adapter {name!r}. "
        f"Customizable variables are: {list(options)}"
    )


def set_options(name, variable, value=None):
    """Set one or more adapter options.

    ``variable`` may be a single option name plus ``value``, or a dictionary of
    option/value pairs. Values set here are read by adapter arguments left as
    ``DEFAULT``.
    """

    get(name)  # validates name
    if isinstance(variable, str):
        variable = {variable: value}
    for key, val in variable.items():
        if key not in _ADAPTER_OPTIONS_IN_USE[name]:
            raise ValueError(
                f"The variable {key!r} is not customizable for adapter {name!r}. "
                f"Customizable variables are: {list(_ADAPTER_OPTIONS_IN_USE[name])}"
            )
        mod_tag = "*" if val != _ADAPTER_OPTIONS[name][key][0] else " "
        _ADAPTER_OPTIONS_IN_USE[name][key] = (
            val,
            _ADAPTER_OPTIONS[name][key][1],
            mod_tag,
        )


def reset_options(name, variable="all"):
    """Reset one, some, or all options for an adapter to shipped defaults."""

    get(name)  # validates name
    if variable == "all":
        _ADAPTER_OPTIONS_IN_USE[name] = dict(_ADAPTER_OPTIONS[name])
    elif isinstance(variable, list):
        for key in variable:
            _ADAPTER_OPTIONS_IN_USE[name][key] = _ADAPTER_OPTIONS[name][key]
    else:
        _ADAPTER_OPTIONS_IN_USE[name][variable] = _ADAPTER_OPTIONS[name][variable]


def default_plot_args(name, adapter_kwargs=None):
    """Return plot arguments an adapter prepares by default."""

    get(name)  # validates name
    adapter_kwargs = adapter_kwargs or {}
    adapter_values = get_options(name, "values") | adapter_kwargs
    defaults = {}
    for plot_arg, adapter_value in _DEFAULT_PLOT_ARGS_FROM_OPTIONS.get(
        name, {}
    ).items():
        if adapter_value in adapter_values:
            defaults[plot_arg] = adapter_values[adapter_value] or "transcript_id"
        else:
            defaults[plot_arg] = adapter_value
    return defaults


def accepted_kwargs(name):
    """Return explicit keyword parameters accepted by an adapter."""

    get(name)  # validates name
    return set(_ADAPTER_OPTIONS[name])


__all__ = [
    "DEFAULT",
    "mRNA",
    "SNP",
    "get",
    "names",
    "describe",
    "get_options",
    "set_options",
    "reset_options",
    "default_plot_args",
    "accepted_kwargs",
]


sys.modules[__name__].__class__ = _AdaptersModule
