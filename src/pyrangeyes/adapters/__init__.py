"""Adapters that prepare domain-specific inputs for :func:`pyrangeyes.plot`."""

from __future__ import annotations

import pandas as pd
import pyranges1 as pr


class _DefaultOption:
    """Sentinel meaning "read this adapter argument from get_options()"."""

    def __repr__(self):
        return "DEFAULT"


DEFAULT = _DefaultOption()

_CDS_FEATURES = {"CDS"}
_EXON_FEATURES = {"exon"}


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

    The adapter returns a PyRanges object enriched with columns suitable for
    ``plot(..., height_col=..., depth_col=...)``. It reproduces the visual
    effect of ``thick_cds=True`` without using that plotting option: exon/UTR
    intervals receive a thinner relative height, CDS intervals receive full
    height, and CDS intervals draw on top.

    Parameters
    ----------
    annotation : pyranges.PyRanges or pandas.DataFrame
        Freshly loaded GTF/GFF-style annotation. GTF inputs usually provide a
        ``transcript_id`` column. GFF3 inputs usually provide ``ID`` and
        ``Parent`` columns; for exon/CDS rows, the first ``Parent`` value is
        used as the transcript identifier when needed.
    id_col : str, None, or DEFAULT, default DEFAULT
        Transcript/mRNA identifier column. ``DEFAULT`` reads the current
        adapter option from ``get_options("mRNA", "id_col")``. The shipped
        option is ``None``, which means ``transcript_id`` (or similar) is used
        when present, otherwise it is derived from GFF3 ``Parent``/``ID``
        columns and stored as ``transcript_id``.
    feature_col : str or DEFAULT, default DEFAULT
        Column containing feature types such as ``exon`` and ``CDS``.
        ``DEFAULT`` reads the current ``feature_col`` adapter option.
    height_col : str or DEFAULT, default DEFAULT
        Output relative-height column for ``plot(height_col=...)``.
        ``DEFAULT`` reads the current ``height_col`` adapter option.
    depth_col : str or DEFAULT, default DEFAULT
        Output draw-order column for ``plot(depth_col=...)``. Exons are ``0``;
        CDS intervals are ``1`` and draw later/on top. ``DEFAULT`` reads the
        current ``depth_col`` adapter option.
    utr_height : float or DEFAULT, default DEFAULT
        Relative height assigned to exon/UTR intervals. Must be in ``[0, 1]``.
        ``DEFAULT`` reads the current ``utr_height`` adapter option.

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

    return pr.PyRanges(df)


_ADAPTERS = {"mRNA": mRNA}
_ADAPTER_DESCRIPTIONS = {
    "mRNA": "Pre-configured mRNA/GTF/GFF visualization: exon/UTR intervals are thin and CDS intervals are thick."
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
            "Output relative-height column used as plot(height_col=...).",
            " ",
        ),
        "depth_col": (
            "__mrna_depth__",
            "Output draw-order column used as plot(depth_col=...). CDS draws on top of exon.",
            " ",
        ),
        "utr_height": (
            0.3,
            "Relative height assigned to exon/UTR intervals; must range from 0 to 1.",
            " ",
        ),
    }
}
_ADAPTER_OPTIONS_IN_USE = {
    adapter_name: dict(options) for adapter_name, options in _ADAPTER_OPTIONS.items()
}
_DEFAULT_PLOT_ARGS_FROM_OPTIONS = {
    "mRNA": {
        "id_col": "id_col",
        "height_col": "height_col",
        "depth_col": "depth_col",
    }
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


def get_options(name, varname="all"):
    """Return adapter options currently in use.

    ``varname="all"`` returns the full ``{name: (value, description,
    modified_tag)}`` mapping. ``varname="values"`` returns only option values
    as a dict. A single option name returns that option value, and a list of
    names returns the corresponding values in order.

    Adapter functions use the public ``DEFAULT`` sentinel in their signatures;
    when an argument is left as ``DEFAULT`` at runtime, the corresponding value
    is fetched through this function. This means ``set_options`` affects both
    ``plot(..., adapter=...)`` and direct adapter calls such as ``mRNA(...)``.
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
    "get",
    "names",
    "describe",
    "get_options",
    "set_options",
    "reset_options",
    "default_plot_args",
    "accepted_kwargs",
]
