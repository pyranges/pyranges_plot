"""Prototype adapters for common annotation inputs."""

import pandas as pd
import pyranges1 as pr

from .plot_main import plot

MRNA_HEIGHT_COL = "__mrna_height__"
MRNA_DEPTH_COL = "__mrna_depth__"

_CDS_FEATURES = {"CDS"}
_EXON_FEATURES = {"exon"}
_TRANSCRIPT_FEATURES = {"mRNA", "transcript"}


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


def prepare_mrna_intervals(
    annotation,
    *,
    transcript_id_col=None,
    feature_col="Feature",
    height_col=MRNA_HEIGHT_COL,
    depth_col=MRNA_DEPTH_COL,
    utr_height=0.3,
):
    """Prepare exon/CDS annotation intervals for mRNA-style plotting.

    Parameters
    ----------
    annotation : pyranges.PyRanges or pandas.DataFrame
        Freshly loaded GTF/GFF-style annotation. GTF inputs usually provide a
        ``transcript_id`` column. GFF3 inputs usually provide ``ID`` and
        ``Parent`` columns; for exon/CDS rows the first ``Parent`` value is used
        as the transcript identifier.
    transcript_id_col : str, optional
        Column identifying the transcript/mRNA. If omitted, the function uses
        ``transcript_id`` when present, otherwise derives it from GFF3
        ``Parent``/``ID`` columns.
    feature_col : str, default "Feature"
        Column containing feature types such as ``exon`` and ``CDS``.
    height_col : str, default "__mrna_height__"
        Output column with relative interval heights for ``height_col`` in
        :func:`pyrangeyes.plot`. UTR/exon intervals are set to ``utr_height``;
        CDS intervals are set to ``1``.
    depth_col : str, default "__mrna_depth__"
        Output column with draw order for ``depth_col`` in :func:`pyrangeyes.plot`.
        Exons are set to ``0`` and CDS intervals to ``1``, so CDS blocks draw on
        top of the thinner exon/UTR blocks.
    utr_height : float, default 0.3
        Relative height assigned to exon/UTR intervals. Must be in ``[0, 1]``.

    Returns
    -------
    pyranges.PyRanges
        A PyRanges object filtered to exon/CDS intervals and enriched with
        transcript ID, height, and depth columns.
    """

    if not 0 <= utr_height <= 1:
        raise ValueError("utr_height must be in the [0, 1] range.")

    if isinstance(annotation, pd.DataFrame):
        df = annotation.copy()
    elif isinstance(annotation, pr.PyRanges):
        df = annotation.copy()
    else:
        raise TypeError("annotation must be a pyranges1.PyRanges or pandas.DataFrame.")

    if feature_col not in df.columns:
        raise ValueError(f"feature_col {feature_col!r} is not present in annotation.")

    df = df[df[feature_col].isin(_EXON_FEATURES | _CDS_FEATURES)].copy()
    if df.empty:
        raise ValueError("annotation does not contain exon or CDS features to plot.")

    if transcript_id_col is None:
        transcript_id_col = _first_existing(
            df.columns, ["transcript_id", "transcript", "transcriptId"]
        )

    if transcript_id_col is None:
        if "Parent" in df.columns:
            transcript_id_col = "transcript_id"
            df[transcript_id_col] = df["Parent"].map(_split_parent)
        elif "ID" in df.columns:
            transcript_id_col = "transcript_id"
            df[transcript_id_col] = df["ID"].map(_split_parent)
        else:
            raise ValueError(
                "Could not infer a transcript identifier column. Provide "
                "transcript_id_col, or load a GTF/GFF with transcript_id or Parent."
            )
    elif transcript_id_col not in df.columns:
        raise ValueError(f"transcript_id_col {transcript_id_col!r} is not present.")

    missing_transcripts = df[transcript_id_col].isna() | (
        df[transcript_id_col].astype(str) == ""
    )
    if missing_transcripts.any() and "Parent" in df.columns:
        df.loc[missing_transcripts, transcript_id_col] = df.loc[
            missing_transcripts, "Parent"
        ].map(_split_parent)
        missing_transcripts = df[transcript_id_col].isna() | (
            df[transcript_id_col].astype(str) == ""
        )

    if missing_transcripts.any():
        raise ValueError(
            f"transcript_id_col {transcript_id_col!r} contains missing values."
        )

    df[height_col] = df[feature_col].map({"exon": utr_height, "CDS": 1.0}).astype(float)
    df[depth_col] = df[feature_col].map({"exon": 0, "CDS": 1}).astype(int)

    return pr.PyRanges(df)


def plot_mrna_annotation(
    annotation,
    *,
    transcript_id_col=None,
    feature_col="Feature",
    height_col=MRNA_HEIGHT_COL,
    depth_col=MRNA_DEPTH_COL,
    utr_height=0.3,
    color_col=None,
    tooltip=None,
    **plot_kwargs,
):
    """Plot a GTF/GFF mRNA annotation using ``height_col`` and ``depth_col``.

    This prototype adapter reproduces the visual effect of ``thick_cds=True``
    without using the ``thick_cds`` option: exon/UTR intervals are rendered as
    thin blocks, CDS intervals as full-height blocks, and CDS intervals draw on
    top of exon intervals.

    Returns whatever :func:`pyrangeyes.plot` returns. Use ``return_plot='fig'``
    to get the backend figure object.
    """

    prepared = prepare_mrna_intervals(
        annotation,
        transcript_id_col=transcript_id_col,
        feature_col=feature_col,
        height_col=height_col,
        depth_col=depth_col,
        utr_height=utr_height,
    )

    if transcript_id_col is None:
        transcript_id_col = _first_existing(
            prepared.columns, ["transcript_id", "transcript", "transcriptId"]
        )

    if color_col is None:
        color_col = transcript_id_col
    if tooltip is None:
        tooltip = "{" + feature_col + "}<br>{" + transcript_id_col + "}"

    return plot(
        prepared,
        id_col=transcript_id_col,
        color_col=color_col,
        height_col=height_col,
        depth_col=depth_col,
        tooltip=tooltip,
        **plot_kwargs,
    )
