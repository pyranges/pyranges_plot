import matplotlib
import pytest
import pyranges1 as pr
import pyrangeyes as pre

matplotlib.use("Agg")


def _gtf_like_annotation():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1", "1"],
            "Start": [10, 20, 25, 50],
            "End": [40, 30, 35, 80],
            "Strand": ["+", "+", "+", "+"],
            "Feature": ["exon", "CDS", "CDS", "exon"],
            "transcript_id": ["tx1", "tx1", "tx1", "tx1"],
        }
    )


def _gff_like_annotation():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1"],
            "Start": [10, 20, 50],
            "End": [40, 30, 80],
            "Strand": ["+", "+", "+"],
            "Feature": ["exon", "CDS", "exon"],
            "ID": ["exon1", "cds1", "exon2"],
            "Parent": ["tx1", "tx1", "tx1"],
        }
    )


def _patches(fig):
    return [patch for patch in fig.axes[0].patches if patch.get_width() > 0]


def test_prepare_mrna_intervals_adds_height_and_depth_columns_for_gtf():
    prepared = pre.prepare_mrna_intervals(_gtf_like_annotation())
    df = prepared.sort_values(["Start", "End"]).reset_index(drop=True)

    assert df["__mrna_height__"].tolist() == [0.3, 1.0, 1.0, 0.3]
    assert df["__mrna_depth__"].tolist() == [0, 1, 1, 0]
    assert df["transcript_id"].tolist() == ["tx1"] * 4


def test_prepare_mrna_intervals_derives_transcript_id_from_gff_parent():
    prepared = pre.prepare_mrna_intervals(_gff_like_annotation())

    assert "transcript_id" in prepared.columns
    assert prepared["transcript_id"].tolist() == ["tx1", "tx1", "tx1"]


def test_plot_mrna_annotation_uses_height_and_depth_columns_matplotlib():
    pre.set_engine("matplotlib")

    fig = pre.plot_mrna_annotation(
        _gtf_like_annotation(),
        color_col="Feature",
        colormap={"exon": "lightgrey", "CDS": "gold"},
        return_plot="fig",
        warnings=False,
    )

    heights = sorted(round(patch.get_height(), 3) for patch in _patches(fig))
    assert heights == [0.18, 0.18, 0.6, 0.6]


def test_plot_mrna_annotation_rejects_missing_transcript_identifier():
    annotation = _gff_like_annotation().drop(columns=["Parent", "ID"])

    with pytest.raises(ValueError, match="Could not infer"):
        pre.prepare_mrna_intervals(annotation)
