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


def _snp_like_variants():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1"],
            "Start": [15, 32, 70],
            "End": [16, 33, 71],
            "ID": ["rs1", ".", "rs3"],
            "REF": ["A", "G", "T"],
            "ALT": ["C", "A", "G"],
        }
    )


def _patches(fig):
    return [patch for patch in fig.axes[0].patches if patch.get_width() > 0]


def test_mrna_adapter_adds_height_and_depth_columns_for_gtf():
    prepared = pre.adapters.mRNA(_gtf_like_annotation())
    df = prepared.sort_values(["Start", "End"]).reset_index(drop=True)

    assert df["__mrna_height__"].tolist() == [0.3, 1.0, 1.0, 0.3]
    assert df["__mrna_depth__"].tolist() == [0, 1, 1, 0]
    assert df["transcript_id"].tolist() == ["tx1"] * 4


def test_mrna_adapter_derives_transcript_id_from_gff_parent():
    prepared = pre.adapters.mRNA(_gff_like_annotation())

    assert "transcript_id" in prepared.columns
    assert prepared["transcript_id"].tolist() == ["tx1", "tx1", "tx1"]


def test_plot_calls_mrna_adapter_before_rendering():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        pre.Track(_gtf_like_annotation(), "mRNA"),
        fill_col="Feature",
        colormap={"exon": "lightgrey", "CDS": "gold"},
        return_plot="fig",
        warnings=False,
    )

    heights = sorted(round(patch.get_height(), 3) for patch in _patches(fig))
    assert heights == [0.18, 0.18, 0.6, 0.6]


def test_plot_passes_adapter_specific_arguments_and_rejects_unknown_arguments():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        pre.Track(_gtf_like_annotation(), "mRNA", utr_height=0.5),
        return_plot="fig",
        warnings=False,
    )
    heights = sorted(round(patch.get_height(), 3) for patch in _patches(fig))
    assert heights == [0.3, 0.3, 0.6, 0.6]

    with pytest.raises(Exception, match="Unknown Track option"):
        pre.plot(
            pre.Track(_gtf_like_annotation(), "mRNA", not_an_adapter_arg=True),
            return_plot="fig",
            warnings=False,
        )


def test_adapter_options_are_printable_settable_and_used_by_plot(capsys):
    pre.set_engine("matplotlib")
    pre.reset_options(adapter="mRNA")

    assert pre.print_options(adapter="mRNA", return_keys=True) == {
        "id_col",
        "feature_col",
        "height_col",
        "depth_col",
        "utr_height",
    }
    pre.print_options(adapter="mRNA")
    printed = capsys.readouterr().out
    assert "utr_height" in printed
    assert "<infer>" in printed

    pre.set_options(adapter="mRNA", variable="utr_height", value=0.5)
    assert pre.get_options("utr_height", adapter="mRNA") == 0.5

    prepared = pre.adapters.mRNA(_gtf_like_annotation())
    assert sorted(prepared["__mrna_height__"].unique()) == [0.5, 1.0]

    fig = pre.plot(
        pre.Track(_gtf_like_annotation(), "mRNA"),
        return_plot="fig",
        warnings=False,
    )
    heights = sorted(round(patch.get_height(), 3) for patch in _patches(fig))
    assert heights == [0.3, 0.3, 0.6, 0.6]

    pre.reset_options(adapter="mRNA")


def test_adapter_describe_lists_available_adapters(capsys):
    descriptions = pre.adapters.describe()

    assert "mRNA" in descriptions
    assert "SNP" in descriptions
    assert "CDS" in descriptions["mRNA"]
    assert pre.adapters.describe("mRNA") == descriptions["mRNA"]

    adapters_repr = repr(pre.adapters)
    assert "Available pyrangeyes adapters" in adapters_repr
    assert "pe.print_options('mRNA')" in adapters_repr

    pre.print_options("mRNA")
    printed = capsys.readouterr().out
    assert "Adapter option" in printed
    assert "utr_height" in printed


def test_snp_adapter_adds_fixed_size_shape_and_generated_ids():
    prepared = pre.adapters.SNP(_snp_like_variants(), shape="triangle-up", size=20)

    assert prepared["__adapter_shape__"].tolist() == ["triangle-up"] * 3
    assert prepared["__adapter_height__"].tolist() == [0.8, 0.8, 0.8]
    assert prepared["__adapter_marker_size__"].tolist() == [20, 20, 20]
    assert prepared["End"].sub(prepared["Start"]).tolist() == [1, 1, 1]


def test_plot_accepts_one_adapter_per_input_object():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        [
            pre.Track(_gtf_like_annotation(), "mRNA"),
            pre.Track(_snp_like_variants(), "SNP"),
        ],
        return_plot="fig",
        warnings=False,
    )

    assert fig.axes


def test_plot_rejects_adapter_list_length_mismatch():
    with pytest.raises(Exception, match="adapter"):
        pre.plot(_gtf_like_annotation(), adapter="mRNA")


def test_mrna_adapter_rejects_missing_transcript_identifier():
    annotation = _gff_like_annotation().drop(columns=["Parent", "ID"])

    with pytest.raises(ValueError, match="Could not infer"):
        pre.adapters.mRNA(annotation)
