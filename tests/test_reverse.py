import pytest
import pyranges1 as pr
import pyrangeyes as pre


def _strand_data():
    return pr.PyRanges(
        {
            "Chromosome": ["chr1", "chr1", "chr2", "chr2"],
            "Start": [10, 40, 100, 130],
            "End": [20, 55, 115, 145],
            "Strand": ["-", "-", "+", "-"],
            "tx": ["neg", "neg", "mixed_a", "mixed_b"],
        }
    )


def test_reverse_true_mirrors_matplotlib_coordinates_and_labels_abs_ticks():
    pre.set_engine("plt")
    fig = pre.plot(_strand_data(), id_col="tx", reverse=True, return_plot="fig")
    assert fig.axes[0].get_xlim()[0] < fig.axes[0].get_xlim()[1]
    assert fig.axes[0].get_xlim()[1] < 0
    tick_labels = [label.get_text() for label in fig.axes[0].get_xticklabels()]
    assert all(not label.startswith("-") for label in tick_labels if label)


def test_reverse_auto_only_mirrors_all_negative_panels():
    pre.set_engine("ply")
    fig = pre.plot(_strand_data(), id_col="tx", reverse="auto", return_plot="fig")
    assert fig.layout.xaxis.range[0] < fig.layout.xaxis.range[1]
    assert fig.layout.xaxis.range[0] < 0
    assert fig.layout.xaxis2.range[0] < fig.layout.xaxis2.range[1]
    assert fig.layout.xaxis2.range[1] > 0
    assert all(int(label) >= 0 for label in fig.layout.xaxis.ticktext)


def test_reverse_selects_display_chromosome_for_region_panels():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        regions=[("chr1", 0, 60), ("chr2", 90, 160)],
        reverse=["chr2"],
        return_plot="fig",
    )
    assert fig.layout.xaxis.range[0] < fig.layout.xaxis.range[1]
    assert fig.layout.xaxis.range[1] > 0
    assert fig.layout.xaxis2.range[0] < fig.layout.xaxis2.range[1]
    assert fig.layout.xaxis2.range[0] < 0


def test_reverse_selects_explicit_region_tuple():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        regions=[("chr1", 0, 60), ("chr2", 90, 160)],
        reverse=[("chr1", 0, 60)],
        return_plot="fig",
    )
    assert fig.layout.xaxis.range[0] < fig.layout.xaxis.range[1]
    assert fig.layout.xaxis.range[0] < 0
    assert fig.layout.xaxis2.range[0] < fig.layout.xaxis2.range[1]
    assert fig.layout.xaxis2.range[1] > 0


def test_reverse_accepts_single_explicit_region_tuple():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        regions=[("chr1", 0, 60), ("chr2", 90, 160)],
        reverse=("chr2", 90, 160),
        return_plot="fig",
    )
    assert fig.layout.xaxis.range[0] < fig.layout.xaxis.range[1]
    assert fig.layout.xaxis.range[1] > 0
    assert fig.layout.xaxis2.range[0] < fig.layout.xaxis2.range[1]
    assert fig.layout.xaxis2.range[0] < 0


def test_reverse_rejects_unknown_selector():
    pre.set_engine("ply")
    with pytest.raises(ValueError, match="did not match any panel"):
        pre.plot(_strand_data(), id_col="tx", reverse=["missing"], return_plot="fig")


def test_reverse_custom_tooltips_use_original_coordinates():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        reverse=True,
        tooltip="coords {Start}-{End} strand {Strand}",
        return_plot="fig",
    )
    hover_text = [str(trace.text) for trace in fig.data if "coords" in str(trace.text)]
    assert any("coords 10-20 strand -" in text for text in hover_text)
    assert all("coords -" not in text for text in hover_text)


def test_reverse_default_tooltips_use_original_strand():
    pre.set_engine("ply")
    fig = pre.plot(_strand_data(), id_col="tx", reverse=True, return_plot="fig")
    hover_text = [str(trace.text) for trace in fig.data if "ID: neg" in str(trace.text)]
    assert any("[-]" in text for text in hover_text)
    assert not any("[+]" in text for text in hover_text)


def test_reverse_shrink_keeps_shrunk_region_highlights_visible_with_original_labels():
    pre.set_engine("ply")
    data = pr.PyRanges(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [10, 100, 200],
            "End": [20, 110, 210],
            "Strand": ["-", "-", "-"],
            "tx": ["neg", "neg", "neg"],
        }
    )
    fig = pre.plot(data, id_col="tx", shrink=True, reverse=True, return_plot="fig")
    highlights = [
        trace for trace in fig.data if getattr(trace, "fill", None) == "toself"
    ]
    assert highlights
    assert any("[20 - 100]" in str(trace.text) for trace in highlights)
    assert all(max(trace.x) < 0 for trace in highlights)


def test_reverse_x_ticks_accept_original_positive_coordinates():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        reverse=["chr1"],
        x_ticks={"chr1": [10, 20, 55]},
        return_plot="fig",
    )
    assert tuple(fig.layout.xaxis.tickvals) == (-10, -20, -55)
    assert tuple(fig.layout.xaxis.ticktext) == (10, 20, 55)


def test_reverse_title_chr_exposes_orientation_and_rev_flag():
    pre.set_engine("ply")
    fig = pre.plot(
        _strand_data(),
        id_col="tx",
        reverse=["chr2"],
        title_chr="{chrom} {orientation}{rev_flag}",
        return_plot="fig",
    )
    titles = [ann.text for ann in fig.layout.annotations[:2]]
    assert titles == ["chr1 fwd", "chr2 rev (rev)"]
