import pytest
import pyranges1 as pr
import pyrangeyes as pre


def _data():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1"],
            "Start": [10, 30],
            "End": [20, 40],
            "id": ["a", "b"],
            "kind": ["exon", "CDS"],
        }
    )


def test_text_template_formats_plotly_annotations():
    pre.set_engine("ply")
    fig = pre.plot(
        _data(),
        id_col="id",
        label="{kind}:{id}",
        label_position="right",
        label_fit=False,
        label_color="red",
        label_size=13,
        return_plot="fig",
    )
    anns = [ann for ann in fig.layout.annotations if ann.text in {"exon:a", "CDS:b"}]
    assert [ann.text for ann in anns] == ["exon:a", "CDS:b"]
    assert all(ann.xanchor == "left" for ann in anns)
    assert all(ann.font.color == "red" for ann in anns)
    assert all(ann.font.size == 13 for ann in anns)


def test_label_position_center():
    pre.set_engine("ply")
    fig = pre.plot(
        _data(),
        id_col="id",
        label="{id}",
        label_position="center",
        return_plot="fig",
    )
    anns = [ann for ann in fig.layout.annotations if ann.text in {"a", "b"}]
    assert all(ann.xanchor == "center" for ann in anns)


def test_text_rejects_dict_options():
    pre.set_engine("ply")
    with pytest.raises(TypeError, match="label must be None, bool, or a format string"):
        pre.plot(_data(), id_col="id", label={"label": "{id}"}, return_plot="fig")


def test_text_uses_group_span_for_positioning():
    pre.set_engine("ply")
    grouped = pr.PyRanges(
        {
            "Chromosome": ["1", "1"],
            "Start": [10, 30],
            "End": [20, 40],
            "id": ["tx1", "tx1"],
        }
    )

    centered = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="center",
        label_pad=0,
        return_plot="fig",
    )
    ann = next(ann for ann in centered.layout.annotations if ann.text == "tx1")
    assert ann.x == 25

    right = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=0,
        return_plot="fig",
    )
    ann = next(ann for ann in right.layout.annotations if ann.text == "tx1")
    assert ann.x == 40


def test_label_pad_moves_group_labels_in_both_engines():
    grouped = pr.PyRanges(
        {
            "Chromosome": ["1", "1"],
            "Start": [10, 30],
            "End": [20, 40],
            "id": ["tx1", "tx1"],
        }
    )

    pre.set_engine("ply")
    fig0 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=0,
        return_plot="fig",
    )
    fig1 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=10,
        return_plot="fig",
    )
    x0 = next(ann.x for ann in fig0.layout.annotations if ann.text == "tx1")
    x1 = next(ann.x for ann in fig1.layout.annotations if ann.text == "tx1")
    assert x1 > x0
    assert x1 == 43

    fig0 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="above",
        label_pad=0,
        return_plot="fig",
    )
    fig1 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="above",
        label_pad=10,
        return_plot="fig",
    )
    y0 = next(ann.y for ann in fig0.layout.annotations if ann.text == "tx1")
    y1 = next(ann.y for ann in fig1.layout.annotations if ann.text == "tx1")
    assert y1 > y0

    pre.set_engine("plt")
    mpl0 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=0,
        return_plot="fig",
    )
    mpl1 = pre.plot(
        grouped,
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=10,
        return_plot="fig",
    )
    x0 = next(t.get_position()[0] for t in mpl0.axes[0].texts if t.get_text() == "tx1")
    x1 = next(t.get_position()[0] for t in mpl1.axes[0].texts if t.get_text() == "tx1")
    assert x1 > x0


def test_label_pad_percentage_is_panel_specific_for_list_inputs():
    pre.set_engine("ply")
    wide = pr.PyRanges(
        {"Chromosome": ["1"], "Start": [0], "End": [1000], "id": ["wide"]}
    )
    narrow = pr.PyRanges(
        {"Chromosome": ["1"], "Start": [0], "End": [100], "id": ["narrow"]}
    )

    fig = pre.plot(
        [wide, narrow],
        id_col="id",
        label="{id}",
        label_position="right",
        label_pad=10,
        return_plot="fig",
    )
    annotations = {ann.text: ann for ann in fig.layout.annotations}

    assert annotations["wide"].x == 1100
    assert annotations["narrow"].x == 110


def test_plotly_above_below_pad_zero_uses_full_row_height():
    pre.set_engine("ply")
    data = pr.PyRanges(
        {
            "Chromosome": ["1"],
            "Start": [0],
            "End": [10],
            "id": ["utr"],
            "height": [0.3],
        }
    )

    above = pre.plot(
        data,
        id_col="id",
        height_col="height",
        interval_height=0.8,
        label="{id}",
        label_position="above",
        label_pad=0,
        return_plot="fig",
    )
    below = pre.plot(
        data,
        id_col="id",
        height_col="height",
        interval_height=0.8,
        label="{id}",
        label_position="below",
        label_pad=0,
        return_plot="fig",
    )

    above_ann = next(ann for ann in above.layout.annotations if ann.text == "utr")
    below_ann = next(ann for ann in below.layout.annotations if ann.text == "utr")
    assert above_ann.y == pytest.approx(1.05)
    assert below_ann.y == pytest.approx(0.25)


def test_matplotlib_above_below_pad_zero_uses_full_row_height():
    pre.set_engine("plt")
    data = pr.PyRanges(
        {
            "Chromosome": ["1"],
            "Start": [0],
            "End": [10],
            "id": ["utr"],
            "height": [0.3],
        }
    )

    fig = pre.plot(
        data,
        id_col="id",
        height_col="height",
        interval_height=0.8,
        label="{id}",
        label_position="above",
        label_pad=0,
        return_plot="fig",
    )
    label = next(t for t in fig.axes[0].texts if t.get_text() == "utr")
    assert label.get_position()[1] == pytest.approx(1.05)


def test_plotly_label_angle_matches_matplotlib_direction():
    pre.set_engine("ply")
    fig = pre.plot(
        _data(),
        id_col="id",
        label="{id}",
        label_position="above",
        label_angle=25,
        return_plot="fig",
    )
    anns = [ann for ann in fig.layout.annotations if ann.text in {"a", "b"}]
    assert all(ann.textangle == -25 for ann in anns)
