import pytest
import pyranges1 as pr
import pyrangeyes as pe


@pytest.fixture(autouse=True)
def _use_plotly_engine():
    pe.set_engine("plotly")


DATA = pr.PyRanges(
    {
        "Chromosome": ["1"] * 6,
        "Start": [10, 30, 12, 45, 70, 90],
        "End": [20, 40, 25, 55, 80, 100],
        "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx3", "tx3"],
        "Feature": ["exon"] * 6,
    }
)


def _filled_rectangle_heights(fig):
    heights = []
    for trace in fig.data:
        if (
            getattr(trace, "fill", None) != "toself"
            or trace.fillcolor == "white"
            or not trace.y
        ):
            continue
        y_values = [y for y in trace.y if y is not None]
        if y_values:
            heights.append(max(y_values) - min(y_values))
    return heights


def _annotation_texts(fig):
    return [ann.text for ann in fig.layout.annotations]


def _filled_rectangle_centers(fig):
    centers = []
    for trace in fig.data:
        if (
            getattr(trace, "fill", None) != "toself"
            or trace.fillcolor == "white"
            or not trace.y
        ):
            continue
        y_values = [y for y in trace.y if y is not None]
        if y_values:
            centers.append((max(y_values) + min(y_values)) / 2)
    return centers


def test_squish_bool_reduces_interval_height_and_hides_default_text():
    fig = pe.plot(
        DATA,
        id_col="transcript_id",
        squish=True,
        squish_factor=0.25,
        return_plot="fig",
    )

    assert _filled_rectangle_heights(fig)
    assert all(
        height == pytest.approx(0.6 * 0.25) for height in _filled_rectangle_heights(fig)
    )
    assert _annotation_texts(fig) == ["Chromosome 1"]


def test_squish_reduces_stacked_row_spacing():
    normal = pe.plot(DATA, id_col="transcript_id", return_plot="fig")
    squished = pe.plot(
        DATA,
        id_col="transcript_id",
        squish=True,
        squish_factor=0.25,
        return_plot="fig",
    )

    normal_centers = _filled_rectangle_centers(normal)
    squished_centers = _filled_rectangle_centers(squished)
    assert max(squished_centers) - min(squished_centers) == pytest.approx(
        (max(normal_centers) - min(normal_centers)) * 0.25
    )


def test_squish_text_true_still_hides_labels():
    fig = pe.plot(
        DATA, id_col="transcript_id", squish=True, text=True, return_plot="fig"
    )

    assert _annotation_texts(fig) == ["Chromosome 1"]


def test_squish_format_text_still_hides_labels():
    fig = pe.plot(
        DATA,
        id_col="transcript_id",
        squish=True,
        text="{transcript_id}",
        return_plot="fig",
    )

    assert _annotation_texts(fig) == ["Chromosome 1"]


def test_squish_works_with_matplotlib():
    pe.set_engine("plt")
    fig = pe.plot(
        DATA, id_col="transcript_id", squish=True, text=True, return_plot="fig"
    )

    assert sum(len(ax.patches) for ax in fig.axes) > 0
    assert not {"tx1", "tx2", "tx3"}.intersection(
        {text.get_text() for ax in fig.axes for text in ax.texts}
    )


def test_squish_list_selects_tracks():
    fig = pe.plot(
        [DATA, DATA],
        id_col="transcript_id",
        squish=[False, True],
        return_plot="fig",
    )

    # Default text labels remain for the normal track and are hidden for the squished track.
    label_texts = [text for text in _annotation_texts(fig) if text != "Chromosome 1"]
    assert sorted(label_texts) == ["tx1", "tx2", "tx3"]


def test_squish_list_must_match_tracks():
    with pytest.raises(ValueError, match="one bool per track"):
        pe.plot([DATA, DATA], id_col="transcript_id", squish=[True], return_plot="fig")


def test_track_labels_must_match_tracks():
    with pytest.raises(ValueError, match="track_labels"):
        pe.plot(
            [DATA, DATA, DATA, DATA],
            id_col="transcript_id",
            track_labels=["one", "two", "three"],
            return_plot="fig",
        )


def test_all_track_labels_render_for_squished_multitrack():
    fig = pe.plot(
        [DATA, DATA, DATA, DATA],
        id_col="transcript_id",
        track_labels=["one", "two", "three", "four"],
        squish=[False, True, False, True],
        return_plot="fig",
    )

    assert list(fig.layout.yaxis.ticktext) == ["one", "two", "three", "four"]


def test_first_track_renders_on_top():
    fig = pe.plot(
        [DATA, DATA, DATA],
        id_col="transcript_id",
        track_labels=["first", "second", "third"],
        return_plot="fig",
    )

    top_to_bottom_labels = [
        label
        for _, label in sorted(
            zip(fig.layout.yaxis.tickvals, fig.layout.yaxis.ticktext), reverse=True
        )
    ]

    assert top_to_bottom_labels == ["first", "second", "third"]


def test_packed_list_must_match_tracks():
    with pytest.raises(ValueError, match="one bool per track"):
        pe.plot([DATA, DATA], id_col="transcript_id", packed=[True], return_plot="fig")


def test_packed_list_entries_must_be_bool():
    with pytest.raises(TypeError, match="list entries"):
        pe.plot(
            [DATA, DATA],
            id_col="transcript_id",
            packed=[True, "no"],
            return_plot="fig",
        )


def test_packed_list_controls_track_layout():
    fig = pe.plot(
        [DATA, DATA],
        id_col="transcript_id",
        track_labels=["packed", "unpacked"],
        packed=[True, False],
        text=False,
        return_plot="fig",
    )

    filled_centers_by_track = {}
    for trace in fig.data:
        if getattr(trace, "fill", None) != "toself" or not trace.y:
            continue
        y_values = [float(y) for y in trace.y if y is not None]
        center = round((min(y_values) + max(y_values)) / 2, 4)
        if trace.fillcolor == "white":
            continue
        filled_centers_by_track.setdefault(center, 0)
        filled_centers_by_track[center] += 1

    divider = min(shape.y0 for shape in fig.layout.shapes if shape.y0 != 0)
    packed_centers = [center for center in filled_centers_by_track if center > divider]
    unpacked_centers = [
        center for center in filled_centers_by_track if center < divider
    ]

    assert len(packed_centers) == 2
    assert len(unpacked_centers) == 3


def test_mixed_squish_compacts_panel_height():
    unsquished_fig = pe.plot(
        [DATA, DATA, DATA, DATA],
        id_col="transcript_id",
        track_labels=["one", "two", "three", "four"],
        squish=False,
        return_plot="fig",
    )
    squished_fig = pe.plot(
        [DATA, DATA, DATA, DATA],
        id_col="transcript_id",
        track_labels=["one", "two", "three", "four"],
        squish=[False, True, False, True],
        return_plot="fig",
    )

    unsquished_span = (
        unsquished_fig.layout.yaxis.range[1] - unsquished_fig.layout.yaxis.range[0]
    )
    squished_span = (
        squished_fig.layout.yaxis.range[1] - squished_fig.layout.yaxis.range[0]
    )

    assert squished_span < unsquished_span


def test_squish_list_entries_must_be_bool():
    with pytest.raises(TypeError, match="list entries"):
        pe.plot(
            [DATA, DATA], id_col="transcript_id", squish=[True, 1], return_plot="fig"
        )


def test_squish_factor_range():
    with pytest.raises(ValueError, match="squish_factor"):
        pe.plot(
            DATA,
            id_col="transcript_id",
            squish=True,
            squish_factor=0,
            return_plot="fig",
        )
