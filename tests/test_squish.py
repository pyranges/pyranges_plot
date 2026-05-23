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
