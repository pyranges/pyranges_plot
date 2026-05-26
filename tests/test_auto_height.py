import pytest
import pyranges1 as pr
import pyrangeyes as pe


DATA = pr.PyRanges(
    {
        "Chromosome": ["1"] * 4,
        "Start": [10, 30, 12, 45],
        "End": [20, 40, 25, 55],
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "Feature": ["exon"] * 4,
    }
)


@pytest.fixture(autouse=True)
def _reset_options():
    pe.reset_options()
    pe.set_engine("plotly")
    yield
    pe.reset_options()


def test_auto_height_is_default_for_plotly_figures():
    fig = pe.plot(DATA, id_col="transcript_id", return_plot="fig")

    assert fig.layout.height != 800
    assert fig.layout.width == 1120


def test_auto_height_px_per_unit_controls_inferred_height():
    default_fig = pe.plot(DATA, id_col="transcript_id", return_plot="fig")
    taller_fig = pe.plot(
        DATA,
        id_col="transcript_id",
        auto_height_px_per_unit=88,
        return_plot="fig",
    )

    assert taller_fig.layout.height > default_fig.layout.height


def test_auto_height_px_per_unit_can_be_set_globally():
    pe.set_options("auto_height_px_per_unit", 88)
    global_fig = pe.plot(DATA, id_col="transcript_id", return_plot="fig")
    local_fig = pe.plot(
        DATA,
        id_col="transcript_id",
        auto_height_px_per_unit=88,
        return_plot="fig",
    )

    assert global_fig.layout.height == local_fig.layout.height


def test_explicit_file_size_overrides_auto_height(tmp_path):
    output = tmp_path / "plot.pdf"
    fig = pe.plot(DATA, id_col="transcript_id", to_file=(str(output), (640, 480)), return_plot="fig")

    assert fig.layout.width == 640
    assert fig.layout.height == 480


def test_auto_height_px_per_unit_must_be_positive():
    with pytest.raises(ValueError, match="auto_height_px_per_unit"):
        pe.plot(DATA, id_col="transcript_id", auto_height_px_per_unit=0, return_plot="fig")
