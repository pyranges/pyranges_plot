import pytest
import pyranges1 as pr
import pyrangeyes as pe

import matplotlib

matplotlib.use("Agg")


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
        auto_height_px_per_unit=120,
        return_plot="fig",
    )

    assert taller_fig.layout.height > default_fig.layout.height


def test_auto_height_px_per_unit_can_be_set_globally():
    pe.set_options("auto_height_px_per_unit", 120)
    global_fig = pe.plot(DATA, id_col="transcript_id", return_plot="fig")
    local_fig = pe.plot(
        DATA,
        id_col="transcript_id",
        auto_height_px_per_unit=120,
        return_plot="fig",
    )

    assert global_fig.layout.height == local_fig.layout.height


def test_auto_size_is_height_only_even_with_legend():
    fig = pe.plot(DATA, id_col="transcript_id", legend=True, return_plot="fig")

    assert fig.layout.width == 1120


def test_explicit_file_size_overrides_auto_height(tmp_path):
    output = tmp_path / "plot.pdf"
    fig = pe.plot(
        DATA,
        id_col="transcript_id",
        to_file=(str(output), (640, 480)),
        return_plot="fig",
    )

    assert fig.layout.width == 640
    assert fig.layout.height == 480


def test_auto_height_keeps_matplotlib_interval_pixel_height_across_tracks():
    pe.set_engine("matplotlib")
    data = pe.example_data.p1
    regions = pr.PyRanges(
        {
            "Chromosome": [1, 1, 2, 2],
            "Start": [15, 75, 35, 130],
            "End": [30, 95, 55, 165],
            "group": ["g1", "g2", "g3", "g4"],
        }
    )

    single = pe.plot(data, id_col="transcript_id", return_plot="fig")
    multi = pe.plot(
        [pe.Track(data), pe.Track(regions, id_col="group", pack=False)],
        id_col="transcript_id",
        return_plot="fig",
    )

    def interval_px(fig):
        fig.canvas.draw()
        ax = fig.axes[0]
        y_min, y_max = ax.get_ylim()
        px_per_unit = ax.get_window_extent().height / abs(y_max - y_min)
        return px_per_unit * pe.get_options("interval_height")

    assert interval_px(multi) == pytest.approx(interval_px(single), rel=0.02)


def test_auto_height_reserves_matplotlib_title_space():
    pe.set_engine("matplotlib")

    fig = pe.plot(DATA, id_col="transcript_id", warnings=False, return_plot="fig")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for ax in fig.axes:
        title_box = ax.title.get_window_extent(renderer)
        axes_box = ax.get_window_extent(renderer)
        assert title_box.y0 >= axes_box.y1 + 6


def test_disabled_track_labels_do_not_reserve_above_label_band():
    unlabeled = pe.plot(
        [
            pe.Track(DATA, id_col="transcript_id", label="{transcript_id}"),
            pe.Track(DATA, id_col="transcript_id", pack=False, label=False),
        ],
        id_col="transcript_id",
        return_plot="fig",
    )
    labeled = pe.plot(
        [
            pe.Track(DATA, id_col="transcript_id", label="{transcript_id}"),
            pe.Track(DATA, id_col="transcript_id", pack=False, label="{transcript_id}"),
        ],
        id_col="transcript_id",
        return_plot="fig",
    )

    assert labeled.layout.height > unlabeled.layout.height


def test_auto_height_px_per_unit_must_be_positive():
    with pytest.raises(ValueError, match="auto_height_px_per_unit"):
        pe.plot(
            DATA, id_col="transcript_id", auto_height_px_per_unit=0, return_plot="fig"
        )
