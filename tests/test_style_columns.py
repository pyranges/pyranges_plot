import matplotlib
import matplotlib.colors as mcolors
import pytest
import pyranges1 as pr
import pyrangeyes as pre

matplotlib.use("Agg")


def _style_data():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1"],
            "Start": [10, 30],
            "End": [20, 40],
            "id": ["a", "b"],
            "fill": ["#ff0000", "#0000ff"],
            "outline": ["#00ff00", "#ffff00"],
            "kind": ["x", "y"],
            "status": ["ok", "warn"],
            "height": [0.5, 1.0],
            "bad_height": [0.5, 1.2],
            "text_height": ["short", "tall"],
            "depth": [2, 1],
            "bad_depth": ["front", "back"],
        }
    )


def _patches(fig):
    return [patch for patch in fig.axes[0].patches if patch.get_width() > 0]


def _hex(color):
    return mcolors.to_hex(color)


def test_matplotlib_direct_color_and_outline_columns():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="fill",
        outline_col="outline",
        colormap="direct",
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#ff0000", "#0000ff"]
    assert [_hex(p.get_edgecolor()) for p in patches] == ["#00ff00", "#ffff00"]


def test_matplotlib_channel_colormaps_for_color_and_outline():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="kind",
        outline_col="status",
        colormap={
            "color": {"x": "#111111", "y": "#222222"},
            "outline": {"ok": "#333333", "warn": "#444444"},
        },
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#111111", "#222222"]
    assert [_hex(p.get_edgecolor()) for p in patches] == ["#333333", "#444444"]


def test_matplotlib_outline_color_option_overrides_outline_col():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="fill",
        outline_col="outline",
        colormap="direct",
        outline_color="black",
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#ff0000", "#0000ff"]
    assert [_hex(p.get_edgecolor()) for p in patches] == ["#000000", "#000000"]


def test_matplotlib_height_col_scales_interval_height():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="fill",
        colormap="direct",
        height_col="height",
        interval_height=0.8,
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [round(p.get_height(), 3) for p in patches] == [0.4, 0.8]


def test_matplotlib_height_col_validates_presence():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="height_col .* not present"):
        pre.plot(
            _style_data(),
            id_col="id",
            height_col="missing_height",
            return_plot="fig",
            warnings=False,
        )


def test_matplotlib_height_col_validates_numeric_values():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="height_col .* numeric"):
        pre.plot(
            _style_data(),
            id_col="id",
            height_col="text_height",
            return_plot="fig",
            warnings=False,
        )


def test_matplotlib_height_col_validates_unit_range():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="height_col .* range from 0 to 1"):
        pre.plot(
            _style_data(),
            id_col="id",
            height_col="bad_height",
            return_plot="fig",
            warnings=False,
        )


def test_matplotlib_depth_col_validates_presence():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="depth_col .* not present"):
        pre.plot(
            _style_data(),
            id_col="id",
            depth_col="missing_depth",
            return_plot="fig",
            warnings=False,
        )


def test_matplotlib_depth_col_validates_numeric_values():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="depth_col .* numeric"):
        pre.plot(
            _style_data(),
            id_col="id",
            depth_col="bad_depth",
            return_plot="fig",
            warnings=False,
        )


def test_matplotlib_depth_col_draws_larger_values_on_top():
    pre.set_engine("matplotlib")
    data = pr.PyRanges(
        {
            "Chromosome": ["1", "1"],
            "Start": [10, 10],
            "End": [30, 30],
            "id": ["same", "same"],
            "fill": ["#ff0000", "#0000ff"],
            "depth": [2, 1],
        }
    )

    fig = pre.plot(
        data,
        id_col="id",
        color_col="fill",
        colormap="direct",
        depth_col="depth",
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#0000ff", "#ff0000"]


def test_plotly_direct_color_and_outline_columns():
    pre.set_engine("plotly")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="fill",
        outline_col="outline",
        colormap="direct",
        return_plot="fig",
        warnings=False,
    )

    box_traces = [
        trace
        for trace in fig.data
        if getattr(trace, "fill", None) == "toself" and trace.fillcolor != "white"
    ]
    assert [trace.fillcolor for trace in box_traces] == ["#ff0000", "#0000ff"]
    assert [trace.line.color for trace in box_traces] == ["#00ff00", "#ffff00"]
