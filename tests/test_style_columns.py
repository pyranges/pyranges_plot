import matplotlib
import matplotlib.colors as mcolors
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
