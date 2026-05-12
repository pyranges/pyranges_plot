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
            "kind": ["x", "y"],
            "status": ["ok", "warn"],
            "score": [0.0, 1.0],
            "outline_score": [1.0, 0.0],
            "bad_score": ["low", "high"],
        }
    )


def _patches(fig):
    return [patch for patch in fig.axes[0].patches if patch.get_width() > 0]


def _hex(color):
    return mcolors.to_hex(color)


def test_matplotlib_quantitative_colormap_auto_range():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="score",
        colormap={"type": "quantitative", "colors": ["blue", "red"]},
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#0000ff", "#ff0000"]


def test_matplotlib_quantitative_colormap_manual_range():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="score",
        colormap={
            "type": "quantitative",
            "colors": ["blue", "white", "red"],
            "range": (0, 2),
        },
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#0000ff", "#fffefe"]


def test_matplotlib_quantitative_colormap_for_fill_and_outline():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="score",
        outline_col="outline_score",
        colormap={
            "color": {"type": "quantitative", "colors": ["blue", "red"]},
            "outline": {"type": "quantitative", "colors": ["black", "white"]},
        },
        return_plot="fig",
        warnings=False,
    )

    patches = _patches(fig)
    assert [_hex(p.get_facecolor()) for p in patches] == ["#0000ff", "#ff0000"]
    assert [_hex(p.get_edgecolor()) for p in patches] == ["#ffffff", "#000000"]


def test_plotly_quantitative_colormap_auto_range():
    pre.set_engine("plotly")

    fig = pre.plot(
        _style_data(),
        id_col="id",
        color_col="score",
        colormap={"type": "quantitative", "colors": ["blue", "red"]},
        return_plot="fig",
        warnings=False,
    )

    box_traces = [
        trace
        for trace in fig.data
        if getattr(trace, "fill", None) == "toself" and trace.fillcolor != "white"
    ]
    assert [trace.fillcolor for trace in box_traces] == ["#0000ff", "#ff0000"]


def test_channel_outline_spec_requires_outline_col():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="requires outline_col"):
        pre.plot(
            _style_data(),
            id_col="id",
            color_col="kind",
            colormap={
                "color": {"x": "#111111", "y": "#222222"},
                "outline": {"ok": "#333333", "warn": "#444444"},
            },
            return_plot="fig",
            warnings=False,
        )


def test_fixed_outline_color_is_not_a_colormap_outline_spec():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="outline_color='black'"):
        pre.plot(
            _style_data(),
            id_col="id",
            color_col="kind",
            outline_col="status",
            colormap={"color": {"x": "#111111", "y": "#222222"}, "outline": "black"},
            return_plot="fig",
            warnings=False,
        )


def test_colormap_outline_and_outline_color_are_mutually_exclusive():
    pre.set_engine("matplotlib")

    with pytest.raises(
        ValueError, match=r"both colormap\['outline'\] and outline_color"
    ):
        pre.plot(
            _style_data(),
            id_col="id",
            color_col="kind",
            outline_col="status",
            outline_color="black",
            colormap={
                "color": {"x": "#111111", "y": "#222222"},
                "outline": {"ok": "#333333", "warn": "#444444"},
            },
            return_plot="fig",
            warnings=False,
        )


@pytest.mark.parametrize(
    "colormap, message",
    [
        ({"type": "quantitative"}, "requires a 'colors' entry"),
        ({"type": "continuous", "colors": "viridis"}, "type must be 'quantitative'"),
        ({"type": "quantitative", "colors": {"low": "blue"}}, "dict mappings"),
        (
            {"type": "quantitative", "colors": "viridis", "range": (1, 0)},
            "min must be smaller",
        ),
        (
            {"type": "quantitative", "colors": [(0, "blue"), (2, "red")]},
            "between 0 and 1",
        ),
        ({"color": "Set3", "outline": "Set3", "extra": "bad"}, "only accepts keys"),
    ],
)
def test_invalid_colormap_specs_raise_informative_errors(colormap, message):
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match=message):
        pre.plot(
            _style_data(),
            id_col="id",
            color_col="score",
            outline_col="outline_score",
            colormap=colormap,
            return_plot="fig",
            warnings=False,
        )


def test_quantitative_colormap_requires_numeric_values():
    pre.set_engine("matplotlib")

    with pytest.raises(ValueError, match="requires numeric values"):
        pre.plot(
            _style_data(),
            id_col="id",
            color_col="bad_score",
            colormap={"type": "quantitative", "colors": "viridis"},
            return_plot="fig",
            warnings=False,
        )
