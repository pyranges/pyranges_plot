import matplotlib
import pyranges1 as pr
import pyrangeyes as pre

matplotlib.use("Agg")


def _top_to_bottom_ylabels(ax):
    labels_and_ticks = zip(
        [tick.get_text() for tick in ax.get_yticklabels()],
        ax.get_yticks(),
    )
    return [
        label for label, _ in sorted(labels_and_ticks, key=lambda x: x[1], reverse=True)
    ]


def _out_of_genomic_order_data():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1", "1"],
            "Start": [5, 13, 18, 4],
            "End": [10, 15, 22, 7],
            "transcript_id": ["first", "second", "second", "third"],
        }
    )


def test_packed_false_preserves_input_row_order_by_default_matplotlib():
    pre.set_engine("matplotlib")

    fig = pre.plot(
        _out_of_genomic_order_data(),
        id_col="transcript_id",
        packed=False,
        return_plot="fig",
        text=True,
        warnings=False,
        title_chr=" ",
    )

    assert _top_to_bottom_ylabels(fig.axes[0]) == ["first", "second", "third"]


def test_text_none_defaults_to_packed_only_matplotlib():
    pre.set_engine("matplotlib")

    packed_fig = pre.plot(
        _out_of_genomic_order_data(),
        id_col="transcript_id",
        packed=True,
        return_plot="fig",
        warnings=False,
        title_chr=" ",
    )
    unpacked_fig = pre.plot(
        _out_of_genomic_order_data(),
        id_col="transcript_id",
        packed=False,
        return_plot="fig",
        warnings=False,
        title_chr=" ",
    )

    packed_text = {text.get_text() for text in packed_fig.axes[0].texts}
    unpacked_text = {text.get_text() for text in unpacked_fig.axes[0].texts}

    assert {"first", "second", "third"}.issubset(packed_text)
    assert unpacked_text == {""}
