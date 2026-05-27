import pyranges1 as pr
import pyrangeyes as pe


DATA = pr.PyRanges(
    {
        "Chromosome": ["1", "1"],
        "Start": [10, 30],
        "End": [20, 40],
        "id": ["a", "b"],
        "kind": ["x", "y"],
    }
)


def test_track_wraps_data_without_copying_and_keeps_adapter_options():
    track = pe.Track(DATA, "mRNA", name="transcripts", squish=True, fill_col="kind")

    assert track.data is DATA
    assert track.adapter == "mRNA"
    assert track.options == {"name": "transcripts", "squish": True, "fill_col": "kind"}


def test_register_methods_adds_plot_and_track_methods():
    pe.register_methods("plotly")

    track = DATA.track("mRNA", name="transcripts", squish=True)

    assert isinstance(track, pe.Track)
    assert track.data is DATA
    assert track.adapter == "mRNA"
    assert track.options == {"name": "transcripts", "squish": True}
    assert DATA.plot(id_col="id", return_plot="fig") is not None


def test_track_option_overrides_plot_default():
    pe.set_engine("matplotlib")

    fig = pe.plot(
        [pe.Track(DATA, name="packed"), pe.Track(DATA, name="unpacked", pack=False)],
        id_col="id",
        pack=True,
        label=False,
        return_plot="fig",
    )

    assert [tick.get_text() for tick in fig.axes[0].get_yticklabels()] == [
        "packed",
        "unpacked",
    ]


def test_track_plot_alias(monkeypatch):
    calls = []

    def fake_plot(obj, **kwargs):
        calls.append((obj, kwargs))
        return "ok"

    import pyrangeyes.plot_main as plot_main

    monkeypatch.setattr(plot_main, "plot", fake_plot)

    track = pe.Track(DATA, name="example")
    result = track.plot(id_col="transcript_id")

    assert result == "ok"
    assert calls == [(track, {"id_col": "transcript_id"})]
