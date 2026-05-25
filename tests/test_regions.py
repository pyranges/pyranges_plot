import matplotlib
import pyranges1 as pr
import pyrangeyes as pre

matplotlib.use("Agg")


def _data():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1", "1", "2"],
            "Start": [1000, 2000, 50000, 55000, 10],
            "End": [1500, 2500, 51000, 56000, 20],
            "Strand": ["+"] * 5,
            "tx": ["a", "a", "b", "b", "c"],
            "region_group": ["left", "left", "right", "right", "chr2"],
            "Feature": ["exon"] * 5,
        }
    )


def test_regions_list_replaces_chromosome_layout_matplotlib():
    pre.set_engine("matplotlib")
    fig = pre.plot(
        _data(),
        id_col="tx",
        regions=[("1", 50_000, 60_000), ("1", 1_000, 5_000)],
        return_plot="fig",
        panel_title="{chrom}:{start}-{end}",
    )

    assert [ax.get_title() for ax in fig.axes] == [
        "1:50,000-60,000",
        "1:1,000-5,000",
    ]


def test_regions_list_accepts_embedded_pyranges_matplotlib():
    pre.set_engine("matplotlib")
    regions_pr = pr.PyRanges({"Chromosome": ["2"], "Start": [0], "End": [100]})
    fig = pre.plot(
        _data(),
        id_col="tx",
        regions=[("1", 1_000, 5_000), regions_pr],
        return_plot="fig",
        panel_title="{chrom}:{start}-{end}",
    )

    assert [ax.get_title() for ax in fig.axes] == ["1:1,000-5,000", "2:0-100"]


def test_regions_column_uses_groups_as_panels_matplotlib():
    pre.set_engine("matplotlib")
    fig = pre.plot(_data(), id_col="tx", regions="region_group", return_plot="fig")

    assert [ax.get_title() for ax in fig.axes] == ["left", "right", "chr2"]


def test_regions_list_replaces_chromosome_layout_plotly():
    pre.set_engine("plotly")
    fig = pre.plot(
        _data(),
        id_col="tx",
        regions=[("1", 50_000, 60_000), ("1", 1_000, 5_000)],
        return_plot="fig",
        panel_title="{chrom}:{start}-{end}",
    )

    title_texts = [annotation.text for annotation in fig.layout.annotations[:2]]
    assert title_texts == ["1:50,000-60,000", "1:1,000-5,000"]
