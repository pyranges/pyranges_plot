import pyrangeyes as pre
import pyranges1 as pr
import pytest
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["savefig.dpi"] = 400

data1 = pr.PyRanges(
    {
        "Chromosome": ["1"] * 9,
        "Strand": ["+", "+", "-", "-", "-", "+", "+", "+", "-"],
        "Start": [i * 100 for i in [5, 35, 3, 13, 35, 45, 49, 56, 60]],
        "End": [i * 100 for i in [15, 37, 6, 17, 39, 47, 51, 57, 67]],
        "transcript_id": ["t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3", "t4"],
        "second_id": ["a"] * 4 + ["b"] * 5,
        "Feature": ["mRNA"] * 9,
    }
)

data2 = pr.PyRanges(
    {
        "Chromosome": ["1"] * 10 + ["2"] * 10 + ["4"],
        "Strand": ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"]
        + ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"]
        + ["+"],
        "Start": [90, 61, 104, 228, 9, 142, 52, 149, 218, 151]
        + [5, 27, 37, 47, 1, 7, 42, 37, 60, 80]
        + [20],
        "End": [92, 64, 113, 229, 12, 147, 57, 155, 224, 153]
        + [8, 32, 40, 50, 5, 10, 46, 40, 70, 90]
        + [50],
        "transcript_id": ["t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t3", "t3"]
        + ["t4", "t4", "t4", "t4", "t5", "t5", "t5", "t5", "t6", "t6", "t7"],
        "Feature": [
            "CDS",
            "exon",
            "exon",
            "exon",
            "exon",
            "exon",
            "exon",
            "exon",
            "exon",
            "exon",
        ]
        + ["exon"] * 8
        + ["exon"] * 3,
    }
)

data3 = pr.PyRanges(
    {
        "Chromosome": ["1", "1", "2", "2", "2", "2", "2", "3", "4", "4", "4", "5"],
        "Strand": ["+", "+", "-", "-", "+", "+", "+", "+", "-", "-", "-", "+"],
        "Start": [1, 40, 10, 70, 85, 110, 150, 140, 5, 170, 240, 100],
        "End": [11, 60, 25, 80, 100, 115, 180, 152, 150, 200, 260, 200],
        "transcript_id": [
            "t1",
            "t1",
            "T2",
            "T2",
            "T3",
            "T3",
            "T3",
            "T4",
            "T5",
            "T5",
            "T5",
            "T6",
        ],
        "Feature": ["exon"] * 12,
    }
)

data4 = pr.PyRanges(
    {
        "Start": [10, 30],
        "End": [40, 60],
        "Chromosome": [1, 1],
        "id": [1, 1],
        "depth": [0, 1],
    }
)

data5 = data4.copy()
data5["depth"] = [1, 0]

data6 = data4.copy()
data6["height"] = [0.4, 1.0]

data7 = pr.PyRanges(
    {
        "Start": [10, 30, 50],
        "End": [25, 45, 65],
        "Chromosome": [1, 1, 1],
        "id": ["a", "b", "c"],
        "fill": ["#ff0000", "#00aa00", "#0000ff"],
        "outline": ["black", "gold", "navy"],
        "score": [0.0, 0.5, 1.0],
        "outline_score": [1.0, 0.5, 0.0],
    }
)

data8 = pr.PyRanges(
    {
        "Chromosome": ["1", "1", "2", "2", "4"],
        "Start": [80, 145, 12, 44, 35],
        "End": [81, 146, 13, 45, 36],
        "ID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        "REF": ["A", "G", "C", "T", "A"],
        "ALT": ["T", "A", "G", "C", "G"],
    }
)


pre.set_engine("plt")
pre.set_id_col("transcript_id")


# test id_col
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test01():
    pre.plot(data1, color_col="transcript_id", outline_color="black", sort_ranges=True)
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test02():
    pre.plot(
        data1,
        id_col="second_id",
        color_col="transcript_id",
        shrink=True,
        outline_color="black",
        sort_ranges=True,
    )
    fig = plt.gcf()
    return fig  # 1 id_col


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test03():
    pre.plot(
        data1,
        id_col=["transcript_id", "second_id"],
        color_col="transcript_id",
        packed=False,
        sort_ranges=True,
        # to_file="tests/img/test03.png",
    )  # +1 id_col, 1 pr packed False
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test04():
    pre.plot(data1[data1["transcript_id"] == "t4"], id_col="transcript_id", shrink=True)
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test05():
    pre.plot(data1[data1["transcript_id"] == "t2"], id_col="transcript_id", shrink=True)
    fig = plt.gcf()
    return fig


# test +1 pr
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test06():
    pre.plot(
        [data2, data3],
        color_col="transcript_id",
        sort_ranges=True,
        # to_file="tests/baseline_mpl/test06.png"
    )  # no id col
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test07():
    pre.plot(
        [data2, data3],
        id_col="Feature",
        color_col="transcript_id",
        text="{Feature}",
        sort_ranges=True,
    )  # 1 id_col, text
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test08():
    pre.plot(
        [data1, data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        y_labels=[1, 2, 3],
        shrink=True,
        sort_ranges=True,
        # to_file="tests/baseline_mpl/test08.png"
    )  # shrink and y_labels
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test09():
    pre.plot(
        [data2, data2],
        id_col="transcript_id",
        packed=False,
        thick_cds=True,
        sort_ranges=True,
    )  # repeated rows in different pr, same chromosome, thick_cds with exon+cds
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test10():
    pre.plot(
        data2, thick_cds=True, limits=(75, 125), text="{Feature}", sort_ranges=True
    )  # thick_cds not all exon+cds, text, limits as tuple
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test11():
    pre.plot(
        data3,
        id_col="transcript_id",
        limits={"1": (None, 1000), "2": (20, 40), "3": None, "4": (-1000, None)},
        # to_file="tests/baseline_mpl/test11.png"
    )  # limits as dict
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test12():
    pre.plot(
        data2,
        id_col="transcript_id",
        limits=data3,
        arrow_size=0.1,
        arrow_color="red",
        sort_ranges=True,
    )  # limit as other pr, arrow_size,colorprp.plot(data2, id_col="transcript_id", limits=data3)
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test13():
    pre.plot(
        [data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        legend=True,
        title_chr="TITLE {chrom}",
        sort_ranges=True,
        # to_file="tests/baseline_mpl/test13.png"
    )  # legend, title string, intron and exon color
    fig = plt.gcf()
    return fig


# depth
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test14():
    pre.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="pastel",
        # to_file="tests/baseline_mpl/test14.png"
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test15():
    pre.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="dark",
        # to_file="tests/baseline_mpl/test15.png"
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test16():
    pre.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        colormap={"0": "#505050", "1": "goldenrod"},
        # to_file="tests/baseline_mpl/test16.png"
    )  # colormap as dict
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test17():
    pre.plot(
        data1,
        id_col=["transcript_id", "second_id"],
        color_col="transcript_id",
        packed=False,
        # to_file="tests/img/test03.png",
    )  # +1 id_col, 1 pr packed False
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test18():
    pre.plot(
        [data2, data2],
        id_col="transcript_id",
        packed=False,
        thick_cds=True,
        # to_file="tests/baseline_mpl/test18.png",
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test19():
    pre.plot(
        data6,
        id_col="id",
        color_col="height",
        height_col="height",
        interval_height=0.8,
        theme="pastel",
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test20():
    pre.plot(
        data7,
        id_col="id",
        color_col="fill",
        outline_col="outline",
        colormap="direct",
        legend=True,
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test21():
    pre.plot(
        data7,
        id_col="id",
        color_col="score",
        outline_col="outline_score",
        colormap={
            "color": {
                "type": "quantitative",
                "colors": ["blue", "white", "red"],
                "range": (0, 1),
            },
            "outline": {"type": "quantitative", "colors": ["black", "gold"]},
        },
        legend=True,
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test22():
    pre.plot(
        data2,
        "mRNA",
        color_col="Feature",
        colormap={"exon": "lightgrey", "CDS": "gold"},
        sort_ranges=True,
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test23():
    pre.plot(
        data2,
        adapter="mRNA",
        utr_height=0.6,
        color_col="Feature",
        colormap={"exon": "lightgrey", "CDS": "gold"},
        sort_ranges=True,
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test24():
    pre.plot(
        data8,
        "SNP",
        color_col="ALT",
        width=8,
        sort_ranges=True,
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test25():
    pre.plot(
        [data2, data8],
        adapter=["mRNA", "SNP"],
        sort_ranges=True,
    )
    fig = plt.gcf()
    return fig
