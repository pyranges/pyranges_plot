import pyranges_plot as prp
import pyranges as pr
import pytest
import matplotlib.pyplot as plt

data1 = pr.PyRanges(
    {
        "Chromosome": ["1"] * 9,
        "Strand": ["+", "+", "-", "-", "-", "+", "+", "+", "-"],
        "Start": [i * 100 for i in [5, 35, 3, 13, 35, 45, 49, 56, 60]],
        "End": [i * 100 for i in [15, 37, 6, 17, 39, 47, 51, 57, 67]],
        "transcript_id": ["t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3", "t4"],
        "second_id": ["a"] * 4 + ["b"] * 5,
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


prp.set_engine("plt")
prp.set_id_col("transcript_id")


# test id_col
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test01():
    prp.plot(data1, color_col="transcript_id", exon_border="black")
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test02():
    prp.plot(
        data1,
        id_col="second_id",
        color_col="transcript_id",
        shrink=True,
        exon_border="black",
    )
    fig = plt.gcf()
    return fig  # 1 id_col


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test03():
    prp.plot(
        data1,
        id_col=["transcript_id", "second_id"],
        color_col="transcript_id",
        packed=False,
        # to_file="tests/img/test03.png",
    )  # +1 id_col, 1 pr packed False
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test04():
    prp.plot(data1[data1["transcript_id"] == "t4"], id_col="transcript_id", shrink=True)
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test05():
    prp.plot(data1[data1["transcript_id"] == "t2"], id_col="transcript_id", shrink=True)
    fig = plt.gcf()
    return fig


# test +1 pr
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test06():
    prp.plot([data2, data3], color_col="transcript_id")  # no id col
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test07():
    prp.plot(
        [data2, data3], id_col="Feature", color_col="transcript_id", text="{Feature}"
    )  # 1 id_col, text
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test08():
    prp.plot(
        [data1, data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        y_labels=[1, 2, 3],
        shrink=True,
    )  # shrink and y_labels
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test09():
    prp.plot(
        [data2, data2], id_col="transcript_id", packed=False, thick_cds=True
    )  # repeated rows in different pr, same chromosome, thick_cds with exon+cds
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test10():
    prp.plot(
        data2, thick_cds=True, limits=(75, 125), text="{Feature}"
    )  # thick_cds not all exon+cds, text, limits as tuple
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test11():
    prp.plot(
        data3,
        id_col="transcript_id",
        limits={"1": (None, 1000), "2": (20, 40), "3": None, "4": (-1000, None)},
    )  # limits as dict
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test12():
    prp.plot(
        data2, id_col="transcript_id", limits=data3, arrow_size=0.1, arrow_color="red"
    )  # limit as other pr, arrow_size,colorprp.plot(data2, id_col="transcript_id", limits=data3)
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test13():
    prp.plot(
        [data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        legend=True,
        title_chr="TITLE {chrom}",
    )  # legend, title string, intron and exon color
    fig = plt.gcf()
    return fig


# depth
@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test14():
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="Mariotti_lab",
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test15():
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="dark",
    )
    fig = plt.gcf()
    return fig


@pytest.mark.mpl_image_compare(baseline_dir="baseline_mpl")
def test16():
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        colormap={"0": "#505050", "1": "goldenrod"},
    )  # colormap as dict
    fig = plt.gcf()
    return fig
