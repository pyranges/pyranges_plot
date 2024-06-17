import pyranges_plot as prp
import pyranges as pr
from PIL import Image
import os

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


def test_plot_calls():
    # test id_col
    prp.plot(
        data1,
        color_col="transcript_id",
        exon_border="black",
        to_file="tests/img/test01.png",
    )  # no id_col
    prp.plot(
        data1,
        id_col="second_id",
        color_col="transcript_id",
        shrink=True,
        exon_border="black",
        to_file="tests/img/test02.png",
    )  # 1 id_col
    prp.plot(
        data1,
        id_col=["transcript_id", "second_id"],
        color_col="transcript_id",
        packed=False,
        to_file="tests/img/test03.png",
    )  # +1 id_col, 1 pr packed False
    prp.plot(
        data1[data1["transcript_id"] == "t4"],
        id_col="transcript_id",
        shrink=True,
        to_file="tests/img/test04.png",
    )
    prp.plot(
        data1[data1["transcript_id"] == "t2"],
        id_col="transcript_id",
        shrink=True,
        to_file="tests/img/test05.png",
    )

    # test +1 pr
    prp.plot(
        [data2, data3],
        color_col="transcript_id",
        to_file="tests/img/test06.png",
    )  # no id col
    prp.plot(
        [data2, data3],
        id_col="Feature",
        color_col="transcript_id",
        text="{Feature}",
        to_file="tests/img/test07.png",
    )  # 1 id_col, text
    prp.plot(
        [data1, data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        y_labels=[1, 2, 3],
        shrink=True,
        to_file="tests/img/test08.png",
    )  # shrink and y_labels
    prp.plot(
        [data2, data2],
        id_col="transcript_id",
        packed=False,
        thick_cds=True,
        to_file="tests/img/test09.png",
    )  # repeated rows in different pr, same chromosome, thick_cds with exon+cds
    prp.plot(
        data2,
        thick_cds=True,
        limits=(75, 125),
        text="{Feature}",
        to_file="tests/img/test10.png",
    )  # thick_cds not all exon+cds, text, limits as tuple

    prp.plot(
        data3,
        id_col="transcript_id",
        limits={"1": (None, 1000), "2": (20, 40), "3": None, "4": (-1000, None)},
        to_file="tests/img/test11.png",
    )  # limits as dict
    prp.plot(
        data2,
        id_col="transcript_id",
        limits=data3,
        arrow_size=0.1,
        arrow_color="red",
        to_file="tests/img/test12.png",
    )  # limit as other pr, arrow_size,colorprp.plot(data2, id_col="transcript_id", limits=data3)

    prp.plot(
        [data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        legend=True,
        title_chr="TITLE {chrom}",
        to_file="tests/img/test13.png",
    )  # legend, title string, intron and exon color

    # depth
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="Mariotti_lab",
        to_file="tests/img/test14.png",
    )
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="dark",
        to_file="tests/img/test15.png",
    )
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        colormap={"0": "#505050", "1": "goldenrod"},
        to_file="tests/img/test16.png",
    )  # colormap as dict


def test_matplotlib():
    # Create the plots
    prp.set_engine("plt")
    test_plot_calls()

    images = [
        Image.open("tests/img/" + f).convert("RGB")
        for f in sorted([i for i in os.listdir("tests/img/")])
    ]
    images[0].save(
        "tests/matplotlib_tests.pdf", "PDF", save_all=True, append_images=images[1:]
    )


def test_plotly():
    # Create the plots
    prp.set_engine("ply")
    test_plot_calls()

    images = [
        Image.open("tests/img/" + f).convert("RGB")
        for f in sorted([i for i in os.listdir("tests/img/")])
    ]
    images[0].save(
        "tests/plotly_tests.pdf", "PDF", save_all=True, append_images=images[1:]
    )
