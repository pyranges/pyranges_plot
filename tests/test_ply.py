import pyranges_plot as prp
import pyranges as pr
from PIL import Image, ImageChops
import os
import pytest
import matplotlib.pyplot as plt
import plotly.graph_objects as go

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

vcf = pr.PyRanges(
    {
        "Chromosome": ["1"] * 9,  # CHROM renamed to Chromosome
        "Start": [500, 3500, 300, 1300, 3500, 4500, 4900, 5600, 6000],  # POS column renamed to Start
        "End": [501, 3501, 301, 1301, 3501, 4501, 4901, 5601, 6001],  # End is calculated as Start + 1
        "ID": ["."] * 9,  # ID from the VCF file
        "REF": ["A"] * 9,  # REF column
        "ALT": ["T"] * 9,  # ALT column
        "QUAL": ["."] * 9,  # QUAL column
        "FILTER": ["PASS"] * 9,  # FILTER column
        "transcript_id": ["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9"],  # Extracted from INFO
        "second_id": ["a", "a", "a", "a", "b", "b", "b", "b", "b"],  # Extracted from INFO
        "Count": [1, 2, 3, 4, 5, 6, 7, 8, 9]
    }
)

# Create a scatter plot trace with the adapted Start positions
aligned_traces = [
    (go.Scatter(
        x=[500, 3500, 300, 1300, 3500, 4500, 4900, 5600, 6000],  # Adapted Start positions
        y=[1, 2, 3, 4, 5, 6, 7, 8, 9],  # Arbitrary y-values for demonstration
        mode='markers'
    ), {"title": "VCF Scatter Plot"})
]

prp.set_engine("ply")
prp.set_id_col("transcript_id")

BASELINE_DIR = "tests/baseline_ply"
RESULTS_DIR = "tests/results_plotly"

@pytest.mark.parametrize("test_name", ["test01"])
# test id_col
def test01(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data1,
        color_col="transcript_id",
        exon_border="black",
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test02"])
def test02(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data1,
        id_col="second_id",
        color_col="transcript_id",
        shrink=True,
        exon_border="black",
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test03"])
def test03(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data1,
        id_col=["transcript_id", "second_id"],
        color_col="transcript_id",
        packed=False,
        to_file=output_path
        #to_file="tests/img/test03.png",
    )  # +1 id_col, 1 pr packed False
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test04"])
def test04(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data1[data1["transcript_id"] == "t4"],
        id_col="transcript_id",
        shrink=True,
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test05"])
def test05(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data1[data1["transcript_id"] == "t2"],
        id_col="transcript_id",
        shrink=True,
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

# test +1 pr
@pytest.mark.parametrize("test_name", ["test06"])
def test06(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data2, data3],
        color_col="transcript_id",
        to_file=output_path
    )  # no id col
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test07"])
def test07(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data2, data3],
        id_col="Feature",
        color_col="transcript_id",
        text="{Feature}",
        to_file=output_path
    )  # 1 id_col, text
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test08"])
def test08(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data1, data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        y_labels=[1, 2, 3],
        shrink=True,
        to_file=output_path
    )  # shrink and y_labels
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test09"])
def test09(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data2, data2],
        id_col="transcript_id",
        packed=False,
        thick_cds=True,
        to_file=output_path
    )  # repeated rows in different pr, same chromosome, thick_cds with exon+cds
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test10"])
def test10(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data2,
        thick_cds=True,
        limits=(75, 125),
        text="{Feature}",
        to_file=output_path
    )  # thick_cds not all exon+cds, text, limits as tuple
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test11"])
def test11(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data3,
        id_col="transcript_id",
        limits={"1": (None, 1000), "2": (20, 40), "3": None, "4": (-1000, None)},
        to_file=output_path
    )  # limits as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test12"])
def test12(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        data2,
        id_col="transcript_id",
        limits=data3,
        arrow_size=0.1,
        arrow_color="red",
        to_file=output_path
    )  # limit as other pr, arrow_size,colorprp.plot(data2, id_col="transcript_id", limits=data3)
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test13"])
def test13(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data2, data3],
        id_col="transcript_id",
        color_col="Feature",
        legend=True,
        title_chr="TITLE {chrom}",
        to_file=output_path
    )  # legend, title string, intron and exon color
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

# depth
@pytest.mark.parametrize("test_name", ["test14"])
def test14(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="Mariotti_lab",
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test15"])
def test15(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        theme="dark",
        to_file=output_path
    )
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test16"])
def test16(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data4, data5],
        id_col="id",
        color_col="depth",
        depth_col="depth",
        tooltip="{depth}",
        colormap={"0": "#505050", "1": "goldenrod"},
        to_file=output_path
    )  # colormap as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

@pytest.mark.parametrize("test_name", ["test17"])
def test17(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data1, vcf],
        id_col="transcript_id",
        to_file=output_path
    )  # colormap as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

aligned = prp.make_scatter(vcf,y='Count')
@pytest.mark.parametrize("test_name", ["test18"])
def test18(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data1, vcf],
        id_col="transcript_id",
        add_aligned_plots = [aligned],
        to_file=output_path,
    )  # colormap as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

aligned1 = prp.make_scatter(vcf,y='Count',color_by="second_id", title="Human Variants", title_color="Magenta",title_size=18)
@pytest.mark.parametrize("test_name", ["test19"])
def test19(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data1, vcf],
        id_col="transcript_id",
        add_aligned_plots = [aligned1],
        to_file=output_path,
    )  # colormap as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")

aligned2 = prp.make_scatter(vcf,y='Count',color_by="second_id", title="Human Variants", title_color="Magenta",title_size=18,height=0.5,y_space=0.5)
@pytest.mark.parametrize("test_name", ["test20"])
def test20(test_name):
    output_path = os.path.join(RESULTS_DIR, f"{test_name}.png")
    prp.plot(
        [data1, vcf],
        id_col="transcript_id",
        add_aligned_plots = [aligned2],
        to_file=output_path,
    )  # colormap as dict
    baseline_path = os.path.join(BASELINE_DIR, f"{test_name}.png")

    # Opening images
    baseline = Image.open(baseline_path).convert("RGB")
    result = Image.open(output_path).convert("RGB")

    # Resizing images to avoid pytest crash
    if result.size != baseline.size:
        result = result.resize(baseline.size)

    diff = ImageChops.difference(baseline, result)

    if diff.getbbox():
        # Save the diff image for inspection
        diff_path = os.path.join(RESULTS_DIR, f"{test_name}_diff.png")
        diff.save(diff_path)
        pytest.fail(f"{test_name} does not match the baseline image. See difference at {diff_path}")