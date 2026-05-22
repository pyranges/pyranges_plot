import pyranges1 as pr
from pyrangeyes.example_data import ncbi_gff

COLORS = ["#5b8ff9", "#61dDAA", "#65789B", "#F6BD16", "#7262fd", "#78D3F8", "#9661BC"]
OUTLINES = ["#1d39c4", "#237804", "#3a3a3a", "#ad6800", "#391085", "#0050b3", "#531dab"]


def overlapping_data():
    n = 7
    return pr.PyRanges(
        {
            "Chromosome": ["chr1"] * n,
            "Start": [100, 106, 112, 118, 124, 130, 136],
            "End": [104, 110, 116, 122, 128, 134, 140],
            "id": [f"tx{i}" for i in range(1, n + 1)],
            "label": [f"long_label_{i}_overlap_demo" for i in range(1, n + 1)],
            "kind": ["exon", "CDS", "UTR", "SNP", "motif", "peak", "repeat"],
            "fill": COLORS,
            "outline": OUTLINES,
            "score": [0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00],
            "height": [0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95],
        }
    )


def transcript_data():
    return ncbi_gff().loci["1", "-", 173903000:173950000]


def transcript_list_data():
    tx = transcript_data()
    return [tx, tx.loci["1", "-", 173930000:173950000]]


SIMPLE_CASES = [
    dict(
        data=overlapping_data,
        kwargs=dict(id_col="id", packed=True, text="{id}"),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            text="{label}",
            text_position="right",
            text_fit=False,
            text_pad=1,
            text_size=8,
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            text="{label}",
            text_position="right",
            text_fit=True,
            text_pad=1,
            text_size=8,
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            text="{kind}",
            text_position="above",
            text_size=9,
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            text="{kind}",
            text_position="below",
            text_size=9,
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            color_col="fill",
            outline_col="outline",
            colormap="direct",
            interval_height=0.7,
            text="{id}:{kind}",
            text_position="center",
            text_color="white",
            text_size=8,
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            color_col="score",
            height_col="height",
            interval_height=0.9,
            colormap={
                "color": {"type": "quantitative", "colors": ["#d9f7be", "#237804"]}
            },
            text="{score}",
            text_position="above",
            text_angle=25,
            text_size=8,
        ),
    ),
]

POSITION_PADS = {
    "left": [0, 1, 5],
    "right": [0, 1, 5],
    "above": [0, 1, 5],
    "below": [0, 1, 5],
    "center": [0, 1, 5],
}


TEXT_COLORMAP_CASES = [
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            color_col="kind",
            colormap={
                "color": {
                    "exon": "#d9e8ff",
                    "CDS": "#f6bd16",
                    "UTR": "#b7eb8f",
                    "SNP": "#ffadd2",
                    "motif": "#d3adf7",
                    "peak": "#87e8de",
                    "repeat": "#ffd591",
                },
                "text": {
                    "exon": "#003a8c",
                    "CDS": "#ad4e00",
                    "UTR": "#237804",
                    "SNP": "#c41d7f",
                    "motif": "#531dab",
                    "peak": "#006d75",
                    "repeat": "#ad6800",
                },
            },
            text="{kind}",
            text_position="above",
            text_size=9,
            text_color_col="kind",
        ),
    ),
    dict(
        data=overlapping_data,
        kwargs=dict(
            id_col="id",
            packed=True,
            color_col="kind",
            colormap={
                "exon": "#d9e8ff",
                "CDS": "#f6bd16",
                "UTR": "#b7eb8f",
                "SNP": "#ffadd2",
                "motif": "#d3adf7",
                "peak": "#87e8de",
                "repeat": "#ffd591",
            },
            text="{kind}",
            text_position="center",
            text_size=8,
            text_color="#262626",
        ),
    ),
]


def iter_showcase_cases(engine):
    plot_no = 1 if engine == "plt" else 8
    for case in SIMPLE_CASES:
        yield plot_no, case["data"](), case["kwargs"]
        plot_no += 1

    plot_no = 15 if engine == "plt" else 30
    for position, pads in POSITION_PADS.items():
        for pad in pads:
            yield (
                plot_no,
                transcript_data(),
                dict(
                    adapter="mRNA",
                    packed=True,
                    text="{transcript_id}",
                    text_position=position,
                    text_pad=pad,
                    text_size=12,
                    color_col="Feature",
                    colormap={"exon": "#d9e8ff", "CDS": "#f6bd16"},
                    sort_ranges=True,
                ),
            )
            plot_no += 1

    plot_no = 45 if engine == "plt" else 47
    for position, pad in (("right", 5), ("above", 5)):
        yield (
            plot_no,
            transcript_list_data(),
            dict(
                adapter=["mRNA", "mRNA"],
                packed=True,
                text="{transcript_id}",
                text_position=position,
                text_pad=pad,
                text_size=12,
                color_col="Feature",
                colormap={"exon": "#d9e8ff", "CDS": "#f6bd16"},
                sort_ranges=True,
            ),
        )
        plot_no += 1

    plot_no = 60 if engine == "plt" else 70
    for case in TEXT_COLORMAP_CASES:
        yield plot_no, case["data"](), case["kwargs"]
        plot_no += 1
