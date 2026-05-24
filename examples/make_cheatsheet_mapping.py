import pyranges1 as pr
import pyrangeyes as pre
import matplotlib.pyplot as plt

# based on ggsci npg, with some modifications
color_list = [
    "#F39B7F",  # salmon / peach
    "#91D1C2",  # pale turquoise
    "#E64B35",  # vivid red
    "#DC0000",  # crimson red
    "#8491B4",  # muted periwinkle
    "#999999",  # gray
    "#8a705a",  # taupe brown
    "#B09C85",  # light khaki / beige
    "#F0E442",  # light yellow
    "#4DBBD5",  # bright sky-blue
    "#00A087",  # teal / sea-green
]

# Load data
gr = pr.PyRanges(
    {
        "Chromosome": ["chr1", "chr1", "chr1", "chr1"],
        "Start": [100, 300, 500, 600],
        "End": [200, 400, 550, 700],
        "Strand": ["+", "+", "-", "-"],
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
    }
)
tr = pr.PyRanges(
    {
        "Chromosome": ["tx1", "tx1", "tx1", "tx2", "tx2"],
        "Start": [0, 120, 160, 0, 80],
        "End": [80, 140, 170, 20, 130],
        "Strand": ["-", "+", "-", "+", "+"],
        "transcript_id": ["a", "b", "c", "d", "e"],
    }
)


w = pr.PyRanges(
    {
        "Chromosome": ["chr1", "chr1"],
        "Start": [180, 320],
        "End": [200, 340],
        "Strand": ["+", "-"],
        "transcript_id": ["w1", "w2"],
    }
)

z = pr.PyRanges(
    {
        "Chromosome": ["chr1", "chr1"],
        "Start": [500, 570],
        "End": [520, 590],
        "Strand": ["+", "+"],
        "transcript_id": ["z1", "z2"],
    }
)


t1 = tr.loci["tx1"]
t2 = tr.loci["tx2"]

mapped_t1 = t1.map_to_global(gr, "transcript_id")
mapped_t2 = t2.map_to_global(gr, "transcript_id")
mapped_w = w.map_to_local(gr, "transcript_id")
mapped_z = z.map_to_local(gr, "transcript_id")


pre.set_engine("plt")
# customize left margin to fit titles
ori_margin = plt.rcParams["figure.subplot.left"]
plt.rcParams["figure.subplot.left"] = 0.4

for ext in ["png"]:  # , 'pdf']:
    # Plot

    pre.plot(
        [gr, w, z, mapped_t1, mapped_t2, t1, t2, mapped_w, mapped_z],
        id_col="transcript_id",
        fill_col="transcript_id",
        panel_title="Sequence: {chrom}",
        warnings=False,
        y_labels=[
            "g",
            "w",
            "z",
            "t1.map_to_global(g,id)",
            "t2.map_to_global(g,id)",
            "t1",
            "t2",
            "w.map_to_local(g,id)",
            "z.map_to_local(g,id)",
        ],
        limits={"chr1": (60, 700), "tx1": (0, 200), "tx2": (0, 150)},
        theme="light",
        to_file=(f"cheatsheet_mapping.{ext}", (700, 450)),
        arrow_line_width=1,
        arrow_color="#4A4A4A",
        outline_color="#4A4A4A",
        colormap=color_list,
    )

# reset rcparams
plt.rcParams["figure.subplot.left"] = ori_margin
