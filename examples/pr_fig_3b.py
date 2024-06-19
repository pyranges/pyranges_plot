import pyranges as pr
import pyranges_plot as prp
import matplotlib.pyplot as plt

# Start r object
r = pr.PyRanges(
    {
        "Chromosome": [1, 1, 1],
        "Start": [3, 11, 16],
        "End": [8, 15, 19],
        "transcript_id": ["r", "r", "r"],
    }
)

# r.window object
r_window = r.window(2)
r_window["transcript_id"] = ["r.window(2)"] * len(r_window)  # rename for plot

# tile_genome object
tile_g = pr.tile_genome(
    pr.PyRanges({"Chromosome": [1], "Start": [0], "End": [max(r["End"])]}), 2
)
tile_g["transcript_id"] = ["pr.tile_genome(2)"] * len(tile_g)  # rename for plot

# r.tile object
r_tile = r.tile(2)
r_tile["transcript_id"] = ["r.tile(2)"] * len(r_tile)  # rename for plot

# Get plot
prp.set_engine("plt")

# customize left margin to fit titles
ori_margin = plt.rcParams["figure.subplot.left"]
plt.rcParams["figure.subplot.left"] = 0.27

prp.plot(
    [r, r_window, tile_g, r_tile],
    packed=False,
    text=False,
    id_col="transcript_id",
    exon_border="black",
    title_chr=" ",
    to_file=("fig_3b.png", (600, 400)),
    theme="Mariotti_lab",
)

# reset rcparams
plt.rcParams["figure.subplot.left"] = ori_margin
