import pyranges as pr
import pyrangeyes as pre
import matplotlib.pyplot as plt

# based on ggsci npg, with some modifications
colors = {
    "red": "#DC0000",
    "redish": "#E64B35",
    "pink": "#F39B7F",
    "green": "#00A087",
    "greenish": "#91D1C2",
    "blue": "#3C5488",
    "blueish": "#8491B4",
    "azure": "#4DBBD5",
    "brown": "#7E6148",
    "brownish": "#B09C85",
}


colormap = {
    "xa": colors["redish"],
    "a": colors["redish"],
    "a2": colors["brownish"],
    "a3": colors["pink"],
    "z": colors["azure"],
}

# Load data
r = pr.PyRanges(
    {
        "Chromosome": [1, 1, 1],
        "Start": [3, 12, 18],
        "End": [9, 14, 22],
        "transcript_id": ["r", "r", "r"],
        "Strand": ["+"] * 3,
        "to_color": ["xa", "xa", "xa"],
    }
)

# extend methods
r_ext = r.extend_ranges(1)
r_ext["transcript_id"] = ["r.extend_ranges(1)"] * len(r_ext)
r_ext_id = r.extend_ranges(2, group_by="transcript_id")
r_ext_id["transcript_id"] = ["r.extend_ranges(2,id)"] * len(r_ext_id)

r_ext_id_out = r.extend_ranges(ext_5=5, group_by="transcript_id")
r_ext_id_out["transcript_id"] = ["k:=r.extend_ranges(ext_5=5,id)"] * len(r_ext_id_out)

r_ext_id_out_clip = r_ext_id_out.clip_ranges()
r_ext_id_out_clip["transcript_id"] = ["k.clip_ranges()"] * len(r_ext_id_out_clip)

for i in [r_ext, r_ext_id, r_ext_id_out, r_ext_id_out_clip]:
    i["to_color"] = ["a3"] * len(i)

# boundaries
b = r.outer_ranges("transcript_id")
b["transcript_id"] = ["r.outer_ranges(id)"] * len(b)
b["to_color"] = ["a2"] * len(b)
# b_subseq = b.subsequence(0, 10)
# b_subseq["transcript_id"] = ["r.boundaries().subsequence(0,10)"] * len(b_subseq)

# upstream, downstream
r_upstream = r.upstream(2, group_by="transcript_id")
r_upstream["transcript_id"] = ["r.upstream(2,id)"] * len(r_upstream)
r_downstream = r.downstream(2, group_by="transcript_id")
r_downstream["transcript_id"] = ["r.downstream(2,id)"] * len(r_downstream)
r_upstream["to_color"] = ["a2"] * len(r_upstream)
r_downstream["to_color"] = ["a2"] * len(r_downstream)

# complement
r_complement = r.complement_ranges("transcript_id")
r_complement["transcript_id"] = ["r.complement_ranges(id)"] * len(r_complement)
r_complement["to_color"] = ["a2"] * len(r_complement)


# subsequence
# g_subseq_id = r.subsequence(0, 10, "transcript_id")
# g_subseq_id["transcript_id"] = ["r.subsequence(0,10,id)"] * len(g_subseq_id)
# g_subseq = r.subsequence(0, 3)
# g_subseq["transcript_id"] = ["r.subsequence(0,3)"] * len(g_subseq)

# spliced subsequence
r_spl_subseq_for = r.slice_ranges(0, 8, "transcript_id")
r_spl_subseq_for["transcript_id"] = ["r.slice_ranges(0,8,id)"] * len(r_spl_subseq_for)
r_spl_subseq_for2 = r.slice_ranges(3, None, "transcript_id")
r_spl_subseq_for2["transcript_id"] = ["r.slice_ranges(3,None,id)"] * len(
    r_spl_subseq_for2
)


r_spl_subseq_rev = r.slice_ranges(-5, -1, "transcript_id")
r_spl_subseq_rev["transcript_id"] = ["r.slice_ranges(-5,-1,id)"] * len(r_spl_subseq_rev)

r_spl_subseq_noid = r.slice_ranges(0, 3)
r_spl_subseq_noid["transcript_id"] = ["r.slice_ranges(0,3)"] * len(r_spl_subseq_noid)

# colors for slices
for i in [
    r_spl_subseq_for,
    r_spl_subseq_for2,
    r_spl_subseq_rev,
    r_spl_subseq_noid,
]:
    i["to_color"] = ["a3"] * len(i)

### add tile thingies
g_window = r.window_ranges(2)
g_window["transcript_id"] = ["r.window_ranges(2)"] * len(g_window)

tile_g = pr.tile_genome(
    pr.PyRanges({"Chromosome": [1], "Start": [0], "End": [max(r["End"]) + 2]}), 2
)
tile_g["transcript_id"] = ["pyranges.tile_genome(2)"] * len(tile_g)  # rename for plot


g_tile = r.tile_ranges(2)
g_tile["transcript_id"] = ["r.tile_ranges(2)"] * len(g_tile)


g_window["to_color"] = ["a3"] * len(g_window)
g_tile["to_color"] = ["a3"] * len(g_tile)
tile_g["to_color"] = ["z"] * len(tile_g)


# Get plot
pre.set_engine("plt")
# customize left margin to fit titles
ori_margin = plt.rcParams["figure.subplot.left"]
plt.rcParams["figure.subplot.left"] = 0.4


for ext in ["png"]:  # , 'pdf']:
    pre.plot(
        [
            r,
            pr.concat([b, r_upstream, r_downstream, r_complement]).reset_index(
                drop=True
            ),
            pr.concat([r_ext_id, r_ext_id_out, r_ext_id_out_clip]).reset_index(
                drop=True
            ),
            # pr.concat([g_subseq_id, g_subseq, ]),
            pr.concat(
                [r_spl_subseq_for, r_spl_subseq_for2, r_spl_subseq_rev]
            ).reset_index(drop=True),
            # pr.concat([r_ext, r_spl_subseq_noid]).reset_index(drop=True),
            pr.concat([g_window, tile_g, g_tile]).reset_index(drop=True),
        ],
        packed=False,
        text=False,
        warnings=False,
        id_col="transcript_id",
        color_col="to_color",
        arrow_line_width=1,
        arrow_color="#4A4A4A",
        outline_color="#4A4A4A",
        theme="light",
        title_chr=" ",
        to_file=(f"cheatsheet_genomic_ops.{ext}", (700, 350)),
        colormap=colormap,
    )


# reset rcparams
plt.rcParams["figure.subplot.left"] = ori_margin
