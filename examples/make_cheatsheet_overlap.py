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
    "xb": colors["blueish"],
    "o": colors["greenish"],
}


# Load data
a = pr.PyRanges(
    {
        "Chromosome": [1] * 8,
        "Start": [2, 12, 17, 22, 27, 24, 32, 33],
        "End": [5, 14, 20, 26, 29, 28, 37, 36],
        "Strand": ["+"] * 2 + ["-"] * 4 + ["+"] * 2,
        "to_color": ["xa"] * 8,
    }
)

b = pr.PyRanges(
    {
        "Chromosome": [1] * 5,
        "Start": [6, 11, 18, 24, 34],
        "End": [8, 13, 19, 28, 36],
        "Strand": ["+"] * 3 + ["-"] * 1 + ["+"] * 1,
        "to_color": ["xb"] * 5,
    }
)

# overlap
a_ov_b = a.overlap(b)
a_ov_b_slack = a.overlap(b, slack=2)
a_ov_b_nostrand = a.overlap(b, strand_behavior="ignore")
a_ov_b_opstrand = a.overlap(b, strand_behavior="opposite")
for pr_obj in [a_ov_b, a_ov_b_slack, a_ov_b_nostrand, a_ov_b_opstrand]:
    print(pr_obj)
    pr_obj["to_color"] = ["a"] * len(pr_obj)

# count overlaps
a_count_b = a.count_overlaps(b)
a_count_b["to_color"] = ["a"] * len(a_count_b)
a_count_b["text"] = a_count_b["Count"]

# intersection
a_inters_b = a.intersect_overlaps(b)
a_inters_b["to_color"] = ["o"] * len(a_inters_b)
a_setinters_b = a.set_intersect_overlaps(b)
a_setinters_b["to_color"] = ["o"] * len(a_setinters_b)

# set union
a_union_b = a.set_union_overlaps(b)
a_union_b["to_color"] = ["o"] * len(a_union_b)


# subtract
a_subt_b = a.subtract_overlaps(b)
a_subt_b["to_color"] = ["o"] * len(a_subt_b)

# merge
a_merge = a.merge_overlaps()
a_merge["to_color"] = ["a2"] * len(a_merge)

# split
a_split = a.split_overlaps()
a_split["to_color"] = ["a2"] * len(a_split)

# max_disjoint
a_max_disjoint = a.max_disjoint_overlaps()
a_max_disjoint["to_color"] = ["a2"] * len(a_max_disjoint)

# cluster
a_cluster = a.cluster_overlaps()
a_cluster["to_color"] = ["a"] * len(a_cluster)
a_cluster["text"] = a_cluster["Cluster"]

# concatenate
a_concat_b = pr.concat([a, b])
# a_concat_b["to_color"] = ["concat"] * len(a_concat_b)


######
# Get plot
pre.set_engine("plt")

# customize left margin to fit titles
ori_margin = plt.rcParams["figure.subplot.left"]
plt.rcParams["figure.subplot.left"] = 0.4

data = [
    a,
    a_merge,
    a_split,
    a_max_disjoint,
    a_cluster,
    b,
    a_ov_b,
    a_ov_b_nostrand,
    a_ov_b_opstrand,
    a_ov_b_slack,
    a_count_b,
    a_inters_b,
    a_setinters_b,
    a_union_b,
    a_subt_b,
    a_concat_b,
]

for i, x in enumerate(data):
    if not "text" in x.columns:
        x["text"] = ""
    # overriding color
    ## x["to_color"] = str(i)

for ext in ["png"]:  # , 'pdf']:
    pre.plot(
        data=data,
        y_labels=[
            "a",
            "a.merge_overlaps()",
            "a.split_overlaps()",
            "a.max_disjoint_overlaps()",
            "a.cluster_overlaps()",
            "b",
            "a.overlap(b)",
            "a.overlap(b,\nstrand_behavior='ignore')",
            "a.overlap(b,\nstrand_behavior='opposite')",
            "a.overlap(b, slack=2)",
            "a.count_overlaps(b)",
            "a.intersect_overlaps(b)",
            "a.set_intersect_overlaps(b)",
            "a.set_union_overlaps(b)",
            "a.subtract_overlaps(b)",
            "pyranges.concat([a, b])",
        ],
        title_chr=" ",
        warnings=False,
        text="{text}",
        text_size=8,
        to_file=(f"cheatsheet_overlap.{ext}", (700, 800)),
        color_col="to_color",
        arrow_line_width=1,
        arrow_color="#4A4A4A",
        outline_color="#4A4A4A",
        theme="light",
        colormap=colormap,
    )

# reset rcparams
plt.rcParams["figure.subplot.left"] = ori_margin
