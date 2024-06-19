import pyranges as pr
import pyranges_plot as prp
import matplotlib.pyplot as plt

# Load data
a = pr.PyRanges(
    {
        "Chromosome": [1] * 7,
        "Start": [2, 12, 17, 22, 27, 32, 33],
        "End": [5, 14, 20, 26, 29, 37, 36],
        "Strand": ["+"] * 2 + ["-"] * 3 + ["+"] * 2,
        "to_color": ["a"] * 7,
    }
)

b = pr.PyRanges(
    {
        "Chromosome": [1] * 5,
        "Start": [6, 11, 18, 24, 34],
        "End": [8, 13, 19, 28, 36],
        "Strand": ["+"] * 3 + ["-"] * 1 + ["+"] * 1,
        "to_color": ["b"] * 5,
    }
)

# overlap
a_ov_b = a.overlap(b)
a_ov_b_slack = a.overlap(b, slack=2)
a_ov_b_nostrand = a.overlap(b, strand_behavior="ignore")
a_ov_b_opstrand = a.overlap(b, strand_behavior="opposite")
for pr_obj in [a_ov_b, a_ov_b_slack, a_ov_b_nostrand, a_ov_b_opstrand]:
    print(pr_obj)
    pr_obj["to_color"] = ["overlap"] * len(pr_obj)

# intersection
a_inters_b = a.intersect(b)
a_inters_b["to_color"] = ["intersect"] * len(a_inters_b)
a_setinters_b = a.set_intersect(b)
a_setinters_b["to_color"] = ["setintersect"] * len(a_setinters_b)


# subtract
a_subt_b = a.subtract_ranges(b)
a_subt_b["to_color"] = ["subtract"] * len(a_subt_b)

# Get plot
prp.set_engine("plt")

# customize left margin to fit titles
ori_margin = plt.rcParams["figure.subplot.left"]
plt.rcParams["figure.subplot.left"] = 0.4

prp.plot(
    [
        a,
        b,
        a_ov_b,
        a_ov_b_slack,
        a_ov_b_nostrand,
        a_ov_b_opstrand,
        a_inters_b,
        a_setinters_b,
        a_subt_b,
    ],
    y_labels=[
        "a",
        "b",
        "a.overlap(b)",
        "a.overlap(b, slack=2)",
        "a.overlap(b, strand_behavior='ignore')",
        "a.overlap(b, strand_behavior='opposite')",
        "a.intersect(b)",
        "a.set_intersect(b)",
        "a.subtract_ranges(b)",
    ],
    title_chr=" ",
    warnings=False,
    text=False,
    to_file=("fig_3c.png", (800, 400)),
    color_col="to_color",
    arrow_color="black",
    arrow_line_width=0.5,
    theme="Mariotti_lab",
)

# reset rcparams
plt.rcParams["figure.subplot.left"] = ori_margin
