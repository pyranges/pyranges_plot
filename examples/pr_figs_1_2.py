import pyranges as pr
import pyranges_plot as prp

# Load data
gr = pr.PyRanges(
    {
        "Chromosome": ["1"] * 10 + ["2"] * 10,
        "Strand": ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"]
        + ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"],
        "Start": [90, 61, 104, 228, 9, 142, 52, 149, 218, 151]
        + [6, 27, 37, 47, 1, 7, 42, 37, 60, 80],
        "End": [92, 64, 113, 229, 12, 147, 57, 155, 224, 153]
        + [8, 32, 40, 50, 5, 10, 46, 40, 70, 90],
        "transcript_id": ["t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t3", "t3"]
        + ["t4", "t4", "t4", "t4", "t5", "t5", "t5", "t5", "t6", "t6"],
    }
)

# Subset data
gr_1 = gr[gr["Chromosome"] == "1"]

gr_2 = pr.example_data.ncbi_gff
gr_2 = gr_2[gr_2.Feature.isin(["CDS", "exon"])]
gr_2 = gr_2[gr_2.Parent.isin(["rna-DGYR_LOCUS12552-2", "rna-DGYR_LOCUS12552"])]
gr_2 = gr_2[["Chromosome", "Feature", "Start", "End", "Strand", "Parent"]]


# Figure 1
# interactive plot to show tooltip (Figure 1.A)
prp.set_engine("plt")
prp.plot(gr, title_size=30, text_size=15)

# show id, color col and cmap and save as png (Figure 1.B)
prp.plot(
    gr_1,
    id_col="transcript_id",
    color_col="Strand",
    colormap={"+": "lightgreen", "-": "lightblue"},
    to_file=("fig_1b.png", (500, 200)),
    limits=(-5, None),
    text_pad=2,
)


# Figure 2
# show thick_cds and save png (Figure 2.A)
prp.plot(
    gr_2,
    id_col="Parent",
    thick_cds=True,
    text=False,
    to_file=("fig_2a.png", (500, 300)),
    arrow_size=0.04,
)

# show introns off and save png (Figure 2.B)
prp.plot(
    gr_1,
    id_col="transcript_id",
    shrink=True,
    to_file=("fig_2b.png", (500, 300)),
    arrow_size=0.04,
    text_pad=2,
    limits=(5, None),
)
