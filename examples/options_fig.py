import pyranges as pr
import pyranges_plot as prp

prp.set_engine("plt")
prp.set_id_col("transcript_id")

p = pr.PyRanges(
    {
        "Chromosome": ["1"] * 9,
        "Strand": ["+", "+", "-", "-", "-", "+", "+", "+", "-"],
        "Start": [i * 100 for i in [5, 35, 3, 13, 35, 45, 49, 56, 60]],
        "End": [i * 100 for i in [15, 37, 6, 17, 39, 47, 51, 57, 67]],
        "transcript_id": ["t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3", "t4"],
    }
)

prp.plot(
    p,
    colormap="Alphabet",
    arrow_size=0.03,
    shrink=True,
    shrink_threshold=500,
    text=True,
    text_pad=100,
    title_size=25,
    exon_border="black",
    text_size=15,
)
