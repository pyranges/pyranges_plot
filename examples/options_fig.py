import pyranges1 as pr
import pyrangeyes as pre

pre.set_engine("plt")
pre.set_id_col("transcript_id")

p = pr.PyRanges(
    {
        "Chromosome": ["1"] * 9,
        "Strand": ["+", "+", "-", "-", "-", "+", "+", "+", "-"],
        "Start": [i * 100 for i in [5, 35, 3, 13, 35, 45, 49, 56, 60]],
        "End": [i * 100 for i in [15, 37, 6, 17, 39, 47, 51, 57, 67]],
        "transcript_id": ["t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3", "t4"],
    }
)

pre.plot(
    p,
    colormap="Alphabet",
    arrow_size=0.03,
    shrink=True,
    shrink_threshold=500,
    label=True,
    label_pad=100,
    title_size=25,
    outline_color="black",
    label_size=15,
)
