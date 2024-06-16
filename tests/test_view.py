# import pyranges_plot

from pyranges_plot.example_data import data1, data2, data3
from pyranges_plot import plot
from pyranges_plot.core import set_engine, set_warnings

set_engine("plt")
set_warnings(False)

# test id_col
plot(data1, color_col="transcript_id")  # no id_col
plot(data1, id_col="second_id", color_col="transcript_id")  # 1 id_col
plot(
    data1,
    id_col=["transcript_id", "second_id"],
    color_col="transcript_id",
    packed=False,
)  # +1 id_col, 1 pr packed False
plot(data1[data1["transcript_id"] == "t4"], id_col="transcript_id", shrink=True)
plot(data1[data1["transcript_id"] == "t2"], id_col="transcript_id", shrink=True)

# test +1 pr
plot([data2, data3], color_col="transcript_id")  # no id col
plot(
    [data2, data3], id_col="Feature", color_col="transcript_id", text="{Feature}"
)  # 1 id_col, text
plot(
    [data1, data2, data3],
    id_col="transcript_id",
    color_col="Feature",
    y_labels=[1, 2, 3],
    shrink=True,
)  # shrink and y_labels
plot(
    [data2, data2], packed=False, thick_cds=True
)  # repeated rows in different pr, same chromosome, thick_cds
