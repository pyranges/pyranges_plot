# Project Notes

## 2026-05-13

- Matplotlib categorical legends should label entries with the mapped column name, e.g. `Strand: +`, not the rendering channel, e.g. `color: +`. The same applies to mapped outline legends.
- Matplotlib legends/colorbars should live in a reserved bottom band outside the plotting axes. This avoids overlap with plotted intervals, titles, aligned plots, and with other legends when categorical and quantitative legends are both present. Use conservative row slots for figure legends because their rendered text/frame height is only known after drawing.
- Plotly legends/colorbars should also be anchored below the plotting area with enough bottom margin. Use horizontal legends and horizontal quantitative colorbars stacked below the panels so Matplotlib and Plotly showcase the same non-overlapping bottom-legend behavior.
