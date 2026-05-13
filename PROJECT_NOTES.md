# Project Notes

## 2026-05-13

- Matplotlib categorical legends should label entries with the mapped column name, e.g. `Strand: +`, not the rendering channel, e.g. `color: +`. The same applies to mapped outline legends.
- Matplotlib legends/colorbars should live in a reserved bottom band outside the plotting axes. This avoids overlap with plotted intervals, titles, aligned plots, and with other legends when categorical and quantitative legends are both present.
