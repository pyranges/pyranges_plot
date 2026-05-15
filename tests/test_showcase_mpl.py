import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

import pyrangeyes as pre

sys.path.insert(0, str(Path(__file__).parent))
from text_options_showcase_cases import iter_showcase_cases


def _make_test(plot_no, data, kwargs):
    @pytest.mark.mpl_image_compare(baseline_dir="baseline_showcase_mpl")
    def test_func():
        pre.set_engine("plt")
        fig = pre.plot(data, return_plot="fig", warnings=False, **kwargs)
        fig.set_size_inches(15, 6.2)
        return fig

    test_func.__name__ = f"test_showcase_plot_{plot_no:02d}"
    return test_func


for _plot_no, _data, _kwargs in iter_showcase_cases("plt"):
    globals()[f"test_showcase_plot_{_plot_no:02d}"] = _make_test(
        _plot_no, _data, _kwargs
    )

del _plot_no, _data, _kwargs
