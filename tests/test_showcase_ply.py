import json
import os
import sys
from pathlib import Path

from deepdiff import DeepDiff
import numpy as np
import pytest

import pyrangeyes as pre

sys.path.insert(0, str(Path(__file__).parent))
from text_options_showcase_cases import iter_showcase_cases

BASELINE_DIR = os.path.join(os.path.dirname(__file__), "baseline_showcase_ply")


def normalize_plotly_json(obj):
    if isinstance(obj, dict):
        return {k: normalize_plotly_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [normalize_plotly_json(i) for i in obj]
    if isinstance(obj, tuple):
        return [normalize_plotly_json(i) for i in obj]
    if isinstance(obj, np.ndarray):
        return normalize_plotly_json(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


CASES = list(iter_showcase_cases("ply"))


@pytest.mark.parametrize(
    "plot_no,data,kwargs",
    CASES,
    ids=lambda value: f"plot{value:02d}" if isinstance(value, int) else None,
)
def test_showcase_plotly_baseline(plot_no, data, kwargs):
    pre.set_engine("ply")
    fig = pre.plot(data, return_plot="fig", warnings=False, **kwargs)
    fig.update_layout(width=1700, height=700, margin=dict(l=25, r=25, t=45, b=35))
    actual_json = normalize_plotly_json(fig.to_plotly_json())
    baseline_path = os.path.join(BASELINE_DIR, f"plot{plot_no:02d}.json")
    with open(baseline_path) as handle:
        expected_json = normalize_plotly_json(json.load(handle))
    diff = DeepDiff(
        expected_json, actual_json, ignore_order=True, report_repetition=True
    )
    assert diff == {}, (
        f"showcase plot {plot_no:02d} does not match baseline:\n{diff.pretty()}"
    )
