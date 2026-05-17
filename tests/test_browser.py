import pyranges1 as pr
import pytest

import pyrangeyes as pre


def _browser_data():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "2"],
            "Start": [10, 30, 5],
            "End": [20, 40, 15],
            "id": ["a", "b", "c"],
        }
    )


def test_browser_is_plotly_only():
    pre.set_engine("plt")
    with pytest.raises(ValueError, match="Plotly-only"):
        pre.browser(_browser_data(), id_col="id")


def test_browser_adds_one_mode_button_per_panel():
    pre.set_engine("ply")
    fig = pre.browser(_browser_data(), id_col="id")

    assert len(fig.layout.updatemenus) == 2
    assert [button.label for button in fig.layout.updatemenus[0].buttons] == [
        "Zip",
        "Squish",
        "Packed",
        "Full",
    ]
    assert all(menu.x <= 1 for menu in fig.layout.updatemenus)
    assert fig.layout.dragmode == "select"
    assert fig.layout.selectdirection == "h"


def test_browser_buttons_update_only_their_panel_traces():
    pre.set_engine("ply")
    fig = pre.browser(_browser_data(), id_col="id")
    first_menu, second_menu = fig.layout.updatemenus
    first_zip = first_menu.buttons[0]
    second_zip = second_menu.buttons[0]

    first_trace_indices = list(first_zip.args[2])
    second_trace_indices = list(second_zip.args[2])
    assert first_trace_indices
    assert second_trace_indices
    assert set(first_trace_indices).isdisjoint(second_trace_indices)
    assert len(first_zip.args[0]["visible"]) == len(first_trace_indices)
    assert len(second_zip.args[0]["visible"]) == len(second_trace_indices)


def test_browser_uses_plot_like_chrom_panels_for_multiple_objects():
    pre.set_engine("ply")
    data = _browser_data()
    fig = pre.browser([data, data], id_col="id")

    assert len(fig.layout.updatemenus) == 2
    assert fig.layout.yaxis2 is not None


def test_browser_default_mode_controls_initial_visibility():
    pre.set_engine("ply")
    fig = pre.browser(_browser_data(), id_col="id", default_mode="zip")

    zip_button = fig.layout.updatemenus[0].buttons[0]
    panel_trace_indices = list(zip_button.args[2])
    zip_visible_mask = list(zip_button.args[0]["visible"])
    actual_visible_mask = [
        fig.data[ix].visible is not False for ix in panel_trace_indices
    ]
    assert actual_visible_mask == zip_visible_mask
