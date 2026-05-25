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
        pre.browse(_browser_data(), id_col="id")


def test_browser_adds_one_mode_button_per_panel():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id")

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
    fig = pre.browse(_browser_data(), id_col="id")
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
    fig = pre.browse([data, data], id_col="id")

    assert len(fig.layout.updatemenus) == 2
    assert fig.layout.yaxis2 is not None


def test_browser_mode_buttons_update_panel_shapes_and_y_range():
    pre.set_engine("ply")
    data = _browser_data()
    fig = pre.browse([data, data], id_col="id")

    assert len(fig.layout.shapes) > 0
    first_menu = fig.layout.updatemenus[0]
    full_button = first_menu.buttons[3]
    layout_update = full_button.args[1]

    assert any(key.startswith("shapes[") for key in layout_update)
    assert "yaxis.range" in layout_update


def test_browser_default_mode_controls_initial_visibility():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id", default_mode="zip")

    zip_button = fig.layout.updatemenus[0].buttons[0]
    panel_trace_indices = list(zip_button.args[2])
    zip_visible_mask = list(zip_button.args[0]["visible"])
    actual_visible_mask = [
        fig.data[ix].visible is not False for ix in panel_trace_indices
    ]
    assert actual_visible_mask == zip_visible_mask


def test_browser_full_button_restores_axis_group_labels():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id")

    full_button = fig.layout.updatemenus[0].buttons[3]
    layout_update = full_button.args[1]

    assert layout_update["yaxis.ticktext"] == ["a", "b"]
    assert layout_update["yaxis.tickvals"] == [1.4, 0.55]


def test_browser_modes_update_axis_geometry_and_clear_stale_labels():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id")

    squish_button = fig.layout.updatemenus[0].buttons[1]
    pack_button = fig.layout.updatemenus[0].buttons[2]
    full_button = fig.layout.updatemenus[0].buttons[3]

    assert "yaxis.domain" in squish_button.args[1]
    assert "yaxis2.domain" in squish_button.args[1]
    assert "yaxis.range" in squish_button.args[1]
    assert squish_button.args[1]["yaxis.range"] != pack_button.args[1]["yaxis.range"]

    squish_height = (
        squish_button.args[1]["yaxis.domain"][1]
        - squish_button.args[1]["yaxis.domain"][0]
    ) * (squish_button.args[1]["height"] - 120)
    pack_height = (
        pack_button.args[1]["yaxis.domain"][1] - pack_button.args[1]["yaxis.domain"][0]
    ) * (pack_button.args[1]["height"] - 120)
    squish_other_height = (
        squish_button.args[1]["yaxis2.domain"][1]
        - squish_button.args[1]["yaxis2.domain"][0]
    ) * (squish_button.args[1]["height"] - 120)
    pack_other_height = (
        pack_button.args[1]["yaxis2.domain"][1]
        - pack_button.args[1]["yaxis2.domain"][0]
    ) * (pack_button.args[1]["height"] - 120)
    assert squish_height < pack_height
    assert abs(squish_other_height - pack_other_height) < 1e-6

    assert full_button.args[1]["yaxis.domain"] != pack_button.args[1]["yaxis.domain"]
    assert pack_button.args[1]["yaxis.ticktext"] == []


def test_browser_zip_hides_axes_and_is_shorter_than_squish_for_single_panel():
    pre.set_engine("ply")
    data = pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1"],
            "Start": [10, 20, 30],
            "End": [15, 25, 35],
            "id": ["a", "b", "c"],
        }
    )
    fig = pre.browse(data, id_col="id")

    zip_button = fig.layout.updatemenus[0].buttons[0]
    squish_button = fig.layout.updatemenus[0].buttons[1]

    assert zip_button.args[1]["xaxis.visible"] is False
    assert zip_button.args[1]["yaxis.visible"] is False
    assert (
        squish_button.args[1]["height"]
        < fig.layout.updatemenus[0].buttons[2].args[1]["height"]
    )


def test_browser_mode_switch_moves_subplot_titles_with_domains():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id")

    squish_button = fig.layout.updatemenus[0].buttons[1]
    layout_update = squish_button.args[1]

    assert "annotations[0].y" in layout_update
    assert "annotations[1].y" in layout_update
    assert layout_update["annotations[1].y"] > layout_update["yaxis2.domain"][1]


def test_browser_titles_have_pixel_gap_above_each_panel_after_resizing():
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id")
    for menu in fig.layout.updatemenus:
        for button in menu.buttons:
            layout_update = button.args[1]
            plot_height = layout_update["height"] - 120
            for panel_ix, axis_name in enumerate(["yaxis", "yaxis2"]):
                title_y = layout_update[f"annotations[{panel_ix}].y"]
                domain_top = layout_update[f"{axis_name}.domain"][1]
                assert (title_y - domain_top) * plot_height >= 16


def test_squish_ignores_label_padding_when_packing_rows():
    pre.set_engine("ply")
    data = pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1", "1"],
            "Start": [100, 110, 120, 130],
            "End": [106, 116, 126, 136],
            "id": ["a", "b", "c", "d"],
        }
    )
    fig = pre.browse(data, id_col="id")

    squish_button = fig.layout.updatemenus[0].buttons[1]

    assert squish_button.args[1]["yaxis.range"] == [0.0, 0.123]
