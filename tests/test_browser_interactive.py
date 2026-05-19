"""Interactive browser regression tests.

These tests intentionally click through the rendered HTML in a real browser and
measure rendered pixels. They are guardrails for the browser UX, not immutable
visual snapshots. If a future intentional visual change improves the browser but
moves these pixel-level invariants, update/xfail the specific assertion only with
an explicit reason and maintainer approval.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyranges1 as pr
import pytest
from plotly.utils import PlotlyJSONEncoder

import pyrangeyes as pre

pytestmark = pytest.mark.browser_interactive

CHROME_CANDIDATES = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
]


@pytest.fixture(scope="module")
def sync_playwright():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright as _sync_playwright

    with _sync_playwright() as playwright:
        yield playwright


@pytest.fixture(scope="module")
def browser(sync_playwright):
    launch_kwargs = {"headless": True}
    executable = next((path for path in CHROME_CANDIDATES if Path(path).exists()), None)
    if executable is not None:
        launch_kwargs["executable_path"] = executable
    try:
        browser = sync_playwright.chromium.launch(**launch_kwargs)
    except Exception as exc:  # pragma: no cover - environment-dependent skip
        pytest.skip(f"Chromium/Playwright browser is not available: {exc}")
    yield browser
    browser.close()


@pytest.fixture(scope="module")
def browser_html(tmp_path_factory):
    pre.set_engine("ply")
    fig = pre.browse(_browser_data(), id_col="id", default_mode="packed")
    fig.update_layout(
        autosize=True,
        width=None,
        margin=dict(l=90, r=16, t=84, b=36),
        dragmode="select",
        selectdirection="h",
    )
    meta = _browser_meta(fig)
    fig.layout.updatemenus = ()
    html = fig.to_html(full_html=True, include_plotlyjs="cdn", config={"scrollZoom": True})
    html = html.replace("<head>", "<head>" + _style())
    html = html.replace("</body>", _enhancement_script(meta) + "</body>")
    path = tmp_path_factory.mktemp("browser-html") / "browser.html"
    path.write_text(html, encoding="utf-8")
    return path


def _browser_data():
    return pr.PyRanges(
        {
            "Chromosome": ["chr1"] * 10 + ["chr2"] * 5 + ["chr3"] * 3,
            "Start": [
                100,
                112,
                124,
                136,
                148,
                160,
                172,
                184,
                196,
                208,
                40,
                52,
                64,
                76,
                88,
                10,
                26,
                42,
            ],
            "End": [
                108,
                120,
                132,
                144,
                156,
                168,
                180,
                192,
                204,
                216,
                48,
                60,
                72,
                84,
                96,
                20,
                36,
                52,
            ],
            "id": [f"chr1_gene_{i}" for i in range(1, 11)]
            + [f"chr2_gene_{i}" for i in range(1, 6)]
            + [f"chr3_gene_{i}" for i in range(1, 4)],
        }
    )


def _style():
    return """
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>
html, body { margin: 0; padding: 0; width: 100%; min-height: 100%; font-family: system-ui, sans-serif; }
.plotly-graph-div { width: 100% !important; }
.pyrangeyes-browser-shell { position: relative; }
.pre-view-control { position: absolute; z-index: 50; display: inline-block; }
.pre-view-button { width: 28px; height: 28px; border-radius: 7px; border: 1px solid #9aa4b2; background: rgba(255,255,255,.92); }
.pre-view-menu { display: none; position: absolute; right: 32px; top: 0; min-width: 98px; padding: 4px; border: 1px solid #d0d7de; border-radius: 8px; background: #fff; box-shadow: 0 8px 24px rgba(140,149,159,.35); }
.pre-view-control.open .pre-view-menu { display: block; }
.pre-view-option { display: block; width: 100%; border: 0; background: transparent; padding: 7px 8px; text-align: left; }
.pre-view-option .check { display: inline-block; width: 16px; color: #0969da; font-weight: 700; }
.pre-zip-block { position: absolute; z-index: 35; display: none; align-items: center; justify-content: center; text-align: center; box-sizing: border-box; padding: 2px 44px 2px 8px; color: #57606a; font-size: 13px; line-height: 1.3; pointer-events: auto; user-select: text; cursor: text; background: #fff; }
</style>
"""


def _xaxis_name(ix):
    return "xaxis" if ix == 0 else f"xaxis{ix + 1}"


def _zip_annotation_indices(menus):
    out = []
    for menu in menus:
        found = None
        for button in menu.get("buttons", []):
            if button.get("label") != "Zip":
                continue
            for key, value in button.get("args", [{}, {}])[1].items():
                if key.startswith("annotations[") and key.endswith("].visible") and value is True:
                    found = int(key.split("[")[1].split("]")[0])
                    break
        out.append(found)
    return out


def _browser_meta(fig):
    menus = [menu.to_plotly_json() for menu in fig.layout.updatemenus]
    bounds = {}
    for ix, _menu in enumerate(menus):
        axis = getattr(fig.layout, _xaxis_name(ix))
        if axis.range is not None:
            bounds[_xaxis_name(ix)] = list(axis.range)
    return {
        "menus": menus,
        "bounds": bounds,
        "titleIndices": list(range(len(menus))),
        "zipAnnotationIndices": _zip_annotation_indices(menus),
    }


def _enhancement_script(meta):
    meta_json = json.dumps(meta, cls=PlotlyJSONEncoder)
    return f"""
<script>
(function() {{
  const META = {meta_json};
  const GAP_PX = 60;
  const yaxisName = ix => ix === 0 ? 'yaxis' : 'yaxis' + (ix + 1);
  const xaxisName = ix => ix === 0 ? 'xaxis' : 'xaxis' + (ix + 1);
  const buttonFor = (ix, mode) => META.menus[ix].buttons.find(b => b.label.toLowerCase() === mode);
  function modeHeight(ix, mode) {{
    if (mode === 'zip') return 58;
    const update = buttonFor(ix, mode).args[1];
    const r = update[yaxisName(ix) + '.range'] || [0, 1];
    const span = Math.max(0.1, Math.abs(Number(r[1]) - Number(r[0])) || 1);
    if (mode === 'squish') return Math.max(20, Math.round(span / 0.25 * 8));
    return Math.max(48, Math.round(span / 0.6 * 18));
  }}
  function domainsForStates(states, gd) {{
    const fl = gd._fullLayout || {{}};
    const margin = fl.margin || {{l: 90, r: 16, t: 84, b: 36}};
    const panelHeights = states.map((mode, ix) => modeHeight(ix, mode));
    const plotH = panelHeights.reduce((a, b) => a + b, 0) + GAP_PX * Math.max(0, panelHeights.length - 1);
    const out = {{height: Math.round(margin.t + plotH + margin.b)}};
    let topPx = plotH;
    panelHeights.forEach((h, ix) => {{
      const top = topPx / plotH;
      out[yaxisName(ix) + '.domain'] = [(topPx - h) / plotH, top];
      out['annotations[' + META.titleIndices[ix] + '].y'] = top + 20 / plotH;
      topPx -= h + GAP_PX;
    }});
    return out;
  }}
  function cleanRelayout(update, ix) {{
    const out = {{}};
    const titleKeys = new Set(META.titleIndices.map(i => 'annotations[' + i + '].y'));
    const zipAnnKeys = new Set(META.zipAnnotationIndices.filter(i => i !== null).map(i => 'annotations[' + i + '].visible'));
    Object.keys(update || {{}}).forEach(k => {{
      if (k === 'height' || k.includes('.domain') || k.startsWith('updatemenus[') || titleKeys.has(k) || zipAnnKeys.has(k)) return;
      if (k.startsWith(yaxisName(ix) + '.') || k.startsWith(xaxisName(ix) + '.') || k.startsWith('annotations[') || k.startsWith('shapes[')) out[k] = update[k];
    }});
    META.zipAnnotationIndices.forEach(annIx => {{ if (annIx !== null) out['annotations[' + annIx + '].visible'] = false; }});
    return out;
  }}
  function install(gd) {{
    const shell = document.createElement('div');
    shell.className = 'pyrangeyes-browser-shell';
    gd.parentNode.insertBefore(shell, gd);
    shell.appendChild(gd);
    const states = META.menus.map(() => 'packed');
    const controls = [];
    const zipBlocks = [];
    META.menus.forEach((menu, ix) => {{
      const control = document.createElement('div');
      control.className = 'pre-view-control';
      const btn = document.createElement('button');
      btn.className = 'pre-view-button';
      btn.type = 'button';
      btn.textContent = '☰';
      const list = document.createElement('div');
      list.className = 'pre-view-menu';
      menu.buttons.forEach(button => {{
        const opt = document.createElement('button');
        opt.className = 'pre-view-option';
        opt.type = 'button';
        opt.dataset.mode = button.label.toLowerCase();
        opt.innerHTML = '<span class="check"></span>' + button.label;
        opt.addEventListener('click', () => applyMode(ix, button));
        list.appendChild(opt);
      }});
      btn.addEventListener('click', ev => {{ ev.stopPropagation(); control.classList.toggle('open'); updatePositions(); }});
      control.appendChild(btn);
      control.appendChild(list);
      shell.appendChild(control);
      controls.push(control);
      const zip = document.createElement('div');
      zip.className = 'pre-zip-block';
      ['wheel', 'mousedown', 'mousemove', 'mouseup', 'touchstart', 'touchmove', 'touchend'].forEach(type => zip.addEventListener(type, ev => ev.stopPropagation(), {{passive: true}}));
      shell.appendChild(zip);
      zipBlocks.push(zip);
    }});
    function applyMode(ix, button) {{
      states[ix] = button.label.toLowerCase();
      const relayout = Object.assign(cleanRelayout(button.args[1], ix), domainsForStates(states, gd));
      Plotly.restyle(gd, button.args[0], button.args[2]).then(() => Plotly.relayout(gd, relayout)).then(() => {{ controls[ix].classList.remove('open'); updateChecks(); updatePositions(); }});
    }}
    function updateChecks() {{
      controls.forEach((control, ix) => control.querySelectorAll('.pre-view-option').forEach(opt => opt.querySelector('.check').textContent = opt.dataset.mode === states[ix] ? '✓' : ''));
    }}
    function updatePositions() {{
      const fl = gd._fullLayout;
      if (!fl) return;
      if (fl.height) gd.style.height = Math.round(fl.height) + 'px';
      const m = fl.margin;
      const plotH = gd.clientHeight - m.t - m.b;
      const plotW = gd.clientWidth - m.l - m.r;
      controls.forEach((control, ix) => {{
        const ax = ix === 0 ? fl.yaxis : fl['yaxis' + (ix + 1)];
        const domain = ax.domain;
        const top = m.t + (1 - domain[1]) * plotH - 40;
        control.style.left = Math.round(gd.clientWidth - 40) + 'px';
        control.style.top = Math.max(4, Math.round(top)) + 'px';
        const z = zipBlocks[ix];
        const annIx = META.zipAnnotationIndices[ix];
        if (annIx !== null && gd.layout.annotations[annIx]) z.innerHTML = gd.layout.annotations[annIx].text || '';
        z.style.left = Math.round(m.l) + 'px';
        z.style.top = Math.round(m.t + (1 - domain[1]) * plotH) + 'px';
        z.style.width = Math.round(plotW) + 'px';
        z.style.height = Math.max(38, Math.round((domain[1] - domain[0]) * plotH)) + 'px';
        z.style.display = states[ix] === 'zip' ? 'flex' : 'none';
      }});
    }}
    updateChecks();
    Plotly.relayout(gd, domainsForStates(states, gd)).then(updatePositions);
    gd.on('plotly_afterplot', updatePositions);
    window.__preBrowserTest = {{states}};
  }}
  const ready = fn => document.readyState === 'loading' ? document.addEventListener('DOMContentLoaded', fn) : fn();
  ready(() => install(document.querySelector('.plotly-graph-div')));
}})();
</script>
"""


def _open_page(browser, browser_html, width=390, height=844):
    page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=2, is_mobile=width <= 500)
    page.goto(browser_html.as_uri(), wait_until="networkidle")
    page.wait_for_selector(".pre-view-button", timeout=15000)
    return page


def _click_mode(page, panel_ix, mode):
    page.locator(".pre-view-button").nth(panel_ix).click()
    option = page.locator(".pre-view-control").nth(panel_ix).locator(".pre-view-option", has_text=mode)
    option.click()
    page.wait_for_timeout(500)


def _layout_metrics(page):
    return page.evaluate(
        """
        () => {
          const gd = document.querySelector('.plotly-graph-div');
          const graph = gd.getBoundingClientRect();
          const fl = gd._fullLayout;
          const m = fl.margin;
          const plotH = fl.height - m.t - m.b;
          const axes = [fl.yaxis, fl.yaxis2, fl.yaxis3].map(ax => ({
            top: graph.y + m.t + (1 - ax.domain[1]) * plotH,
            bottom: graph.y + m.t + (1 - ax.domain[0]) * plotH,
          }));
          const titles = [...document.querySelectorAll('text.annotation-text')]
            .filter(t => t.textContent.includes('Chromosome'))
            .map(t => { const b = t.getBoundingClientRect(); return {top: b.y, bottom: b.bottom, text: t.textContent}; });
          const controls = [...document.querySelectorAll('.pre-view-control')]
            .map(c => { const b = c.getBoundingClientRect(); return {top: b.y, bottom: b.bottom, right: b.right}; });
          return {graphRight: graph.right, axes, titles, controls};
        }
        """
    )


def _assert_layout_clean(page):
    metrics = _layout_metrics(page)
    assert len(metrics["titles"]) == 3
    for ix, title in enumerate(metrics["titles"]):
        assert title["bottom"] <= metrics["axes"][ix]["top"] - 8, (ix, title, metrics["axes"][ix])
        if ix:
            assert title["top"] >= metrics["axes"][ix - 1]["bottom"] + 8, (ix, title, metrics["axes"][ix - 1])
        control = metrics["controls"][ix]
        assert abs(control["top"] - title["top"]) <= 4, (ix, control, title)
        assert metrics["graphRight"] - control["right"] <= 14, (ix, control, metrics["graphRight"])


def _panel_brick_heights(page, panel_ix):
    return page.evaluate(
        """
        (panelIx) => {
          const gd = document.querySelector('.plotly-graph-div');
          const graph = gd.getBoundingClientRect();
          const fl = gd._fullLayout;
          const m = fl.margin;
          const plotH = fl.height - m.t - m.b;
          const ax = panelIx === 0 ? fl.yaxis : fl['yaxis' + (panelIx + 1)];
          const top = graph.y + m.t + (1 - ax.domain[1]) * plotH;
          const bottom = graph.y + m.t + (1 - ax.domain[0]) * plotH;
          return [...gd.querySelectorAll('.scatterlayer .trace path')].map(p => {
            const b = p.getBoundingClientRect();
            const st = getComputedStyle(p);
            return {h: b.height, cy: b.y + b.height / 2, fill: st.fill, w: b.width};
          }).filter(o => o.w > 1 && o.h > 1 && o.cy >= top && o.cy <= bottom && o.fill !== 'rgb(255, 255, 255)' && o.fill !== 'none')
            .map(o => Math.round(o.h * 10) / 10);
        }
        """,
        panel_ix,
    )


def test_interactive_mode_cycling_keeps_titles_and_menus_outside_plots(browser, browser_html):
    page = _open_page(browser, browser_html, width=390)
    assert page.locator("g.updatemenu").count() == 0
    _assert_layout_clean(page)

    for panel_ix in range(3):
        for mode in ["Zip", "Squish", "Packed", "Full", "Packed"]:
            _click_mode(page, panel_ix, mode)
            _assert_layout_clean(page)
    page.close()


@pytest.mark.parametrize("width,height", [(390, 844), (1200, 900)])
def test_interactive_brick_heights_are_mode_consistent(browser, browser_html, width, height):
    page = _open_page(browser, browser_html, width=width, height=height)
    packed = _panel_brick_heights(page, 0)
    assert packed and max(packed) - min(packed) <= 1, packed

    _click_mode(page, 0, "Squish")
    squish = _panel_brick_heights(page, 0)
    assert squish and max(squish) < min(packed), (packed, squish)

    _click_mode(page, 0, "Full")
    full = _panel_brick_heights(page, 0)
    assert full and abs(max(full) - max(packed)) <= 1, (packed, full)
    page.close()


def test_interactive_panel_resize_is_independent(browser, browser_html):
    page = _open_page(browser, browser_html, width=390)

    before = page.evaluate(
        """
        () => {
          const fl = document.querySelector('.plotly-graph-div')._fullLayout;
          const plotH = fl.height - fl.margin.t - fl.margin.b;
          const d = fl.yaxis2.domain;
          return (d[1] - d[0]) * plotH;
        }
        """
    )
    _click_mode(page, 0, "Squish")
    after = page.evaluate(
        """
        () => {
          const fl = document.querySelector('.plotly-graph-div')._fullLayout;
          const plotH = fl.height - fl.margin.t - fl.margin.b;
          const d = fl.yaxis2.domain;
          return (d[1] - d[0]) * plotH;
        }
        """
    )
    assert abs(after - before) < 2, (before, after)
    page.close()


def test_interactive_zip_is_text_layer_not_selectable_plot(browser, browser_html):
    page = _open_page(browser, browser_html, width=390)
    _click_mode(page, 1, "Zip")
    zip_block = page.locator(".pre-zip-block").nth(1)
    from playwright.sync_api import expect

    expect(zip_block).to_be_visible()
    before = page.evaluate("document.querySelector('.plotly-graph-div')._fullLayout.xaxis2.range.slice()")
    box = zip_block.bounding_box()
    page.mouse.move(box["x"] + 8, box["y"] + box["height"] - 8)
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] - 8, box["y"] + 8, steps=8)
    page.mouse.up()
    page.wait_for_timeout(500)
    after = page.evaluate("document.querySelector('.plotly-graph-div')._fullLayout.xaxis2.range.slice()")
    assert after == before, (before, after)
    page.close()
