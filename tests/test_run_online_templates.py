from pathlib import Path

import pytest
import tornado.web


pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:pkg_resources is deprecated as an API:UserWarning:pandexo\\.engine\\.utils"
    ),
    pytest.mark.filterwarnings(
        "ignore:Deprecated call to `pkg_resources\\.declare_namespace.*:DeprecationWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:The notation 'flux\\([HJ]\\)' is deprecated.*:DeprecationWarning"
    ),
]


def test_base_handler_render_supplies_bokeh_resources(monkeypatch):
    from pandexo.engine.run_online import BaseHandler

    captured = {}

    def fake_render(self, template_name, **kwargs):
        captured["template_name"] = template_name
        captured["kwargs"] = kwargs
        return "rendered"

    monkeypatch.setattr(tornado.web.RequestHandler, "render", fake_render)

    handler = object.__new__(BaseHandler)

    assert handler.render("new.html") == "rendered"
    assert captured["template_name"] == "new.html"
    assert "bokeh_css" in captured["kwargs"]
    assert "bokeh_js" in captured["kwargs"]


def test_base_handler_render_preserves_explicit_bokeh_resources(monkeypatch):
    from pandexo.engine.run_online import BaseHandler

    captured = {}

    def fake_render(self, template_name, **kwargs):
        captured["kwargs"] = kwargs

    monkeypatch.setattr(tornado.web.RequestHandler, "render", fake_render)

    handler = object.__new__(BaseHandler)
    handler.render("new.html", bokeh_css="custom-css", bokeh_js="custom-js")

    assert captured["kwargs"]["bokeh_css"] == "custom-css"
    assert captured["kwargs"]["bokeh_js"] == "custom-js"


def test_new_calculation_template_lists_new_miri_lrs_subarrays():
    template = Path("pandexo/engine/templates/new.html").read_text()

    assert 'value="slitlessprism_ip"' in template
    assert 'value="slitlessprism_ips"' in template
    assert 'SLITLESSPRISM (Not Recommended' in template
    assert 'value="full"' in template
    assert "FULL (tframe=2.775)" in template
    assert 'value="subslit"' in template
    assert "miriSubarraysByMode" in template
    assert 'lrsslit: ["full", "subslit"]' in template


def test_new_calculation_template_lists_nirspec_prism_multistripe_subarrays():
    template = Path("pandexo/engine/templates/new.html").read_text()

    assert 'value="s256m2_prm"' in template
    assert 'value="s128m4_prm"' in template
    assert 'value="s64m8_prm"' in template
    assert 'value="s32m16_prm"' in template


def test_miri_reference_uses_supported_fastr1_readout():
    reference = Path("pandexo/engine/reference/miri_input.json").read_text()

    assert '"readout_pattern":"fastr1"' in reference
    assert '"readmode": "fastr1"' in reference
    assert '"readout_pattern":"fast"' not in reference
    assert '"readmode": "fast"' not in reference
