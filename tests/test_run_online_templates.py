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
