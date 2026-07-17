import tornado.testing

from pandexo.engine.run_online import Application


class TestCalculationRedirect(tornado.testing.AsyncHTTPTestCase):
    def get_app(self):
        return Application()

    def test_calculation_redirects_to_dashboard(self):
        response = self.fetch("/calculation", follow_redirects=False)

        assert response.code == 301
        assert response.headers["Location"] == "/dashboard"
