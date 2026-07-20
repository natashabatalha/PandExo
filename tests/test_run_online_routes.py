import tornado.testing
from urllib.parse import urljoin

from pandexo.engine.run_online import Application


class TestCalculationRedirect(tornado.testing.AsyncHTTPTestCase):
    def get_app(self):
        return Application()

    def test_calculation_redirects_to_dashboard(self):
        for path, location in (
            ("/calculation", "dashboard"),
            ("/calculation/", "../dashboard"),
        ):
            response = self.fetch(path, follow_redirects=False)

            assert response.code == 301
            assert response.headers["Location"] == location
            assert (
                urljoin(
                    f"https://exoctk-dev.stsci.edu/pandexo{path}", location
                )
                == "https://exoctk-dev.stsci.edu/pandexo/dashboard"
            )
