import numpy as np
import pytest

from pandexo.engine.utils.plot_helpers import add_wavelength_gap_breaks


def test_add_wavelength_gap_breaks_inserts_nan_after_large_jump():
    x, y = add_wavelength_gap_breaks(
        np.array([3.70, 3.71, 3.72, 3.84, 3.85]),
        np.array([10.0, 11.0, 12.0, 20.0, 21.0]),
    )

    assert np.isnan(y[3])
    assert x[3] == pytest.approx(3.78)
    np.testing.assert_allclose(x, [3.70, 3.71, 3.72, 3.78, 3.84, 3.85])


def test_add_wavelength_gap_breaks_leaves_regular_grid_unchanged():
    x = np.array([3.70, 3.71, 3.72, 3.73])
    y = np.array([10.0, 11.0, 12.0, 13.0])

    xout, yout = add_wavelength_gap_breaks(x, y)

    np.testing.assert_array_equal(xout, x)
    np.testing.assert_array_equal(yout, y)


def test_add_wavelength_gap_breaks_allows_smooth_dispersion_change():
    steps = np.array([
        0.010,
        0.011,
        0.013,
        0.017,
        0.025,
        0.040,
        0.055,
        0.070,
        0.085,
    ])
    x = np.concatenate(([5.0], 5.0 + np.cumsum(steps)))
    y = np.arange(x.size, dtype=float)

    xout, yout = add_wavelength_gap_breaks(x, y)

    np.testing.assert_array_equal(xout, x)
    np.testing.assert_array_equal(yout, y)
