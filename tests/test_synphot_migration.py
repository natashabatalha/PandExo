import os

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import BlackBody

from pandexo.engine.create_input import bothTrans, outTrans
from pandexo.engine.synphot_compat import (
    calculate_bin_edges,
    load_phoenix_spectrum,
    make_array_spectrum,
    sample_spectrum_micron_mjy,
)


def _has_j_bandpass():
    refdata = os.environ.get("PYSYN_CDBS")
    if refdata is None:
        return False
    path = os.path.join(refdata, "comp", "nonhst", "bessell_j_003_syn.fits")
    return _is_readable_file(path)


def _has_phoenix_grid():
    refdata = os.environ.get("PYSYN_CDBS")
    if refdata is None:
        return False
    path = os.path.join(refdata, "grid", "phoenix", "catalog.fits")
    return _is_readable_file(path)


def _is_readable_file(path):
    if path is None or not os.path.exists(path):
        return False
    try:
        with open(path, "rb") as handle:
            handle.read(1)
    except OSError:
        return False
    return True


def test_calculate_bin_edges_linear_grid():
    edges = calculate_bin_edges([1.0, 2.0, 3.0])
    np.testing.assert_allclose(edges, [0.5, 1.5, 2.5, 3.5])


def test_calculate_bin_edges_nonuniform_grid():
    edges = calculate_bin_edges([1.0, 2.0, 4.0])
    np.testing.assert_allclose(edges, [0.5, 1.5, 3.0, 5.0])


def test_calculate_bin_edges_requires_at_least_two_centers():
    with pytest.raises(ValueError, match="at least two"):
        calculate_bin_edges([1.0])


def test_array_spectrum_samples_unsorted_um_jy_as_sorted_micron_mjy():
    spectrum = make_array_spectrum(
        wave=[1.5, 1.0, 1.25],
        flux=[3.0, 1.0, 2.0],
        wave_unit="um",
        flux_unit="Jy",
    )

    wave, flux = sample_spectrum_micron_mjy(spectrum)

    np.testing.assert_allclose(wave, [1.0, 1.25, 1.5])
    np.testing.assert_allclose(flux, [1000.0, 2000.0, 3000.0])


def test_array_spectrum_preserves_historical_fnu_to_jy_conversion():
    spectrum = make_array_spectrum(
        wave=[1.0, 2.0],
        flux=[1e-23, 2e-23],
        wave_unit="um",
        flux_unit="erg/cm2/s/Hz",
    )

    wave, flux = sample_spectrum_micron_mjy(spectrum)

    np.testing.assert_allclose(wave, [1.0, 2.0])
    np.testing.assert_allclose(flux, [1000.0, 2000.0])


def test_array_spectrum_rejects_invalid_wavelength_units():
    with pytest.raises(Exception, match="Units are not correct"):
        make_array_spectrum([1.0, 2.0], [1.0, 2.0], "meterish", "Jy")


def test_array_spectrum_rejects_invalid_flux_units():
    with pytest.raises(Exception, match="Units are not correct"):
        make_array_spectrum([1.0, 2.0], [1.0, 2.0], "um", "watts")


@pytest.mark.skipif(
    not _has_j_bandpass(),
    reason="PYSYN_CDBS J-band normalization bandpass is unavailable",
)
def test_outtrans_user_spectrum_sorts_and_returns_micron_mjy():
    star = {
        "type": "user",
        "starpath": {
            "w": np.array([1.5, 1.0, 1.25, 1.1]),
            "f": np.array([1.2, 0.9, 1.0, 1.1]),
        },
        "w_unit": "um",
        "f_unit": "Jy",
        "ref_wave": 1.25,
        "mag": 8.0,
    }

    out = outTrans(star)

    assert np.all(np.diff(out["wave"]) > 0)
    np.testing.assert_allclose(out["wave"], [1.0, 1.1, 1.25, 1.5])
    assert np.all(np.isfinite(out["flux_out_trans"]))
    assert np.nanmedian(out["flux_out_trans"]) > 1.0
    assert "stellar_flux" in out
    assert np.all(np.isfinite(out["stellar_flux"]))


@pytest.mark.skipif(
    not _has_j_bandpass(),
    reason="PYSYN_CDBS J-band normalization bandpass is unavailable",
)
def test_outtrans_phoenix_spectrum_when_reference_grid_is_available():
    star = {
        "type": "phoenix",
        "temp": 5000,
        "metal": 0.0,
        "logg": 4.5,
        "ref_wave": 1.25,
        "mag": 8.0,
    }

    try:
        out = outTrans(star)
    except Exception as exc:
        pytest.skip(f"PHOENIX reference grid is unavailable: {exc}")

    assert np.all(np.diff(out["wave"]) > 0)
    assert np.all(np.isfinite(out["flux_out_trans"]))
    assert np.all(np.isfinite(out["stellar_flux"]))


@pytest.mark.skipif(
    not _has_phoenix_grid(),
    reason="PYSYN_CDBS PHOENIX reference grid is unavailable",
)
def test_phoenix_invalid_grid_point_raises_clear_error():
    with pytest.raises(ValueError, match="PHOENIX.*no valid flux data"):
        load_phoenix_spectrum(4780.0, 0.31, 4.66)


def test_bothtrans_constant_fpfs_uses_sampled_stellar_flux():
    out_trans = {
        "wave": np.array([0.4, 0.8, 1.0, 1.5, 2.0, 2.5, 16.0]),
        "flux_out_trans": np.array([9e6, 8e6, 7e6, 6e6, 5e6, 4e6, 3e6]),
        "stellar_flux": np.array([900.0, 800.0, 700.0, 600.0, 500.0, 400.0, 300.0]),
        "phoenix": None,
    }
    planet = {
        "type": "constant",
        "f_unit": "fp/f*",
        "radius": 1.0,
        "r_unit": "R_jup",
        "temp": 1000,
    }
    star_radius = {"radius": 1.0, "r_unit": "R_sun"}

    result = bothTrans(out_trans, planet, star=star_radius)

    wave = np.array([1.0, 1.5, 2.0])
    rplan = (planet["radius"] * u.Unit(planet["r_unit"])).to(u.km)
    rstar = (star_radius["radius"] * u.Unit(star_radius["r_unit"])).to(u.km)
    bb = BlackBody(temperature=planet["temp"] * u.K)
    planet_flux = (bb(wave * u.micron) * np.pi * u.sr).to(u.mJy)
    expected = np.array(
        (planet_flux / (np.array([700.0, 600.0, 500.0]) * u.mJy))
        * (rplan / rstar) ** 2.0
    )
    wrong_flux_out_trans = np.array(
        (planet_flux / (np.array([7e6, 6e6, 5e6]) * u.mJy))
        * (rplan / rstar) ** 2.0
    )

    np.testing.assert_allclose(result["model_spec"], expected)
    assert not np.allclose(result["model_spec"], wrong_flux_out_trans)
