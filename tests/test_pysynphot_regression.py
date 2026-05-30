import copy
import contextlib
import io
import json
import os
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import BlackBody

from pandexo.engine.create_input import bothTrans, outTrans

pytestmark = pytest.mark.filterwarnings("ignore:.*alltrue.*:DeprecationWarning")


USER_CASES = {
    "user_j_um_jy": {
        "type": "user",
        "starpath": {"w": [1.5, 1.0, 1.25, 1.1], "f": [1.2, 0.9, 1.0, 1.1]},
        "w_unit": "um",
        "f_unit": "Jy",
        "ref_wave": 1.25,
        "mag": 8.0,
    },
    "user_h_um_jy": {
        "type": "user",
        "starpath": {"w": [1.9, 1.45, 1.7, 1.6], "f": [1.8, 1.1, 1.5, 1.35]},
        "w_unit": "um",
        "f_unit": "Jy",
        "ref_wave": 1.65,
        "mag": 8.0,
    },
    "user_k_um_jy": {
        "type": "user",
        "starpath": {"w": [2.6, 1.9, 2.3, 2.1], "f": [2.4, 1.5, 2.0, 1.8]},
        "w_unit": "um",
        "f_unit": "Jy",
        "ref_wave": 2.2,
        "mag": 8.0,
    },
    "user_j_nm_fnu": {
        "type": "user",
        "starpath": {"w": [1500.0, 1000.0, 1250.0, 1100.0], "f": [1.2e-23, 0.9e-23, 1.0e-23, 1.1e-23]},
        "w_unit": "nm",
        "f_unit": "erg/cm2/s/Hz",
        "ref_wave": 1.25,
        "mag": 8.0,
    },
    "user_h_angs_flam": {
        "type": "user",
        "starpath": {"w": [19000.0, 14500.0, 17000.0, 16000.0], "f": [8.0e-12, 5.0e-12, 7.0e-12, 6.0e-12]},
        "w_unit": "Angs",
        "f_unit": "FLAM",
        "ref_wave": 1.65,
        "mag": 8.0,
    },
}

PHOENIX_CASE = {
    "type": "phoenix",
    "temp": 5000,
    "metal": 0.0,
    "logg": 4.5,
    "ref_wave": 1.25,
    "mag": 8.0,
}

PHOENIX_GRID_CASES = {
    "phoenix_j_5000_m0_g45": PHOENIX_CASE,
    "phoenix_h_3500_m0_g45": {
        "type": "phoenix",
        "temp": 3500,
        "metal": 0.0,
        "logg": 4.5,
        "ref_wave": 1.65,
        "mag": 8.0,
    },
    "phoenix_k_6000_m0_g45": {
        "type": "phoenix",
        "temp": 6000,
        "metal": 0.0,
        "logg": 4.5,
        "ref_wave": 2.2,
        "mag": 8.0,
    },
}

FPFS_PLANET = {
    "type": "constant",
    "f_unit": "fp/f*",
    "radius": 1.0,
    "r_unit": "R_jup",
    "temp": 1000,
}
FPFS_STAR = {"radius": 1.0, "r_unit": "R_sun"}


def _require_valid_pandeia_refdata():
    try:
        import pandeia.engine
    except Exception as exc:
        pytest.skip(f"Pandeia cannot be imported without configured refdata: {exc}")

    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        pandeia.engine.pandeia_version()
    status = output.getvalue()
    invalid_pandeia = (
        "Pandeia RefData version: INVALID INSTALLATION" in status
        or "Pandeia RefData version: ENVIRONMENT VARIABLE UNSET" in status
        or "Pandeia PSFs version:    INVALID INSTALLATION" in status
        or "Pandeia PSFs version:    ENVIRONMENT VARIABLE UNSET" in status
    )
    if invalid_pandeia:
        pytest.skip(
            "Pandeia reference data is invalid for the installed pandeia.engine.\n"
            + status
        )


def _default_smoke_exo_dict(star):
    try:
        import pandexo.engine.justdoit as jdi
    except Exception as exc:
        pytest.skip(f"Pandeia-dependent run machinery is unavailable: {exc}")

    exo_dict = jdi.load_exo_dict()
    exo_dict["observation"]["sat_level"] = 80
    exo_dict["observation"]["sat_unit"] = "%"
    exo_dict["observation"]["noccultations"] = 2
    exo_dict["observation"]["R"] = None
    exo_dict["observation"]["baseline"] = 1.0
    exo_dict["observation"]["baseline_unit"] = "frac"
    exo_dict["observation"]["noise_floor"] = 0
    exo_dict["star"].update(copy.deepcopy(star))
    exo_dict["star"]["radius"] = 1
    exo_dict["star"]["r_unit"] = "R_sun"
    exo_dict["planet"]["type"] = "constant"
    exo_dict["planet"]["radius"] = 1
    exo_dict["planet"]["r_unit"] = "R_jup"
    exo_dict["planet"]["transit_duration"] = 2.0 * 60.0 * 60.0
    exo_dict["planet"]["td_unit"] = "s"
    exo_dict["planet"]["f_unit"] = "rp^2/r*^2"
    return exo_dict


def _case_to_numpy(case):
    case = copy.deepcopy(case)
    if case.get("type") == "user":
        case["starpath"] = {
            "w": np.asarray(case["starpath"]["w"], dtype=float),
            "f": np.asarray(case["starpath"]["f"], dtype=float),
        }
    return case


def _filter_name(ref_wave):
    ref_wave = float(ref_wave)
    if 1.2 <= ref_wave <= 1.3:
        return "J"
    if 1.6 <= ref_wave <= 1.7:
        return "H"
    if 2.1 <= ref_wave <= 2.3:
        return "K"
    raise ValueError("Only J H and K zeropoints are included")


def _bandpass_path(case):
    refdata = os.environ.get("PYSYN_CDBS")
    if refdata is None:
        return None
    filenames = {
        "J": "bessell_j_003_syn.fits",
        "H": "bessell_h_004_syn.fits",
        "K": "bessell_k_003_syn.fits",
    }
    return os.path.join(refdata, "comp", "nonhst", filenames[_filter_name(case["ref_wave"])])


def _skip_without_bandpass(case):
    path = _bandpass_path(case)
    if path is None or not os.path.exists(path):
        pytest.skip(f"required normalization bandpass is unavailable: {path}")


def _legacy_outtrans(case):
    _skip_without_bandpass(case)
    psyn = pytest.importorskip("pysynphot")
    case = _case_to_numpy(case)

    if case["type"] == "user":
        wave = np.asarray(case["starpath"]["w"], dtype=float)
        flux = np.asarray(case["starpath"]["f"], dtype=float)
        order = np.argsort(wave)
        wave = wave[order]
        flux = flux[order]

        waveunits = {
            "um": "um",
            "nm": "nm",
            "cm": "cm",
            "Angs": "angstrom",
            "Hz": "Hz",
        }[case["w_unit"]]

        if case["f_unit"] == "Jy":
            fluxunits = "jy"
        elif case["f_unit"] == "FLAM":
            fluxunits = "FLAM"
        elif case["f_unit"] == "erg/cm2/s/Hz":
            flux = flux * 1e23
            fluxunits = "jy"
        else:
            raise ValueError(case["f_unit"])

        spectrum = psyn.ArraySpectrum(
            wave,
            flux,
            waveunits=waveunits,
            fluxunits=fluxunits,
        )
        spectrum.convert("nm")
        spectrum.convert("jy")
    elif case["type"] == "phoenix":
        try:
            spectrum = psyn.Icat(
                "phoenix",
                case["temp"],
                min(case["metal"], 0.5),
                case["logg"],
            )
        except Exception as exc:
            pytest.skip(f"PHOENIX reference grid is unavailable: {exc}")
        spectrum.convert("nm")
        spectrum.convert("jy")
    else:
        raise ValueError(case["type"])

    bandpass = psyn.FileBandpass(_bandpass_path(case))
    spectrum.convert("angstroms")
    bandpass.convert("angstroms")
    normalized = spectrum.renorm(float(case["mag"]), "vegamag", bandpass)
    normalized.convert("microns")
    normalized.convert("mjy")

    return {
        "flux_out_trans": np.asarray(normalized.flux, dtype=float),
        "wave": np.asarray(normalized.wave, dtype=float),
        "phoenix": spectrum,
        "stellar_flux": (np.asarray(spectrum.flux, dtype=float) * u.Jy).to(u.mJy).value,
        "stellar_wave": np.asarray(normalized.wave, dtype=float),
    }


def _legacy_bothtrans_constant_fpfs(out_trans, planet, star):
    planet = copy.deepcopy(planet)
    rplan = (planet["radius"] * u.Unit(planet["r_unit"])).to(u.km)
    rstar = (star["radius"] * u.Unit(star["r_unit"])).to(u.km)
    mask = (out_trans["wave"] > 0.5) & (out_trans["wave"] < 15)
    wave_planet = out_trans["wave"][mask]
    flux_star = (out_trans["phoenix"].flux * u.Jy).to(u.mJy)[mask]

    bb = BlackBody(temperature=planet["temp"] * u.K)
    flux_planet = (bb(wave_planet * u.micron) * np.pi * u.sr).to(u.mJy)
    flux_planet = np.asarray((flux_planet / flux_star) * (rplan / rstar) ** 2.0)

    wave_star = out_trans["wave"]
    flux_star_out = out_trans["flux_out_trans"]
    wavemin = max([min(wave_planet), min(wave_star), 0.5])
    wavemax = min([max(wave_planet), max(wave_star), 15])
    trim = (wave_planet > wavemin) & (wave_planet < wavemax)
    wave_planet = wave_planet[trim]
    flux_planet = flux_planet[trim]
    flux_out_trans = np.interp(wave_planet, wave_star, flux_star_out)
    depth_fraction = 1.0 + flux_planet
    flux_in_trans = flux_out_trans * depth_fraction
    return {
        "wave": wave_planet,
        "flux_in_trans": flux_in_trans,
        "flux_out_trans": flux_out_trans,
        "model_wave": wave_planet,
        "model_spec": flux_planet,
        "frac": depth_fraction,
    }


def _max_relative_difference(new, old):
    new = np.asarray(new, dtype=float)
    old = np.asarray(old, dtype=float)
    denom = np.maximum(np.abs(old), 1e-300)
    return float(np.max(np.abs(new - old) / denom))


@pytest.mark.parametrize("name,case", USER_CASES.items())
def test_live_user_outtrans_matches_legacy_pysynphot(name, case):
    legacy = _legacy_outtrans(case)
    new = outTrans(_case_to_numpy(case))

    np.testing.assert_allclose(new["wave"], legacy["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(new["flux_out_trans"], legacy["flux_out_trans"], rtol=1e-4, atol=1e-8)
    np.testing.assert_allclose(new["stellar_flux"], legacy["stellar_flux"], rtol=1e-10, atol=1e-8)


@pytest.mark.parametrize("name,case", PHOENIX_GRID_CASES.items())
def test_live_phoenix_outtrans_matches_legacy_pysynphot(name, case):
    legacy = _legacy_outtrans(case)
    new = outTrans(copy.deepcopy(case))

    np.testing.assert_allclose(new["wave"], legacy["wave"], rtol=0, atol=1e-10)
    # stsynphot and pysynphot load the same CDBS grid, but their normalization
    # and unit conversion paths differ slightly across the very large PHOENIX
    # wavelength range.
    np.testing.assert_allclose(new["flux_out_trans"], legacy["flux_out_trans"], rtol=1e-3, atol=1e-18)
    np.testing.assert_allclose(new["stellar_flux"], legacy["stellar_flux"], rtol=1e-3, atol=1e-18)


def test_live_constant_fpfs_matches_legacy_pysynphot():
    case = USER_CASES["user_j_um_jy"]
    legacy_out = _legacy_outtrans(case)
    new_out = outTrans(_case_to_numpy(case))

    legacy = _legacy_bothtrans_constant_fpfs(legacy_out, FPFS_PLANET, FPFS_STAR)
    new = bothTrans(new_out, copy.deepcopy(FPFS_PLANET), star=copy.deepcopy(FPFS_STAR))

    np.testing.assert_allclose(new["wave"], legacy["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(new["model_spec"], legacy["model_spec"], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(new["flux_in_trans"], legacy["flux_in_trans"], rtol=1e-4, atol=1e-8)


def test_live_nirspec_precision_matches_legacy_pysynphot_outtrans(monkeypatch):
    _require_valid_pandeia_refdata()
    legacy_out = _legacy_outtrans(PHOENIX_CASE)

    try:
        import pandexo.engine.justdoit as jdi
        import pandexo.engine.jwst as jwst
    except Exception as exc:
        pytest.skip(f"Pandeia-dependent run machinery is unavailable: {exc}")

    exo_dict = _default_smoke_exo_dict(PHOENIX_CASE)
    new = jdi.run_pandexo(copy.deepcopy(exo_dict), ["NIRSpec G140H"], save_file=False)

    def legacy_outtrans_for_jwst(star):
        assert star["type"] == PHOENIX_CASE["type"]
        assert star["temp"] == PHOENIX_CASE["temp"]
        assert star["metal"] == PHOENIX_CASE["metal"]
        assert star["logg"] == PHOENIX_CASE["logg"]
        assert star["ref_wave"] == PHOENIX_CASE["ref_wave"]
        assert star["mag"] == PHOENIX_CASE["mag"]
        return legacy_out

    monkeypatch.setattr(jwst.create, "outTrans", legacy_outtrans_for_jwst)
    legacy = jdi.run_pandexo(copy.deepcopy(exo_dict), ["NIRSpec G140H"], save_file=False)

    # This is the closest regression to the proposer-facing quantity: only the
    # stellar-spectrum backend changes, while the full PandExo/Pandeia
    # simulation and precision calculation are exercised.
    np.testing.assert_allclose(
        new["FinalSpectrum"]["wave"],
        legacy["FinalSpectrum"]["wave"],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        new["FinalSpectrum"]["spectrum"],
        legacy["FinalSpectrum"]["spectrum"],
        rtol=1e-4,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        new["FinalSpectrum"]["error_w_floor"],
        legacy["FinalSpectrum"]["error_w_floor"],
        rtol=1e-4,
        atol=1e-12,
    )


def _load_frozen_reference():
    fixture_path = Path(__file__).parent / "data" / "pysynphot_legacy_reference.json"
    with fixture_path.open() as handle:
        return json.load(handle)


@pytest.mark.parametrize("name", USER_CASES)
def test_frozen_user_outtrans_regression(name):
    reference = _load_frozen_reference()["outtrans"][name]
    case = _case_to_numpy(reference["input"])
    _skip_without_bandpass(case)

    new = outTrans(case)

    np.testing.assert_allclose(new["wave"], reference["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(new["flux_out_trans"], reference["flux_out_trans"], rtol=1e-4, atol=1e-8)
    np.testing.assert_allclose(new["stellar_flux"], reference["stellar_flux"], rtol=1e-10, atol=1e-8)


def test_frozen_constant_fpfs_regression():
    reference = _load_frozen_reference()["constant_fpfs"]
    case = _case_to_numpy(reference["outtrans_input"])
    _skip_without_bandpass(case)

    new = bothTrans(
        outTrans(case),
        copy.deepcopy(reference["planet"]),
        star=copy.deepcopy(reference["star"]),
    )

    np.testing.assert_allclose(new["wave"], reference["wave"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(new["model_spec"], reference["model_spec"], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(new["flux_in_trans"], reference["flux_in_trans"], rtol=1e-4, atol=1e-8)


def test_frozen_reference_records_observed_legacy_differences():
    reference = _load_frozen_reference()
    assert reference["metadata"]["generator"]
    for case in reference["outtrans"].values():
        assert case["max_relative_difference"]["flux_out_trans"] < 1e-4
        assert case["max_relative_difference"]["stellar_flux"] < 1e-10
