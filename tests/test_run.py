import contextlib
import io

import numpy as np
import pytest


def _import_compute_timing():
    try:
        from pandexo.engine.jwst import compute_timing
    except Exception as exc:
        pytest.skip(f"Pandeia-dependent JWST module is unavailable: {exc}")
    return compute_timing


def _import_justdoit():
    try:
        import pandexo.engine.justdoit as jdi
    except Exception as exc:
        pytest.skip(f"Pandeia-dependent run machinery is unavailable: {exc}")
    return jdi


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
            "Pandeia reference data is invalid for the installed pandeia.engine. "
            "pandeia.engine 2026.2 requires matching 2026.2 refdata with detector "
            "subarray nsuperstripe entries.\n"
            + status
        )


def _default_smoke_exo_dict(jdi):
    exo_dict = jdi.load_exo_dict()
    exo_dict["observation"]["sat_level"] = 80
    exo_dict["observation"]["sat_unit"] = "%"
    exo_dict["observation"]["noccultations"] = 2
    exo_dict["observation"]["R"] = None
    exo_dict["observation"]["baseline"] = 1.0
    exo_dict["observation"]["baseline_unit"] = "frac"
    exo_dict["observation"]["noise_floor"] = 0
    exo_dict["star"]["type"] = "phoenix"
    exo_dict["star"]["mag"] = 8.0
    exo_dict["star"]["ref_wave"] = 1.25
    exo_dict["star"]["temp"] = 5500
    exo_dict["star"]["metal"] = 0.0
    exo_dict["star"]["logg"] = 4.0
    exo_dict["star"]["radius"] = 1
    exo_dict["star"]["r_unit"] = "R_sun"
    exo_dict["planet"]["type"] = "constant"
    exo_dict["planet"]["radius"] = 1
    exo_dict["planet"]["r_unit"] = "R_jup"
    exo_dict["planet"]["transit_duration"] = 2.0 * 60.0 * 60.0
    exo_dict["planet"]["td_unit"] = "s"
    exo_dict["planet"]["f_unit"] = "rp^2/r*^2"
    return exo_dict


def test_compute_timing_keeps_first_minus_last_on_source_time_positive():
    compute_timing = _import_compute_timing()

    timing, flags = compute_timing(
        {
            "maxexptime_per_int": 0.1,
            "tframe": 1.0,
            "nframe": 1,
            "mingroups": 1,
            "nskip": 0,
        },
        transit_duration=10.0,
        expfact_out=1.0,
        noccultations=1,
    )

    assert timing["APT: Num Groups per Integration"] == 2
    assert timing["Num Integrations In Transit"] > 0
    assert (
        timing["Seconds per Frame"]
        * (timing["APT: Num Groups per Integration"] - 1)
        * timing["Num Integrations In Transit"]
        > 0
    )
    assert flags["flag_default"] == "NGROUPS<2SET TO NGROUPS=2"


@pytest.mark.parametrize("instrument", ["NIRSpec G140H", "MIRI LRS"])
def test_run_pandexo_smoke_has_sorted_wavelengths(instrument):
    _require_valid_pandeia_refdata()
    jdi = _import_justdoit()

    result = jdi.run_pandexo(
        _default_smoke_exo_dict(jdi),
        [instrument],
        save_file=False,
    )

    wave = np.asarray(result["FinalSpectrum"]["wave"])
    assert np.all(np.diff(wave) >= 0), f"{instrument} produced unsorted wavelengths"
