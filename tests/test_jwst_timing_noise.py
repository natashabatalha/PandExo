"""JWST timing, calculation-routing, and noise-scaling regression tests."""

import contextlib
import io
import json
from copy import deepcopy

import numpy as np
import pytest

from pandexo.engine.compute_noise import ExtractSpec
from pandexo.engine.jwst import (
    build_timing_display_div,
    compute_timing,
    select_calculation,
    update_timing_measurement_time,
)


def _timing(nsuperstripe, ngroup=3, mingroups=2):
    timing, _ = compute_timing(
        {
            "ngroup": ngroup,
            "tframe": 2.0,
            "nframe": 1,
            "mingroups": mingroups,
            "nskip": 0,
            "nsuperstripe": nsuperstripe,
        },
        transit_duration=24.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )
    return timing


def _timing_with_pandeia_cycle(nsuperstripe, exposure_time_per_int, ngroup=3):
    timing, _ = compute_timing(
        {
            "ngroup": ngroup,
            "tframe": 2.0,
            "nframe": 1,
            "mingroups": 1,
            "nskip": 0,
            "nsuperstripe": nsuperstripe,
            "exposure_time_per_int": exposure_time_per_int,
            "exposure_time_ngroup": ngroup,
        },
        transit_duration=48.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )
    return timing


def _pandeia_out(instrument=None, detector=None):
    return {
        "input": {
            "configuration": {
                "instrument": instrument or {
                    "instrument": "niriss",
                    "mode": "soss",
                    "filter": "clear",
                    "aperture": "soss",
                    "disperser": "gr700xd",
                },
                "detector": detector or {
                    "subarray": "substrip96",
                    "readout_pattern": "nisrapid",
                },
            }
        }
    }


def _fml_result(timing):
    wave = np.array([1.0, 2.0])
    flux = np.array([100.0, 100.0])
    bg = np.array([0.0, 0.0])
    out = {
        "scalar": {"measurement_time": 10.0},
        "1d": {
            "extracted_flux": [wave, flux],
            "extracted_bg_only": [wave, bg],
        },
    }
    inn = {
        "scalar": {"measurement_time": 10.0},
        "1d": {
            "extracted_flux": [wave, flux],
            "extracted_bg_only": [wave, bg],
        },
    }
    return ExtractSpec(inn, out, rn=0.0, extraction_area=1.0, timing=timing).run_f_minus_l()


def _slope_result(timing):
    wave = np.array([1.0, 2.0])
    flux = np.array([100.0, 100.0])
    noise = np.array([10.0, 10.0])
    bg = np.array([0.0, 0.0])
    out = {
        "scalar": {"measurement_time": 10.0},
        "1d": {
            "extracted_flux": [wave, flux],
            "extracted_noise": [wave, noise],
            "extracted_bg_only": [wave, bg],
        },
    }
    inn = deepcopy(out)
    return ExtractSpec(inn, out, rn=0.0, extraction_area=1.0, timing=timing).run_slope_method()


def _phase_result(timing, method):
    wave = np.array([1.0, 2.0])
    flux = np.array([100.0, 200.0])
    noise = np.array([3.0, 4.0])
    bg = np.array([1.0, 2.0])
    out = {
        "scalar": {"measurement_time": 12.0},
        "1d": {
            "extracted_flux": [wave, flux],
            "extracted_noise": [wave, noise],
            "extracted_bg_only": [wave, bg],
        },
    }
    inn = {
        "time": np.array([0.0, 24.0]),
        "planet_phase": np.array([0.0, 0.01]),
    }
    extract = ExtractSpec(inn, out, rn=2.0, extraction_area=3.0, timing=timing)
    return getattr(extract, method)()


def _fractional_uncertainty(result):
    photon_out = result["photon_out_1d"]
    photon_in = result["photon_in_1d"]
    var_out = result["var_out_1d"]
    var_in = result["var_in_1d"]
    to = result["on_source_out"]
    ti = result["on_source_in"]
    var_tot = (
        (to / ti / photon_out) ** 2.0 * var_in
        + (photon_in * to / ti / photon_out ** 2.0) ** 2.0 * var_out
    )
    return np.sqrt(var_tot)


@pytest.mark.parametrize(
    ("planet_wave_unit", "nsuperstripe", "is_dhs", "expected"),
    [
        ("um", 1, False, "fml"),
        ("um", 3, False, "slope method"),
        ("um", 1, True, "slope method"),
        ("sec", 1, False, "phase_spec_fml"),
        ("sec", 3, False, "phase_spec_slope"),
        ("sec", 1, True, "phase_spec_slope"),
    ],
)
def test_calculation_depends_on_detector_readout(
    planet_wave_unit, nsuperstripe, is_dhs, expected
):
    assert select_calculation(planet_wave_unit, nsuperstripe, is_dhs) == expected


def test_nsuperstripe_one_preserves_effective_timing():
    timing = _timing(nsuperstripe=1)

    assert timing["Num Superstripes"] == 1
    assert timing["Effective Integrations In Transit"] == timing["Num Integrations In Transit"]
    assert timing["Effective Integrations Out of Transit"] == timing["Num Integrations Out of Transit"]
    assert timing["Effective On Source Time In Transit"] == timing["On Source Time In Transit"]
    assert timing["Effective On Source Time Out of Transit"] == timing["On Source Time Out of Transit"]


def test_single_group_timing_uses_one_positive_measurement_frame():
    timing = _timing(nsuperstripe=1, ngroup=1, mingroups=1)

    assert timing["APT: Num Groups per Integration"] == 1
    assert timing["Zero Frame Efficiency Loss"] == 0
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(2.0)
    assert timing["On Source Time In Transit"] > 0
    assert timing["Effective On Source Time In Transit"] > 0


def test_compute_timing_keeps_single_group_on_source_time_positive():
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

    assert timing["APT: Num Groups per Integration"] == 1
    assert timing["Num Integrations In Transit"] > 0
    assert (
        timing["Seconds per Frame"]
        * (
            timing["APT: Num Groups per Integration"]
            + timing["Zero Frame Efficiency Loss"]
        )
        * timing["Num Integrations In Transit"]
        > 0
    )
    assert flags["flag_default"] == "Optimized NGROUPS below minimum (1). SET TO NGROUPS=1"


def test_multigroup_measurement_time_preserves_first_minus_last_interval():
    timing = _timing(nsuperstripe=1, ngroup=3)

    assert timing["Zero Frame Efficiency Loss"] == -1
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(4.0)


def test_multistripe_effective_time_is_divided_by_nsuperstripe():
    timing = _timing(nsuperstripe=4)

    assert timing["Effective Integrations In Transit"] == pytest.approx(
        timing["Num Integrations In Transit"] / 4.0
    )
    assert timing["Effective Integrations Out of Transit"] == pytest.approx(
        timing["Num Integrations Out of Transit"] / 4.0
    )
    assert timing["On Source Time In Transit"] == pytest.approx(
        timing["Effective Integrations In Transit"]
        * timing["Measurement Time per Integration (sec)"]
    )
    assert timing["Effective On Source Time In Transit"] == pytest.approx(
        timing["On Source Time In Transit"]
    )
    assert timing["Effective On Source Time Out of Transit"] == pytest.approx(
        timing["On Source Time Out of Transit"]
    )


def test_pandeia_measurement_time_update_controls_on_source_metadata():
    timing = _timing(nsuperstripe=4)
    update_timing_measurement_time(timing, measurement_time_per_int=32.0)

    assert timing["Measurement Time per Integration (sec)"] == 32.0
    assert timing["On Source Time In Transit"] == pytest.approx(
        timing["Effective Integrations In Transit"] * 32.0
    )
    assert timing["Effective On Source Time In Transit"] == pytest.approx(
        timing["Effective Integrations In Transit"] * 32.0
    )


def test_multistripe_timing_uses_pandeia_full_cycle_clock_for_real_integrations():
    timing = _timing_with_pandeia_cycle(
        nsuperstripe=4,
        exposure_time_per_int=40.0,
        ngroup=3,
    )

    assert timing["Time/Integration incl reset (sec)"] == pytest.approx(40.0)
    assert timing["Num Integrations In Transit"] == 2
    assert timing["Effective Integrations In Transit"] == pytest.approx(0.5)
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(16.0)
    assert timing["On Source Time In Transit"] == pytest.approx(8.0)


def test_dhs_clock_time_includes_pandeia_fixed_overhead():
    timing, _ = compute_timing(
        {
            "ngroup": 30,
            "tframe": 1.36765,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 1,
            "exposure_time_per_int": 42.41763,
            "exposure_time_ngroup": 30,
        },
        transit_duration=2.8032 * 3600.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=101,
    )

    assert timing["Time/Integration incl reset (sec)"] == pytest.approx(42.41763)
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(39.66185)


def test_soss_sub17stripe_clock_time_matches_pandeia_full_cycle():
    timing, _ = compute_timing(
        {
            "ngroup": 8,
            "tframe": 0.06164,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 120,
            "exposure_time_per_int": 69.0288,
            "exposure_time_ngroup": 8,
        },
        transit_duration=2.8032 * 3600.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )

    assert timing["Time/Integration incl reset (sec)"] == pytest.approx(69.0288)
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(51.7776)
    assert timing["Measurement Time per Integration (sec)"] / timing[
        "Num Superstripes"
    ] == pytest.approx(0.43148)
    assert timing["Num Integrations In Transit"] == 147
    assert timing["Effective Integrations In Transit"] == pytest.approx(147 / 120.0)


def test_multistripe_timing_display_uses_apt_and_calculation_tables():
    timing = _timing_with_pandeia_cycle(
        nsuperstripe=4,
        exposure_time_per_int=40.0,
        ngroup=3,
    )
    html = build_timing_display_div(_pandeia_out(), timing).decode()

    assert "APT Inputs" in html
    assert "Calculation Details" in html
    assert "Groups per Integration" in html
    assert "Integrations per Occultation" in html
    assert "Number of Stripes" in html
    assert "Elapsed Time per APT Integration incl. Reset (sec)" in html
    assert "Science Time per Full Multistripe Cycle excl. Reset (sec)" in html
    assert "Science Time per Stripe excl. Reset (sec)" in html
    assert "Effective Per-Wavelength" not in html
    assert "Time/Integration incl reset" not in html
    assert "Measurement Time per Integration" not in html


def test_non_multistripe_timing_display_omits_stripe_rows():
    timing = _timing(nsuperstripe=1)
    html = build_timing_display_div(_pandeia_out(), timing).decode()

    assert "Elapsed Time per Integration incl. Reset (sec)" in html
    assert "Science Time per Integration excl. Reset (sec)" in html
    assert "Number of Stripes" not in html
    assert "Full Multistripe Cycle" not in html
    assert "Science Time per Stripe" not in html


def test_nircam_timing_display_shows_channel_and_pupil_rows():
    timing = _timing(nsuperstripe=1)
    html = build_timing_display_div(
        _pandeia_out(
            instrument={
                "instrument": "nircam",
                "mode": "sw_tsgrism",
                "filter": "f150w2",
                "pandexofilterpair": "f322w2",
                "aperture": "dhs0spec8",
                "disperser": "dhs0",
            },
            detector={
                "subarray": "sub260s4_8-spectra",
                "readout_pattern": "rapid",
            },
        ),
        timing,
    ).decode()

    assert "SW Channel Mode" in html
    assert "GRISM" in html
    assert "SUB260S4_8-SPECTRA" in html
    assert "No. of Output Channels" in html
    assert "<td>4</td>" in html
    assert "Short Pupil+Filter" in html
    assert "GDHS0+F150W2" in html
    assert "Long Pupil+Filter" in html
    assert "GRISMR+F322W2" in html


def test_nircam_lw_only_timing_display_shows_imaging_sw_channel():
    timing = _timing(nsuperstripe=1)
    html = build_timing_display_div(
        _pandeia_out(
            instrument={
                "instrument": "nircam",
                "mode": "lw_tsgrism",
                "filter": "f322w2",
                "aperture": "lw",
                "disperser": "grismr",
            },
            detector={
                "subarray": "subgrism64",
                "readout_pattern": "rapid",
            },
        ),
        timing,
    ).decode()

    assert "SW Channel Mode" in html
    assert "IMAGING" in html
    assert "SUBGRISM64" in html
    assert "No. of Output Channels" in html
    assert "<td>4</td>" in html
    assert "Short Pupil+Filter" in html
    assert "CHOOSE THIS USING ETC" in html
    assert "Long Pupil+Filter" in html
    assert "GRISMR+F322W2" in html


def test_slope_uncertainty_increases_by_sqrt_nsuperstripe():
    nsuperstripe = 9
    normal = _slope_result(_timing(nsuperstripe=1))
    multistripe = _slope_result(_timing(nsuperstripe=nsuperstripe))

    ratio = _fractional_uncertainty(multistripe) / _fractional_uncertainty(normal)

    assert ratio == pytest.approx(np.sqrt(nsuperstripe))
    assert multistripe["on_source_in"] == pytest.approx(normal["on_source_in"] / nsuperstripe)
    assert multistripe["nint_in"] == pytest.approx(normal["nint_in"] / nsuperstripe)


def test_single_group_fml_path_stays_finite():
    result = _fml_result(_timing(nsuperstripe=1, ngroup=1, mingroups=1))
    uncertainty = _fractional_uncertainty(result)

    assert result["on_source_in"] > 0
    assert np.all(np.isfinite(uncertainty))


def test_phase_spec_fml_uses_out_measurement_time_without_in_scalar():
    result = _phase_result(_timing(nsuperstripe=1), "run_phase_spec_fml")

    assert np.diff(result["time"]) == pytest.approx([4.0, 4.0, 4.0, 4.0, 4.0])
    assert result["on_source_in"] == pytest.approx(4.0)
    assert np.all(np.isfinite(_fractional_uncertainty(result)))


def test_phase_spec_slope_uses_full_multistripe_measurement_time_and_pandeia_noise():
    timing = _timing(nsuperstripe=3)
    update_timing_measurement_time(timing, measurement_time_per_int=12.0)
    result = _phase_result(timing, "run_phase_spec_slope")

    assert np.diff(result["time"]) == pytest.approx([12.0])
    assert result["on_source_in"] == pytest.approx(12.0)
    assert result["photon_out_1d"] == pytest.approx([3600.0, 3600.0])
    assert result["var_out_1d"] == pytest.approx([3600.0, 3600.0])
    assert np.all(np.isfinite(_fractional_uncertainty(result)))


def _skip_if_pandeia_refdata_invalid():
    from pandeia.engine import pandeia_version

    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        pandeia_version()
    version_report = stream.getvalue()
    if "INVALID INSTALLATION" in version_report or "ENVIRONMENT VARIABLE UNSET" in version_report:
        pytest.skip(version_report)


def _niriss_config(subarray):
    with open("pandexo/engine/reference/niriss_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf = deepcopy(conf)
    conf["detector"]["subarray"] = subarray
    conf["detector"]["ngroup"] = 2
    return conf


@pytest.mark.parametrize(
    ("subarray", "expected_nsuperstripe"),
    [
        ("sub17stripe_soss", 120),
        ("sub60stripe_soss", 34),
        ("sub204stripe_soss", 10),
        ("sub680stripe_soss", 3),
    ],
)
def test_pandeia_soss_multistripe_exposes_superstripes(subarray, expected_nsuperstripe):
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    exp_pars = InstrumentFactory(config=_niriss_config(subarray)).the_detector.exposure_spec
    timing = _timing(nsuperstripe=int(getattr(exp_pars, "nsuperstripe", 1) or 1))

    assert exp_pars.nsuperstripe == expected_nsuperstripe
    assert select_calculation("um", exp_pars.nsuperstripe) == "slope method"
    assert timing["Num Superstripes"] > 1
    assert "Num Integrations In Transit" in timing
    assert "Effective Integrations In Transit" in timing
    assert "On Source Time In Transit" in timing
    assert "Effective On Source Time In Transit" in timing


@pytest.mark.parametrize("subarray", ["substrip96", "substrip256"])
def test_pandeia_standard_soss_substrips_have_one_superstripe(subarray):
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    exp_pars = InstrumentFactory(config=_niriss_config(subarray)).the_detector.exposure_spec
    timing = _timing(nsuperstripe=int(getattr(exp_pars, "nsuperstripe", 1) or 1))

    assert exp_pars.nsuperstripe == 1
    assert select_calculation("um", exp_pars.nsuperstripe) == "fml"
    assert timing["Num Superstripes"] == 1
    assert timing["Effective Integrations In Transit"] == timing["Num Integrations In Transit"]


def test_pandeia_dhs_readout_uses_slope_without_superstripes():
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    with open("pandexo/engine/reference/nircam_dhs_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf["detector"]["ngroup"] = 2
    exp_pars = InstrumentFactory(config=conf).the_detector.exposure_spec

    assert exp_pars.nsuperstripe == 1
    assert select_calculation("um", exp_pars.nsuperstripe, is_dhs=True) == "slope method"
