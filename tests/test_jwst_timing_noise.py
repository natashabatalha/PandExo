"""JWST timing, calculation-routing, and noise-scaling regression tests."""

import contextlib
import io
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from pandexo.engine.compute_noise import ExtractSpec
from pandexo.engine.jwst import (
    DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB,
    DHS_READOUT_PATTERNS,
    NIRCAM_READOUT_PATTERNS,
    _nircam_dhs_optimization_configs,
    _table_html,
    add_warnings,
    apt_exposure_parameters,
    build_timing_display_div,
    compute_timing,
    estimate_dhs_data_excess,
    estimate_nircam_data_excess,
    nircam_dhs_no_ta_overhead,
    nircam_no_ta_overhead,
    nirspec_prism_apt_parameters,
    nirspec_prism_exposure_warning,
    select_calculation,
    update_nirspec_prism_apt_timing,
    update_apt_timing,
    update_timing_measurement_time,
    validate_miri_lrs_subarray,
)


@pytest.mark.parametrize(
    ("subarray", "readout", "ngroup", "expected_rate"),
    [
        ("sub41s1_2-spectra", "rapid", 2, 6.02),
        ("sub260s4_8-spectra", "rapid", 2, 6.21),
        ("sub260s4_8-spectra", "bright1", 3, 3.88),
        ("sub260s4_8-spectra", "dhs4", 3, 1.07),
        ("sub260s4_8-spectra", "dhs6", 3, -0.13),
    ],
)
def test_dhs_data_excess_estimate_matches_stsci_table(
    subarray, readout, ngroup, expected_rate
):
    rate, total = estimate_dhs_data_excess(
        subarray, readout, ngroup, exposure_hours=2.0
    )

    assert rate == pytest.approx(expected_rate, abs=0.01)
    assert total == pytest.approx(max(0.0, 2.0 * expected_rate), abs=0.02)


def test_dhs_readout_order_reaches_first_recommended_six_hour_setup():
    selected = None
    for readout in DHS_READOUT_PATTERNS:
        _, total = estimate_dhs_data_excess(
            "sub260s4_8-spectra", readout, ngroup=2, exposure_hours=6.0
        )
        if total <= DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB:
            selected = readout
            break

    assert selected == "dhs3"


def test_dhs_no_ta_overhead_includes_standard_initial_slew():
    overhead = nircam_dhs_no_ta_overhead(tframe=1.36765)

    assert overhead == pytest.approx(2794.683825)
    assert overhead == nircam_no_ta_overhead(tframe=1.36765)


def test_dhs_data_excess_uses_no_ta_allocation_overhead():
    _, without_overhead = estimate_dhs_data_excess(
        "sub260s4_8-spectra", "dhs3", 100, exposure_hours=5.7931374583
    )
    _, with_overhead = estimate_dhs_data_excess(
        "sub260s4_8-spectra",
        "dhs3",
        100,
        exposure_hours=5.7931374583,
        allocation_overhead_seconds=nircam_dhs_no_ta_overhead(1.36765),
    )

    assert without_overhead == pytest.approx(9.006139368, abs=1e-6)
    assert with_overhead == pytest.approx(6.574764707, abs=1e-6)


def test_standard_nircam_optimizer_reaches_first_recommended_readout():
    selected = None
    for readout in NIRCAM_READOUT_PATTERNS:
        _, total = estimate_nircam_data_excess(
            "subgrism64", readout, ngroup=100, exposure_hours=6.0
        )
        if total <= DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB:
            selected = readout
            break

    assert selected == "bright1"


def test_standard_nircam_one_output_subarray_can_retain_rapid():
    rate, total = estimate_nircam_data_excess(
        "subgrism64 (noutputs=1)", "rapid", ngroup=100,
        exposure_hours=6.0,
    )

    assert rate < 0
    assert total == 0


def test_standard_nircam_multiframe_readout_includes_frame_zero():
    bright1_rate, _ = estimate_nircam_data_excess(
        "subgrism64", "bright1", ngroup=2, exposure_hours=1.0
    )
    bright2_rate, _ = estimate_nircam_data_excess(
        "subgrism64", "bright2", ngroup=2, exposure_hours=1.0
    )

    assert bright2_rate > bright1_rate


def test_standard_nircam_data_excess_matches_apt_bright1_setup():
    _, total = estimate_nircam_data_excess(
        "subgrism64",
        "bright1",
        ngroup=14,
        exposure_hours=21622.625 / 3600.0,
        allocation_overhead_seconds=nircam_no_ta_overhead(0.34061),
    )

    assert total == pytest.approx(4.5875, abs=0.002)


def test_standard_nircam_data_excess_warning_is_exposed():
    timing = _timing(nsuperstripe=1)
    timing["Estimated NIRCam Data Excess (GB)"] = 10.0
    warnings = add_warnings(
        {"warnings": {}},
        timing,
        sat_level=0.8,
        flags={
            "flag_default": "All good",
            "flag_high": "All good",
            "flag_min_nint": "All good",
            "flag_nircam_readout": "Selected BRIGHT1 after RAPID.",
        },
        instrument="nircam",
    )

    assert warnings["NIRCam Readout Optimization"].startswith(
        "Selected BRIGHT1"
    )
    assert "above the 5 GB lower threshold" in warnings[
        "NIRCam Data Excess?"
    ]
    assert "10.0 GB" in warnings["NIRCam Data Excess?"]
    assert "10.00 GB" not in warnings["NIRCam Data Excess?"]


def test_dhs_optimization_configs_are_independent_of_displayed_channel():
    base_conf = {
        'instrument': {
            'instrument': 'nircam',
            'mode': 'sw_tsgrism',
            'filter': 'f150w2',
            'pandexofilterpair': 'f322w2',
            'aperture': 'dhs0spec8',
            'disperser': 'dhs0',
        },
        'detector': {
            'readout_pattern': 'rapid',
            'subarray': 'sub260s4_8-spectra',
            'ngroup': 'optimize',
        },
    }
    long_wave_conf = deepcopy(base_conf)
    long_wave_conf['instrument'].update(
        mode='lw_tsgrism',
        filter='f322w2',
        pandexofilterpair='f150w2',
        aperture='lw',
        disperser='grismr',
    )

    short_display_configs = _nircam_dhs_optimization_configs(base_conf)
    long_display_configs = _nircam_dhs_optimization_configs(long_wave_conf)

    for short_display, long_display in zip(
        short_display_configs, long_display_configs
    ):
        for key in (
            'filter', 'pandexofilterpair', 'mode', 'aperture', 'disperser'
        ):
            assert short_display['instrument'][key] == (
                long_display['instrument'][key]
            )
        assert short_display['detector'] == long_display['detector']
    assert short_display_configs[0]['instrument']['filter'] == 'f150w2'
    assert short_display_configs[1]['instrument']['filter'] == 'f322w2'
    assert base_conf['instrument']['aperture'] == 'dhs0spec8'
    assert long_wave_conf['instrument']['aperture'] == 'lw'


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


def _timing_with_pandeia_cycle(
        nsuperstripe, exposure_time_per_int, ngroup=3,
        transit_duration=48.0):
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
        transit_duration=transit_duration,
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


def test_skipped_frames_contribute_to_elapsed_integration_time():
    timing, _ = compute_timing(
        {
            "ngroup": 3,
            "tframe": 2.0,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 2,
            "nsuperstripe": 1,
        },
        transit_duration=24.0,
        expfact_out=1.0,
        noccultations=1,
    )

    assert timing["Time/Integration incl reset (sec)"] == pytest.approx(16.0)
    assert timing["Measurement Time per Integration (sec)"] == pytest.approx(12.0)


def test_multistripe_cycles_are_not_divided_by_nsuperstripe():
    timing = _timing(nsuperstripe=4)

    assert (
        timing["Effective Integrations In Transit"]
        == timing["Num Integrations In Transit"]
    )
    assert (
        timing["Effective Integrations Out of Transit"]
        == timing["Num Integrations Out of Transit"]
    )
    assert timing["On Source Time In Transit"] == pytest.approx(
        timing["Num Integrations In Transit"]
        * timing["Measurement Time per Integration (sec)"]
        / timing["Num Superstripes"]
    )
    assert timing["Effective On Source Time In Transit"] == pytest.approx(
        timing["On Source Time In Transit"]
    )
    assert timing["Effective On Source Time Out of Transit"] == pytest.approx(
        timing["On Source Time Out of Transit"]
    )


def test_slope_method_reports_per_stripe_electrons_per_real_integration():
    result = _slope_result(_timing(nsuperstripe=4))

    assert result["photon_out_1d_per_int"] == pytest.approx([250.0, 250.0])


def test_fml_method_reports_electrons_per_real_integration():
    result = _fml_result(_timing(nsuperstripe=1))

    assert result["photon_out_1d_per_int"] == pytest.approx([1000.0, 1000.0])


def test_pandeia_measurement_time_update_controls_on_source_metadata():
    timing = _timing(nsuperstripe=4)
    update_timing_measurement_time(timing, measurement_time_per_int=32.0)

    assert timing["Measurement Time per Integration (sec)"] == 32.0
    assert timing["On Source Time In Transit"] == pytest.approx(
        timing["Num Integrations In Transit"] * 32.0 / 4.0
    )
    assert timing["Effective On Source Time In Transit"] == pytest.approx(
        timing["On Source Time In Transit"]
    )


def test_multistripe_timing_uses_pandeia_full_cycle_clock_for_real_integrations():
    timing = _timing_with_pandeia_cycle(
        nsuperstripe=4,
        exposure_time_per_int=40.0,
        ngroup=3,
    )

    assert timing["Time/Integration incl reset (sec)"] == pytest.approx(40.0)
    assert timing["Num Integrations In Transit"] == 2
    assert timing["Effective Integrations In Transit"] == 2
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
    assert timing["Effective Integrations In Transit"] == 147
    assert timing["On Source Time In Transit"] == pytest.approx(147 * 0.43148)


def test_soss_multistripe_efficiency_uses_per_stripe_science_time():
    timing, _ = compute_timing(
        {
            "ngroup": 1111,
            "tframe": 0.06164,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 120,
            "exposure_time_per_int": 8227.6992,
            "exposure_time_ngroup": 1111,
        },
        transit_duration=2.8032 * 3600.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )

    assert timing["Observing Efficiency (%)"] == pytest.approx(0.8315860648)


def test_optimized_multistripe_timing_reduces_ngroups_for_minimum_integrations():
    timing, flags = compute_timing(
        {
            "maxexptime_per_int": 1109 * 0.06164,
            "tframe": 0.06164,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 120,
            "tfffr": 0.02048,
            "nreset1": 1,
            "ndrop1": 0,
            "ndrop3": 0,
        },
        transit_duration=2.8032 * 3600.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )

    assert timing["APT: Num Groups per Integration"] == 369
    assert timing["Num Integrations In Transit"] >= 3
    assert "Reduced NGROUPS from 1109 to 369" in flags["flag_min_nint"]


def test_user_multistripe_timing_warns_without_changing_ngroups_below_minimum():
    timing, flags = compute_timing(
        {
            "ngroup": 1111,
            "tframe": 0.06164,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 120,
            "exposure_time_per_int": 8227.6992,
            "exposure_time_ngroup": 1111,
        },
        transit_duration=2.8032 * 3600.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=65536,
    )

    warnings = add_warnings(
        {"warnings": {}},
        timing,
        sat_level=0.8,
        flags=flags,
        instrument="niriss",
    )

    assert timing["APT: Num Groups per Integration"] == 1111
    assert timing["Num Integrations In Transit"] == 2
    assert "User-specified NGROUPS produces 2" in flags["flag_min_nint"]
    assert warnings["Minimum Integrations?"] == flags["flag_min_nint"]


def test_dhs_lower_threshold_warning_is_exposed():
    timing = _timing(nsuperstripe=1)
    timing["Estimated DHS Data Excess (GB)"] = 10.0
    warnings = add_warnings(
        {"warnings": {}},
        timing,
        sat_level=0.8,
        flags={
            "flag_default": "All good",
            "flag_high": "All good",
            "flag_min_nint": "All good",
            "flag_dhs_readout": "Selected DHS3 after RAPID, BRIGHT1 exceeded the limit.",
        },
        instrument="nircam",
    )

    assert warnings["DHS Readout Optimization"].startswith("Selected DHS3")
    assert "above the 5 GB lower threshold" in warnings["DHS Data Excess?"]
    assert "acceptable for DHS" in warnings["DHS Data Excess?"]
    assert "10.0 GB" in warnings["DHS Data Excess?"]
    assert "10.00 GB" not in warnings["DHS Data Excess?"]


def test_multistripe_timing_display_uses_apt_and_calculation_tables():
    timing = _timing_with_pandeia_cycle(
        nsuperstripe=4,
        exposure_time_per_int=40.0,
        ngroup=3,
    )
    apt_div, calculation_div = build_timing_display_div(_pandeia_out(), timing)
    html = (apt_div + calculation_div).decode()

    assert "APT Inputs" not in html
    assert "Calculation Details" not in html
    assert html.count("<table") == 2
    assert "Groups per Integration" in html
    assert "Exposures/Dith" in html
    assert "Integrations/Exp" in html
    assert "Number of Stripes" in html
    assert "Elapsed Time per Full Multistripe Cycle incl. Reset (sec)" in html
    assert "Elapsed Time per APT Stripe Integration" not in html
    assert "Science Time per Full Multistripe Cycle excl. Reset (sec)" in html
    assert "Science Time per Stripe excl. Reset (sec)" in html
    assert "Effective Per-Wavelength" not in html
    assert "Time/Integration incl reset" not in html
    assert "Measurement Time per Integration" not in html


def test_timing_display_includes_estimated_dhs_data_excess():
    timing = _timing(nsuperstripe=1)
    timing["Estimated DHS Data Excess (GB)"] = 4.26
    timing["Assumed DHS Allocation Overhead (sec)"] = 2794.68
    _, calculation_div = build_timing_display_div(_pandeia_out(), timing)

    html = calculation_div.decode()

    assert "Estimated DHS Data Excess (GB)" in html
    assert "4.3 (Verify using APT)" in html
    assert "4.26" not in html
    assert "Assumed No-TA Scheduling + Slew Overhead (sec)" in html


def test_total_observing_time_is_quantized_once_for_multistripe_timing():
    cycle_time = 1.09312
    requested_cycles = 30810
    total_observing_time = requested_cycles * cycle_time
    transit_duration = 3.816 * 3600.0

    timing, _ = compute_timing(
        {
            "ngroup": 3,
            "tframe": 0.02904,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 8,
            "exposure_time_per_int": cycle_time,
            "exposure_time_ngroup": 3,
        },
        transit_duration=transit_duration,
        expfact_out=transit_duration / (
            total_observing_time - transit_duration
        ),
        noccultations=1,
        total_observing_time=total_observing_time,
    )

    assert timing["Num Integrations In Transit"] == 12568
    assert timing["Num Integrations Out of Transit"] == 18242
    assert timing["APT: Num Integrations per Occultation"] == requested_cycles
    assert timing["Transit+Baseline, no overhead (hrs)"] == pytest.approx(
        total_observing_time / 3600.0
    )


def test_nirspec_prism_apt_parameters_pad_published_example_to_full_cycles():
    timing = {
        "Num Superstripes": 8,
        "Num Integrations In Transit": 15405,
        "Num Integrations Out of Transit": 15405,
        "Time/Integration incl reset (sec)": 1.09312,
    }

    parameters = nirspec_prism_apt_parameters(timing)

    assert parameters["cycles"] == 30810
    assert parameters["required_integrations"] == 246480
    assert parameters["exposures_per_dither"] == 4
    assert parameters["integrations_per_exposure"] == 61624
    assert parameters["scheduled_integrations"] == 246496
    assert parameters["integrations_per_exposure"] % 8 == 0
    assert parameters["stripe_elapsed_time"] == pytest.approx(0.13664)
    assert parameters["exposure_duration"] == pytest.approx(33681.21344)


def test_nirspec_prism_apt_values_are_available_to_offline_outputs():
    timing = {
        "Num Superstripes": 8,
        "Num Integrations In Transit": 15405,
        "Num Integrations Out of Transit": 15405,
        "Time/Integration incl reset (sec)": 1.09312,
        "APT: Num Integrations per Occultation": 30810,
    }

    update_nirspec_prism_apt_timing(timing)

    assert timing["Num Multistripe Cycles per Occultation"] == 30810
    assert timing["APT: Exposures/Dith"] == 4
    assert timing["APT: Num Integrations per Exposure"] == 61624
    assert timing["APT: Num Integrations per Occultation"] == 246496


def test_nirspec_prism_apt_parameters_round_up_without_exceeding_limit():
    timing = {
        "Num Superstripes": 8,
        "Num Integrations In Transit": 8192,
        "Num Integrations Out of Transit": 8192,
        "Time/Integration incl reset (sec)": 1.09312,
    }

    parameters = nirspec_prism_apt_parameters(timing)

    assert parameters["required_integrations"] == 131072
    assert parameters["exposures_per_dither"] == 3
    assert parameters["integrations_per_exposure"] == 43696
    assert parameters["integrations_per_exposure"] <= 65535
    assert parameters["integrations_per_exposure"] % 8 == 0
    assert parameters["scheduled_integrations"] == 131088


def test_apt_integration_limit_is_applied_to_standard_modes():
    timing = {
        "Num Integrations In Transit": 40000,
        "Num Integrations Out of Transit": 40000,
        "Time/Integration incl reset (sec)": 2.0,
    }

    update_apt_timing(timing)

    assert timing["APT: Exposures/Dith"] == 2
    assert timing["APT: Num Integrations per Exposure"] == 40000
    assert timing["APT: Num Integrations per Occultation"] == 80000


def test_frame_limit_can_require_additional_exposures():
    timing = {
        "Num Integrations In Transit": 25000,
        "Num Integrations Out of Transit": 25000,
        "Time/Integration incl reset (sec)": 10.0,
    }

    parameters = apt_exposure_parameters(
        timing,
        max_frames_per_exposure=196608,
        frames_per_integration=10,
    )

    assert parameters["max_integrations_per_exposure"] == 19660
    assert parameters["exposures_per_dither"] == 3
    assert parameters["integrations_per_exposure"] == 16667
    assert parameters["integrations_per_exposure"] * 10 <= 196608


def test_exposure_limit_must_fit_one_complete_multistripe_cycle():
    timing = {
        "Num Integrations In Transit": 1,
        "Num Integrations Out of Transit": 1,
        "Time/Integration incl reset (sec)": 8.0,
    }

    with pytest.raises(
            ValueError, match="one complete multistripe cycle"):
        apt_exposure_parameters(
            timing,
            integration_multiplier=8,
            max_integrations_per_exposure=7,
        )


def test_niriss_soss_groups_are_limited_to_30():
    from pandexo.engine.jwst import max_ngroup

    assert max_ngroup["niriss"] == 30


def test_user_specified_groups_are_capped_at_instrument_limit():
    timing, flags = compute_timing(
        {
            "ngroup": 31,
            "tframe": 2.0,
            "nframe": 1,
            "mingroups": 2,
            "nskip": 0,
            "nsuperstripe": 1,
        },
        transit_duration=240.0,
        expfact_out=1.0,
        noccultations=1,
        max_ngroup_instrument=30,
    )

    assert timing["APT: Num Groups per Integration"] == 30
    assert "User-specified NGROUPS (31) exceeds the maximum (30)" in flags["flag_high"]


def test_nirspec_prism_timing_display_uses_apt_stripe_integrations():
    timing = _timing_with_pandeia_cycle(
        nsuperstripe=8,
        exposure_time_per_int=1.09312,
        ngroup=3,
    )
    timing.update({
        "Num Integrations In Transit": 15405,
        "Num Integrations Out of Transit": 15405,
        "Effective Integrations In Transit": 15405,
        "Effective Integrations Out of Transit": 15405,
        "APT: Num Integrations per Occultation": 30810,
    })
    apt_div, calculation_div = build_timing_display_div(
        _pandeia_out(
            instrument={
                "instrument": "nirspec",
                "mode": "bots",
                "filter": "clear",
                "aperture": "s1600a1",
                "disperser": "prism",
            },
            detector={
                "subarray": "s64m8_prm",
                "readout_pattern": "nrsrapid",
            },
        ),
        timing,
    )
    apt_html = apt_div.decode()
    calculation_html = calculation_div.decode()

    assert "Exposures/Dith" in apt_html
    assert "<td>4</td>" in apt_html
    assert "Integrations/Exp" in apt_html
    assert "<td>61624</td>" in apt_html
    assert "Integrations per Occultation" not in apt_html
    assert "Elapsed Time per Full Multistripe Cycle incl. Reset (sec)" in calculation_html
    assert "<td>1.093120</td>" in calculation_html
    assert "Elapsed Time per APT Stripe Integration incl. Reset (sec)" in calculation_html
    assert "<td>0.136640</td>" in calculation_html


def test_nirspec_multistripe_hga_warning_is_included_when_relevant():
    long_timing = {
        "Num Superstripes": 8,
        "Num Integrations In Transit": 15405,
        "Num Integrations Out of Transit": 15405,
        "Time/Integration incl reset (sec)": 1.09312,
    }
    warning = nirspec_prism_exposure_warning(long_timing)
    warnings = add_warnings(
        {"warnings": {}},
        _timing(nsuperstripe=8),
        sat_level=0.8,
        flags={
            "flag_default": "All good",
            "flag_high": "All good",
            "flag_min_nint": "All good",
            "flag_hga_repoint": warning,
        },
        instrument="nirspec",
    )

    assert warnings["Exposure Duration?"].startswith("APT will warn")
    assert "NIRSpec BOTS permits longer exposures" in warnings["Exposure Duration?"]


def test_view_template_owns_timing_table_headings():
    template = Path("pandexo/engine/templates/view.html").read_text()

    assert "<h3>APT Inputs</h3>" in template
    assert "<h3>Calculation Details</h3>" in template
    assert "{% raw div['apt_div']  %}" in template
    assert "{% raw div['calculation_div']  %}" in template
    assert "{% raw div['timing_div']  %}" in template


def test_summary_tables_use_shared_fixed_column_layout():
    template = Path("pandexo/engine/templates/view.html").read_text()
    html = _table_html([('Parameter', 'Value')])

    assert 'pandexo-summary-table' in html
    assert '.pandexo-summary-table {' in template
    assert 'table-layout: fixed;' in template
    assert 'width: 50%;' in template


def test_non_multistripe_timing_display_omits_stripe_rows():
    timing = _timing(nsuperstripe=1)
    apt_div, calculation_div = build_timing_display_div(_pandeia_out(), timing)
    html = (apt_div + calculation_div).decode()

    assert "Elapsed Time per Integration incl. Reset (sec)" in html
    assert "Science Time per Integration excl. Reset (sec)" in html
    assert "Number of Stripes" not in html
    assert "Full Multistripe Cycle" not in html
    assert "Science Time per Stripe" not in html


def test_timing_display_formats_transit_and_integration_counts_as_integers():
    timing = _timing(nsuperstripe=1)
    timing['Number of Transits'] = 1.0
    timing['Num Integrations In Transit'] = 12.0
    timing['Num Integrations Out of Transit'] = 15.0

    apt_div, calculation_div = build_timing_display_div(_pandeia_out(), timing)
    html = (apt_div + calculation_div).decode()

    assert "<td>1</td>" in html
    assert "<td>12</td>" in html
    assert "<td>15</td>" in html
    assert "<td>1.0</td>" not in html
    assert "<td>12.0</td>" not in html
    assert "<td>15.0</td>" not in html


def test_niriss_timing_display_omits_filter_row():
    timing = _timing(nsuperstripe=1)
    apt_div, _ = build_timing_display_div(_pandeia_out(), timing)
    html = apt_div.decode()

    assert "<th>Filter</th>" not in html
    assert "<td>clear</td>" not in html


def test_nircam_timing_display_shows_channel_and_pupil_rows():
    timing = _timing(nsuperstripe=1)
    apt_div, calculation_div = build_timing_display_div(
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
    )
    html = (apt_div + calculation_div).decode()

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
    apt_div, calculation_div = build_timing_display_div(
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
    )
    html = (apt_div + calculation_div).decode()

    assert "SW Channel Mode" in html
    assert "IMAGING" in html
    assert "SUBGRISM64" in html
    assert "No. of Output Channels" in html
    assert "<td>4</td>" in html
    assert "Short Pupil+Filter" in html
    assert "CHOOSE THIS USING ETC" in html
    assert "Long Pupil+Filter" in html
    assert "GRISMR+F322W2" in html


def test_nircam_timing_display_formats_estimated_data_excess():
    timing = _timing(nsuperstripe=1)
    timing['Estimated NIRCam Data Excess (GB)'] = 4.5875
    timing['Assumed NIRCam Allocation Overhead (sec)'] = 2794.0
    _, calculation_div = build_timing_display_div(
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
                "readout_pattern": "bright1",
            },
        ),
        timing,
    )
    html = calculation_div.decode()

    assert "Estimated NIRCam Data Excess (GB)" in html
    assert "4.6 (Verify using APT)" in html
    assert "4.6 GB" not in html
    assert "4.5875" not in html


@pytest.mark.parametrize(
    ('disperser', 'filter_name', 'expected'),
    [
        ('prism', 'clear', 'PRISM/CLEAR'),
        ('g140h', 'f070lp', 'G140H/F070LP'),
    ],
)
def test_nirspec_timing_display_combines_grating_and_filter(
    disperser, filter_name, expected
):
    timing = _timing(nsuperstripe=1)
    apt_div, calculation_div = build_timing_display_div(
        _pandeia_out(
            instrument={
                'instrument': 'nirspec',
                'mode': 'bots',
                'filter': filter_name,
                'aperture': 's1600a1',
                'disperser': disperser,
            }
        ),
        timing,
    )
    html = (apt_div + calculation_div).decode()

    assert '<th>Grating/Filter</th>' in html
    assert expected in html
    assert '<th>Filter</th>' not in html


@pytest.mark.parametrize(
    ('subarray', 'expected'),
    [
        ('s256m2_prm', 'S256M2_PRISM'),
        ('s128m4_prm', 'S128M4_PRISM'),
        ('s64m8_prm', 'S64M8_PRISM'),
        ('s32m16_prm', 'S32M16_PRISM'),
    ],
)
def test_nirspec_multistripe_subarray_display_uses_prism_suffix(
    subarray, expected
):
    timing = _timing(nsuperstripe=1)
    apt_div, _ = build_timing_display_div(
        _pandeia_out(
            instrument={
                'instrument': 'nirspec',
                'mode': 'bots',
                'filter': 'clear',
                'aperture': 's1600a1',
                'disperser': 'prism',
            },
            detector={
                'subarray': subarray,
                'readout_pattern': 'nrsrapid',
            },
        ),
        timing,
    )
    html = apt_div.decode()

    assert expected in html
    assert '_PRM' not in html


def test_miri_lrs_timing_display_uses_dither_not_filter_and_uppercase_readout():
    timing = _timing(nsuperstripe=1)
    apt_div, calculation_div = build_timing_display_div(
        _pandeia_out(
            instrument={
                "instrument": "miri",
                "mode": "lrsslitless",
                "filter": None,
                "aperture": "imager",
                "disperser": "p750l",
            },
            detector={
                "subarray": "slitlessprism_ip",
                "readout_pattern": "fastr1",
            },
        ),
        timing,
    )
    html = (apt_div + calculation_div).decode()

    assert "<th>Filter</th>" not in html
    assert "<th>Dither</th>" in html
    assert "<td>None</td>" in html
    assert "<td>FASTR1</td>" in html
    assert "<td>fastr1</td>" not in html


def test_slope_noise_uses_real_multistripe_cycle_count():
    nsuperstripe = 9
    normal = _slope_result(_timing(nsuperstripe=1))
    multistripe = _slope_result(_timing(nsuperstripe=nsuperstripe))

    ratio = _fractional_uncertainty(multistripe) / _fractional_uncertainty(normal)

    assert ratio == pytest.approx(1.0)
    assert multistripe["on_source_in"] == pytest.approx(normal["on_source_in"] / nsuperstripe)
    assert multistripe["nint_in"] == normal["nint_in"]


def test_fixed_wall_time_retains_multistripe_cycle_penalty():
    normal = _slope_result(_timing_with_pandeia_cycle(
        1, 10.0, transit_duration=270.0
    ))
    multistripe = _slope_result(_timing_with_pandeia_cycle(
        9, 90.0, transit_duration=270.0
    ))

    ratio = _fractional_uncertainty(multistripe) / _fractional_uncertainty(normal)

    assert normal["nint_in"] == 27
    assert multistripe["nint_in"] == 3
    assert ratio == pytest.approx(3.0)


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


def _nirspec_prism_config(subarray):
    with open("pandexo/engine/reference/nirspec_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf = deepcopy(conf)
    conf["instrument"]["disperser"] = "prism"
    conf["instrument"]["filter"] = "clear"
    conf["detector"]["subarray"] = subarray
    conf["detector"]["ngroup"] = 2
    return conf


def _miri_lrs_slit_config(subarray):
    with open("pandexo/engine/reference/miri_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf = deepcopy(conf)
    conf["instrument"]["mode"] = "lrsslit"
    conf["instrument"]["aperture"] = "lrsslit"
    conf["detector"]["subarray"] = subarray
    conf["detector"]["ngroup"] = 2
    return conf


@pytest.mark.parametrize(
    ("mode", "subarray"),
    [
        ("lrsslitless", "slitlessprism"),
        ("lrsslitless", "slitlessprism_ip"),
        ("lrsslitless", "slitlessprism_ips"),
        ("lrsslit", "full"),
        ("lrsslit", "subslit"),
    ],
)
def test_validate_miri_lrs_subarray_accepts_supported_pairs(mode, subarray):
    with open("pandexo/engine/reference/miri_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf = deepcopy(conf)
    conf["instrument"]["mode"] = mode
    conf["detector"]["subarray"] = subarray

    validate_miri_lrs_subarray(conf)


@pytest.mark.parametrize(
    ("mode", "subarray"),
    [
        ("lrsslitless", "full"),
        ("lrsslitless", "subslit"),
        ("lrsslit", "slitlessprism"),
        ("lrsslit", "slitlessprism_ip"),
        ("lrsslit", "slitlessprism_ips"),
    ],
)
def test_validate_miri_lrs_subarray_rejects_unsupported_pairs(mode, subarray):
    with open("pandexo/engine/reference/miri_input.json") as handle:
        conf = json.load(handle)["configuration"]
    conf = deepcopy(conf)
    conf["instrument"]["mode"] = mode
    conf["detector"]["subarray"] = subarray

    with pytest.raises(ValueError, match="MIRI LRS"):
        validate_miri_lrs_subarray(conf)


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


@pytest.mark.parametrize(
    ("subarray", "expected_nsuperstripe"),
    [
        ("s256m2_prm", 2),
        ("s128m4_prm", 4),
        ("s64m8_prm", 8),
        ("s32m16_prm", 16),
    ],
)
def test_pandeia_nirspec_prism_multistripe_exposes_superstripes(
    subarray, expected_nsuperstripe
):
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    try:
        exp_pars = InstrumentFactory(
            config=_nirspec_prism_config(subarray)
        ).the_detector.exposure_spec
    except Exception as exc:
        pytest.skip(f"NIRSpec PRISM multistripe subarray is unavailable: {exc}")

    timing = _timing(nsuperstripe=int(getattr(exp_pars, "nsuperstripe", 1) or 1))

    assert exp_pars.nsuperstripe == expected_nsuperstripe
    assert select_calculation("um", exp_pars.nsuperstripe) == "slope method"
    assert timing["Num Superstripes"] == expected_nsuperstripe
    assert timing["Observing Efficiency (%)"] < 100.0


def test_pandeia_nirspec_multistripe_nint_counts_complete_cycles():
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    one_cycle_config = _nirspec_prism_config("s64m8_prm")
    one_cycle_config["detector"]["nint"] = 1
    one_cycle_config["detector"]["nexp"] = 1
    eight_cycle_config = deepcopy(one_cycle_config)
    eight_cycle_config["detector"]["nint"] = 8

    one_cycle = InstrumentFactory(
        config=one_cycle_config
    ).the_detector.exposure_spec
    eight_cycles = InstrumentFactory(
        config=eight_cycle_config
    ).the_detector.exposure_spec

    assert one_cycle.nsuperstripe == 8
    assert one_cycle.nramps == 1
    assert eight_cycles.nramps == 8
    assert eight_cycles.measurement_time == pytest.approx(
        8 * one_cycle.measurement_time
    )


def test_pandeia_miri_lrs_slit_supports_subslit():
    _skip_if_pandeia_refdata_invalid()
    from pandeia.engine.instrument_factory import InstrumentFactory

    try:
        exp_pars = InstrumentFactory(
            config=_miri_lrs_slit_config("subslit")
        ).the_detector.exposure_spec
    except Exception as exc:
        pytest.skip(f"MIRI LRS SUBSLIT subarray is unavailable: {exc}")

    assert exp_pars.tframe == pytest.approx(0.27904)
    assert exp_pars.nsuperstripe == 1


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
    assert conf["detector"]["readout_pattern"] == "optimize"
    conf["detector"]["readout_pattern"] = "rapid"
    conf["detector"]["ngroup"] = 2
    exp_pars = InstrumentFactory(config=conf).the_detector.exposure_spec

    assert exp_pars.nsuperstripe == 1
    assert select_calculation("um", exp_pars.nsuperstripe, is_dhs=True) == "slope method"


def test_dhs_throughput_resolves_optimized_readout_before_pandeia(monkeypatch):
    from pandexo.engine import justdoit

    captured = {}

    class FakeInstrument:
        def __init__(self, config):
            captured.update(config)

        def get_wave_range(self):
            return {"wmin": 1.0, "wmax": 2.0}

        def get_total_eff(self, wave):
            return np.ones_like(wave)

    monkeypatch.setattr(justdoit, "InstrumentFactory", FakeInstrument)

    justdoit.get_thruput("NIRCam DHS")

    assert captured["detector"]["readout_pattern"] == "rapid"
