"""CI-friendly tests for PandExo saturation warnings and errors."""

import numpy as np
import pytest

from pandexo.engine.jwst import (
    mask_fully_saturated_final_spectrum,
    nirspec_valid_channel_mask,
    validate_saturation_state,
)


OBSERVING_MODE_CONFIGS = [
    (
        "miri_lrs",
        {
            "instrument": {
                "instrument": "miri",
                "mode": "lrsslitless",
                "disperser": "prism",
                "filter": "clear",
            },
            "detector": {
                "subarray": "slitlessprism",
                "readout_pattern": "fastr1",
            },
        },
    ),
    (
        "nirspec_prism",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "prism",
                "filter": "clear",
            },
            "detector": {
                "subarray": "sub512",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g140m",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g140m",
                "filter": "f100lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g140h",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g140h",
                "filter": "f100lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g235m",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g235m",
                "filter": "f170lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g235h",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g235h",
                "filter": "f170lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g395m",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g395m",
                "filter": "f290lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "nirspec_g395h",
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "bots",
                "disperser": "g395h",
                "filter": "f290lp",
            },
            "detector": {
                "subarray": "sub2048",
                "readout_pattern": "nrsrapid",
            },
        },
    ),
    (
        "niriss_soss_legacy",
        {
            "instrument": {
                "instrument": "niriss",
                "mode": "soss",
                "disperser": "gr700xd",
                "filter": "clear",
            },
            "detector": {
                "subarray": "substrip96",
                "readout_pattern": "nisrapid",
            },
        },
    ),
    (
        "niriss_soss_multistripe",
        {
            "instrument": {
                "instrument": "niriss",
                "mode": "soss",
                "disperser": "gr700xd",
                "filter": "clear",
            },
            "detector": {
                "subarray": "sub680stripe_soss",
                "readout_pattern": "nisrapid",
            },
        },
    ),
    (
        "nircam_f322w2",
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "lw_tsgrism",
                "aperture": "lw",
                "disperser": "grismr",
                "filter": "f322w2",
            },
            "detector": {
                "subarray": "subgrism64",
                "readout_pattern": "rapid",
            },
        },
    ),
    (
        "nircam_f444w",
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "lw_tsgrism",
                "aperture": "lw",
                "disperser": "grismr",
                "filter": "f444w",
            },
            "detector": {
                "subarray": "subgrism64",
                "readout_pattern": "rapid",
            },
        },
    ),
    (
        "nircam_dhs",
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "dhs",
                "aperture": "dhs0bright",
                "filter": "f150w2",
            },
            "detector": {
                "subarray": "sub260s4_8-spectra",
                "readout_pattern": "rapid",
            },
        },
    ),
]


def _pandeia_output_with_warning(warning_key, warning_text, **scalar):
    return {
        "warnings": {
            warning_key: warning_text,
        },
        "scalar": scalar,
    }


@pytest.mark.parametrize(
    "mode_name, conf",
    OBSERVING_MODE_CONFIGS,
    ids=[mode_name for mode_name, _ in OBSERVING_MODE_CONFIGS],
)
def test_partial_saturation_warns_but_keeps_valid_channels(mode_name, conf):
    """Treat partially saturated flux levels as warnings, not hard failures.

    This synthetic Pandeia output mimics a bright-but-not-hopeless target:
    at least one spectral channel still has positive ``extracted_noise``, while
    Pandeia reports ``partial_saturated``.  The test is parameterized across the
    JWST observing-mode families used by the regression notebook without
    requiring CI to have Pandeia reference data installed.
    """
    pandeia_output = _pandeia_output_with_warning(
        "partial_saturated",
        "partial saturation in detector image",
        fraction_saturation=0.7,
        sat_ngroups=2,
    )
    extracted_noise = np.array([12.0, 0.0, np.nan, 3.0])

    with pytest.warns(UserWarning, match="partial saturation") as warnings:
        validate_saturation_state(conf, pandeia_output, extracted_noise)

    message = str(warnings[0].message)
    assert mode_name.split("_")[0] in message
    assert "fraction_saturation=0.7" in message
    assert "sat_ngroups=2" in message


@pytest.mark.parametrize(
    "mode_name, conf",
    OBSERVING_MODE_CONFIGS,
    ids=[mode_name for mode_name, _ in OBSERVING_MODE_CONFIGS],
)
def test_sparse_full_saturation_warns_but_keeps_valid_channels(mode_name, conf):
    """Warn when Pandeia reports full saturation in only part of the image.

    A few fully saturated pixels should be visible to users, but they should not
    abort the calculation while Pandeia still returns usable spectral channels.
    The hard error is reserved for the case where every channel has invalid
    ``extracted_noise``.
    """
    pandeia_output = _pandeia_output_with_warning(
        "full_saturated",
        "some detector pixels are fully saturated",
        fraction_saturation=1.2,
        sat_ngroups=1,
    )
    extracted_noise = np.array([12.0, 0.0, np.nan, 3.0])

    with pytest.warns(UserWarning, match="full saturation") as warnings:
        validate_saturation_state(conf, pandeia_output, extracted_noise)

    message = str(warnings[0].message)
    assert mode_name.split("_")[0] in message
    assert "fraction_saturation=1.2" in message
    assert "sat_ngroups=1" in message


@pytest.mark.parametrize(
    "mode_name, conf",
    OBSERVING_MODE_CONFIGS,
    ids=[mode_name for mode_name, _ in OBSERVING_MODE_CONFIGS],
)
def test_full_spectrum_saturation_raises_clear_error(mode_name, conf):
    """Reject flux levels where every spectral channel is unusable.

    This represents the fully saturated regression case: Pandeia may still
    return wavelength samples, but every ``extracted_noise`` value is zero,
    negative, or non-finite.  PandExo should raise before downstream code
    collapses into an empty spectrum with a cryptic NumPy reduction error.
    """
    pandeia_output = _pandeia_output_with_warning(
        "full_saturated",
        "all spectral channels saturated",
        fraction_saturation=15.4,
        sat_ngroups=0,
    )
    extracted_noise = np.array([0.0, np.nan, -1.0])

    with pytest.raises(ValueError, match="All spectral channels"):
        validate_saturation_state(conf, pandeia_output, extracted_noise)


def test_nirspec_gap_mask_keeps_fully_saturated_channels_for_nan_masking():
    """Keep saturated NIRSpec channels even when their noise is invalid.

    The NIRSpec detector-gap mask uses non-positive or non-finite Pandeia
    ``extracted_noise`` values to remove unobserved wavelengths.  Fully
    saturated wavelengths can have the same invalid-noise signature, but they
    should survive this gap-removal step so the final spectrum can mark them as
    NaN rather than silently dropping them.
    """
    conf = dict(OBSERVING_MODE_CONFIGS[1][1])
    extracted_noise = np.array([10.0, 0.0, np.nan, 5.0])
    full_saturation = np.array([0.0, 2.0, 0.0, 0.0])

    valid_channel = nirspec_valid_channel_mask(
        conf, extracted_noise, full_saturation
    )

    np.testing.assert_array_equal(
        valid_channel,
        np.array([True, True, False, True]),
    )


def test_full_saturation_mask_sets_final_spectrum_values_to_nan():
    """Mark final spectrum values as NaN for hard-saturated wavelengths.

    Pandeia may return finite flux and noise estimates for channels with
    ``n_full_saturated > 0``.  PandExo should keep the wavelength samples but
    make the invalid precision/products obvious in the public ``FinalSpectrum``
    arrays.
    """
    finalspec = {
        "wave": np.array([1.0, 2.0, 3.0]),
        "spectrum": np.array([0.01, 0.02, 0.03]),
        "spectrum_w_rand": np.array([0.011, 0.022, 0.033]),
        "error_w_floor": np.array([100.0, 200.0, 300.0]),
    }
    full_saturation_mask = np.array([False, True, False])

    masked = mask_fully_saturated_final_spectrum(
        finalspec, full_saturation_mask
    )

    assert np.isfinite(masked["wave"]).all()
    assert np.isnan(masked["spectrum"][1])
    assert np.isnan(masked["spectrum_w_rand"][1])
    assert np.isnan(masked["error_w_floor"][1])
    assert masked["spectrum"][0] == 0.01
    np.testing.assert_array_equal(
        masked["full_saturation_mask"], full_saturation_mask
    )
