import numpy as np
import pytest

from pandexo.engine.jwst import (
    dhs_f150w_wavelength_mask,
    nirspec_valid_channel_mask,
    nirspec_prism_multistripe_wavelength_mask,
)


def test_nirspec_mask_rejects_unobserved_channels():
    """Keep only NIRSpec channels with finite positive extracted noise.

    The test array represents Pandeia's per-channel ``extracted_noise`` values:
    ``12.0`` and ``3.0`` are usable observed channels, while ``0.0``, ``nan``,
    and ``-1.0`` are placeholders for channels with no valid extraction.
    """
    conf = {
        "instrument": {
            "instrument": "nirspec",
        }
    }

    mask = nirspec_valid_channel_mask(
        conf, np.array([12.0, 0.0, np.nan, -1.0, 3.0])
    )

    assert mask.tolist() == [True, False, False, False, True]


@pytest.mark.parametrize(
    ("instrument", "disperser"),
    [
        ("miri", "lrs"),
        ("niriss", "soss"),
        ("nircam", "grismr"),
    ],
)
def test_nirspec_mask_is_only_for_nirspec(instrument, disperser):
    """Do not apply the NIRSpec extracted-noise rule to other instruments.

    The ``0.0`` value would be rejected for NIRSpec, but this helper returns
    ``None`` for other instruments so their existing wavelength/flux filtering
    remains unchanged.
    """
    conf = {
        "instrument": {
            "instrument": instrument,
            "disperser": disperser,
        }
    }

    assert nirspec_valid_channel_mask(conf, np.array([1.0, 0.0])) is None


def test_dhs_f150w_mask_truncates_below_minimum_wavelength():
    conf = {
        "instrument": {
            "instrument": "nircam",
            "mode": "dhs",
            "aperture": "dhs0spec8",
            "filter": "f150w2",
        }
    }

    mask = dhs_f150w_wavelength_mask(
        conf, np.array([0.94, 0.96, 1.0, 1.5])
    )

    assert mask.tolist() == [False, True, True, True]


@pytest.mark.parametrize(
    "conf",
    [
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "lw_tsgrism",
                "aperture": "lw",
                "filter": "f322w2",
            }
        },
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "ssgrism",
                "aperture": "lw",
                "filter": "f150w2",
            }
        },
        {
            "instrument": {
                "instrument": "nirspec",
                "mode": "fixed_slit",
                "aperture": "s1600a1",
                "filter": "f100lp",
            }
        },
    ],
)
def test_dhs_f150w_mask_only_applies_to_dhs_f150w(conf):
    assert dhs_f150w_wavelength_mask(conf, np.array([0.94, 1.0])) is None


@pytest.mark.parametrize(
    ("subarray", "max_wavelength"),
    [
        ("s256m2_prm", 5.22),
        ("s128m4_prm", 5.14),
        ("s64m8_prm", 4.98),
        ("s32m16_prm", 4.66),
    ],
)
def test_nirspec_prism_multistripe_mask_truncates_red_edge(
    subarray, max_wavelength
):
    conf = {
        "instrument": {
            "instrument": "nirspec",
            "disperser": "prism",
            "filter": "clear",
        },
        "detector": {"subarray": subarray},
    }

    mask = nirspec_prism_multistripe_wavelength_mask(
        conf, np.array([max_wavelength - 0.01, max_wavelength, max_wavelength + 0.01])
    )

    assert mask.tolist() == [True, True, False]


@pytest.mark.parametrize(
    "conf",
    [
        {
            "instrument": {
                "instrument": "nirspec",
                "disperser": "prism",
                "filter": "clear",
            },
            "detector": {"subarray": "sub512"},
        },
        {
            "instrument": {
                "instrument": "nirspec",
                "disperser": "g395m",
                "filter": "f290lp",
            },
            "detector": {"subarray": "s32m16_prm"},
        },
        {
            "instrument": {
                "instrument": "nirspec",
                "disperser": "prism",
                "filter": "f100lp",
            },
            "detector": {"subarray": "s32m16_prm"},
        },
    ],
)
def test_nirspec_prism_multistripe_mask_only_applies_to_supported_modes(conf):
    assert nirspec_prism_multistripe_wavelength_mask(
        conf, np.array([4.0, 5.0])
    ) is None
