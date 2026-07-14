import numpy as np
import pytest

from pandexo.engine.jwst import miri_lrs_wavelength_mask


def test_miri_lrs_ips_mask_truncates_red_edge():
    conf = {
        "instrument": {
            "instrument": "miri",
            "mode": "lrsslitless",
        },
        "detector": {
            "subarray": "slitlessprism_ips",
        },
    }

    mask = miri_lrs_wavelength_mask(
        conf, np.array([12.4, 12.5, 12.51])
    )

    assert mask.tolist() == [True, True, False]


@pytest.mark.parametrize(
    "conf",
    [
        {
            "instrument": {
                "instrument": "miri",
                "mode": "lrsslitless",
            },
            "detector": {"subarray": "slitlessprism"},
        },
        {
            "instrument": {
                "instrument": "miri",
                "mode": "lrsslitless",
            },
            "detector": {"subarray": "slitlessprism_ip"},
        },
        {
            "instrument": {
                "instrument": "miri",
                "mode": "lrsslit",
            },
            "detector": {"subarray": "slitlessprism_ips"},
        },
        {
            "instrument": {
                "instrument": "nircam",
                "mode": "dhs",
            },
            "detector": {"subarray": "slitlessprism_ips"},
        },
    ],
)
def test_miri_lrs_ips_mask_only_applies_to_ips_slitless(conf):
    assert miri_lrs_wavelength_mask(conf, np.array([11.0, 13.0])) is None
