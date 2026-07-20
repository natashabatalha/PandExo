import itertools

import pytest

from pandexo.engine.run_online import (
    MIRI_LRS_ALLOWED_SUBARRAYS,
    NIRCAM_DHS_WEB_FILTERS,
    NIRCAM_DHS_WEB_READOUTS,
    NIRCAM_DHS_WEB_SUBARRAYS,
    NIRCAM_WEB_FILTERS,
    NIRCAM_WEB_READOUTS,
    NIRCAM_WEB_SUBARRAYS,
    NIRISS_WEB_SUBARRAYS,
    NIRSPEC_PRISM_MULTISTRIPE_SUBARRAYS,
    NIRSPEC_STANDARD_SUBARRAYS,
    NIRSPEC_WEB_MODES,
    validate_online_instrument_configuration,
)


def _configuration(instrument, detector, **instrument_values):
    return {
        "instrument": {"instrument": instrument, **instrument_values},
        "detector": detector,
    }


def test_all_website_instrument_configurations_are_valid():
    configurations = []
    for mode, subarrays in MIRI_LRS_ALLOWED_SUBARRAYS.items():
        configurations.extend(
            _configuration("miri", {"subarray": subarray}, mode=mode)
            for subarray in subarrays
        )

    for mode in NIRSPEC_WEB_MODES:
        subarrays = NIRSPEC_STANDARD_SUBARRAYS
        if mode == "prismclear":
            subarrays = tuple(
                subarray
                for subarray in NIRSPEC_STANDARD_SUBARRAYS
                if subarray != "sub1024a"
            ) + NIRSPEC_PRISM_MULTISTRIPE_SUBARRAYS
        configurations.extend(
            _configuration(
                "nirspec",
                {"subarray": subarray},
                disperser=mode[:5],
                filter=mode[5:],
            )
            for subarray in subarrays
        )

    configurations.extend(
        _configuration(
            "nircam",
            {"subarray": subarray, "readout_pattern": readout},
            mode="lw_tsgrism",
            filter=filt,
        )
        for filt, subarray, readout in itertools.product(
            NIRCAM_WEB_FILTERS, NIRCAM_WEB_SUBARRAYS, NIRCAM_WEB_READOUTS
        )
    )

    configurations.extend(
        _configuration("niriss", {"subarray": subarray})
        for subarray in NIRISS_WEB_SUBARRAYS
    )

    for display_mode in ("sw_tsgrism", "lw_tsgrism"):
        for sw_filter, lw_filter, subarray, readout in itertools.product(
            NIRCAM_DHS_WEB_FILTERS,
            NIRCAM_WEB_FILTERS,
            NIRCAM_DHS_WEB_SUBARRAYS,
            NIRCAM_DHS_WEB_READOUTS,
        ):
            if display_mode == "sw_tsgrism":
                filt, pair = sw_filter, lw_filter
            else:
                filt, pair = lw_filter, sw_filter
            configurations.append(
                _configuration(
                    "nircam",
                    {"subarray": subarray, "readout_pattern": readout},
                    mode=display_mode,
                    filter=filt,
                    pandexofilterpair=pair,
                )
            )

    for conf in configurations:
        assert validate_online_instrument_configuration(conf) is None


@pytest.mark.parametrize(
    "conf",
    [
        _configuration("miri", {"subarray": "full"}, mode="lrsslitless"),
        _configuration(
            "nirspec",
            {"subarray": "s64m8_prm"},
            disperser="g395h",
            filter="f290lp",
        ),
        _configuration(
            "nircam",
            {"subarray": "", "readout_pattern": "rapid"},
            mode="lw_tsgrism",
            filter="f444w",
        ),
        _configuration("niriss", {"subarray": ""}),
    ],
)
def test_online_instrument_validation_rejects_invalid_choices(conf):
    with pytest.raises(ValueError):
        validate_online_instrument_configuration(conf)
