import contextlib
import io
import os

import numpy as np
import pytest

KNOWN_PANDEIA_ZERO_NOISE_WARNINGS = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in divide:RuntimeWarning:pandeia\\.engine\\.report"
    ),
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in divide:RuntimeWarning:pandeia\\.engine\\.projection"
    ),
]


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
    _require_readable_phoenix_grid()


def _require_readable_phoenix_grid():
    refdata = os.environ.get("PYSYN_CDBS")
    path = None
    if refdata is not None:
        path = os.path.join(refdata, "grid", "phoenix", "catalog.fits")
    if path is None or not os.path.exists(path):
        pytest.skip(f"PYSYN_CDBS PHOENIX reference grid is unavailable: {path}")
    try:
        with open(path, "rb") as handle:
            handle.read(1)
    except OSError as exc:
        pytest.skip(f"PYSYN_CDBS PHOENIX reference grid is unreadable: {exc}")


def _require_fortney_grid():
    path = os.environ.get("FORTGRID_DIR")
    if path is None or not os.path.isfile(path):
        pytest.skip(f"FORTGRID_DIR is unavailable: {path}")
    try:
        with open(path, "rb") as handle:
            handle.read(1)
    except OSError as exc:
        pytest.skip(f"FORTGRID_DIR is unreadable: {exc}")


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


@pytest.mark.parametrize(
    "instrument",
    [
        pytest.param(
            "NIRSpec G140H",
            marks=KNOWN_PANDEIA_ZERO_NOISE_WARNINGS,
        ),
        "MIRI LRS",
    ],
)
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


def test_run_pandexo_wasp12b_grid_nirspec_prism_multistripe():
    """Exercise the Fortney-grid path with the NIRSpec PRISM multistripe mode."""
    _require_valid_pandeia_refdata()
    _require_fortney_grid()
    jdi = _import_justdoit()

    exo_dict = _default_smoke_exo_dict(jdi)
    exo_dict["star"].update(
        {
            "mag": 10.19,
            "temp": 6300,
            "metal": 0.3,
            "logg": 4.2,
            "radius": 1.63,
        }
    )
    exo_dict["planet"].update(
        {
            "type": "grid",
            "temp": 2500,
            "chem": "noTiO",
            "cloud": "ray10",
            "mass": 1.4,
            "m_unit": "M_jup",
            "radius": 1.9,
            "transit_duration": 3.0 * 60.0 * 60.0,
        }
    )
    inst_dict = jdi.load_mode_dict("NIRSpec Prism")
    inst_dict["configuration"]["detector"].update(
        {
            "subarray": "s64m8_prm",
            "readout_pattern": "nrsrapid",
            "ngroup": "optimize",
            "nint": 1,
        }
    )

    result = jdi.run_pandexo(
        exo_dict, inst_dict, save_file=False, verbose=False
    )

    final_spectrum = result["FinalSpectrum"]
    wave = np.asarray(final_spectrum["wave"])
    assert len(wave) > 0
    assert np.all(np.diff(wave) >= 0)
    assert wave[-1] <= 4.98
    assert len(final_spectrum["spectrum"]) == len(wave)
    assert len(final_spectrum["error_w_floor"]) == len(wave)
