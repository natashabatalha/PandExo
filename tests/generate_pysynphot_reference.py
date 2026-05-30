"""Generate small legacy regression fixtures from a pysynphot environment.

This script is for maintainers with a working legacy ``pysynphot`` install and
the required ``PYSYN_CDBS`` reference data. It is not used by normal test runs.
The generated fixture is intentionally tiny: it protects against future drift
in representative migration cases, not every possible stellar SED workflow.
"""

import copy
import json
from pathlib import Path

import numpy as np

from pandexo.engine.create_input import bothTrans, outTrans
from test_pysynphot_regression import (
    FPFS_PLANET,
    FPFS_STAR,
    USER_CASES,
    _case_to_numpy,
    _legacy_bothtrans_constant_fpfs,
    _legacy_outtrans,
    _max_relative_difference,
)


def _json_array(values):
    return np.asarray(values, dtype=float).tolist()


def main():
    outtrans = {}
    for name, case in USER_CASES.items():
        legacy = _legacy_outtrans(case)
        new = outTrans(_case_to_numpy(case))
        outtrans[name] = {
            "input": copy.deepcopy(case),
            "wave": _json_array(legacy["wave"]),
            "flux_out_trans": _json_array(legacy["flux_out_trans"]),
            "stellar_flux": _json_array(legacy["stellar_flux"]),
            "max_relative_difference": {
                "flux_out_trans": _max_relative_difference(
                    new["flux_out_trans"], legacy["flux_out_trans"]
                ),
                "stellar_flux": _max_relative_difference(
                    new["stellar_flux"], legacy["stellar_flux"]
                ),
            },
        }

    fpfs_case = USER_CASES["user_j_um_jy"]
    legacy_fpfs = _legacy_bothtrans_constant_fpfs(
        _legacy_outtrans(fpfs_case),
        copy.deepcopy(FPFS_PLANET),
        copy.deepcopy(FPFS_STAR),
    )
    new_fpfs = bothTrans(
        outTrans(_case_to_numpy(fpfs_case)),
        copy.deepcopy(FPFS_PLANET),
        star=copy.deepcopy(FPFS_STAR),
    )

    fixture = {
        "metadata": {
            "generator": "tests/generate_pysynphot_reference.py",
            "description": "Small legacy pysynphot outputs for synphot/stsynphot migration regression tests.",
        },
        "outtrans": outtrans,
        "constant_fpfs": {
            "outtrans_input": copy.deepcopy(fpfs_case),
            "planet": copy.deepcopy(FPFS_PLANET),
            "star": copy.deepcopy(FPFS_STAR),
            "wave": _json_array(legacy_fpfs["wave"]),
            "model_spec": _json_array(legacy_fpfs["model_spec"]),
            "flux_in_trans": _json_array(legacy_fpfs["flux_in_trans"]),
            "max_relative_difference": {
                "model_spec": _max_relative_difference(
                    new_fpfs["model_spec"], legacy_fpfs["model_spec"]
                ),
                "flux_in_trans": _max_relative_difference(
                    new_fpfs["flux_in_trans"], legacy_fpfs["flux_in_trans"]
                ),
            },
        },
    }

    output_path = Path(__file__).parent / "data" / "pysynphot_legacy_reference.json"
    output_path.parent.mkdir(exist_ok=True)
    output_path.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
