import inspect

import pandexo.engine.justdoit as jdi


class _ImmediateParallel:
    """Execute joblib delayed calls synchronously for branch tests."""

    def __init__(self, n_jobs):
        self.n_jobs = n_jobs

    def __call__(self, tasks):
        return [function(*args, **kwargs) for function, args, kwargs in tasks]


def test_load_exo_dict_contains_required_sections():
    exo_dict = jdi.load_exo_dict()

    assert "observation" in exo_dict
    assert "star" in exo_dict
    assert "planet" in exo_dict


def test_public_defaults_are_resolved_at_call_time():
    load_exo_signature = inspect.signature(jdi.load_exo_dict)
    run_signature = inspect.signature(jdi.run_pandexo)

    assert load_exo_signature.parameters["pl_kwargs"].default is None
    assert run_signature.parameters["output_path"].default is None
    assert run_signature.parameters["num_cores"].default is None


def test_run_pandexo_uses_the_call_time_working_directory(
    monkeypatch, tmp_path
):
    """The public default must not capture the directory at import time."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(jdi, "wrapper", lambda *args, **kwargs: {"ok": True})

    result = jdi.run_pandexo({}, {}, verbose=False)

    assert result == {"ok": True}
    assert (tmp_path / "singlerun.p").is_file()


def test_run_pandexo_honors_an_explicit_output_directory(
    monkeypatch, tmp_path
):
    output_path = tmp_path / "results"
    output_path.mkdir()
    monkeypatch.setattr(jdi, "wrapper", lambda *args, **kwargs: {"ok": True})

    jdi.run_pandexo({}, {}, output_path=output_path, verbose=False)

    assert (output_path / "singlerun.p").is_file()


def test_run_pandexo_returns_a_dict_for_single_named_instrument(
    monkeypatch,
):
    monkeypatch.setattr(
        jdi, "load_mode_dict", lambda instrument: {"name": instrument}
    )
    monkeypatch.setattr(
        jdi, "wrapper", lambda *args, **kwargs: {"single": True}
    )

    result = jdi.run_pandexo({}, ["MIRI LRS"], save_file=False, verbose=False)

    assert result == {"single": True}


def test_run_pandexo_returns_keyed_lists_for_multi_run_branches(
    monkeypatch,
):
    monkeypatch.setattr(jdi, "Parallel", _ImmediateParallel)
    monkeypatch.setattr(
        jdi,
        "run_inst_space",
        lambda instrument, exo, verbose=False: {
            instrument: {"instrument": instrument}
        },
    )
    monkeypatch.setattr(
        jdi,
        "run_param_space",
        lambda value, exo, inst, param_space, verbose=False: {
            str(value): {"parameter": value}
        },
    )
    monkeypatch.setattr(jdi, "ALL", {"First": False, "Second": False})

    selected = jdi.run_pandexo(
        {},
        ["NIRSpec G140M", "NIRSpec G235M"],
        save_file=False,
        verbose=False,
    )
    all_modes = jdi.run_pandexo(
        {}, ["RUN ALL"], save_file=False, verbose=False
    )
    parameter_space = jdi.run_pandexo(
        {},
        ["NIRSpec G140M"],
        param_space="star+temp",
        param_range=[5000, 5500],
        save_file=False,
        verbose=False,
    )

    assert selected == [
        {"NIRSpec G140M": {"instrument": "NIRSpec G140M"}},
        {"NIRSpec G235M": {"instrument": "NIRSpec G235M"}},
    ]
    assert all_modes == [
        {"First": {"instrument": "First"}},
        {"Second": {"instrument": "Second"}},
    ]
    assert parameter_space == [
        {"5000": {"parameter": 5000}},
        {"5500": {"parameter": 5500}},
    ]
