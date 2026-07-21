import inspect

import pandexo.engine.justdoit as jdi


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
