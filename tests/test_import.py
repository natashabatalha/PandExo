import pandexo.engine.justdoit as jdi


def test_load_exo_dict_contains_required_sections():
    exo_dict = jdi.load_exo_dict()

    assert "observation" in exo_dict
    assert "star" in exo_dict
    assert "planet" in exo_dict
