"""Уровни доверия I/II и условия прецедента (решение 0007, R11)."""

from chem_agent.evidence import EvidenceIndex


def _index(tmp_path):
    corpus = tmp_path / "corpus.csv"
    corpus.write_text(
        "id,rxn_smiles,source\nAB-1,CCCCCCCCCCCCO.C1CO1.[Na+].[OH-]>>CCCCCCCCCCCCOCCO,US1\n",
        encoding="utf-8",
    )
    meta = tmp_path / "meta.csv"
    meta.write_text(
        "id,temperature_c_min,temperature_c_max,solvent\nAB-1,140,160,\n", encoding="utf-8"
    )
    manual = tmp_path / "manual.csv"
    manual.write_text(
        "id,rxn_smiles,source,note\nMAN-S01,S1SSSSSSS1.O=O>>O=S=O,Сжигание серы,\n",
        encoding="utf-8",
    )
    return EvidenceIndex.build(str(corpus), str(meta), str(manual))


def test_known_reaction_is_level_one(tmp_path):
    idx = _index(tmp_path)
    assert idx.level(["CCCCCCCCCCCCO", "C1CO1"], "CCCCCCCCCCCCOCCO") == ("I", "AB-1", "US1")
    assert idx.level(["S1SSSSSSS1", "O=O"], "O=S=O")[0] == "I"  # ручная библиотека


def test_analog_is_level_two(tmp_path):
    idx = _index(tmp_path)
    assert idx.level(["CCCCCCCCCCCCO", "CC1CO1"], "CCCCCCCCCCCCOCC(C)O")[0] == "II"


def test_precedent_conditions(tmp_path):
    idx = _index(tmp_path)
    assert idx.conditions["AB-1"] == {"T мин, °C": "140", "T макс, °C": "160"}
