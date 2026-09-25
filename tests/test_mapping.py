"""Медленные тесты: требуют rxnmapper. Запуск: make test-all."""

import csv

import pytest

from chem_agent.template_engine import main


@pytest.mark.slow
def test_map_produces_atom_maps(tmp_path):
    pytest.importorskip("rxnmapper")
    corpus = tmp_path / "c.csv"
    corpus.write_text("id,rxn_smiles\nX1,CCCCCCCCCCCCO.C1CO1>>CCCCCCCCCCCCOCCO\n", encoding="utf-8")
    out = tmp_path / "m.csv"
    main(["map", "--corpus", str(corpus), "--out", str(out)])
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    assert len(rows) == 1
    assert ":1]" in rows[0]["mapped_rxn"]
