"""Тесты подготовки корпуса all_balanced_reactions (docs/decisions/0008)."""

import csv

from chem_agent.corpus_prep import ELEMENTAL_SULFUR, NO_SOURCE, main, prepare


def raw(record_id, rxn, patent="", subset="ELEMENTAL_SULFUR_by_SMILES_deep", **extra):
    row = {
        "record_id": record_id,
        "source_subset": subset,
        "reaction_smiles": rxn,
        "reaction_smiles_balanced": "",
        "agents_smiles": "",
        "patent": patent,
    }
    row.update(extra)
    return row


def by_id(rows):
    return {r["id"]: r for r in rows}


def test_single_arrow_format_becomes_engine_format():
    corpus, _, _ = prepare([raw("1", "CC#N.[S]>CC(N)=S", agents_smiles="CCO")])
    assert corpus[0]["rxn_smiles"].count(">") == 2
    reac, agents, prod = corpus[0]["rxn_smiles"].split(">")
    assert agents == "CCO"
    assert prod == "CC(N)=S"


def test_elemental_sulfur_atom_becomes_s8():
    """[S] в патентах — элементарная сера; движок и пользователь пишут S8."""
    corpus, meta, _ = prepare([raw("1", "CC#N.[S]>CC(N)=S")])
    assert ELEMENTAL_SULFUR in corpus[0]["rxn_smiles"].split(">")[0].split(".")
    assert meta[0]["sulfur_normalized"] == 1


def test_hydrogen_sulfide_is_not_turned_into_sulfur():
    """'S' — это H2S, а не сера. Заменять его нельзя."""
    corpus, meta, _ = prepare([raw("1", "CC#N.S>CC(N)=S")])
    assert "S" in corpus[0]["rxn_smiles"].split(">")[0].split(".")
    assert ELEMENTAL_SULFUR not in corpus[0]["rxn_smiles"]
    assert meta[0]["sulfur_normalized"] == 0


def test_duplicates_merged_and_patent_preferred_as_source():
    rows = [
        raw("1", "CCO.C1CO1>>CCOCCO", subset="2naoh_dataset"),  # без патента
        raw("2", "C1CO1.OCC>OCCOCC", patent="US123"),  # та же реакция
    ]
    corpus, meta, stats = prepare(rows)
    assert len(corpus) == 1
    assert stats["duplicates_merged"] == 1
    m = meta[0]
    assert m["n_records"] == 2
    assert m["record_ids"] == "1|2"
    assert m["source"] == "US123"
    assert corpus[0]["id"] == "AB-2"  # представитель — запись с патентом


def test_reaction_without_patent_marked_no_source():
    corpus, meta, _ = prepare([raw("7", "CCO.C1CO1>>CCOCCO", subset="2naoh_dataset")])
    assert corpus[0]["source"] == NO_SOURCE
    assert meta[0]["has_source"] == 0


def test_unparsable_and_bad_format_are_counted_not_silently_lost():
    rows = [raw("1", "C1CC>>CC"), raw("2", "CCO"), raw("3", "CCO.C1CO1>>CCOCCO")]
    corpus, _, stats = prepare(rows)
    assert len(corpus) == 1
    assert stats["skip_unparsable"] == 1
    assert stats["skip_bad_format"] == 1


def test_balanced_column_mismatch_is_flagged():
    ok = raw("1", "CCO.C1CO1>CCOCCO", reaction_smiles_balanced="CCO.C1CO1>>CCOCCO")
    bad = raw("2", "CCCO.C1CO1>CCCOCCO", reaction_smiles_balanced="CC#N.S>>CC(N)=S")
    _, meta, _ = prepare([ok, bad])
    flags = {m["record_ids"]: m["balanced_consistent"] for m in meta}
    assert flags == {"1": 1, "2": 0}


def test_subset_filter():
    rows = [
        raw("1", "CCO.C1CO1>>CCOCCO", subset="NAOH_reactions_reactant_strict"),
        raw("2", "CCCO.C1CO1>>CCCOCCO", subset="ETHYLENE_OXIDE_reactions_reactant_strict"),
    ]
    corpus, _, stats = prepare(rows, exclude=["NAOH", "2naoh"])
    assert [r["id"] for r in corpus] == ["AB-2"]
    assert stats["skip_subset"] == 1


def test_cli_writes_corpus_meta_and_report(tmp_path):
    src = tmp_path / "raw.csv"
    fields = ["record_id", "source_subset", "reaction_smiles", "patent"]
    with open(src, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(
            {
                "record_id": "1",
                "source_subset": "O2",
                "reaction_smiles": "[S].O=O>O=S=O",
                "patent": "US1",
            }
        )
    out = tmp_path / "out"
    main(["--raw", str(src), "--out-dir", str(out)])
    corpus = list(csv.DictReader(open(out / "corpus.csv", encoding="utf-8")))
    assert corpus[0]["rxn_smiles"] == f"{ELEMENTAL_SULFUR}.O=O>>O=S=O"
    assert (out / "meta.csv").exists()
    assert "Уникальных реакций" in (out / "prep_report.md").read_text(encoding="utf-8")
