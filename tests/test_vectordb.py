"""Векторная БД реакций (ТЗ 2.2): уникальность, поиск похожих, хранение."""

from chem_agent.vectordb import VectorDB

REACTIONS = {
    "eo_c12": ("CCCCCCCCCCCCO.C1CO1>>CCCCCCCCCCCCOCCO", "алкоксилирование"),
    "po_c12": ("CCCCCCCCCCCCO.CC1CO1>>CCCCCCCCCCCCOCC(C)O", "алкоксилирование"),
    "eo_pg": ("CC(O)CO.C1CO1>>CC(O)COCCO", "алкоксилирование"),
    "s_burn": ("S1SSSSSSS1.O=O>>O=S=O", "окисление серы"),
    "contact": ("O=S=O.O=O>>O=S(=O)=O", "окисление серы"),
    "sulf_c12": ("CCCCCCCCCCCCO.O=S(=O)=O>>CCCCCCCCCCCCOS(=O)(=O)O", "сульфатирование"),
    "sulf_c12e1": ("CCCCCCCCCCCCOCCO.O=S(=O)=O>>CCCCCCCCCCCCOCCOS(=O)(=O)O", "сульфатирование"),
    "hydr_eo": ("C1CO1.O>>OCCO", "гидратация"),
    "hydr_po": ("CC1CO1.O>>CC(O)CO", "гидратация"),
    "isom_po": ("CC1CO1>>C=CCO", "изомеризация"),
    "bisulfite": ("O=S=O.[Na+].[OH-]>>O=S([O-])O", "абсорбция SO2"),
}


def _db():
    db = VectorDB()
    for key, (rxn, kind) in REACTIONS.items():
        assert db.add(key, rxn, {"type": kind})
    return db


def test_at_least_ten_unique_vectors():
    """Критерий 2.2: ≥10 реакций, у каждой уникальный вектор."""
    db = _db()
    assert len(db) >= 10
    assert db.unique_vectors() == len(db)


def test_similar_reactions_have_same_type():
    db = _db()
    hits = db.search("CCCCCCCCCCO.C1CO1>>CCCCCCCCCCOCCO", k=2)  # деканол + ЭО
    assert {db.meta[db.keys.index(k)]["type"] for _, k, _ in hits} == {"алкоксилирование"}
    hits = db.search("CCCCCCCCCCO.O=S(=O)=O>>CCCCCCCCCCOS(=O)(=O)O", k=1)
    assert hits[0][1] == "sulf_c12"


def test_duplicate_key_not_added():
    db = _db()
    assert not db.add("eo_c12", REACTIONS["eo_c12"][0])


def test_save_and_load(tmp_path):
    db = _db()
    db.save(str(tmp_path / "db"))
    db2 = VectorDB.load(str(tmp_path / "db"))
    assert db2.keys == db.keys
    assert db2.search(REACTIONS["s_burn"][0], k=1)[0][1] == "s_burn"


def test_charge_state_gives_distinct_vectors():
    """Кислота и её анион, SO3 + H2O и SO3 + OH- — разные реакции, разные векторы."""
    db = VectorDB()
    db.add("acid", "O=C(O)CO.O=S(=O)=O>>O=C(O)COS(=O)(=O)O")
    db.add("anion", "O=C([O-])CO.O=S(=O)=O>>O=C([O-])COS(=O)(=O)O")
    db.add("water", "O.O=S(=O)=O>>O=S(=O)(O)O")
    db.add("hydroxide", "O=S(=O)=O.[Na+].[OH-]>>O=S(=O)([O-])O")
    assert db.unique_vectors() == 4
