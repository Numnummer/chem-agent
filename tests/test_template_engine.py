"""
Регрессионные тесты химического движка.

Каждый тест фиксирует решение или найденный дефект (см. docs/decisions/).
Если тест падает после изменения кода — это химическая регрессия, а не
«устаревший тест»: разберитесь в причине, прежде чем править ожидания.
"""

import csv
import json

from conftest import INPUTS, PURCHASABLE, write_inputs

from chem_agent.template_engine import canon, frag_set, type_key

with open(INPUTS, encoding="utf-8") as f:
    USER_SET = [r["smiles"] for r in csv.DictReader(f) if r["smiles"] != "resource"]


def products(rows):
    return {r["product"] for r in rows}


def load_templates(path):
    return [json.loads(line) for line in open(path, encoding="utf-8")]


# --- представление веществ -------------------------------------------------


def test_sulfur_smiles_is_not_hydrogen_sulfide():
    """SMILES 'S' — это H2S; элементарную серу пишем как S8."""
    assert canon("S") != canon("S1SSSSSSS1")


def test_counter_ions_ignored_in_fragment_sets():
    assert frag_set("[Na+].[OH-]") == frozenset({"[OH-]"})


# --- извлечение шаблонов ---------------------------------------------------


def test_all_templates_reproduce_their_precedents(templates_path):
    for t in load_templates(templates_path):
        assert t["selfcheck_rate"] == 1.0, t["id"]


def test_manual_templates_are_trusted(templates_path):
    trusted = [t for t in load_templates(templates_path) if t.get("trusted")]
    assert len(trusted) >= 7


def test_type_key_merges_same_transformation():
    """Спирт + эпоксид: общий и специфичный (метанол) шаблоны — один тип."""
    general = (
        "[C:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1>>"
        "[C:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]"
    )
    methyl = (
        "[C;D1;H3:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1>>"
        "[C;D1;H3:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]"
    )
    assert type_key(general) == type_key(methyl)


# --- применение ------------------------------------------------------------


def test_propylene_oxide_opens_at_ch2(run_apply):
    """Региохимия переносится из прецедентов: атака по CH2, вторичный спирт."""
    rows = run_apply("--purchasable", str(PURCHASABLE))
    prods = products(r for r in rows if r["input_name"] == "Оксид пропилена")
    assert canon("COCC(C)O") in prods
    assert canon("COC(C)CO") not in prods


def test_roundtrip_rejects_epoxide_deoxygenation(run_apply):
    """Слишком общий шаблон этерификации не должен «раскрывать» эпоксид с потерей O."""
    rows = run_apply("--purchasable", str(PURCHASABLE))
    bad = canon("CCCOC(C)=O")  # пропилацетат из ПО + AcOH — химическая ошибка
    good = canon("CC(=O)OCC(C)O")  # гидроксипропилацетат — правильный продукт
    prods = products(r for r in rows if r["input_name"] == "Оксид пропилена")
    assert bad not in prods
    assert good in prods


def test_roundtrip_off_lets_false_positive_through(run_apply):
    """Контроль: без обратной проверки ложный кандидат появляется — проверка нужна."""
    rows = run_apply("--purchasable", str(PURCHASABLE), "--no-roundtrip")
    assert canon("CCCOC(C)=O") in products(rows)


def test_naoh_participates_as_reagent(run_apply):
    rows = run_apply("--purchasable", str(PURCHASABLE))
    naoh = [r for r in rows if r["input_name"] == "Гидроксид натрия"]
    assert any(r["role"] == "reagent" for r in naoh)


def test_neutralization_requires_base_in_set(run_apply, tmp_path):
    """Без щёлочи в наборе кислота не должна «нейтрализоваться сама»."""
    no_base = write_inputs(
        tmp_path / "no_base.csv", [("Додеканол", "CCCCCCCCCCCCO"), ("Кислород", "O=O")]
    )
    rows = run_apply("--depth", "2", "--internal-only", inputs=no_base)
    assert canon("CCCCCCCCCCCC(=O)O") in products(rows)  # окисление есть
    assert canon("CCCCCCCCCCCC(=O)[O-]") not in products(rows)  # нейтрализации нет


def test_internal_only_uses_only_the_set(run_apply):
    rows = run_apply("--internal-only")
    assert rows
    assert {r["partner_source"] for r in rows} <= {"internal", "internal+intermediate"}


def test_internal_label_means_all_participants_from_the_set(run_apply):
    """На стадиях ≥2 входное вещество — промежуточный продукт. Реакция с ним
    не может считаться «все участники из исходного набора» (критерий 2.1)."""
    rows = run_apply("--depth", "3", "--internal-only")
    user = frozenset().union(*(frag_set(s) for s in USER_SET))
    internal = [r for r in rows if r["partner_source"] == "internal"]
    assert internal
    for r in internal:
        reac, reag, _ = r["reaction_smiles"].split(">")
        assert frag_set(reac + ("." + reag if reag else "")) <= user, r["reaction_smiles"]
    # додеканол + SO3: SO3 получен из серы за две стадии, это не «внутри набора»
    sulfation = [r for r in rows if r["product"] == canon("CCCCCCCCCCCCOS(=O)(=O)O")]
    assert sulfation
    assert all(r["partner_source"] == "internal+intermediate" for r in sulfation)


def test_sles_route_found_from_priority_set(run_apply):
    """Из перечня ТЗ за 5 стадий собирается лауретсульфат натрия (SLES)."""
    rows = run_apply("--depth", "5", "--internal-only")
    assert canon("CCCCCCCCCCCCOCCOS(=O)(=O)[O-]") in products(rows)


def test_every_candidate_has_precedent(run_apply):
    rows = run_apply("--purchasable", str(PURCHASABLE))
    assert all(r["precedent_id"] for r in rows)
