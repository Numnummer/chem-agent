"""
Регрессия по ревью прогона full-v1 (outputs/2026-09-27-full-v1/review_sample.csv,
docs/findings.md): роль «сопутствующий реагент» (R12) и дубли по средней части (R13).

Шаблоны — реальные, из прогона full-v1 (tests/regression/r1_templates.jsonl):
T00136 спирт → альдегид (в прецедентах окислитель Десса–Мартина), T02167 и
T06628 эпоксид → спирт (в прецедентах гидриды), T00336 спирт + эпоксид,
T03405 спирт + O2 → карбоксилат (в прецедентах NaOH), T01782 и T00791 — с S8
в прецедентах, T00317 ацилирование, T15666 сжигание серы.
"""

import csv
import re
from pathlib import Path

import pytest
from conftest import INPUTS

from chem_agent.template_engine import canon, main

TEMPLATES = Path(__file__).parent / "regression" / "r1_templates.jsonl"
NAOH, S8 = "[Na+].[OH-]", "S1SSSSSSS1"
C12OH, EO, PO = "CCCCCCCCCCCCO", "C1CO1", "CC1CO1"


def run(tmp_path, *extra):
    out, report = tmp_path / "c.csv", tmp_path / "r.md"
    main(
        ["apply", "--templates", str(TEMPLATES), "--inputs", str(INPUTS)]
        + ["--min-count", "1", "--jobs", "1", "--out", str(out), "--report", str(report)]
        + list(extra)
    )
    with open(out, encoding="utf-8") as f:
        return list(csv.DictReader(f)), report.read_text(encoding="utf-8")


def key(reactants, reagent, product):
    return (frozenset(canon(s) for s in reactants), reagent, canon(product))


def keys(rows):
    out = set()
    for r in rows:
        reac, reag, prod = r["reaction_smiles"].split(">")
        out.add(key(reac.split("."), reag, prod))
    return out


@pytest.fixture
def internal(tmp_path):
    return run(tmp_path, "--internal-only")


# --- R12: вещество набора — сопутствующий реагент только по делу -------------


@pytest.mark.parametrize(
    "reactants, reagent, product, why",
    [
        ([C12OH], S8, "CCCCCCCCCCCC=O", "строка 1: S8 не окисляет спирт"),
        ([C12OH], NAOH, "CCCCCCCCCCCC=O", "строка 16: NaOH не окислитель"),
        ([PO], NAOH, "CC(C)O", "строка 9: восстановление без восстановителя"),
        ([PO], NAOH, "CCCO", "строка 10: то же"),
        ([EO], NAOH, "CCO", "строка 15: то же"),
    ],
)
def test_reagent_without_needed_function_rejected(internal, reactants, reagent, product, why):
    rows, _ = internal
    assert key(reactants, reagent, product) not in keys(rows), why


def test_competing_epoxide_is_not_auxiliary_reagent(internal):
    """Строки 6, 7, 12, 13: ЭО/ПО сами подходят в слот «спирт + эпоксид» —
    это конкурирующий алкилирующий агент или избыток, а не вспомогательное
    вещество. Эпоксид как акцептор кислоты (строка 31) не моделируем."""
    rows, _ = internal
    assert not [r for r in rows if r["reagent"] in (canon(EO), canon(PO))]


@pytest.mark.parametrize(
    "reactants, reagent, product",
    [
        ([C12OH, PO], NAOH, "CCCCCCCCCCCCOCC(C)O"),  # строка 17: щелочное пропоксилирование
        ([C12OH, EO], NAOH, "CCCCCCCCCCCCOCCO"),  # строка 18: оксиэтилирование
        ([C12OH, "O=O"], NAOH, "CCCCCCCCCCCC(=O)[O-]"),  # строка 19: окисление в щёлочи
        ([S8, "O=O"], "", "O=S=O"),  # строка 2: сжигание серы
        ([C12OH, PO], "", "CCCCCCCCCCCCOCC(C)O"),  # строка 5
    ],
)
def test_correct_reactions_kept(internal, reactants, reagent, product):
    rows, _ = internal
    assert key(reactants, reagent, product) in keys(rows)


def test_sulfur_as_oxidant_kept_where_precedents_use_it(tmp_path):
    """Открытый режим. Строка 21: S8 в прецедентах T01782 (синтез имидазолинов) —
    реагент по делу. Строка 25: T00791 — восстановительное аминирование, S8 там
    случайна (в прецедентах борогидрид)."""
    rows, _ = run(tmp_path)
    s8_reagent = {r["template_id"] for r in rows if r["reagent"] == S8}
    assert "T01782" in s8_reagent
    assert "T00791" not in s8_reagent


# --- R13: одна реакция с разной средней частью — одна реакция ----------------


def test_reagent_variants_counted_once(internal):
    rows, report = internal
    cores = {r["reaction_core"] for r in rows}
    variants = {r["reaction_key"] for r in rows}
    assert len(variants) > len(cores)  # есть варианты с NaOH и без
    reported = int(re.search(r"Уникальных реакций всего: \*\*(\d+)\*\*", report).group(1))
    assert reported == len(cores)


def test_reaction_core_ignores_reagent(internal):
    rows, _ = internal
    for r in rows:
        reac, _, prod = r["reaction_smiles"].split(">")
        assert r["reaction_core"] == ".".join(sorted(reac.split("."))) + ">>" + prod
