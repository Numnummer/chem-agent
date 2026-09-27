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
from conftest import INPUTS, write_inputs

from chem_agent.template_engine import canon, intrinsic_functions, main

TEMPLATES = Path(__file__).parent / "regression" / "r1_templates.jsonl"
NAOH, S8 = "[Na+].[OH-]", "S1SSSSSSS1"
C12OH, EO, PO = "CCCCCCCCCCCCO", "C1CO1", "CC1CO1"


def run(tmp_path, *extra, inputs=INPUTS):
    out, report = tmp_path / "c.csv", tmp_path / "r.md"
    main(
        ["apply", "--templates", str(TEMPLATES), "--inputs", str(inputs)]
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
    ],
)
def test_correct_reactions_kept(internal, reactants, reagent, product):
    rows, _ = internal
    assert key(reactants, reagent, product) in keys(rows)


def test_alkoxylation_proposed_only_with_base(internal):
    """Строки 5 и 17 — одна реакция (R13). В 13 из 20 прецедентов «спирт +
    эпоксид» есть основание (KOtBu, KOH, NaOH, NaH): без катализатора
    алкоксилирование не идёт, поэтому реакция выдаётся только с NaOH
    (решение 0010)."""
    rows, _ = internal
    variants = {
        r["reagent"]
        for r in rows
        if r["reaction_core"] == ".".join(sorted([C12OH, PO])) + ">>CCCCCCCCCCCCOCC(C)O"
    }
    assert variants == {NAOH}


def test_sulfur_as_oxidant_kept_where_precedents_use_it(tmp_path):
    """Открытый режим. Строка 21: S8 в прецедентах T01782 (синтез имидазолинов) —
    реагент по делу. Строка 25: T00791 — восстановительное аминирование, S8 там
    случайна (в прецедентах борогидрид)."""
    rows, _ = run(tmp_path)
    s8_reagent = {r["template_id"] for r in rows if r["reagent"] == S8}
    assert "T01782" in s8_reagent
    assert "T00791" not in s8_reagent


# --- R13: одна реакция с разной средней частью — одна реакция ----------------


def test_reagent_variants_counted_once(tmp_path):
    """Два основания из прецедентов алкоксилирования (NaOH, Et3N) — два
    варианта одной реакции, в отчёте — одна реакция."""
    inputs = write_inputs(
        tmp_path / "in.csv",
        [("Додеканол", C12OH), ("ПО", PO), ("NaOH", NAOH), ("Et3N", "CCN(CC)CC")],
    )
    rows, report = run(tmp_path, "--internal-only", inputs=inputs)
    cores = {r["reaction_core"] for r in rows}
    variants = {r["reaction_key"] for r in rows}
    assert len(variants) > len(cores)
    reported = int(re.search(r"Уникальных реакций всего: \*\*(\d+)\*\*", report).group(1))
    assert reported == len(cores)


def test_reaction_core_ignores_reagent(internal):
    rows, _ = internal
    for r in rows:
        reac, _, prod = r["reaction_smiles"].split(">")
        assert r["reaction_core"] == ".".join(sorted(reac.split("."))) + ">>" + prod


# --- R15: какой реагент нужен самому превращению ------------------------------

INTRINSIC = [
    # строка 15 ревью full-v2: вторичный спирт → кетон, NaOH «окислял»
    (
        "[C:1]-[C@H;D3;+0:2](-[OH;D1;+0:3])-[C:4]>>[C:1]-[C;H0;D3;+0:2](=[O;H0;D1;+0:3])-[C:4]",
        {"oxidant"},
    ),
    # строка 24: карбоксилат → кислота, NaOH «протонировал»
    ("[O-;H0;D1:1]-[C:2]=[O;D1;H0:3]>>[O;D1;H0:3]=[C:2]-[OH;D1;+0:1]", {"acid"}),
    # строка 31: дегидроксилирование фенола — восстановление
    ("O-[c;H0;D3;+0:1](:[c:2]):[c:3]>>[c:2]:[cH;D2;+0:1]:[c:3]", {"reductant"}),
    # строка 42: гидрирование C=C, O2 выписан «восстановителем»
    (
        "[C:1]/[C;H0;D3;+0:2](-[C;D1;H3:3])=[CH;D2;+0:4]/[c:5]"
        ">>[C:1]-[C@@H;D3;+0:2](-[C;D1;H3:3])-[CH2;D2;+0:4]-[c:5]",
        {"reductant"},
    ),
    # строка 30: метилкетон → кислота с потерей C — окисление (галоформ)
    (
        "C-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH-;D0:4]"
        ">>[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[OH;D1;+0:4]",
        {"oxidant"},
    ),
    # спирт + O2 → карбоксилат: окислитель в шаблоне, но продукт — анион
    (
        "O=[O;H0;D1;+0:1].[C:2]-[CH2;D2;+0:3]-[OH;D1;+0:4]"
        ">>[C:2]-[C;H0;D3;+0:3](-[O-;H0;D1:4])=[O;H0;D1;+0:1]",
        {"base"},
    ),
    # без реагента: алкоксилирование, гидратация эпоксида, сжигание серы
    (
        "[C:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1"
        ">>[C:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]",
        set(),
    ),
    (
        "[C:1]1-[CH2;D2;+0:2]-[O;H0;D2;+0:3]-1.[OH-;D0:4]"
        ">>[OH;D1;+0:3]-[C:1]-[CH2;D2;+0:2]-[OH;D1;+0:4]",
        set(),
    ),
    (
        "S1-S-S-S-[S;H0;D2;+0:1]-S-S-S-1.[O;H0;D1;+0:2]=[O;H0;D1;+0:3]"
        ">>[O;H0;D1;+0:2]=[S;H0;D2;+0:1]=[O;H0;D1;+0:3]",
        set(),
    ),
]


@pytest.mark.parametrize("forward, expected", INTRINSIC)
def test_intrinsic_reagent_function(forward, expected):
    assert intrinsic_functions(forward) == expected


def test_reagent_must_cover_intrinsic_function():
    """NaOH (основание) не может быть реагентом окисления, даже если был
    в прецедентах; O2 (окислитель) не может быть реагентом гидрирования."""
    from rdkit import Chem

    from chem_agent.template_engine import Template

    oxidation, _ = INTRINSIC[0]
    hydrogenation, _ = INTRINSIC[3]
    ex = {"id": "X", "reactants": [], "product": "", "spectators": ["[Na+]", "[OH-]", "O=O"]}
    naoh, o2 = Chem.MolFromSmiles(NAOH), Chem.MolFromSmiles("O=O")
    t_ox = Template("T", oxidation, oxidation, 2, 1.0, [ex], "k").build()
    t_red = Template("T", hydrogenation, hydrogenation, 2, 1.0, [ex], "k").build()
    assert not t_ox.can_be_reagent(NAOH, naoh)
    assert not t_red.can_be_reagent("O=O", o2)
