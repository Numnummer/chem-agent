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
    # внутренний перенос протона: HSO3- + ЭО -> изэтионат (OH бисульфита отдаёт
    # H кислороду раскрытого эпоксида) — основание не нужно
    (
        "[C:1]1-[CH2;D2;+0:2]-[O;H0;D2;+0:3]-1.[O-;H0;D1:4]-[S;H0;D3;+0:5](=[O;D1;H0:6])-[OH;D1;+0:7]"
        ">>[O-;H0;D1:7]-[S;H0;D4;+0:5](=[O;D1;H0:6])(=[O;H0;D1;+0:4])-[CH2;D2;+0:2]-[C:1]-[OH;D1;+0:3]",
        set(),
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


# --- R17: хемоселективность ---------------------------------------------------


@pytest.mark.parametrize(
    "smiles, atom_symbol, expected",
    [
        ("CC(=O)O", "O", {"acid_O", "carbonyl_O"}),
        ("CCO", "O", {"alcohol"}),
        ("Oc1ccccc1", "O", {"phenol"}),
        ("CC(N)=O", "N", {"amide_N"}),
        ("CN", "N", {"amine"}),
        ("Nc1ccccc1", "N", {"aniline"}),
        ("CCS", "S", {"thiol"}),
        ("NC(N)=S", "N", {"amide_N"}),  # тиомочевина: N не аминный
    ],
)
def test_atom_class(smiles, atom_symbol, expected):
    from rdkit import Chem

    from chem_agent.template_engine import atom_class

    m = Chem.MolFromSmiles(smiles)
    classes = {atom_class(m, a.GetIdx()) for a in m.GetAtoms() if a.GetSymbol() == atom_symbol}
    assert classes == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Clc1ccc(O)cc1", "aryl_X"),  # строка 48 ревью full-v1: SNAr без активации
        ("O=[N+]([O-])c1ccc(Cl)c([N+](=O)[O-])c1", "aryl_X_activated"),
        ("Clc1ccccn1", "aryl_X_activated"),
    ],
)
def test_aryl_halide_activation(smiles, expected):
    from rdkit import Chem

    from chem_agent.template_engine import atom_class

    m = Chem.MolFromSmiles(smiles)
    c = next(a for a in m.GetAtoms() if any(n.GetSymbol() == "Cl" for n in a.GetNeighbors()))
    assert atom_class(m, c.GetIdx()) == expected


ALKOXYLATION = (
    "[C:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1"
    ">>[C:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]"
)


@pytest.mark.parametrize(
    "alcohol, allowed, why",
    [
        ("CCCCCCCCCCCCO", True, "спирт, как в прецедентах"),
        ("OCCO", True, "вторая OH — равный конкурент, не сильнее"),
        ("CC(=O)O", False, "OH кислоты — не спирт (строка 40 ревью full-v1)"),
        ("OCCS", False, "SH сильнее OH (строка 32 ревью full-v2)"),
        ("NCCO", False, "алифатический NH2 сильнее OH"),
    ],
)
def test_alkoxylation_selectivity(alcohol, allowed, why):
    from rdkit import Chem

    from chem_agent.template_engine import Template, run_products_traced

    ex = {
        "id": "X",
        "reactants": ["CCCCO", "C1CO1"],
        "product": "CCCCOCCO",
        "spectators": [],
        "center": {"2": ["alcohol", 0], "5": ["ether", 0]},
    }
    t = Template("T", ALKOXYLATION, ALKOXYLATION, 5, 1.0, [ex], "k").build()
    mols = [Chem.MolFromSmiles(alcohol), Chem.MolFromSmiles("C1CO1")]
    traced = run_products_traced(t.rxn, mols)
    assert traced, "шаблон должен совпасть структурно"
    assert any(t.selective(p, mols) for p in traced.values()) == allowed, why


def test_nitro_oxygen_is_not_a_nucleophile():
    """[O-] нитрогруппы — не алкоксид: иначе любое нитросоединение «содержит
    сильный нуклеофил» и верные реакции с ним отклоняются."""
    from rdkit import Chem

    from chem_agent.template_engine import atom_class, competitor_rank

    m = Chem.MolFromSmiles("O=[N+]([O-])c1ccc(CO)cc1")
    o_minus = next(a for a in m.GetAtoms() if a.GetFormalCharge() == -1)
    assert atom_class(m, o_minus.GetIdx()) != "alkoxide"
    alcohol_o = next(a for a in m.GetAtoms() if a.GetSymbol() == "O" and a.GetTotalNumHs() == 1)
    assert competitor_rank(m, alcohol_o.GetIdx()) == 0


def test_substitution_center_is_changed():
    """SNAr: у атома кольца меняется сосед (Cl -> O), а не H и степень.
    Такой атом — центр, его класс (активирован ли арилгалогенид) проверяется."""
    from rdkit.Chem import AllChem

    from chem_agent.template_engine import changed_mapnos

    snar = (
        "Cl-[c;H0;D3;+0:1](:[c:2]):[c:3].[C:4]-[OH;D1;+0:5]"
        ">>[C:4]-[O;H0;D2;+0:5]-[c;H0;D3;+0:1](:[c:2]):[c:3]"
    )
    assert changed_mapnos(AllChem.ReactionFromSmarts(snar)) == {1, 5}


# --- R16: миграция гетероатома — признак ошибочной записи корпуса -------------


@pytest.mark.parametrize(
    "forward, migrates, why",
    [
        (
            "O-[CH;D3;+0:1](-[C:2])-[CH3;D1;+0:3].[O;D1;H0:4]=[C:5]-[OH;D1;+0:6]"
            ">>[C:2]-[CH2;D2;+0:1]-[CH2;D2;+0:3]-[O;H0;D2;+0:6]-[C:5]=[O;D1;H0:4]",
            True,
            "T01991: 2-додеканол -> эфир 1-додеканола (строки 22, 28, 29 ревью full-v2)",
        ),
        (ALKOXYLATION, False, "раскрытие эпоксида: связь с O меняется у одного углерода"),
        (
            "Br-[CH2;D2;+0:1]-[C:2].[OH-;D0:3]>>[C:2]-[CH2;D2;+0:1]-[OH;D1;+0:3]",
            False,
            "SN2",
        ),
        (
            "[C:1]=[CH2;D1;+0:2].[OH2;D0;+0:3]>>[C:1](-[OH;D1;+0:3])-[CH3;D1;+0:2]",
            False,
            "гидратация алкена: присоединение по двойной связи",
        ),
    ],
)
def test_heteroatom_migration(forward, migrates, why):
    from rdkit.Chem import AllChem

    from chem_agent.template_engine import heteroatom_migration

    assert heteroatom_migration(AllChem.ReactionFromSmarts(forward)) == migrates, why


# --- Функции реагента: кислоты Льюиса, катализаторы, окислитель по делу -------


@pytest.mark.parametrize(
    "spectators, expected",
    [
        (["Cl[Al](Cl)Cl"], {"acid"}),
        (["[Al+3]", "[Cl-]", "[Cl-]", "[Cl-]"], {"acid"}),
        (["[Al+3]", "[H-]", "[Li+]"], {"reductant"}),  # LiAlH4, не кислота Льюиса
        (["FB(F)F"], {"acid"}),
        (["[Pd]"], {"metal_catalyst"}),
        (["[Pt]", "[Na+]", "[OH-]"], {"metal_catalyst", "base"}),
    ],
)
def test_lewis_acids_and_catalysts(spectators, expected):
    from chem_agent.template_engine import reagent_functions

    assert reagent_functions(spectators) == expected


FRIEDEL_CRAFTS = (
    "Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[cH;D2;+0:4](:[c:5]):[c:6]"
    ">>[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[c;H0;D3;+0:4](:[c:5]):[c:6]"
)
ESTERIFICATION = (
    "[C:1]-[C;H0;D3;+0:2](=[O;D1;H0:3])-[OH;D1;+0:4].[c:5]-[OH;D1;+0:6]"
    ">>[C:1]-[C;H0;D3;+0:2](=[O;D1;H0:3])-[O;H0;D2;+0:6]-[c:5]"
)


def _template(forward, spectator_lists):
    from chem_agent.template_engine import Template

    exs = [
        {"id": f"X{i}", "reactants": [], "product": "", "spectators": sp}
        for i, sp in enumerate(spectator_lists)
    ]
    return Template("T", forward, forward, len(exs), 1.0, exs, "k").build()


def test_naoh_not_reagent_where_lewis_acid_needed():
    """Строка 47 ревью full-v3: ацилирование по Фриделю–Крафтсу — нужен
    AlCl3; NaOH в прецедентах — только обработка."""
    from rdkit import Chem

    t = _template(FRIEDEL_CRAFTS, [["Cl[Al](Cl)Cl", "[Na+]", "[OH-]"]] * 3)
    assert not t.can_be_reagent(NAOH, Chem.MolFromSmiles(NAOH))


def test_oxidant_not_reagent_in_non_redox_reaction():
    """Строка 54 ревью full-v3: этерификация не окисление — O2, записанный в
    одном прецеденте из трёх (атмосфера), реагентом не становится."""
    from rdkit import Chem

    t = _template(ESTERIFICATION, [["O=O"], [], []])
    assert not t.can_be_reagent("O=O", Chem.MolFromSmiles("O=O"))


def test_consistent_oxidant_reagent_kept():
    """Строка 21 ревью full-v1: S8 в 3 из 4 прецедентов (синтез имидазолинов) —
    реагент по делу, хотя превращение формально не окислительное."""
    from rdkit import Chem

    t = _template(ESTERIFICATION, [[S8], [S8], [S8], []])
    assert t.can_be_reagent(S8, Chem.MolFromSmiles(S8))


def test_catalyst_does_not_block_and_is_reported():
    """Катализатор — условие процесса, не сырьё: реакция не блокируется, а
    катализатор из прецедентов указывается."""
    t = _template(ESTERIFICATION, [["[Pt]"], ["[Pt]"], []])
    assert not t.needs_reagent
    assert t.catalysts == ["[Pt]"]


def test_new_stereocenter_not_invented():
    """Строки 33, 52 ревью full-v3: стереоцентр из шаблона не переносится на
    атом, который в исходном веществе стереоцентром не был."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from chem_agent.template_engine import run_products_traced

    reduction = AllChem.ReactionFromSmarts(
        "[C:1]-[C;H0;D3;+0:2](-[C:3])=[O;H0;D1;+0:4]>>[C:1]-[C@H;D3;+0:2](-[C:3])-[OH;D1;+0:4]"
    )
    prods = run_products_traced(reduction, [Chem.MolFromSmiles("CCC(C)=O")], keep_new_stereo=False)
    assert prods and not any("@" in p for p in prods)


def test_new_stereocenter_kept_on_chiral_substrate():
    """Строка 40 ревью full-v2 (OK): ене-реакция Шенка на стероиде — новый
    стереоцентр задаётся хиральным каркасом, его не снимаем."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from chem_agent.template_engine import run_products_traced

    reduction = AllChem.ReactionFromSmarts(
        "[C:1]-[C;H0;D3;+0:2](-[C:3])=[O;H0;D1;+0:4]>>[C:1]-[C@H;D3;+0:2](-[C:3])-[OH;D1;+0:4]"
    )
    prods = run_products_traced(
        reduction, [Chem.MolFromSmiles("C[C@H](O)CC(C)=O")], keep_new_stereo=False
    )
    assert any(p.count("@") == 2 for p in prods)


# --- Активация OH как уходящей группы -----------------------------------------

DIOL_TO_EPOXIDE = "O-[CH2;D2;+0:1]-[C:2]-[OH;D1;+0:3]>>[C:2]1-[CH2;D2;+0:1]-[O;H0;D2;+0:3]-1"


@pytest.mark.parametrize(
    "spectators",
    [["CS(=O)(=O)Cl"], ["O=S(Cl)Cl"], ["ClP(Cl)(Cl)=O"], ["c1ccc(P(c2ccccc2)c2ccccc2)cc1"]],
)
def test_activators_classified(spectators):
    from chem_agent.template_engine import reagent_functions

    assert "activator" in reagent_functions(spectators)


def test_hydroxyl_leaving_sp3_needs_activation():
    """Строки 14, 20 ревью full-v3: диол + NaOH -> эпоксид (T01070). OH не
    уходит с sp3-углерода без активации (MsCl в прецеденте); NaOH её не даёт."""
    from rdkit import Chem

    assert "activation" in intrinsic_functions(DIOL_TO_EPOXIDE)
    t = _template(DIOL_TO_EPOXIDE, [["CS(=O)(=O)Cl", "[Na+]", "[OH-]"], ["[Na+]", "[OH-]"]])
    assert not t.can_be_reagent(NAOH, Chem.MolFromSmiles(NAOH))
    # в прецедентах и активатор (MsCl), и основание: одно вещество их не закрывает
    from chem_agent.template_engine import covers

    assert not covers(t.required_functions, {"activator"})
    assert covers(t.required_functions, {"activator", "base"})


def test_epoxide_opening_and_esterification_need_no_activation():
    assert "activation" not in intrinsic_functions(ALKOXYLATION)
    assert "activation" not in intrinsic_functions(ESTERIFICATION)
