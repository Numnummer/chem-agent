"""Балансировщик уравнений (roadmap R4, решение 0015)."""

import csv

import pytest

from chem_agent.balancer import balance, main

NAOH = "[Na+].[OH-]"


def coeffs(result):
    """{SMILES: коэффициент}; реагенты — отрицательные."""
    return {s: -c for s, c in result.reactants} | dict(result.products)


@pytest.mark.parametrize(
    "reactants, reagent, product, expected, why",
    [
        (
            ["CCCCCCCCCCCCO", "C1CO1"],
            "",
            "CCCCCCCCCCCCOCCO",
            {"CCCCCCCCCCCCO": -1, "C1CO1": -1, "CCCCCCCCCCCCOCCO": 1},
            "этоксилирование: присоединение без побочных продуктов",
        ),
        (
            ["CCCCCCCCCCCCO", "O=S(=O)=O"],
            "",
            "CCCCCCCCCCCCOS(=O)(=O)O",
            {"CCCCCCCCCCCCO": -1, "O=S(=O)=O": -1, "CCCCCCCCCCCCOS(=O)(=O)O": 1},
            "сульфатирование SO3",
        ),
        (
            ["S1SSSSSSS1", "O=O"],
            "",
            "O=S=O",
            {"S1SSSSSSS1": -1, "O=O": -8, "O=S=O": 8},
            "сжигание серы: S8 + 8 O2 -> 8 SO2",
        ),
        (
            ["CCCCCCCCCCCCO", "O=O"],
            "",
            "CCCCCCCCCCCC(=O)O",
            {"CCCCCCCCCCCCO": -1, "O=O": -1, "CCCCCCCCCCCC(=O)O": 1, "O": 1},
            "окисление спирта в кислоту: + вода",
        ),
        (
            ["CC(=O)O", "CCO"],
            "",
            "CCOC(C)=O",
            {"CC(=O)O": -1, "CCO": -1, "CCOC(C)=O": 1, "O": 1},
            "этерификация: + вода",
        ),
    ],
)
def test_balanced_without_reagent(reactants, reagent, product, expected, why):
    r = balance(reactants, reagent, product)
    assert r.status == "ok", (why, r.reason)
    assert coeffs(r) == expected, why


def test_neutralization_consumes_naoh():
    """Кислота + NaOH (сопутствующий реагент) -> соль + вода: NaOH расходуется,
    продукт-анион записывается солью натрия."""
    r = balance(["CCCCCCCCCCCC(=O)O"], NAOH, "CCCCCCCCCCCC(=O)[O-]")
    assert r.status == "reagent_consumed", r.reason
    c = coeffs(r)
    assert c[NAOH] == -1 and c["O"] == 1
    assert c["CCCCCCCCCCCC(=O)[O-].[Na+]"] == 1


def test_oxidation_to_carboxylate_with_base():
    """Додеканол + O2 (NaOH) -> лаурат натрия + вода (сеть full-v8, стадия 1)."""
    r = balance(["CCCCCCCCCCCCO", "O=O"], NAOH, "CCCCCCCCCCCC(=O)[O-]")
    assert r.status == "reagent_consumed", r.reason
    c = coeffs(r)
    assert c["O=O"] == -1 and c[NAOH] == -1 and c["O"] == 2


def test_sn2_gives_salt():
    """R-Br + NaOH -> R-OH + NaBr: побочный продукт — соль, а не ионы по отдельности."""
    r = balance(["CCCCBr", NAOH], "", "CCCCO")
    assert r.status == "ok", r.reason
    assert coeffs(r)["[Br-].[Na+]"] == 1


def test_catalytic_reagent_not_consumed():
    """Алкоксилирование в NaOH: основание — катализатор, в уравнение не входит."""
    r = balance(["CCCCCCCCCCCCO", "CC1CO1"], NAOH, "CCCCCCCCCCCCOCC(C)O")
    assert r.status == "ok"
    assert NAOH not in coeffs(r)


def test_boc_removal_gives_isobutylene_and_co2():
    r = balance(["CC(C)(C)OC(=O)NCCO"], "", "NCCO")
    assert r.status == "ok", r.reason
    assert set(r.byproducts) == {"C=C(C)C", "O=C=O"}


def test_unbalanceable_reports_remainder():
    """Уходящий фрагмент не из пула (трифенилметильная защита) — явный отказ с
    остатком, а не молчаливая ошибка."""
    r = balance(["OCCOC(c1ccccc1)(c1ccccc1)c1ccccc1"], "", "OCCO")
    assert r.status == "failed"
    assert "C19H14" in r.reason


def test_ester_hydrolysis_consumes_water():
    """Метиловый эфир + NaOH -> кислота: реагирует вода (ресурс), уходит метанол
    (ревью балансировщика: 10% отказов открытого режима — гидролиз эфиров)."""
    r = balance(["COC(=O)c1ccccc1"], NAOH, "O=C(O)c1ccccc1")
    assert r.status == "hydrolysis", r.reason
    assert r.equation_formula == "C8H8O2 + H2O -> C7H6O2 + CH4O"


def test_acetate_leaving_group():
    r = balance(["CC(=O)OCCCCCCCCCCCC"], NAOH, "CCCCCCCCCCCCO")
    assert r.status == "reagent_consumed", r.reason
    assert "CC(=O)[O-].[Na+]" in r.byproducts


def test_equation_strings():
    r = balance(["S1SSSSSSS1", "O=O"], "", "O=S=O")
    assert r.equation == "S1SSSSSSS1 + 8 O=O -> 8 O=S=O"
    assert r.equation_formula == "S8 + 8 O2 -> 8 O2S"


def test_cli_adds_columns(tmp_path):
    src = tmp_path / "c.csv"
    with open(src, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_smiles", "product"])
        w.writeheader()
        w.writerow({"reaction_smiles": "S1SSSSSSS1.O=O>>O=S=O", "product": "O=S=O"})
    out = tmp_path / "b.csv"
    main(["--in", str(src), "--out", str(out)])
    row = next(csv.DictReader(open(out, encoding="utf-8")))
    assert row["balance_status"] == "ok"
    assert row["equation"] == "S1SSSSSSS1 + 8 O=O -> 8 O=S=O"
    assert row["byproducts"] == ""


def test_anions_are_sodium_salts():
    """Движок хранит промежуточные анионы без противоиона (HSO3-); в уравнении это
    натриевые соли: NaHSO3 + NaOH -> Na2SO3 + H2O (сеть full-v8, стадия 3)."""
    r = balance(["O=S([O-])O"], NAOH, "O=S([O-])[O-]")
    assert r.status == "reagent_consumed", r.reason
    assert r.equation_formula == "HNaO3S + HNaO -> Na2O3S + H2O"


def test_isethionate_from_sodium_bisulfite():
    r = balance(["C1CO1", "O=S([O-])O"], "", "O=S(=O)([O-])CCO")
    assert r.status == "ok"
    assert r.equation_formula == "C2H4O + HNaO3S -> C2H5NaO4S"


def test_oxidant_reagent_consumed_instead_of_h2():
    """Вторичный спирт >O2> кетон: O2 расходуется (+ вода), а не выделяется H2."""
    r = balance(["CC(C)O"], "O=O", "CC(C)=O")
    assert r.status == "reagent_consumed"
    assert r.equation_formula == "2 C3H8O + O2 -> 2 C3H6O + 2 H2O"


def test_hydroxide_acting_as_water():
    """ЭО + NaOH -> ЭГ: OH- здесь катализатор, реагирует вода (замечание ревью
    full-v3, строки 4, 6)."""
    r = balance(["C1CO1", NAOH], "", "OCCO")
    assert r.status == "naoh_catalyst", r.reason
    assert r.equation_formula == "C2H4O + H2O -> C2H6O2"


@pytest.mark.parametrize(
    "reactants, product",
    [
        (["C=CCO", "CC1CO1"], "C=CCOCC(C)O"),  # изомеры C3H6O
        (["CC(=O)COCC(C)O", "C1CO1"], "CC(=O)COCC(C)OCCO"),  # C6H12O3 = 3 x C2H4O
    ],
)
def test_isomeric_or_proportional_reactants(reactants, product):
    """Состав одного реагента кратен составу другого: решение неоднозначно,
    берётся стехиометрия шаблона — по одной молекуле на слот."""
    r = balance(reactants, NAOH, product)
    assert r.status == "ok", r.reason
    assert [c for _, c in r.reactants] == [1, 1]


def test_bromination_with_atmospheric_o2_gives_hbr():
    """O2 как «атмосфера» из прецедента не расходуется, если без него всё
    сходится чисто: RH + Br2 -> RBr + HBr."""
    r = balance(["ClCC(=O)C", "BrBr"], "O=O", "ClCC(=O)CBr")
    assert r.status == "ok", r.reason
    assert r.equation_formula == "C3H5ClO + Br2 -> C3H4BrClO + HBr"


def test_hydrogen_release_is_flagged():
    """Выделение H2 почти всегда значит, что окислитель не указан — помечаем."""
    r = balance(["CC(C)O"], "", "CC(C)=O")
    assert r.status == "ok"
    assert "H2" in r.reason
