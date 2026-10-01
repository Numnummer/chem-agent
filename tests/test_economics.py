"""Экономика на 1 т продукта (решение 0016)."""

from pathlib import Path

import pytest

from chem_agent.balancer import balance
from chem_agent.economics import Price, PriceBook, evaluate

ROOT = Path(__file__).parent.parent
BOOK = PriceBook(str(ROOT / "data/prices/prices.csv"), str(ROOT / "data/prices/classes.csv"))


def test_consumption_follows_stoichiometry():
    """S8 + 8 O2 -> 8 SO2: на 1 т SO2 — 32.06/64.07 = 0.500 т серы."""
    b = balance(["S1SSSSSSS1", "O=O"], "", "O=S=O")
    eco = evaluate(b.reactants, b.products, BOOK)
    sulfur = next(line for line in eco.raw if line.smiles == "S1SSSSSSS1")
    assert sulfur.t_per_t == pytest.approx(0.5005, abs=1e-3)


def test_yield_increases_consumption():
    b = balance(["S1SSSSSSS1", "O=O"], "", "O=S=O")
    full = evaluate(b.reactants, b.products, BOOK, yield_=1.0)
    ninety = evaluate(b.reactants, b.products, BOOK, yield_=0.9)
    assert ninety.raw_cost == pytest.approx(full.raw_cost / 0.9)


def test_margin_and_product_price():
    b = balance(["C1CO1", "O"], "", "OCCO")
    eco = evaluate(b.reactants, b.products, BOOK)
    assert eco.product_price.basis == "price"  # этиленгликоль в прайсе
    assert eco.margin == pytest.approx(eco.revenue - eco.raw_cost)


def test_class_price_for_unlisted_product():
    """C12E2 нет в прайсе — цена по классу «алкоксилат жирного спирта»."""
    b = balance(["CCCCCCCCCCCCOCCO", "C1CO1"], "[Na+].[OH-]", "CCCCCCCCCCCCOCCOCCO")
    eco = evaluate(b.reactants, b.products, BOOK)
    assert eco.product_price.basis == "class"


def test_route_cost_used_for_intermediate():
    """SO3 нет в прайсе: его цена — себестоимость маршрута из сети."""
    b = balance(["CCCCCCCCCCCCO", "O=S(=O)=O"], "", "CCCCCCCCCCCCOS(=O)(=O)O")
    no_route = evaluate(b.reactants, b.products, BOOK)
    assert no_route.raw_cost is None and "O=S(=O)=O" in no_route.unknown
    eco = evaluate(
        b.reactants,
        b.products,
        BOOK,
        route_cost=lambda s: Price(200.0, "route") if s == "O=S(=O)=O" else None,
    )
    assert eco.raw_cost is not None
    assert next(line for line in eco.raw if line.smiles == "O=S(=O)=O").price.basis == "route"


def test_byproduct_credit():
    """NaHSO4 + NaOH -> Na2SO4 + H2O: вода не кредитуется."""
    b = balance(["CC(=O)OC"], "[Na+].[OH-]", "CO")
    eco = evaluate(b.reactants, b.products, BOOK, credit_share=0.5)
    assert all(line.smiles != "O" for line in eco.byproducts)


def test_make_or_buy_takes_cheaper():
    """SO2 есть в прайсе (350 $/т), но из серы в сети дешевле — берётся меньшая
    (найдено метрикой точности: SO3 292 $/т против эталона 96)."""
    b = balance(["O=S=O", "O=O"], "", "O=S(=O)=O")
    eco = evaluate(
        b.reactants,
        b.products,
        BOOK,
        route_cost=lambda s: Price(105.0, "route") if s == "O=S=O" else None,
    )
    so2 = next(line for line in eco.raw if line.smiles == "O=S=O")
    assert so2.price.basis == "route" and so2.price.usd_per_t == 105.0


def test_glycol_sulfate_is_not_priced_as_polyglycol():
    """Сульфат диэтиленгликоля — не полигликоль: класс только для C, H, O."""
    assert BOOK.market("O=S(=O)(O)OCCOCCO") is None
    assert BOOK.market("CC(O)COC(C)CO").basis == "class"  # ДПГ
