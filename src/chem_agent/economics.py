"""
economics.py — экономическая оценка реакции на 1 т основного продукта
(ТЗ 2.3, решение 0016).

По уравнению балансировщика (0015):
  - расход сырья, т/т продукта = коэф. × M(сырьё) / (коэф. × M(продукт)) / выход;
  - цена сырья: меньшая из рыночной (прайс) и себестоимости маршрута, если
    вещество получается в сети из набора (купить или сделать: SO2 из серы
    дешевле рынка); иначе цена класса продукта; иначе неизвестна;
  - выручка — цена продукта (точная или по классу);
  - кредит за побочные продукты с ценой (Na2SO4, метанол...) × доля;
  - маржа = выручка + кредит − сырьё.

Цены — data/prices/prices.csv и classes.csv; параметры — data/config/mvp.json.
Энергозатраты в MVP не оцениваются (нет данных по процессам).
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

from chem_agent.template_engine import canon


@dataclass
class Price:
    usd_per_t: float
    basis: str  # price | class | route
    note: str = ""


@dataclass
class Line:
    smiles: str
    t_per_t: float
    price: Price | None

    @property
    def usd(self) -> float | None:
        return None if self.price is None else self.t_per_t * self.price.usd_per_t


@dataclass
class Economics:
    raw: list[Line] = field(default_factory=list)
    byproducts: list[Line] = field(default_factory=list)
    product_price: Price | None = None
    yield_: float = 1.0
    credit_share: float = 0.0

    @property
    def raw_cost(self) -> float | None:
        costs = [line.usd for line in self.raw]
        return None if any(c is None for c in costs) else sum(costs)

    @property
    def credit(self) -> float:
        return self.credit_share * sum(line.usd or 0 for line in self.byproducts)

    @property
    def revenue(self) -> float | None:
        return None if self.product_price is None else self.product_price.usd_per_t

    @property
    def margin(self) -> float | None:
        if self.raw_cost is None or self.revenue is None:
            return None
        return self.revenue + self.credit - self.raw_cost

    @property
    def unknown(self) -> list[str]:
        out = [line.smiles for line in self.raw if line.price is None]
        if self.product_price is None:
            out.append("продукт")
        return out


def mol_weight(smiles: str) -> float:
    return MolWt(Chem.MolFromSmiles(smiles))


class PriceBook:
    def __init__(self, prices_csv: str, classes_csv: str | None = None):
        self.exact: dict[str, Price] = {}
        self.names: dict[str, str] = {}
        with open(prices_csv, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                c = canon(r["smiles"]) or r["smiles"]
                self.exact[c] = Price(float(r["price_usd_per_t"]), "price", r.get("source", ""))
                self.names[c] = r["name"]
        self.classes: list[tuple[str, str, object, Price, set]] = []
        if classes_csv:
            with open(classes_csv, encoding="utf-8") as f:
                for r in csv.DictReader(f):
                    patt = Chem.MolFromSmarts(r["smarts"])
                    self.classes.append(
                        (
                            r["class"],
                            r["name"],
                            patt,
                            Price(float(r["price_usd_per_t"]), "class", r["name"]),
                            {e for e in r.get("only_elements", "").split(",") if e},
                        )
                    )

    def market(self, smiles: str) -> Price | None:
        """Рыночная цена: точное вещество, иначе класс продукта."""
        c = canon(smiles) or smiles
        if c in self.exact:
            return self.exact[c]
        mol = Chem.MolFromSmiles(c)
        if mol is None:
            return None
        elements = {a.GetSymbol() for a in mol.GetAtoms()}
        for _, _, patt, price, only in self.classes:
            # класс «полигликоли» — только C, H, O: сульфат гликоля — уже не полигликоль
            if only and not elements <= only:
                continue
            if mol.HasSubstructMatch(patt):
                return price
        return None


def evaluate(
    reactants: list[tuple[str, int]],
    products: list[tuple[str, int]],
    book: PriceBook,
    route_cost=None,
    yield_: float = 1.0,
    credit_share: float = 0.0,
) -> Economics:
    """
    Экономика уравнения на 1 т основного продукта (products[0]).
    route_cost(smiles) -> Price | None — себестоимость вещества, получаемого в
    сети; используется для сырья без цены в прайсе.
    """
    main, p_coef = products[0]
    p_mass = p_coef * mol_weight(main)
    eco = Economics(yield_=yield_, credit_share=credit_share)
    for s, c in reactants:
        t = c * mol_weight(s) / p_mass / yield_
        # купить или сделать: из рыночной цены и себестоимости маршрута в сети
        # берётся меньшая (SO2 по рынку 350 $/т, из серы — ~105 $/т)
        options = [
            p for p in (book.exact.get(canon(s) or s), route_cost(s) if route_cost else None) if p
        ]
        price = min(options, key=lambda p: p.usd_per_t) if options else book.market(s)
        eco.raw.append(Line(s, t, price))
    for s, c in products[1:]:
        price = book.exact.get(canon(s) or s)
        if price is not None and price.usd_per_t > 5:  # воду не кредитуем
            eco.byproducts.append(Line(s, c * mol_weight(s) / p_mass, price))
    eco.product_price = book.market(main)
    return eco


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)
