#!/usr/bin/env python3
"""
balancer.py — полное уравнение реакции с коэффициентами (roadmap R4, решение 0015).

Генератор выдаёт «реагенты > сопутствующий реагент > основной продукт»:
вода, соли и коэффициенты не выписаны. Экономике нужен расход сырья на
тонну продукта, поэтому здесь уравнение достраивается:

  - каждое вещество — вектор «атомы каждого элемента + заряд»;
  - неизвестные — коэффициенты реагентов, основного продукта и побочных
    продуктов из пула (вода, соли, HX, CO2, H2, O2, NH3, N2, MeOH, EtOH, ионы);
  - наборы побочных продуктов перебираются от меньшего к большему (порядок
    пула — порядок предпочтения); принимается первый, при котором решение
    единственное и все коэффициенты положительны; арифметика точная (дроби);
  - сначала без сопутствующего реагента (катализатор не расходуется), потом
    с ним (нейтрализация: NaOH расходуется);
  - продукт-анион с Na+ среди побочных записывается солью натрия;
  - не сошлось — статус failed и остаток (что не удалось распределить).

Пример:
  python -m chem_agent.balancer --in outputs/<прогон>/network.csv --out outputs/<прогон>/network_balanced.csv
"""

from __future__ import annotations

import argparse
import collections
import csv
import itertools
import multiprocessing
import os
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from math import gcd

from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

# Побочные продукты в порядке предпочтения.
BYPRODUCT_POOL = [
    "O",  # вода
    "[Na+]",
    "[Cl-].[Na+]",
    "[Br-].[Na+]",
    "[I-].[Na+]",
    "Cl",
    "Br",
    "I",
    "O=C=O",
    "[H][H]",
    "O=O",
    "N",
    "N#N",
    "CO",
    "CCO",
    # типовые уходящие группы (снятие защит, гидролиз, кросс-сочетания)
    "CC(=O)[O-].[Na+]",  # ацетат натрия
    "CC(=O)O",
    "OB(O)O",  # борная кислота (Сузуки, Чан–Лэм)
    "C=C(C)C",  # изобутилен (Boc, трет-бутиловые эфиры)
    "Cc1ccccc1",  # толуол (бензильная защита)
    "O=C(O)c1ccccc1",  # бензойная кислота (бензоильная защита)
    "S",  # сероводород
    "[H+]",
    "[OH-]",
    "[K+]",
    "[Cl-]",
    "[Br-]",
    "[I-]",
]
MAX_BYPRODUCTS = 3
NAOH = "[Na+].[OH-]"


def canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else smiles


@lru_cache(maxsize=100_000)
def composition(smiles: str) -> tuple[tuple[str, int], ...]:
    """Атомы каждого элемента (с водородами) и заряд: (('C', 2), ..., ('+', 0))."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"не разбирается SMILES: {smiles}")
    c = collections.Counter()
    charge = 0
    for a in mol.GetAtoms():
        c[a.GetSymbol()] += 1
        c["H"] += a.GetTotalNumHs()
        charge += a.GetFormalCharge()
    c["+"] = charge
    return tuple(sorted(c.items()))


def formula(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return CalcMolFormula(mol) if mol else smiles


def _elements(smiles: str) -> set[str]:
    return {k for k, v in composition(smiles) if k != "+" and v}


def _nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """Базис ядра матрицы (точная арифметика, метод Гаусса)."""
    m = [row[:] for row in matrix]
    n_cols = len(m[0]) if m else 0
    pivots, r = [], 0
    for c in range(n_cols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        piv = m[r][c]
        m[r] = [x / piv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r], strict=True)]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    free = [c for c in range(n_cols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * n_cols
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -m[i][fc]
        basis.append(v)
    return basis


def _solve(left: list[str], right: list[str]) -> list[int] | None:
    """Единственное решение с положительными целыми коэффициентами или None."""
    species = left + right
    keys = sorted({k for s in species for k, _ in composition(s)})
    comp = [dict(composition(s)) for s in species]
    sign = [-1] * len(left) + [1] * len(right)
    matrix = [[Fraction(sign[j] * comp[j].get(k, 0)) for j in range(len(species))] for k in keys]
    # Стехиометрия шаблона: молекула продукта строится из одной молекулы каждого
    # органического реагента — их коэффициенты равны коэффициенту продукта.
    # Неорганические (O2, S8, вода, NaOH) свободны: S8 -> 8 SO2, 2 спирта + O2.
    # Без этого «сходятся» бессмыслицы: 4 метилбензоата -> 4 кислоты + изобутилен.
    # Органический побочный продукт (MeOH, ацетат, изобутилен, CO2) — уходящая
    # группа: одна молекула на молекулу продукта.
    n = len(left)
    if comp[n].get("C"):
        for j in [*range(n), *range(n + 1, len(species))]:
            if comp[j].get("C"):
                matrix.append(
                    [Fraction(1 if i == j else -1 if i == n else 0) for i in range(len(species))]
                )
    basis = _nullspace(matrix)
    if len(basis) > 1:
        # неоднозначно (реагенты-изомеры или кратные составы): стехиометрия
        # шаблона — по одной молекуле каждого реагента
        for j in range(1, len(left)):
            matrix.append(
                [Fraction(1 if i == 0 else -1 if i == j else 0) for i in range(len(species))]
            )
        basis = _nullspace(matrix)
    if len(basis) != 1:
        return None
    v = basis[0]
    if v[len(left)] < 0:  # основной продукт — положительный
        v = [-x for x in v]
    if any(x <= 0 for x in v):
        return None
    lcm = 1
    for x in v:
        lcm = lcm * x.denominator // gcd(lcm, x.denominator)
    ints = [int(x * lcm) for x in v]
    g = 0
    for x in ints:
        g = gcd(g, x)
    return [x // g for x in ints]


@dataclass
class Balanced:
    status: str  # ok | reagent_consumed | failed
    reason: str = ""
    reactants: list[tuple[str, int]] = field(default_factory=list)
    products: list[tuple[str, int]] = field(default_factory=list)

    @property
    def byproducts(self) -> list[str]:
        return [s for s, _ in self.products[1:]]

    def _side(self, items, fmt) -> str:
        return " + ".join((f"{c} " if c != 1 else "") + fmt(s) for s, c in items)

    @property
    def equation(self) -> str:
        if self.status == "failed":
            return ""
        return f"{self._side(self.reactants, str)} -> {self._side(self.products, str)}"

    @property
    def equation_formula(self) -> str:
        if self.status == "failed":
            return ""
        return f"{self._side(self.reactants, formula)} -> {self._side(self.products, formula)}"


def _as_salt(product: str, byproducts: list[tuple[str, int]], p_coef: int):
    """Анион-продукт + Na+ среди побочных -> соль натрия."""
    charge = dict(composition(product))["+"]
    na = dict(byproducts).get("[Na+]")
    if charge < 0 and na == p_coef * -charge:
        salt = canon(product + ".[Na+]" * -charge)
        return salt, [(s, c) for s, c in byproducts if s != "[Na+]"]
    return product, byproducts


def sodium_salt(smiles: str) -> str:
    """Анион без противоиона -> натриевая соль (движок хранит HSO3- без Na+)."""
    charge = dict(composition(smiles))["+"]
    return canon(smiles + ".[Na+]" * -charge) if charge < 0 else smiles


FREE_IONS = {"[Na+]", "[H+]", "[OH-]", "[K+]", "[Cl-]", "[Br-]", "[I-]"}
OXIDANTS = {"O=O", "OO", "O=[O+][O-]"}


def _is_clean(b: Balanced, reagent: str) -> bool:
    """Нет свободных ионов; при окислителе-реагенте не выделяется H2."""
    if any(s in FREE_IONS for s in b.byproducts):
        return False
    return not (reagent in OXIDANTS and "[H][H]" in b.byproducts)


def _search(left: list[str], product: str, status: str) -> Balanced | None:
    """Наименьший набор побочных продуктов из пула, дающий единственное решение."""
    # Побочные продукты — только из элементов, которых не хватает в балансе
    # «реагенты по одному минус продукт» (для окисления это H и O: вода, H2,
    # O2), плюс ионы для заряда. Иначе перебор по всему пулу слишком долгий.
    available = set().union(*(_elements(s) for s in left))
    rest = collections.Counter()
    for s in left:
        rest.update(dict(composition(s)))
    rest.subtract(dict(composition(product)))
    missing = {k for k, v in rest.items() if k != "+" and v} or available

    def useful(s: str) -> bool:
        el = _elements(s)
        return el <= available and (bool(el & missing) or s in FREE_IONS)

    pool = [s for s in BYPRODUCT_POOL if s not in left and s != product and useful(s)]
    # Быстрый отказ: элемента остатка нет ни в одном побочном продукте пула, а
    # содержащие его реагенты органические (коэффициент = продукту) — не сойдётся.
    if dict(composition(product)).get("C"):
        pool_el = set().union(*(_elements(s) for s in pool)) if pool else set()
        for e in missing - pool_el:
            carriers = [s for s in left if dict(composition(s)).get(e)]
            if carriers and all(dict(composition(s)).get("C") for s in carriers):
                return None
    for k in range(MAX_BYPRODUCTS + 1):
        for subset in itertools.combinations(pool, k):
            sol = _solve(left, [product, *subset])
            if sol is None:
                continue
            n = len(left)
            prod, byp = _as_salt(product, list(zip(subset, sol[n + 1 :], strict=True)), sol[n])
            return Balanced(
                status,
                reactants=list(zip(left, sol[:n], strict=True)),
                products=[(prod, sol[n]), *byp],
            )
    return None


def balance(reactants: list[str], reagent: str, product: str) -> Balanced:
    """
    Уравнение для «реагенты (+ сопутствующий реагент) -> основной продукт».
    Попытки по порядку, берётся первое «чистое» решение (без свободных ионов):
      1. анионы как натриевые соли, без реагента (катализатор не расходуется);
      2. то же с реагентом (NaOH при нейтрализации; окислитель-реагент
         расходуется, если без него выделялся бы H2);
      3. NaOH среди реагентов как катализатор, реагирует вода (OH- на месте воды);
      4. то же без солевой записи (свободные ионы) — запасной вариант.
    """
    reactants = [canon(s) for s in reactants]
    product = canon(product)
    reagent = canon(reagent) if reagent else ""
    salts = [sodium_salt(s) for s in reactants]
    salt_product = sodium_salt(product)
    attempts = [(salts, salt_product, "ok")]
    if reagent and reagent not in reactants:
        with_reagent = (salts + [sodium_salt(reagent)], salt_product, "reagent_consumed")
        attempts.append(with_reagent)
    if NAOH in salts:
        water = [s for s in salts if s != NAOH] + ["O"]
        attempts.append((water, salt_product, "naoh_catalyst"))
    if "O" not in salts:
        # гидролиз: реагирует вода (всегда доступный ресурс, решение 0011)
        attempts.append((salts + ["O"], salt_product, "hydrolysis"))
        attempts.append((salts + ["O"], product, "hydrolysis"))
    attempts += [(reactants, product, st) for _, _, st in attempts[:2]]
    if product != salt_product:
        attempts.insert(1, (salts, product, "ok"))
    fallback = None
    for left, prod, status in attempts:
        found = _search(left, prod, status)
        if found and _is_clean(found, reagent):
            if status == "naoh_catalyst":
                found.reason = "OH- — катализатор, реагирует вода"
            if "[H][H]" in found.byproducts:
                found.reason = "выделяется H2: вероятно, не указан окислитель"
            return found
        fallback = fallback or found
    if fallback:
        fallback.reason = "в уравнении свободные ионы: противоион не определён"
        return fallback
    # остаток: реагенты по одному минус продукт
    rest = collections.Counter()
    for s in reactants:
        rest.update(dict(composition(s)))
    rest.subtract(dict(composition(product)))
    parts = [f"{k}{v if v != 1 else ''}" for k, v in sorted(rest.items()) if k != "+" and v > 0]
    hill = sorted(parts, key=lambda p: (p[0] != "C", p[0] != "H", p))
    return Balanced(
        "failed",
        reason="нет решения с побочными продуктами из пула; остаток "
        + ("".join(hill) or "—")
        + (f", заряд {rest['+']:+d}" if rest["+"] else ""),
    )


def _reaction_parts(rxn: str) -> tuple[list[str], str, str]:
    """Реагенты (Na+ и OH- собираются в NaOH), сопутствующий реагент, продукт."""
    reac, reag, prod = rxn.split(">")
    frags = [f for f in reac.split(".") if f]
    if "[Na+]" in frags and "[OH-]" in frags:
        frags.remove("[Na+]")
        frags.remove("[OH-]")
        frags.append(NAOH)
    return frags, reag, prod


def _balance_key(rxn: str) -> Balanced:
    return balance(*_reaction_parts(rxn))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--in", dest="inp", required=True, help="candidates.csv / network.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=0, help="процессов; 0 — все ядра")
    args = ap.parse_args(argv)
    with open(args.inp, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    keys = list(dict.fromkeys(r["reaction_smiles"] for r in rows))
    jobs = args.jobs if args.jobs > 0 else (os.cpu_count() or 1)
    if jobs > 1 and len(keys) > 200 and "fork" in multiprocessing.get_all_start_methods():
        with multiprocessing.get_context("fork").Pool(jobs) as pool:
            results = pool.map(_balance_key, keys, chunksize=max(1, len(keys) // (jobs * 16)))
    else:
        results = [_balance_key(k) for k in keys]
    cache = dict(zip(keys, results, strict=True))
    stats = collections.Counter()
    for r in rows:
        b = cache[r["reaction_smiles"]]
        r["balance_status"] = b.status
        r["equation"] = b.equation
        r["equation_formula"] = b.equation_formula
        r["byproducts"] = ".".join(b.byproducts)
        r["balance_reason"] = b.reason
        stats[b.status] += 1
    fields = list(rows[0]) if rows else []
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    total = sum(stats.values()) or 1
    print(
        "[balance] "
        + ", ".join(f"{k}: {v} ({v / total:.0%})" for k, v in stats.most_common())
        + f" -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
