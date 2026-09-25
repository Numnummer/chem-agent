#!/usr/bin/env python3
"""
corpus_prep.py — приведение выгрузки all_balanced_reactions к формату движка.

Вход: CSV из data/raw/ (колонки record_id, source_subset, reaction_smiles,
reaction_smiles_balanced, patent, условия, коэффициенты ...).
Выход в --out-dir:

  corpus.csv       id, rxn_smiles (реагенты>агенты>продукты), source — для шага map
  meta.csv         по id: патенты, номера записей, условия, уравненная реакция
  prep_report.md   сколько записей отброшено и почему, сколько схлопнуто

Что делается (docs/decisions/0008):
  - исходная запись реакции из патента (reaction_smiles), а не уравненная:
    в уравненной растворители и катализаторы стоят с обеих сторон;
  - формат 'A>P' приводится к 'A>агенты>P', агенты — из agents_smiles;
  - одиночный атом серы [S] (элементарная сера) -> S1SSSSSSS1; сероводород
    'S' и сульфиды не трогаются;
  - дубликаты схлопываются по канонической записи реакции, все номера
    записей и патенты сохраняются в meta.csv;
  - запись без патента получает source = «без источника».

Пример:
  python -m chem_agent.corpus_prep --raw data/raw/all_balanced_reactions_training.csv \\
                                   --out-dir outputs/2026-09-25-full-v1
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

from rdkit import Chem, RDLogger

from chem_agent.template_engine import NO_SOURCE

RDLogger.DisableLog("rdApp.*")

ELEMENTAL_SULFUR_ATOM = "[S]"  # так в патентных выгрузках пишут серу
ELEMENTAL_SULFUR = "S1SSSSSSS1"  # соглашение проекта (chem-conventions)

_canon_cache: dict[str, str | None] = {}


def canon(smiles: str) -> str | None:
    if smiles not in _canon_cache:
        mol = Chem.MolFromSmiles(smiles)
        _canon_cache[smiles] = Chem.MolToSmiles(mol) if mol else None
    return _canon_cache[smiles]


def split_sides(rxn: str) -> tuple[list[str], list[str]] | None:
    """'A.B>>P', 'A.B>X>P' и 'A.B>P' -> ([A, B], [P]); агенты отбрасываются."""
    parts = rxn.strip().split(">")
    if len(parts) not in (2, 3):
        return None
    reac, prod = parts[0], parts[-1]
    reac_f = [f for f in reac.split(".") if f]
    prod_f = [f for f in prod.split(".") if f]
    if not reac_f or not prod_f:
        return None
    return reac_f, prod_f


def normalize_fragments(frags: list[str]) -> tuple[list[str], int] | None:
    """Канонизирует фрагменты и заменяет [S] на S8. None, если что-то не парсится."""
    out, n_sulfur = [], 0
    for f in frags:
        c = canon(f)
        if c is None:
            return None
        if c == ELEMENTAL_SULFUR_ATOM:
            c, n_sulfur = ELEMENTAL_SULFUR, n_sulfur + 1
        out.append(c)
    return out, n_sulfur


def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def balanced_consistent(products: list[str], balanced_rxn: str) -> bool:
    """Есть ли основной продукт исходной записи среди продуктов уравненной."""
    sides = split_sides(balanced_rxn) if balanced_rxn else None
    if sides is None:
        return False
    norm = normalize_fragments(sides[1])
    if norm is None:
        return False
    return max(products, key=heavy_atoms) in set(norm[0])


def subset_selected(subset: str, include: list[str], exclude: list[str]) -> bool:
    if include and not any(subset.startswith(p) for p in include):
        return False
    return not any(subset.startswith(p) for p in exclude)


META_FIELDS = [
    "id",
    "source",
    "has_source",
    "n_records",
    "record_ids",
    "patents",
    "subsets",
    "year",
    "sulfur_normalized",
    "temperature_c_min",
    "temperature_c_max",
    "pressure",
    "solvent",
    "atmosphere",
    "time_minutes_stages",
    "mode",
    "ph",
    "rxn_balanced",
    "reactant_coefficients",
    "product_coefficients",
    "synrbl_method",
    "synrbl_confidence",
    "balanced_consistent",
]


def prepare(rows, include=(), exclude=()):
    """Возвращает (corpus_rows, meta_rows, stats)."""
    stats = collections.Counter()
    groups: dict[str, dict] = {}  # каноническая реакция -> запись

    for row in rows:
        stats["records"] += 1
        subset = row.get("source_subset", "")
        if not subset_selected(subset, list(include), list(exclude)):
            stats["skip_subset"] += 1
            continue
        rxn = row.get("reaction_smiles", "")
        sides = split_sides(rxn)
        if sides is None:
            stats["skip_bad_format"] += 1
            continue
        stats["format_single_arrow"] += rxn.count(">") == 1
        reac = normalize_fragments(sides[0])
        prod = normalize_fragments(sides[1])
        if reac is None or prod is None:
            stats["skip_unparsable"] += 1
            continue
        (reac_f, s_r), (prod_f, s_p) = reac, prod
        agents = [c for a in row.get("agents_smiles", "").split(".") if a and (c := canon(a))]
        key = ".".join(sorted(reac_f)) + ">>" + ".".join(sorted(prod_f))
        patent = row.get("patent", "").strip()

        g = groups.get(key)
        if g is None:
            g = groups[key] = {
                "row": row,
                "reactants": reac_f,
                "agents": agents,
                "products": prod_f,
                "sulfur": s_r + s_p,
                "record_ids": [],
                "patents": [],
                "subsets": [],
            }
        elif patent and not g["row"].get("patent", "").strip():
            # представитель группы — запись с патентом: у неё есть источник
            g.update(row=row, agents=agents)
        g["record_ids"].append(row.get("record_id", ""))
        if patent and patent not in g["patents"]:
            g["patents"].append(patent)
        if subset not in g["subsets"]:
            g["subsets"].append(subset)

    corpus_rows, meta_rows = [], []
    for g in groups.values():
        r = g["row"]
        rid = "AB-" + r.get("record_id", "")
        source = g["patents"][0] if g["patents"] else NO_SOURCE
        corpus_rows.append(
            {
                "id": rid,
                "rxn_smiles": ".".join(g["reactants"])
                + ">"
                + ".".join(g["agents"])
                + ">"
                + ".".join(g["products"]),
                "source": source,
            }
        )
        consistent = balanced_consistent(g["products"], r.get("reaction_smiles_balanced", ""))
        meta_rows.append(
            {
                "id": rid,
                "source": source,
                "has_source": int(bool(g["patents"])),
                "n_records": len(g["record_ids"]),
                "record_ids": "|".join(g["record_ids"]),
                "patents": "|".join(g["patents"]),
                "subsets": "|".join(g["subsets"]),
                "year": r.get("year", ""),
                "sulfur_normalized": int(g["sulfur"] > 0),
                "temperature_c_min": r.get("condition_temperature_c_min", ""),
                "temperature_c_max": r.get("condition_temperature_c_max", ""),
                "pressure": r.get("condition_pressure", ""),
                "solvent": r.get("condition_solvent", ""),
                "atmosphere": r.get("condition_atmosphere", ""),
                "time_minutes_stages": r.get("condition_time_minutes_stages", ""),
                "mode": r.get("condition_mode", ""),
                "ph": r.get("condition_ph", ""),
                "rxn_balanced": r.get("reaction_smiles_balanced", ""),
                "reactant_coefficients": r.get("reactant_stoichiometry_coefficients", ""),
                "product_coefficients": r.get("product_stoichiometry_coefficients", ""),
                "synrbl_method": r.get("synrbl_method", ""),
                "synrbl_confidence": r.get("synrbl_confidence", ""),
                "balanced_consistent": int(consistent),
            }
        )
        stats["unique"] += 1
        stats["unique_with_source"] += bool(g["patents"])
        stats["unique_sulfur_normalized"] += g["sulfur"] > 0
        stats["unique_balanced_inconsistent"] += not consistent
    stats["duplicates_merged"] = (
        stats["records"]
        - stats["skip_subset"]
        - stats["skip_bad_format"]
        - stats["skip_unparsable"]
        - stats["unique"]
    )
    return corpus_rows, meta_rows, stats


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def report(stats, meta_rows, args) -> str:
    by_subset = collections.Counter()
    for m in meta_rows:
        for s in m["subsets"].split("|"):
            by_subset[s] += 1
    L = [
        "# Подготовка корпуса\n",
        f"Команда: `python -m chem_agent.corpus_prep --raw {args.raw} --out-dir {args.out_dir}"
        + "".join(f" --include-subset {s}" for s in args.include_subset)
        + "".join(f" --exclude-subset {s}" for s in args.exclude_subset)
        + "`\n",
        "| Шаг | Записей |",
        "|---|---|",
        f"| Всего записей | {stats['records']} |",
        f"| Отброшено: выборка не выбрана | {stats['skip_subset']} |",
        f"| Отброшено: не разбирается формат | {stats['skip_bad_format']} |",
        f"| Отброшено: SMILES не парсится | {stats['skip_unparsable']} |",
        f"| Схлопнуто дубликатов | {stats['duplicates_merged']} |",
        f"| **Уникальных реакций** | **{stats['unique']}** |",
        f"| — с патентом | {stats['unique_with_source']} |",
        f"| — {NO_SOURCE} | {stats['unique'] - stats['unique_with_source']} |",
        f"| — сера [S] заменена на S8 | {stats['unique_sulfur_normalized']} |",
        f"| — уравненная колонка не совпадает с исходной | "
        f"{stats['unique_balanced_inconsistent']} |",
        "",
        f"Записей в формате `A>P` (одна стрелка): {stats['format_single_arrow']}.\n",
        "## Уникальные реакции по выборкам\n",
        "Реакция, попавшая в несколько выборок, считается в каждой.\n",
        "| Выборка | Реакций |",
        "|---|---|",
    ]
    L += [f"| {s} | {n} |" for s, n in by_subset.most_common()]
    return "\n".join(L) + "\n"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--raw", required=True, help="CSV выгрузки all_balanced_reactions")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--include-subset",
        nargs="*",
        default=[],
        help="брать только выборки с этими префиксами source_subset",
    )
    ap.add_argument(
        "--exclude-subset",
        nargs="*",
        default=[],
        help="исключить выборки с этими префиксами (например, NAOH 2naoh)",
    )
    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    csv.field_size_limit(sys.maxsize)
    with open(args.raw, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    corpus_rows, meta_rows, stats = prepare(rows, args.include_subset, args.exclude_subset)
    os.makedirs(args.out_dir, exist_ok=True)
    write_csv(os.path.join(args.out_dir, "corpus.csv"), corpus_rows, ["id", "rxn_smiles", "source"])
    write_csv(os.path.join(args.out_dir, "meta.csv"), meta_rows, META_FIELDS)
    text = report(stats, meta_rows, args)
    with open(os.path.join(args.out_dir, "prep_report.md"), "w", encoding="utf-8") as f:
        f.write(text)
    print(
        f"[prep] записей: {stats['records']}, уникальных реакций: {stats['unique']} "
        f"(с патентом: {stats['unique_with_source']}) -> {args.out_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
