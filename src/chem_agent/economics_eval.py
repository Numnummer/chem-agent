"""
economics_eval.py — точность экономической оценки (ТЗ 2.3, решение 0016).

Эталон (data/economics/reference.csv): расход сырья на 1 т продукта по
суммарному уравнению всего маршрута от исходного набора, посчитанный
независимо от пайплайна (стандартные атомные массы). Затраты эталона —
нормы × цены из того же прайса.

Прогноз: анализ набора (выход 100%, без кредита за побочные), затраты на
сырьё лучшего маршрута к продукту — проверяется вся цепочка: генерация
маршрута, балансировка стадий, пересчёт себестоимости промежуточных.

Метрика: точность = среднее по продуктам (1 − |прогноз − эталон| / эталон);
порог ТЗ — 85%. Бейзлайн — константный предиктор (средние затраты эталона).
Когда заказчик даст свои расчёты, они заменяют теоретические нормы.

  python -m chem_agent.economics_eval --templates outputs/<прогон>/templates.jsonl --out outputs/<дата>-eval-economics
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

from chem_agent.economics import PriceBook, load_config
from chem_agent.pipeline import ROOT, D, build_parser, run
from chem_agent.template_engine import canon

THRESHOLD = 0.85


def evaluate_costs(results: list[dict], reference_csv: str, book: PriceBook) -> dict:
    ref = collections.defaultdict(list)
    names = {}
    with open(reference_csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            ref[canon(r["product_smiles"])].append((r["raw_smiles"], float(r["reference_t_per_t"])))
            names[canon(r["product_smiles"])] = r["product_name"]
    pred = {}
    for r in results:
        e = r.get("economics_usd_per_t")
        p = canon(r.get("product_balanced") or "")
        if e and e["raw_cost"] is not None and p in ref:
            pred[p] = min(pred.get(p, float("inf")), e["raw_cost"])
    ref_cost = {
        p: sum(t * book.exact[canon(s)].usd_per_t for s, t in raws) for p, raws in ref.items()
    }
    const = sum(ref_cost.values()) / len(ref_cost)
    rows = []
    for p, rc in ref_cost.items():
        pc = pred.get(p)
        acc = None if pc is None else max(0.0, 1 - abs(pc - rc) / rc)
        rows.append(
            {
                "product": names[p],
                "reference_usd_t": round(rc, 1),
                "predicted_usd_t": None if pc is None else round(pc, 1),
                "accuracy": None if acc is None else round(acc, 4),
                "baseline_accuracy": round(max(0.0, 1 - abs(const - rc) / rc), 4),
            }
        )
    found = [r for r in rows if r["accuracy"] is not None]
    # продукт, к которому система не построила маршрут, — точность 0
    accuracy = sum(r["accuracy"] or 0 for r in rows) / len(rows)
    baseline = sum(r["baseline_accuracy"] for r in rows) / len(rows)
    return {
        "accuracy": round(accuracy, 4),
        "baseline_accuracy": round(baseline, 4),
        "threshold": THRESHOLD,
        "passed": accuracy >= THRESHOLD,
        "products": len(rows),
        "routes_found": len(found),
        "rows": rows,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--templates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--inputs", default=os.path.join(ROOT, "data/inputs/priority_set.csv"))
    ap.add_argument("--reference", default=os.path.join(ROOT, "data/economics/reference.csv"))
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=0)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    cfg = load_config(os.path.join(ROOT, D["config"]))
    cfg.update(default_yield=1.0, byproduct_credit=0.0)
    cfg_path = os.path.join(args.out, "config_eval.json")
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False)
    run(
        build_parser().parse_args(
            ["--inputs", args.inputs, "--templates", args.templates, "--out", args.out]
            + ["--depth", str(args.depth), "--jobs", str(args.jobs), "--config", cfg_path]
        )
    )
    with open(os.path.join(args.out, "results.json"), encoding="utf-8") as f:
        results = json.load(f)["reactions"]
    book = PriceBook(os.path.join(ROOT, D["prices"]), os.path.join(ROOT, D["classes"]))
    res = evaluate_costs(results, args.reference, book)
    with open(os.path.join(args.out, "economics_eval.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    lines = [
        "# Точность экономической оценки",
        "",
        f"Точность: **{res['accuracy']:.1%}** (порог {THRESHOLD:.0%}); бейзлайн (константа): "
        f"{res['baseline_accuracy']:.1%}; маршрут найден для {res['routes_found']} из {res['products']}.",
        "",
        "| Продукт | Эталон, $/т | Прогноз, $/т | Точность | Бейзлайн |",
        "|---|---|---|---|---|",
    ]
    for r in res["rows"]:
        lines.append(
            f"| {r['product']} | {r['reference_usd_t']} | {r['predicted_usd_t'] or '—'} | "
            + ("—" if r["accuracy"] is None else f"{r['accuracy']:.1%}")
            + f" | {r['baseline_accuracy']:.1%} |"
        )
    with open(os.path.join(args.out, "economics_eval.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines), file=sys.stderr)


if __name__ == "__main__":
    main()
