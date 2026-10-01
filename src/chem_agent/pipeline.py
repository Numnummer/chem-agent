#!/usr/bin/env python3
"""
pipeline.py — MVP: набор веществ -> отранжированные реакции с уравнениями,
экономикой и доказательной базой (ТЗ 3.1 в составе этапа 2; решение 0016).

  1. генерация: шаблоны внутри набора (+ вода и воздух), до N стадий;
  2. уравнения: балансировщик (0015);
  3. доказательства: уровень I/II, прецедент, источник, условия (0007);
  4. экономика на 1 т продукта: сырьё по маршруту от исходного набора,
     выручка, маржа (2.3);
  5. векторная БД: уникальность векторов, похожие известные процессы (2.2);
  6. ранжирование: маржа, доказательность, длина маршрута, качество
     уравнения, «внутри набора» — веса в data/config/mvp.json;
  7. отчёт: results.json, results.csv, summary.md, report.html.

Пример:
  python -m chem_agent.pipeline --inputs data/inputs/priority_set.csv \\
      --templates outputs/<прогон>/templates.jsonl --out outputs/<дата>-mvp
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import sys
import time

from chem_agent.balancer import balance, formula, sodium_salt
from chem_agent.economics import Price, PriceBook, evaluate, load_config
from chem_agent.evidence import EvidenceIndex
from chem_agent.template_engine import NO_SOURCE, canon
from chem_agent.template_engine import main as engine_main
from chem_agent.vectordb import VectorDB

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
D = {
    "utilities": "data/inputs/utilities.csv",
    "manual": "data/manual/inorganic.csv",
    "prices": "data/prices/prices.csv",
    "classes": "data/prices/classes.csv",
    "config": "data/config/mvp.json",
}

BALANCE_SCORE = {"ok": 1.0, "reagent_consumed": 1.0, "hydrolysis": 1.0, "naoh_catalyst": 0.8}


def _read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


class Names:
    """Человекочитаемые имена: набор пользователя, ресурсы, прайс; иначе формула."""

    def __init__(self, *tables: dict[str, str]):
        self.map: dict[str, str] = {}
        for t in tables:
            for s, n in t.items():
                self.map.setdefault(canon(s) or s, n)

    def __call__(self, smiles: str) -> str:
        c = canon(smiles) or smiles
        # движок хранит анионы без противоиона, прайс — натриевые соли
        return self.map.get(c) or self.map.get(sodium_salt(c)) or formula(smiles)


def _parts(rxn: str):
    reac, reag, prod = rxn.split(">")
    frags = [f for f in reac.split(".") if f]
    if "[Na+]" in frags and "[OH-]" in frags:
        frags.remove("[Na+]")
        frags.remove("[OH-]")
        frags.append("[Na+].[OH-]")
    return frags, reag, prod


def _evidence_score(r) -> float:
    if r["evidence_level"] == "I":
        return 1.0
    n = int(r["template_count"] or 1)
    has_source = r["precedent_source"] not in ("", NO_SOURCE)
    return 0.4 + 0.4 * min(1.0, math.log10(n + 1) / 2) + 0.2 * has_source


def _margin_score(r) -> float:
    if r["margin_usd_t"] is None or not r["revenue_usd_t"]:
        return 0.0
    pct = max(-1.0, min(1.0, r["margin_usd_t"] / r["revenue_usd_t"]))
    return (pct + 1) / 2


def run(args) -> dict:
    t0 = time.time()
    os.makedirs(args.out, exist_ok=True)
    cfg = load_config(args.config)
    book = PriceBook(args.prices, args.classes)
    inputs = _read(args.inputs)
    utilities = _read(args.utilities) if args.utilities else []
    names = Names(
        {r["smiles"]: r["name"] for r in inputs if r["smiles"] != "resource"},
        {r["smiles"]: r["name"] for r in utilities},
        book.names,
        {"[Na+].[OH-]": "Гидроксид натрия"},
    )

    # 1. генерация
    network = os.path.join(args.out, "network.csv")
    engine_main(
        ["apply", "--templates", args.templates, "--inputs", args.inputs, "--internal-only"]
        + (["--utilities", args.utilities] if args.utilities else [])
        + ["--depth", str(args.depth), "--jobs", str(args.jobs)]
        + ["--out", network, "--report", os.path.join(args.out, "network.md")]
    )
    rows = _read(network)
    by_core: dict[str, dict] = {}
    for r in rows:
        rec = by_core.setdefault(r["reaction_core"], {**r, "reagents": set()})
        if r["reagent"]:
            rec["reagents"].add(r["reagent"])
    results = list(by_core.values())

    # 2. уравнения
    for r in results:
        reac, reag, prod = _parts(r["reaction_smiles"])
        b = balance(reac, reag, prod)
        r["_balanced"] = b
        r["balance_status"] = b.status
        r["balance_reason"] = b.reason
        r["equation"] = b.equation
        r["equation_formula"] = b.equation_formula

    # 3. доказательства
    evid = EvidenceIndex.load(args.corpus, args.meta, args.manual)
    for r in results:
        reac, _, prod = _parts(r["reaction_smiles"])
        level, rid, src = evid.level(reac, prod)
        r["evidence_level"] = level
        r["evidence_id"] = rid
        r["evidence_source"] = src
        r["conditions"] = evid.conditions.get(rid or r["precedent_id"], {})

    # 4. экономика: себестоимость веществ, получаемых в сети, — по маршруту
    producers: dict[str, list[dict]] = {}
    for r in results:
        b = r["_balanced"]
        if b.status != "failed":
            producers.setdefault(b.products[0][0], []).append(r)
    yield_ = cfg.get("default_yield", 1.0)
    credit = cfg.get("byproduct_credit", 0.0)
    memo: dict[str, tuple[Price | None, dict | None]] = {}

    def route_cost(smiles: str, visiting=frozenset()):
        key = canon(smiles) or smiles
        if key in memo:
            return memo[key][0]
        if key in visiting or key not in producers:
            return None
        best = (None, None)
        for p in producers[key]:
            b = p["_balanced"]
            eco = evaluate(
                b.reactants,
                b.products,
                book,
                lambda s, v=visiting | {key}: route_cost(s, v),
                yield_,
                credit,
            )
            cost = eco.raw_cost
            if cost is not None:
                cost -= eco.credit
                if best[0] is None or cost < best[0].usd_per_t:
                    best = (Price(cost, "route", f"маршрут: {names(key)}"), p)
        memo[key] = best
        return best[0]

    def made_in_network(smiles: str) -> bool:
        """Вещество посчитано по себестоимости маршрута («сделать» дешевле «купить»)."""
        key = canon(smiles) or smiles
        route = memo.get(key, (None, None))[0]
        market = book.exact.get(key)
        return route is not None and (market is None or route.usd_per_t < market.usd_per_t)

    def route_steps(smiles: str, seen=None) -> list[dict]:
        seen = set() if seen is None else seen
        key = canon(smiles) or smiles
        p = memo.get(key, (None, None))[1]
        if p is None or key in seen or not made_in_network(key):
            return []
        seen.add(key)
        steps = []
        for s, _ in p["_balanced"].reactants:
            steps += route_steps(s, seen)
        return steps + [p]

    for r in results:
        b = r["_balanced"]
        if b.status == "failed":
            r["_eco"] = None
            continue
        eco = evaluate(b.reactants, b.products, book, route_cost, yield_, credit)
        r["_eco"] = eco
        r["_route"] = []
        for s, _ in b.reactants:
            r["_route"] += [x for x in route_steps(s) if x not in r["_route"]]
    for r in results:
        eco = r.get("_eco")
        r["raw_cost_usd_t"] = None if eco is None else eco.raw_cost
        r["revenue_usd_t"] = None if eco is None else eco.revenue
        r["margin_usd_t"] = None if eco is None else eco.margin
        r["product_price_basis"] = (
            "" if eco is None or not eco.product_price else eco.product_price.basis
        )

    # 5. векторная БД
    db = VectorDB()
    for r in results:
        db.add(r["reaction_core"], r["reaction_smiles"], {"type_key": r["type_key"]})
    db.save(os.path.join(args.out, "vectordb"))
    known = VectorDB()
    for m in _read(args.manual):
        known.add(m["id"], m["rxn_smiles"].replace(">>", ">>"), {"source": m["source"]})
    if args.known:
        for c in _read(args.known):
            known.add(c["id"], c["rxn_smiles"], {"source": c.get("source", "")})
    for r in results:
        r["similar_known"] = [
            {"id": k, "similarity": round(s, 3), "source": m.get("source", ""), "rxn": m["rxn"]}
            for s, k, m in known.search(r["reaction_smiles"], k=3)
        ]

    # 6. ранжирование
    w = cfg["ranking_weights"]
    for r in results:
        bal = BALANCE_SCORE.get(r["balance_status"], 0.0) * (
            0.5 if "H2" in r["balance_reason"] else 1
        )
        comp = {
            "margin": _margin_score(r),
            "evidence": _evidence_score(r),
            "route": 1 / int(r["level"]),
            "balance": bal,
            "internal": 1.0 if r["partner_source"] == "internal" else 0.7,
        }
        r["score_parts"] = {k: round(v, 3) for k, v in comp.items()}
        r["score"] = round(sum(w[k] * comp[k] for k in w), 4)
    results.sort(key=lambda r: -r["score"])
    for i, r in enumerate(results, 1):
        r["rank"] = i

    # 7. отчёт
    summary = _summary(results, db, cfg, args, time.time() - t0)
    _write_outputs(results, summary, names, cfg, args)
    return summary


def _summary(results, db, cfg, args, elapsed) -> dict:
    valid = [r for r in results if r["balance_status"] != "failed"]
    series = {r["series"] for r in valid}
    types = {r["type_key"] for r in valid}
    priced = [r for r in valid if r["margin_usd_t"] is not None]
    return {
        "inputs": args.inputs,
        "templates": args.templates,
        "depth": args.depth,
        "elapsed_s": round(elapsed, 1),
        "reactions": len(results),
        "criterion_2_1": {
            "valid_reactions": len(valid),
            "series": len(series),
            "types": len(types),
            "level_I": sum(r["evidence_level"] == "I" for r in valid),
            "passed": len(types) >= 10,
        },
        "criterion_2_2": {
            "vectors": len(db),
            "unique_vectors": db.unique_vectors(),
            "passed": len(db) >= 10 and db.unique_vectors() == len(db),
        },
        "criterion_2_3": {
            "with_full_economics": len(priced),
            "share": round(len(priced) / len(valid), 3) if valid else 0,
            "price_note": cfg.get("price_note", ""),
            "accuracy": "см. make eval-economics (методика — решение 0016)",
        },
    }


def _row(r, names, fx) -> dict:
    eco = r.get("_eco")
    return {
        "rank": r["rank"],
        "score": r["score"],
        "score_parts": r["score_parts"],
        "stage": int(r["level"]),
        "reaction": r["reaction_smiles"],
        "reaction_names": " + ".join(names(s) for s in _parts(r["reaction_smiles"])[0])
        + " → "
        + names(r["product"]),
        "product": r["product"],
        "product_name": names(r["product"]),
        "type_key": r["type_key"],
        "series": r["series"],
        "equation": r["equation"],
        "equation_formula": r["equation_formula"],
        "product_balanced": r["_balanced"].products[0][0] if r["_balanced"].products else "",
        "balance_status": r["balance_status"],
        "balance_reason": r["balance_reason"],
        "reagents": sorted(r["reagents"]),
        "catalyst": r.get("catalyst", ""),
        "economics_usd_per_t": None
        if eco is None
        else {
            "raw_cost": _round(eco.raw_cost),
            "revenue": _round(eco.revenue),
            "byproduct_credit": _round(eco.credit),
            "margin": _round(eco.margin),
            "margin_rub": _round(None if eco.margin is None else eco.margin * fx),
            "product_price_basis": r["product_price_basis"],
            "raw_materials": [
                {
                    "name": names(line.smiles),
                    "t_per_t": round(line.t_per_t, 4),
                    "price": None if line.price is None else round(line.price.usd_per_t, 1),
                    "basis": "" if line.price is None else line.price.basis,
                }
                for line in eco.raw
            ],
            "unknown_prices": [names(s) if s != "продукт" else s for s in eco.unknown],
        },
        "route": [x["equation_formula"] or x["reaction_smiles"] for x in r.get("_route", [])],
        "evidence": {
            "level": r["evidence_level"],
            "in_corpus": r["evidence_id"],
            "precedent_id": r["precedent_id"],
            "precedent_reaction": r["precedent_rxn"],
            "precedent_source": r["precedent_source"],
            "template_id": r["template_id"],
            "template_precedents": int(r["template_count"] or 0),
            "conditions": r["conditions"],
        },
        "similar_known": r["similar_known"],
    }


def _round(x):
    return None if x is None else round(x, 1)


def _write_outputs(results, summary, names, cfg, args):
    fx = cfg.get("fx_rub_per_usd", 1)
    rows = [_row(r, names, fx) for r in results]
    with open(os.path.join(args.out, "results.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": summary, "reactions": rows}, f, ensure_ascii=False, indent=1)
    flat = [
        "rank",
        "score",
        "stage",
        "reaction_names",
        "equation_formula",
        "balance_status",
        "catalyst",
        "type_key",
    ]
    with open(os.path.join(args.out, "results.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            flat
            + ["raw_cost_usd_t", "revenue_usd_t", "margin_usd_t", "evidence_level", "precedent_id"]
            + ["precedent_source", "reaction"]
        )
        for row in rows:
            e = row["economics_usd_per_t"] or {}
            w.writerow(
                [row[k] for k in flat]
                + [e.get("raw_cost"), e.get("revenue"), e.get("margin"), row["evidence"]["level"]]
                + [
                    row["evidence"]["precedent_id"],
                    row["evidence"]["precedent_source"],
                    row["reaction"],
                ]
            )
    with open(os.path.join(args.out, "summary.md"), "w", encoding="utf-8") as f:
        f.write(_summary_md(summary, rows))
    from chem_agent.report import render_html

    with open(os.path.join(args.out, "report.html"), "w", encoding="utf-8") as f:
        f.write(render_html(summary, rows, cfg, _read(args.inputs)))


def _summary_md(s, rows) -> str:
    c1, c2, c3 = s["criterion_2_1"], s["criterion_2_2"], s["criterion_2_3"]
    ok = lambda b: "выполнен" if b else "не выполнен"  # noqa: E731
    lines = [
        "# Сводка анализа",
        "",
        f"Вход: `{s['inputs']}`; стадий: {s['depth']}; время: {s['elapsed_s']} с.",
        "",
        f"- **2.1** — валидных (уравненных) реакций: {c1['valid_reactions']}, серий: "
        f"{c1['series']}, **типов превращений: {c1['types']}**, из них в корпусе (уровень I): "
        f"{c1['level_I']}. Критерий ≥10 по типам — {ok(c1['passed'])}.",
        f"- **2.2** — векторов: {c2['vectors']}, уникальных: {c2['unique_vectors']} — "
        f"{ok(c2['passed'])}.",
        f"- **2.3** — полная экономика (цены сырья и продукта): {c3['with_full_economics']} "
        f"({c3['share']:.0%}). {c3['price_note']}",
        "",
        "## Топ-10",
        "",
        "| # | Реакция | Маржа, $/т | Доказательность |",
        "|---|---|---|---|",
    ]
    for r in rows[:10]:
        e = r["economics_usd_per_t"] or {}
        lines.append(
            f"| {r['rank']} | {html.escape(r['reaction_names'])} | {e.get('margin', '—')} | "
            f"{r['evidence']['level']} ({r['evidence']['precedent_id']}) |"
        )
    return "\n".join(lines) + "\n"


def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--inputs", required=True, help="CSV name,smiles — набор пользователя")
    ap.add_argument("--templates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--utilities", default=os.path.join(ROOT, D["utilities"]))
    ap.add_argument("--manual", default=os.path.join(ROOT, D["manual"]))
    ap.add_argument("--corpus", help="corpus.csv (corpus_prep) — для уровня I")
    ap.add_argument("--meta", help="meta.csv (corpus_prep) — условия прецедентов")
    ap.add_argument(
        "--known", help="CSV id,rxn_smiles,source — известные процессы для поиска похожих"
    )
    ap.add_argument("--prices", default=os.path.join(ROOT, D["prices"]))
    ap.add_argument("--classes", default=os.path.join(ROOT, D["classes"]))
    ap.add_argument("--config", default=os.path.join(ROOT, D["config"]))
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    s = run(args)
    print(
        f"[analyze] реакций {s['reactions']}, типов {s['criterion_2_1']['types']}, "
        f"с экономикой {s['criterion_2_3']['with_full_economics']} -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
