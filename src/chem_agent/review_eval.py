#!/usr/bin/env python3
"""
review_eval.py — измерение качества генерации на реакциях с вердиктом химика.

Размеченный набор (tests/regression/review_verdicts.csv) — реакции из ревью
chemistry-reviewer с вердиктом OK / doubt / error и причиной ошибки
(reagent — неверный сопутствующий реагент, template — неверный шаблон).

Все вещества этих реакций подаются движку как набор пользователя
(apply --internal-only, одна стадия), и проверяется, какие размеченные
реакции система предлагает сейчас:
  - OK: реакция должна остаться (по реагентам и продукту, реагент любой);
  - error/reagent: не должно остаться именно этого варианта с этим реагентом;
  - error/template: не должно остаться реакции с этими реагентами и продуктом.
Противоионы при сравнении не учитываются.

Пример:
  python -m chem_agent.review_eval --templates outputs/<прогон>/templates.jsonl \\
                                   --out-dir outputs/<прогон>/eval
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

from chem_agent.template_engine import COUNTER_IONS, canon, frag_set
from chem_agent.template_engine import main as engine_main

DEFAULT_VERDICTS = os.path.join(
    os.path.dirname(__file__), "..", "..", "tests", "regression", "review_verdicts.csv"
)


def split_reaction(rxn: str):
    reac, reag, prod = rxn.split(">")
    return reac, reag, prod


def core_key(reac: str, prod: str):
    return frag_set(reac), canon(prod)


def full_key(reac: str, reag: str, prod: str):
    return frag_set(reac), frag_set(reag) if reag else frozenset(), canon(prod)


def input_substances(rows) -> list[str]:
    """Вещества размеченных реакций; NaOH — одним веществом, противоионы отдельно не нужны."""
    out = []
    for r in rows:
        reac, reag, _ = split_reaction(r["reaction_smiles"])
        frags = [canon(f) for f in (reac + "." + reag).split(".") if f]
        for f in frags:
            if f is None or f in COUNTER_IONS:
                continue
            out.append("[Na+].[OH-]" if f == "[OH-]" else f)
    return list(dict.fromkeys(out))


def evaluate(rows, candidates) -> list[dict]:
    cores = {core_key(*split_reaction(c)[::2]) for c in candidates}
    fulls = {full_key(*split_reaction(c)) for c in candidates}
    out = []
    for r in rows:
        reac, reag, prod = split_reaction(r["reaction_smiles"])
        if r["verdict"] == "error" and r["cause"] == "reagent":
            present = full_key(reac, reag, prod) in fulls
        else:
            present = core_key(reac, prod) in cores
        out.append({**r, "generated": int(present)})
    return out


def summary(results) -> list[str]:
    c = collections.Counter()
    for r in results:
        g = r["verdict"] + (f"/{r['cause']}" if r["cause"] else "")
        c[(g, r["generated"])] += 1
    lines = ["| Вердикт | Всего | Генерируется сейчас |", "|---|---|---|"]
    for g in ("OK", "doubt", "error/reagent", "error/template"):
        n1, n0 = c[(g, 1)], c[(g, 0)]
        if n1 + n0:
            lines.append(f"| {g} | {n1 + n0} | {n1} ({n1 / (n1 + n0):.0%}) |")
    ok = [r for r in results if r["verdict"] == "OK"]
    err = [r for r in results if r["verdict"] == "error"]
    kept = sum(r["generated"] for r in ok)
    leaked = sum(r["generated"] for r in err)
    lines += [
        "",
        f"- Верные реакции сохранены: **{kept}/{len(ok)}**",
        f"- Ошибочные всё ещё выдаются: **{leaked}/{len(err)}**",
    ]
    return lines


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--verdicts", default=os.path.normpath(DEFAULT_VERDICTS))
    ap.add_argument("--templates", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-count", default="2")
    ap.add_argument("--jobs", default="0")
    args = ap.parse_args(argv)

    with open(args.verdicts, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    os.makedirs(args.out_dir, exist_ok=True)
    inputs = os.path.join(args.out_dir, "inputs.csv")
    with open(inputs, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "smiles"])
        for i, s in enumerate(input_substances(rows)):
            w.writerow([f"S{i:03d}", s])
    cand = os.path.join(args.out_dir, "candidates.csv")
    engine_main(
        ["apply", "--templates", args.templates, "--inputs", inputs, "--internal-only"]
        + ["--min-count", args.min_count, "--jobs", args.jobs]
        + ["--out", cand, "--report", os.path.join(args.out_dir, "coverage.md")]
    )
    with open(cand, encoding="utf-8") as f:
        candidates = [r["reaction_smiles"] for r in csv.DictReader(f)]
    results = evaluate(rows, candidates)
    with open(os.path.join(args.out_dir, "eval.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    lines = ["# Качество на реакциях с вердиктом химика\n", *summary(results)]
    with open(os.path.join(args.out_dir, "eval.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines), file=sys.stderr)


if __name__ == "__main__":
    main()
