#!/usr/bin/env python3
"""
template_engine.py — шаблонный генератор реакций для химического агента.

Идея: вместо того чтобы просить LLM «придумать» реакцию, извлекаем из корпуса
известные превращения в виде шаблонов (reaction templates) и применяем их
к входным соединениям пользователя. Продукт строится по атомной разметке
шаблона, поэтому скелетный баланс атомов гарантирован по конструкции.

Три шага (каждый — отдельная подкоманда, результаты пишутся в файлы):

  map      атомный маппинг корпуса (RXNMapper). Уже размеченные реакции
           (с :номерами в SMILES) пропускаются без изменений.
  extract  извлечение шаблонов (RDChiral), канонизация, подсчёт частот,
           самопроверка: воспроизводит ли шаблон свою же исходную реакцию.
  apply    прямое применение шаблонов к входным соединениям,
           подбор партнёров, дедупликация и отчёт о покрытии.

Пример полного прогона:
  python -m chem_agent.template_engine map     --corpus corpus.csv --out mapped.csv
  python -m chem_agent.template_engine extract --mapped mapped.csv manual_mapped.csv \\
                                               --trusted manual_mapped.csv --out templates.jsonl
  python -m chem_agent.template_engine apply   --templates templates.jsonl --inputs inputs.csv \\
                                               --out candidates.csv --report coverage.md

Зависимости: rdkit, rdchiral, rxnmapper (только для шага map; ему нужен
setuptools<81 из-за pkg_resources).
"""

from __future__ import annotations

import argparse
import collections
import csv
import itertools
import json
import multiprocessing
import os
import re
import signal
import sys
import time
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

MAP_NUM_RE = re.compile(r":\d+\]")

# Прецедент без ссылки на первоисточник (решение 0007, 0008).
NO_SOURCE = "без источника"

# Растворители: если их атомы не вошли в продукт, это агенты, а не
# сопутствующие реагенты (решение 0008). Вода и спирты, вошедшие в продукт
# (гидролиз, переэтерификация), по разметке остаются реагирующими.
SOLVENTS = {
    "O",  # вода
    "CO",  # метанол
    "CCO",  # этанол
    "CC(C)O",  # изопропанол
    "C1CCOC1",  # THF
    "C1COCCO1",  # диоксан
    "COCCOC",  # DME
    "CCOCC",  # диэтиловый эфир
    "ClCCl",  # DCM
    "ClC(Cl)Cl",  # хлороформ
    "CN(C)C=O",  # DMF
    "CS(C)=O",  # DMSO
    "CN1CCCC1=O",  # NMP
    "CC#N",  # ацетонитрил
    "CC(C)=O",  # ацетон
    "CCOC(C)=O",  # этилацетат
    "Cc1ccccc1",  # толуол
    "c1ccccc1",  # бензол
    "CCCCCC",  # гексан
}


# ---------------------------------------------------------------------------
# Общие утилиты
# ---------------------------------------------------------------------------


def canon(smiles: str) -> str | None:
    """Канонический SMILES без атомной разметки; None, если не парсится."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def split_rxn(rxn: str) -> tuple[str, str, str] | None:
    """'A.B>C>D' -> ('A.B', 'C', 'D'); 'A.B>>D' -> ('A.B', '', 'D')."""
    parts = rxn.strip().split(">")
    if len(parts) != 3 or not parts[0] or not parts[2]:
        return None
    return parts[0], parts[1], parts[2]


def mapnums(smiles: str) -> set[int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    return {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}


def read_csv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def type_key(forward: str) -> str:
    """
    «Тип превращения» — скелет шаблона без атомной разметки и без уточнений
    по H/степени атомов. Склеивает шаблоны, которые отличаются только
    порядком атомов в строке или степенью конкретизации (метанол vs любой
    спирт), и различает реально разные превращения. Используется для
    подсчёта разнообразия, а не для применения.
    """
    rxn = AllChem.ReactionFromSmarts(forward)

    def skel(m):
        m = Chem.Mol(m)
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        return Chem.MolToSmiles(m)

    r = sorted(skel(rxn.GetReactantTemplate(i)) for i in range(rxn.GetNumReactantTemplates()))
    p = sorted(skel(rxn.GetProductTemplate(i)) for i in range(rxn.GetNumProductTemplates()))
    return ".".join(r) + ">>" + ".".join(p)


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def run_products_traced(
    rxn, reactant_mols, max_products=200, keep_new_stereo: bool = True
) -> dict[str, dict]:
    """
    Прямой прогон шаблона: {канонический SMILES продукта: {номер атома в
    шаблоне: (номер реагента, индекс атома в нём)}} — откуда пришёл каждый
    атом шаблона (для продукта берётся первое совпадение).
    """
    out: dict[str, dict] = {}
    try:
        product_sets = rxn.RunReactants(tuple(reactant_mols), max_products)
    except Exception:
        return out
    for ps in product_sets:
        for p in ps:
            if not keep_new_stereo:
                # стереоцентр шаблона не переносим на ахиральный субстрат; у
                # хирального субстрата новый центр задаётся каркасом — оставляем
                for a in p.GetAtoms():
                    if a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED or not a.HasProp(
                        "react_atom_idx"
                    ):
                        continue
                    src_mol = reactant_mols[a.GetIntProp("react_idx")]
                    src = src_mol.GetAtomWithIdx(a.GetIntProp("react_atom_idx"))
                    achiral = all(
                        x.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
                        for x in src_mol.GetAtoms()
                    )
                    if src.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED and achiral:
                        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
            except Exception:
                continue
            c = canon(smi)
            if c and c not in out:
                out[c] = {
                    a.GetIntProp("old_mapno"): (
                        a.GetIntProp("react_idx"),
                        a.GetIntProp("react_atom_idx"),
                    )
                    for a in p.GetAtoms()
                    if a.HasProp("old_mapno")
                }
    return out


def run_products(rxn, reactant_mols, max_products=200) -> set[str]:
    """Прямой прогон шаблона; возвращает канонические SMILES продуктов."""
    return set(run_products_traced(rxn, reactant_mols, max_products))


# --- Хемоселективность (R17, решение 0010) ----------------------------------
# Класс реагирующего атома в реальной молекуле и «сила» нуклеофила. Шаблон
# радиуса 1 видит только ближайших соседей: [OH] кислоты для него — тот же
# [OH], что у спирта, а амидный N — тот же N, что у амина.

NUC_RANK = {
    "thiol": 5,
    "thiocarbonyl": 5,
    "alkoxide": 4,
    "amine": 4,
    "aniline": 3,
    "alcohol": 2,
    "phenol": 2,
    "aromatic_NH": 2,
}
_HALOGENS = {"F", "Cl", "Br", "I"}


def _is_acyl_like(atom) -> bool:
    """C(=O), C(=S), C(=N), S(=O), P(=O): соседний атом делает N амидом, O — кислотным."""
    return atom.GetSymbol() in ("C", "S", "P") and any(
        b.GetBondType() == Chem.BondType.DOUBLE
        and b.GetOtherAtom(atom).GetSymbol() in ("O", "S", "N")
        for b in atom.GetBonds()
    )


def _activated_aryl(mol, atom) -> bool:
    """Арилгалогенид активирован для SNAr: гетероароматика или акцептор в кольце."""
    rings = [r for r in mol.GetRingInfo().AtomRings() if atom.GetIdx() in r]
    for ring in rings:
        for i in ring:
            a = mol.GetAtomWithIdx(i)
            if a.GetSymbol() == "N" and a.GetIsAromatic():
                return True
            for n in a.GetNeighbors():
                if n.GetIdx() in ring:
                    continue
                if n.GetSymbol() == "N" and n.GetFormalCharge() > 0:  # нитро
                    return True
                if n.GetSymbol() == "C" and any(
                    b.GetBondType() == Chem.BondType.TRIPLE for b in n.GetBonds()
                ):  # циано
                    return True
                if _is_acyl_like(n) or (
                    n.GetSymbol() == "C"
                    and sum(x.GetSymbol() == "F" for x in n.GetNeighbors()) == 3
                ):
                    return True
    return False


def atom_class(mol, idx: int) -> str | None:
    """Класс атома реакционного центра: спирт, фенол, амин, амид, тиол, арилгалогенид..."""
    a = mol.GetAtomWithIdx(idx)
    sym, q, h, nbrs = a.GetSymbol(), a.GetFormalCharge(), a.GetTotalNumHs(), a.GetNeighbors()
    double = any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds())
    if sym == "S":
        if double:
            return "thiocarbonyl" if a.GetDegree() == 1 else "sulfonyl"
        return "thiol" if (h or q < 0) else "sulfide"
    if sym == "O":
        if double:
            return "carbonyl_O"
        if any(_is_acyl_like(n) for n in nbrs):
            return "acid_O" if (h or q < 0) else "ester_O"
        if q < 0 and any(n.GetFormalCharge() > 0 for n in nbrs):
            return "nitro_O"  # [O-] при [N+], N-оксид: не нуклеофил
        if q < 0:
            return "alkoxide"
        if h:
            if not nbrs:
                return "water"
            return "phenol" if any(n.GetIsAromatic() for n in nbrs) else "alcohol"
        return "ether"
    if sym == "N":
        if a.GetIsAromatic():
            return "aromatic_NH" if h else "aromatic_N"
        if any(_is_acyl_like(n) for n in nbrs):
            return "amide_N"
        if q > 0:
            return "ammonium_N"
        if double or any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
            return "unsaturated_N"
        if any(n.GetIsAromatic() for n in nbrs):
            return "aniline" if h else "aryl_amine3"
        return "amine" if h else "amine3"
    if sym == "C" and a.GetIsAromatic() and any(n.GetSymbol() in _HALOGENS for n in nbrs):
        return "aryl_X_activated" if _activated_aryl(mol, a) else "aryl_X"
    return None


def competitor_rank(mol, idx: int) -> int:
    """Самый сильный нуклеофил молекулы, кроме атома idx (0 — нет)."""
    best = 0
    for a in mol.GetAtoms():
        if a.GetIdx() != idx and a.GetSymbol() in ("N", "O", "S"):
            best = max(best, NUC_RANK.get(atom_class(mol, a.GetIdx()), 0))
    return best


def changed_mapnos(rxn) -> set[int]:
    """Номера атомов шаблона, у которых меняются H, заряд, степень или связи."""

    def smarts(n, get):
        # атом и его соседи: при замещении (Cl -> O у кольца) меняются соседи
        return {
            a.GetAtomMapNum(): (
                re.sub(r":\d+\]$", "]", a.GetSmarts()),
                sorted(
                    (
                        str(b.GetOtherAtom(a).GetAtomMapNum() or b.GetOtherAtom(a).GetSymbol()),
                        b.GetSmarts(),
                    )
                    for b in a.GetBonds()
                ),
            )
            for i in range(n)
            for a in get(i).GetAtoms()
            if a.GetAtomMapNum()
        }

    reac = smarts(rxn.GetNumReactantTemplates(), rxn.GetReactantTemplate)
    prod = smarts(rxn.GetNumProductTemplates(), rxn.GetProductTemplate)
    return {k for k, v in reac.items() if prod.get(k) != v}


def heteroatom_migration(rxn) -> bool:
    """
    R16: в одной молекуле один углерод теряет связь с гетероатомом, а другой
    получает, без изменения кратности связи между ними. Так выглядит ошибочная
    запись корпуса (2-додеканол -> эфир 1-додеканола). Замещение (SN2, раскрытие
    эпоксида) меняет гетероатом у того же углерода, присоединение по кратной
    связи — только добавляет.
    """

    def carbons(n, get):
        count, slot, bonds = {}, {}, {}
        for i in range(n):
            for a in get(i).GetAtoms():
                k = a.GetAtomMapNum()
                if not k or a.GetAtomicNum() != 6:
                    continue
                slot[k] = i
                count[k] = sum(nb.GetAtomicNum() in _HETERO for nb in a.GetNeighbors())
                for b in a.GetBonds():
                    j = b.GetOtherAtom(a).GetAtomMapNum()
                    if j:
                        bonds[frozenset((k, j))] = b.GetBondTypeAsDouble()
        return count, slot, bonds

    rc, slot, rb = carbons(rxn.GetNumReactantTemplates(), rxn.GetReactantTemplate)
    pc, _, pb = carbons(rxn.GetNumProductTemplates(), rxn.GetProductTemplate)
    lose = [k for k in rc if k in pc and pc[k] < rc[k]]
    gain = [k for k in rc if k in pc and pc[k] > rc[k]]
    for a in lose:
        for b in gain:
            if slot[a] == slot[b] and rb.get(frozenset((a, b))) == pb.get(frozenset((a, b))):
                return True
    return False


def center_info(rxn, reactant_mols, provenance: dict) -> dict:
    """{номер атома: [класс, сила конкурента]} для атомов центра прецедента."""
    out = {}
    for k in changed_mapnos(rxn):
        if k not in provenance:
            continue
        ri, ai = provenance[k]
        cls = atom_class(reactant_mols[ri], ai)
        if cls is not None:
            out[str(k)] = [cls, competitor_rank(reactant_mols[ri], ai)]
    return out


# ---------------------------------------------------------------------------
# Шаг 1. Атомный маппинг
# ---------------------------------------------------------------------------


def load_mapper(device: str):
    """
    RXNMapper на выбранном устройстве. RXNMapper сам берёт GPU, если он есть;
    'cpu' — принудительно процессор; 'auto' — GPU, а если на нём нет памяти
    (занят другой моделью), то процессор.
    """
    import torch
    from rxnmapper import RXNMapper

    def on_cpu():
        available = torch.cuda.is_available
        torch.cuda.is_available = lambda: False  # RXNMapper выбирает устройство так
        try:
            return RXNMapper()
        finally:
            torch.cuda.is_available = available

    if device == "cpu" or not torch.cuda.is_available():
        mapper = on_cpu()
    else:
        try:
            mapper = RXNMapper()
        except Exception as e:  # torch.AcceleratorError / RuntimeError
            if device == "cuda" or "out of memory" not in str(e):
                raise
            print("[map] на GPU не хватает памяти — разметка на CPU", file=sys.stderr)
            mapper = on_cpu()
    print(f"[map] устройство: {mapper.device}", file=sys.stderr)
    return mapper


def cmd_map(args):
    rows = read_csv(args.corpus)
    out_rows, todo = [], []

    for i, row in enumerate(rows):
        rid = row.get(args.id_col) or f"R{i:07d}"
        parts = split_rxn(row.get(args.smiles_col, ""))
        if parts is None:
            continue
        reac, agents, prod = parts
        rec = {
            "id": rid,
            "rxn_original": row[args.smiles_col],
            "agents": agents,
            "mapped_rxn": "",
            "map_confidence": "",
            "source": row.get(args.source_col, ""),
        }
        if MAP_NUM_RE.search(reac + ">>" + prod):
            rec["mapped_rxn"] = reac + ">>" + prod  # уже размечено в корпусе
            out_rows.append(rec)
        else:
            # агенты (растворители, катализаторы) в маппинг не передаём:
            # в шаблон они не входят, а для QUARC сохраняются в колонке agents
            todo.append((rec, reac + ">>" + prod))

    print(
        f"[map] всего записей: {len(rows)}, уже размечено: {len(out_rows)}, "
        f"к разметке: {len(todo)}",
        file=sys.stderr,
    )

    if todo:
        mapper = load_mapper(args.device)
        t0 = time.time()
        for start in range(0, len(todo), args.batch):
            chunk = todo[start : start + args.batch]
            try:
                results = mapper.get_attention_guided_atom_maps([r for _, r in chunk])
            except Exception:
                # длинные реакции (>512 токенов) роняют батч — размечаем по одной
                results = []
                for _, r in chunk:
                    try:
                        results.append(mapper.get_attention_guided_atom_maps([r])[0])
                    except Exception:
                        results.append(None)
            for (rec, _), res in zip(chunk, results, strict=True):
                if res is None:
                    continue
                rec["mapped_rxn"] = res["mapped_rxn"]
                rec["map_confidence"] = f"{res['confidence']:.3f}"
                out_rows.append(rec)
            done = min(start + args.batch, len(todo))
            print(
                f"[map] {done}/{len(todo)} ({done / (time.time() - t0):.1f} р/с)", file=sys.stderr
            )

    write_csv(
        args.out,
        out_rows,
        ["id", "rxn_original", "agents", "mapped_rxn", "map_confidence", "source"],
    )
    print(f"[map] записано {len(out_rows)} -> {args.out}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Шаг 2. Извлечение шаблонов
# ---------------------------------------------------------------------------


def prepare_for_extraction(mapped_rxn: str):
    """
    Оставляем один основной продукт (самый тяжёлый фрагмент) — генератор
    нужен для главного продукта, побочные (вода, соли) достраивает балансировщик.
    Возвращает (реагирующие_mapped, основной_продукт_mapped,
                реагирующие_без_разметки, продукт_без_разметки, сопутствующие,
                растворители).
    Сопутствующие — реагенты, ни один атом которых не попал в основной продукт
    (основания, кислоты, окислители). Они не входят в шаблон, но по ним
    превращение привязывается к веществу: NaOH «ведёт» нейтрализацию, хотя
    его атомов в продукте нет. Неучаствующие растворители (SOLVENTS) —
    не сопутствующие реагенты, а агенты: их возвращаем отдельно.
    """
    reac, _, prod = mapped_rxn.split(">")
    prod_frags = [p for p in prod.split(".") if p]
    main = max(prod_frags, key=heavy_atoms)
    main_maps = mapnums(main)
    if not main_maps:
        return None
    reac_frags = [r for r in reac.split(".") if r]
    reacting = [r for r in reac_frags if mapnums(r) & main_maps]
    if not reacting:
        return None
    others = [c for r in reac_frags if not (mapnums(r) & main_maps) and (c := canon(r))]
    return (
        ".".join(reacting),
        main,
        [canon(r) for r in reacting],
        canon(main),
        [s for s in others if s not in SOLVENTS],
        [s for s in others if s in SOLVENTS],
    )


def extract_one(row: dict, timeout: int, no_intra: bool, max_reactants: int):
    """
    Шаблон из одной размеченной реакции. Чистая функция: не трогает общего
    состояния, поэтому реакции можно обрабатывать в разных процессах.
    Возвращает ("extracted", данные) или (причина_пропуска, None).
    """
    from rdchiral.template_extractor import extract_from_reaction

    prepared = prepare_for_extraction(row["mapped_rxn"]) if row.get("mapped_rxn") else None
    if prepared is None:
        return "skip_unparsable", None
    reac_m, prod_m, reac_plain, prod_plain, spectators, solvents = prepared
    if None in reac_plain or prod_plain is None:
        return "skip_unparsable", None

    signal.alarm(timeout)
    try:
        t = extract_from_reaction({"reactants": reac_m, "products": prod_m, "_id": row["id"]})
    except Timeout:
        return "skip_timeout", None
    except Exception:
        return "skip_error", None
    finally:
        signal.alarm(0)

    if not t or "reactants" not in t or not t.get("products"):
        return "skip_no_template", None
    if t.get("intra_only") and no_intra:
        return "skip_intra", None

    forward = f"{t['reactants']}>>{t['products']}"
    rxn = AllChem.ReactionFromSmarts(forward)
    n_slots = rxn.GetNumReactantTemplates()
    if n_slots > max_reactants:
        return "skip_too_many_reactants", None
    if heteroatom_migration(rxn):
        return "skip_heteroatom_migration", None

    # Самопроверка: прямой шаблон на исходных реагентах должен дать
    # записанный основной продукт (перебираем порядок реагентов по слотам).
    mols = [Chem.MolFromSmiles(s) for s in reac_plain]
    ok, center = False, {}
    for combo in itertools.permutations(mols, n_slots) if len(mols) >= n_slots else []:
        traced = run_products_traced(rxn, combo)
        if prod_plain in traced:
            ok = True
            center = center_info(rxn, combo, traced[prod_plain])
            break

    return "extracted", {
        "forward": forward,
        "retro": t["reaction_smarts"],
        "n_reactants": n_slots,
        "necessary_reagent": t.get("necessary_reagent", ""),
        "trusted": row["_trusted"],
        "ok": ok,
        "example": {
            "id": row["id"],
            "reactants": reac_plain,
            "product": prod_plain,
            "spectators": spectators,
            "agents": ".".join(a for a in [row.get("agents", "")] + solvents if a),
            "source": row.get("source", ""),
            "selfcheck": ok,
            "center": center,  # классы атомов центра (R17)
        },
    }


# Параллельная обработка: процессы создаются через fork и получают состояние
# (аргументы, шаблоны, пул веществ) из глобальной переменной, а не через pickle —
# объекты RDKit/rdchiral не сериализуются.
_WORKER_STATE: dict = {}


def n_jobs(requested: int) -> int:
    """0 — все ядра. Без fork (Windows) — один процесс: там и SIGALRM нет."""
    if "fork" not in multiprocessing.get_all_start_methods():
        return 1
    return requested if requested > 0 else (os.cpu_count() or 1)


def parallel_map(func, n_items: int, jobs: int):
    """func(i) для i in range(n_items); результаты — в исходном порядке."""
    if jobs <= 1 or n_items < 2:
        return map(func, range(n_items))
    pool = multiprocessing.get_context("fork").Pool(jobs)
    chunk = max(1, n_items // (jobs * 16))
    return _closing_imap(pool, func, n_items, chunk)


def _closing_imap(pool, func, n_items, chunk):
    try:
        yield from pool.imap(func, range(n_items), chunksize=chunk)
    finally:
        pool.terminate()


def _extract_worker(i: int):
    s = _WORKER_STATE
    return extract_one(s["rows"][i], s["timeout"], s["no_intra"], s["max_reactants"])


def _apply_worker(i: int):
    s = _WORKER_STATE
    nt = len(s["templates"])
    return s["unit"](s["active"][i // nt], s["templates"][i % nt], s["level"])


def cmd_extract(args):
    trusted_paths = set(args.trusted or [])
    rows = []
    for path in args.mapped:
        for row in read_csv(path):
            row["_trusted"] = path in trusted_paths
            rows.append(row)
    signal.signal(signal.SIGALRM, _alarm)  # обработчик наследуется процессами при fork

    templates: dict[str, dict] = {}
    stats = collections.Counter()
    jobs = n_jobs(args.jobs)
    print(f"[extract] процессов: {jobs}", file=sys.stderr)
    _WORKER_STATE.update(
        rows=rows, timeout=args.timeout, no_intra=args.no_intra, max_reactants=args.max_reactants
    )

    # Результаты приходят в порядке строк: шаблоны, частоты и примеры
    # собираются так же, как при последовательной обработке.
    for n, (status, x) in enumerate(parallel_map(_extract_worker, len(rows), jobs), 1):
        if n % 5000 == 0:
            print(f"[extract] {n}/{len(rows)}; шаблонов: {len(templates)}", file=sys.stderr)
        if status != "extracted":
            stats[status] += 1
            continue
        rec = templates.get(x["forward"])
        if rec is None:
            rec = templates[x["forward"]] = {
                "forward": x["forward"],
                "retro": x["retro"],
                "n_reactants": x["n_reactants"],
                "necessary_reagent": x["necessary_reagent"],
                "count": 0,
                "selfcheck_pass": 0,
                "examples": [],
                "trusted": False,
                "_seen": set(),
            }
        rec["trusted"] = rec["trusted"] or x["trusted"]
        stats["extracted"] += 1
        # Частота — число разных прецедентов (реагирующие вещества + продукт),
        # а не записей: одна реакция патента и её копия с другими растворителями
        # (2naoh_dataset) — одно наблюдение (R14, решение 0009).
        ex = x["example"]
        precedent = (tuple(sorted(ex["reactants"])), ex["product"])
        if precedent in rec["_seen"]:
            stats["duplicate_precedent"] += 1
            continue
        rec["_seen"].add(precedent)
        rec["count"] += 1
        rec["selfcheck_pass"] += int(x["ok"])
        if len(rec["examples"]) < args.max_examples:
            rec["examples"].append(ex)
    _WORKER_STATE.clear()

    ordered = sorted(templates.values(), key=lambda r: -r["count"])
    with open(args.out, "w", encoding="utf-8") as f:
        for i, rec in enumerate(ordered, 1):
            del rec["_seen"]
            rec["id"] = f"T{i:05d}"
            rec["type_key"] = type_key(rec["forward"])
            rec["selfcheck_rate"] = round(rec["selfcheck_pass"] / rec["count"], 3)
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    total_ok = sum(r["selfcheck_pass"] for r in ordered)
    print(
        f"[extract] типов превращений (скелетный ключ): {len({r['type_key'] for r in ordered})}",
        file=sys.stderr,
    )
    print(f"[extract] реакций на входе: {len(rows)}", file=sys.stderr)
    for k, v in sorted(stats.items()):
        print(f"  {k}: {v}", file=sys.stderr)
    print(
        f"[extract] уникальных шаблонов: {len(ordered)}; "
        f"с частотой >=2: {sum(r['count'] >= 2 for r in ordered)}",
        file=sys.stderr,
    )
    n_precedents = sum(r["count"] for r in ordered)
    if n_precedents:
        print(
            f"[extract] разных прецедентов: {n_precedents}; самопроверка пройдена: "
            f"{total_ok}/{n_precedents} ({total_ok / n_precedents:.0%})",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Шаг 3. Применение шаблонов к входным соединениям
# ---------------------------------------------------------------------------

# Противоионы не считаем «реагентами» при сравнении наборов веществ.
COUNTER_IONS = {"[Na+]", "[K+]", "[Li+]", "[Cl-]", "[Br-]", "[I-]", "[H+]"}

# Функции сопутствующего реагента (решение 0009). Вещество набора может быть
# сопутствующим реагентом шаблона, только если закрывает функции, которые
# нужны шаблону по его прецедентам: NaOH не окисляет и не восстанавливает.
BASE, OXIDANT, REDUCTANT, ACID = "base", "oxidant", "reductant", "acid"
METAL_CATALYST = "metal_catalyst"
ACTIVATOR = "activator"  # MsCl, TsCl, SOCl2, PCl3/POCl3, PPh3, DEAD, DCC
FUNCTIONS = (BASE, OXIDANT, REDUCTANT, ACID, METAL_CATALYST, ACTIVATOR)
# Требование превращения «активировать OH как уходящую группу»: закрывается
# кислотой или активатором.
ACTIVATION = "activation"
# Условия процесса, а не сырьё: не блокируют реакцию, указываются в выдаче.
CONDITIONS = {METAL_CATALYST}
REQUIRED_SHARE = 0.5  # функция нужна шаблону, если есть в >50% прецедентов
MIN_FUNCTION_SHARE = 0.25  # функция вещества должна встречаться в ≥25% прецедентов

_BASE_SMILES = {"[OH-]", "O=C([O-])[O-]", "O=C([O-])O", "[NH2-]", "O=P([O-])([O-])[O-]"}
_OXIDANT_SMILES = {"O=O", "S1SSSSSSS1", "[O-]Cl", "O=Cl[O-]", "OO"}
_OXIDANT_METALS = {"Mn", "Cr", "Os", "Se"}
_CATALYST_METALS = {"Pd", "Pt", "Ni", "Rh", "Ru", "Ir", "Au", "Ag", "Co", "Cu", "V"}
_LEWIS_ACID_ATOMS = {"Al", "B", "Ti", "Zn", "Sn", "Fe"}
_LEWIS_ACID_IONS = {"[Al+3]", "[Ti+4]", "[Zn+2]", "[Fe+3]", "[Sn+4]"}
_HALIDE_IONS = {"[F-]", "[Cl-]", "[Br-]", "[I-]"}
_ACID_SMILES = {
    "Cl",  # HCl
    "Br",  # HBr
    "O=S(=O)(O)O",  # H2SO4
    "CS(=O)(=O)O",  # MsOH
    "Cc1ccc(S(=O)(=O)O)cc1",  # TsOH
    "O=C(O)C(F)(F)F",  # TFA
    "CC(=O)O",  # AcOH
    "O=P(O)(O)O",  # H3PO4
    "[H+]",
}


def covers(required: set[str], functions: set[str]) -> bool:
    """Закрывает ли набор функций требования шаблона (катализатор — условие)."""
    needed = required - CONDITIONS
    if ACTIVATION in needed:
        if not functions & {ACID, ACTIVATOR}:
            return False
        needed = needed - {ACTIVATION}
    return needed <= functions


def _fragment_functions(smiles: str, has_hydride_carrier: bool) -> set[str]:
    if smiles in _BASE_SMILES:
        return {BASE}
    if smiles in _OXIDANT_SMILES:
        return {OXIDANT}
    if smiles in _ACID_SMILES:
        return {ACID}
    if smiles == "[H][H]":
        return {REDUCTANT}
    if smiles == "[H-]":  # NaH — основание; LiAlH4, NaBH4 — восстановитель
        return {REDUCTANT} if has_hydride_carrier else {BASE}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    out = set()
    for a in mol.GetAtoms():
        sym, q = a.GetSymbol(), a.GetFormalCharge()
        halo_nbr = any(n.GetSymbol() in ("Cl", "Br", "I") for n in a.GetNeighbors())
        dbl = [b for b in a.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        if sym in _CATALYST_METALS:
            out.add(METAL_CATALYST)  # Pd, Pt, Ni, V2O5...
        elif sym == "S" and halo_nbr and dbl:
            out.add(ACTIVATOR)  # MsCl, TsCl, SOCl2
        elif sym == "P" and (halo_nbr or (a.GetDegree() == 3 and not dbl)):
            out.add(ACTIVATOR)  # PCl3, POCl3, PBr3, PPh3
        elif sym == "N" and any(b.GetOtherAtom(a).GetSymbol() == "N" for b in dbl):
            out.add(ACTIVATOR)  # азодикарбоксилаты (DEAD, DIAD)
        elif sym == "C" and sum(b.GetOtherAtom(a).GetSymbol() == "N" for b in dbl) == 2:
            out.add(ACTIVATOR)  # карбодиимиды (DCC, EDC)
        elif sym in ("B", "Al") and (a.GetTotalNumHs() > 0 or q < 0):
            out.add(REDUCTANT)  # борогидриды, алюмогидриды
        elif (
            sym in _LEWIS_ACID_ATOMS
            and q == 0
            and any(n.GetSymbol() in ("F", "Cl", "Br", "I") for n in a.GetNeighbors())
        ):
            out.add(ACID)  # кислоты Льюиса: AlCl3, BF3, TiCl4, ZnCl2
        elif sym in _OXIDANT_METALS or (sym == "I" and a.GetDegree() > 1):
            out.add(OXIDANT)  # перманганат, хромовые, OsO4, Десс–Мартин, IBX
        elif sym == "O" and any(n.GetSymbol() == "O" for n in a.GetNeighbors()):
            out.add(OXIDANT)  # пероксиды, mCPBA
        elif sym == "O" and q < 0 and a.GetDegree() <= 1:
            nbr = a.GetNeighbors()
            # алкоксид, но не карбоксилат и не сульфонат
            if not nbr or all(b.GetBondType() != Chem.BondType.DOUBLE for b in nbr[0].GetBonds()):
                out.add(BASE)
        elif sym == "N" and q == 0 and not a.GetIsAromatic() and a.GetTotalNumHs() == 0:
            if a.GetDegree() == 3 and all(n.GetSymbol() == "C" for n in a.GetNeighbors()):
                if not any(
                    b.GetBondType() == Chem.BondType.DOUBLE
                    for n in a.GetNeighbors()
                    for b in n.GetBonds()
                ):
                    out.add(BASE)  # третичный амин (Et3N, DIPEA), не амид
        elif sym == "N" and a.GetIsAromatic() and a.GetTotalNumHs() == 0 and a.GetDegree() == 2:
            if mol.GetNumHeavyAtoms() <= 8:
                out.add(BASE)  # пиридин, лутидин
    return out


def reagent_functions(fragments) -> set[str]:
    """Функции набора сопутствующих веществ одного прецедента (или одного вещества)."""
    frags = [f for f in fragments if f]
    mols = [m for f in frags if (m := Chem.MolFromSmiles(f)) is not None]
    carrier = any(a.GetSymbol() in ("B", "Al") for m in mols for a in m.GetAtoms())
    out = set()
    for f in frags:
        out |= _fragment_functions(f, carrier)
    # кислота Льюиса в ионной записи: [Al+3] + 3 [Cl-] (без гидрида — не LiAlH4)
    if set(frags) & _LEWIS_ACID_IONS and set(frags) & _HALIDE_IONS and "[H-]" not in frags:
        out.add(ACID)
    return out


# Электроотрицательность (Полинг) для степеней окисления атомов шаблона.
_EN = {
    1: 2.20,
    3: 0.98,
    5: 2.04,
    6: 2.55,
    7: 3.04,
    8: 3.44,
    9: 3.98,
    11: 0.93,
    12: 1.31,
    13: 1.61,
    14: 1.90,
    15: 2.19,
    16: 2.58,
    17: 3.16,
    19: 0.82,
    29: 1.90,
    30: 1.65,
    34: 2.55,
    35: 2.96,
    50: 1.96,
    53: 2.66,
}
_HETERO = {7, 8, 9, 15, 16, 17, 34, 35, 53}


def _query_h_charge(atom) -> tuple[int | None, int]:
    """Число H и заряд атома шаблона RDChiral ('[C&H2&D2&+0:4]'); H=None, если не задано."""
    s = atom.GetSmarts()
    if not s.startswith("["):
        return None, 0
    body = re.sub(r":\d+\]$", "", s).strip("[]")
    h, q = None, 0
    for tok in re.split(r"[&;]", body):
        if re.fullmatch(r"H\d*", tok):
            h = int(tok[1:] or 1)
        elif m := re.fullmatch(r"([+-])(\d*)", tok):
            q = (1 if m.group(1) == "+" else -1) * int(m.group(2) or 1)
    return h, q


def _ox_state(atom, h: int, q: int) -> float:
    en = _EN.get(atom.GetAtomicNum())
    if en is None:
        return 0.0
    # каждый H: -1 у атома электроотрицательнее водорода, +1 у менее электроотрицательного
    ox = q - h * ((en > _EN[1]) - (en < _EN[1]))
    for b in atom.GetBonds():
        en_nb = _EN.get(b.GetOtherAtom(atom).GetAtomicNum(), en)
        if en_nb > en:
            ox += b.GetBondTypeAsDouble()
        elif en_nb < en:
            ox -= b.GetBondTypeAsDouble()
    return ox


def intrinsic_functions(forward: str) -> set[str]:
    """
    Какой реагент нужен самому превращению (R15, решение 0010), независимо от
    того, что записано в прецедентах:
      - размеченный фрагмент окисляется (сумма степеней окисления атомов
        центра растёт), а окислителя в шаблоне нет — нужен окислитель;
        окислителем в шаблоне считается связь гетероатом–гетероатом
        (O2, SO3, S8, Br2, HNO3);
      - восстанавливается, а восстановителя (B, Al, Si, металл, H2) нет —
        нужен восстановитель;
      - атом становится менее отрицательным и получает H — нужна кислота;
        наоборот — основание.
    """
    rxn = AllChem.ReactionFromSmarts(forward)
    slots_all = [rxn.GetReactantTemplate(i) for i in range(rxn.GetNumReactantTemplates())]

    def mapped(n, get):
        return {
            a.GetAtomMapNum(): a for i in range(n) for a in get(i).GetAtoms() if a.GetAtomMapNum()
        }

    reac = mapped(rxn.GetNumReactantTemplates(), rxn.GetReactantTemplate)
    prod = mapped(rxn.GetNumProductTemplates(), rxn.GetProductTemplate)
    delta, need = 0.0, set()
    protonated = deprotonated = gains_h = loses_h = False
    for k, ar in reac.items():
        ap = prod.get(k)
        if ap is None:
            continue
        (hr, qr), (hp, qp) = _query_h_charge(ar), _query_h_charge(ap)
        if hr is None or hp is None:
            continue
        delta += _ox_state(ap, hp, qp) - _ox_state(ar, hr, qr)
        if ar.GetAtomicNum() in _HETERO:
            protonated |= qp > qr and hp > hr
            deprotonated |= qp < qr and hp < hr
            gains_h |= hp > hr
            loses_h |= hp < hr
    # Протон, перешедший внутри реакции (OH бисульфита -> O раскрытого
    # эпоксида), внешней кислоты или основания не требует.
    # OH/OR не уходит с sp3-углерода без активации (кислота или MsCl, PPh3...)
    for m in slots_all:
        for a in m.GetAtoms():
            if a.GetAtomicNum() != 8 or a.GetAtomMapNum():
                continue
            for b in a.GetBonds():
                c = b.GetOtherAtom(a)
                if (
                    c.GetAtomicNum() == 6
                    and c.GetAtomMapNum()
                    and not c.GetIsAromatic()
                    and b.GetBondType() == Chem.BondType.SINGLE
                    and not any(x.GetBondType() != Chem.BondType.SINGLE for x in c.GetBonds())
                ):
                    need.add(ACTIVATION)
    if protonated and not loses_h:
        need.add(ACID)
    if deprotonated and not gains_h:
        need.add(BASE)
    slots = [rxn.GetReactantTemplate(i) for i in range(rxn.GetNumReactantTemplates())]
    has_oxidant = any(
        b.GetBeginAtom().GetAtomicNum() in _HETERO and b.GetEndAtom().GetAtomicNum() in _HETERO
        for m in slots
        for b in m.GetBonds()
    )
    has_reductant = any(
        a.GetAtomicNum() == 1 or (_EN.get(a.GetAtomicNum(), 9) < 2.1)
        for m in slots
        for a in m.GetAtoms()
    )
    if delta >= 1 and not has_oxidant:
        need.add(OXIDANT)
    elif delta <= -1 and not has_reductant:
        need.add(REDUCTANT)
    return need


def frag_set(smiles: str) -> frozenset[str]:
    """Набор канонических фрагментов без противоионов: '[Na+].[OH-]' -> {'[OH-]'}."""
    out = set()
    for f in smiles.split("."):
        c = canon(f)
        if c and c not in COUNTER_IONS:
            out.add(c)
    return frozenset(out)


@dataclass
class Template:
    id: str
    forward: str
    retro: str
    count: int
    selfcheck_rate: float
    examples: list
    type_key: str
    rxn: object = None
    retro_rxn: object = None
    slots: list = field(default_factory=list)
    slot_keys: list = field(default_factory=list)
    spectator_sets: list = field(default_factory=list)

    def build(self):
        from rdchiral.main import rdchiralReaction

        self.rxn = AllChem.ReactionFromSmarts(self.forward)
        self.rxn.Initialize()
        self.slots = [
            self.rxn.GetReactantTemplate(i) for i in range(self.rxn.GetNumReactantTemplates())
        ]
        self.slot_keys = [Chem.MolToSmarts(s) for s in self.slots]
        try:
            self.retro_rxn = rdchiralReaction(self.retro)
        except Exception:
            self.retro_rxn = None
        self.spectator_sets = [frag_set(".".join(ex.get("spectators", []))) for ex in self.examples]
        # Какие функции (основание, окислитель, восстановитель) выполняли
        # сопутствующие вещества в прецедентах и какие из них шаблону нужны.
        ex_functions = [reagent_functions(ex.get("spectators", [])) for ex in self.examples]
        n_ex = max(1, len(ex_functions))
        self.function_share = {f: sum(f in fs for fs in ex_functions) / n_ex for f in FUNCTIONS}
        # Нужные функции: по прецедентам (>50%) и по самому превращению (R15).
        self.intrinsic_functions = intrinsic_functions(self.forward)
        self.required_functions = {
            f for f, share in self.function_share.items() if share > REQUIRED_SHARE
        } | self.intrinsic_functions
        self.ex_functions = ex_functions
        # Катализатор прецедентов — для выдачи (условие процесса, не сырьё).
        # В корпусе он обычно в агентах; считается, только если есть хотя бы
        # в половине прецедентов (разовый металл — не катализатор превращения).
        cat = collections.Counter(
            f
            for ex in self.examples
            for f in set(ex.get("spectators", [])) | set(ex.get("agents", "").split("."))
            if f and METAL_CATALYST in reagent_functions([f])
        )
        n_ex = max(1, len(self.examples))
        self.catalysts = [f for f, c in cat.most_common(2) if c / n_ex >= REQUIRED_SHARE]
        # Классы атомов центра в прецедентах и самая сильная конкурирующая
        # группа, при которой прецедент всё же шёл по этому центру (R17).
        self.center_classes: dict[int, set] = collections.defaultdict(set)
        self.max_competitor: dict[int, int] = collections.defaultdict(int)
        for ex in self.examples:
            for k, (cls, comp) in (ex.get("center") or {}).items():
                self.center_classes[int(k)].add(cls)
                self.max_competitor[int(k)] = max(self.max_competitor[int(k)], comp)
        # Прецедент для выдачи: первый пример со ссылкой на источник, иначе первый.
        self.precedent = next(
            (ex for ex in self.examples if ex.get("source") not in ("", None, NO_SOURCE)),
            self.examples[0],
        )
        # Одномолекулярный шаблон, у которого во всех прецедентах был
        # сопутствующий реагент (нейтрализация, омыление, гидролиз): без этого
        # реагента превращение не идёт, и применять его «в пустоту» нельзя.
        # Шаблон, которому нужна функция реагента (окислитель, основание...),
        # без такого реагента тоже не применяется.
        self.needs_reagent = bool(self.required_functions - CONDITIONS) or (
            len(self.slots) == 1 and bool(self.spectator_sets) and all(self.spectator_sets)
        )
        return self

    def can_be_reagent(self, smiles: str, mol) -> bool:
        """
        Может ли вещество быть сопутствующим реагентом этого шаблона (решение
        0009). Нужны все пять условий:
        1. у вещества есть функция реагента (основание, окислитель, восстановитель);
        2. оно само не подходит в слот шаблона — иначе это конкурирующий
           реагент или избыток (ЭО рядом с ПО в «спирт + эпоксид»);
        3. оно было сопутствующим хотя бы в одном прецеденте;
        4. оно закрывает все функции, нужные шаблону (NaOH не заменит гидрид);
        5. его функция встречается хотя бы в MIN_FUNCTION_SHARE прецедентов.
        """
        functions = reagent_functions(smiles.split("."))
        if not functions:
            return False
        if any(mol.HasSubstructMatch(s) for s in self.slots):
            return False
        fs = frag_set(smiles)
        if not fs or not any(fs <= sp for sp in self.spectator_sets):
            return False
        if not covers(self.required_functions, functions):
            return False
        # окислитель/восстановитель — только в окислительно-восстановительном
        # превращении, либо если прецеденты устойчиво (>=50%) его используют
        redox = functions & {OXIDANT, REDUCTANT}
        if (
            redox
            and not redox & self.intrinsic_functions
            and max(self.function_share[f] for f in redox) < REQUIRED_SHARE
        ):
            return False
        return max(self.function_share[f] for f in functions) >= MIN_FUNCTION_SHARE

    def selective(self, provenance: dict, mols) -> bool:
        """
        Хемоселективность (R17): реагирующий атом кандидата того же класса,
        что в прецедентах (кислота — не спирт, амид — не амин, неактивированный
        арилгалогенид — не активированный), и в молекуле нет более сильного
        нуклеофила, чем реагирующий, если прецеденты такого не показывали.
        """
        for k, classes in self.center_classes.items():
            if k not in provenance:
                continue
            ri, ai = provenance[k]
            cls = atom_class(mols[ri], ai)
            if cls not in classes:
                return False
            rank = NUC_RANK.get(cls)
            if rank is not None:
                comp = competitor_rank(mols[ri], ai)
                if comp > rank and comp > self.max_competitor[k]:
                    return False
        return True

    def roundtrip(self, reactants: list[str], product: str) -> bool:
        """
        Обратная проверка: ретро-шаблон, применённый к полученному продукту,
        должен вернуть ровно те реагенты, из которых продукт построен.
        Отсекает ложные срабатывания слишком общих шаблонов (например,
        «этерификация», которая раскрывает эпоксид, выбрасывая кислород).
        """
        if self.retro_rxn is None:
            return False
        from rdchiral.main import rdchiralReactants, rdchiralRun

        try:
            outcomes = rdchiralRun(self.retro_rxn, rdchiralReactants(product))
        except Exception:
            return False
        target = frag_set(".".join(reactants))
        return any(frag_set(o) == target for o in outcomes)


def load_templates(path, min_count, min_selfcheck) -> list[Template]:
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            # Курируемые шаблоны (ручная библиотека) не режем по частоте:
            # каждая реакция там встречается один раз, но проверена человеком.
            if r["selfcheck_rate"] < min_selfcheck:
                continue
            if r["count"] < min_count and not r.get("trusted"):
                continue
            out.append(
                Template(
                    r["id"],
                    r["forward"],
                    r["retro"],
                    r["count"],
                    r["selfcheck_rate"],
                    r["examples"],
                    r.get("type_key") or type_key(r["forward"]),
                ).build()
            )
    return out


def load_inputs(path) -> list[dict]:
    out = []
    for r in read_csv(path):
        smi = r["smiles"].strip()
        if smi.lower() in ("", "none", "resource"):
            out.append({"name": r["name"], "smiles": None, "mol": None, "canon": None})
            continue
        c = canon(smi)
        out.append(
            {
                "name": r["name"],
                "smiles": smi,
                "canon": c,
                "mol": Chem.MolFromSmiles(c) if c else None,
            }
        )
    return out


def cmd_apply(args):
    t_start = time.time()
    templates = load_templates(args.templates, args.min_count, args.min_selfcheck)
    inputs = load_inputs(args.inputs)
    # Всегда доступные ресурсы (вода, воздух; решение 0011): партнёры и
    # реагенты наравне с набором, но сами стадии не запускают.
    utilities = (
        [u for u in load_inputs(args.utilities) if u["mol"] is not None] if args.utilities else []
    )
    utility_canon = {u["canon"] for u in utilities}
    user_frags = {frag_set(i["canon"]) for i in inputs + utilities if i["mol"] is not None}

    # Пул известных веществ: исходный набор пользователя + продукты прошлых
    # стадий. В режиме --internal-only партнёры берутся только отсюда.
    pool = {i["canon"]: i["mol"] for i in inputs + utilities if i["mol"] is not None}
    pool_level = {c: 0 for c in pool}
    pool_frags = {frag_set(c) for c in pool}
    # Каким шаблоном получено вещество пула: (шаблон, повтор, серия). Шаблон
    # наращивает свой продукт не больше --max-repeat раз (олигомеры, решение 0011).
    made_template: dict[str, tuple[str, int, str]] = {}

    # Словарь внешних партнёров: реагенты из прецедентов корпуса с частотами,
    # опционально — пересечение с покупаемыми (база экономического агента).
    vocab_freq = collections.Counter()
    for t in templates:
        for ex in t.examples:
            for s in ex["reactants"]:
                vocab_freq[s] += 1
    allowed = None
    if args.purchasable:
        allowed = {canon(r["smiles"]) for r in read_csv(args.purchasable)} - {None}
        vocab_freq = collections.Counter({s: c for s, c in vocab_freq.items() if s in allowed})
        for s in allowed:  # покупаемое, но без прецедентов — тоже партнёр
            vocab_freq.setdefault(s, 0)
    vocab_mols = {s: Chem.MolFromSmiles(s) for s in vocab_freq}

    slot_cache: dict[str, list[str]] = {}

    def vocab_matches(slot, key):
        if key not in slot_cache:
            hits = [s for s, m in vocab_mols.items() if m is not None and m.HasSubstructMatch(slot)]
            slot_cache[key] = sorted(hits, key=lambda s: -vocab_freq[s])
        return slot_cache[key]

    def partners_for(t: Template, j: int) -> list[tuple[str, str]]:
        """Кандидаты в слот j по приоритету: пул (свой набор + промежуточные)
        -> прецеденты -> словарь. В --internal-only только пул."""
        slot, seen, out = t.slots[j], set(), []

        def add(s, src):
            if s not in seen:
                seen.add(s)
                out.append((s, src))

        for s, m in pool.items():
            if m.HasSubstructMatch(slot):
                src = "utility" if s in utility_canon else "user"
                add(s, src if pool_level[s] == 0 else "intermediate")
        if not args.internal_only:
            for ex in t.examples:
                for s in ex["reactants"]:
                    if allowed is not None and s not in allowed:
                        continue
                    m = vocab_mols.get(s) or Chem.MolFromSmiles(s)
                    if m is not None and m.HasSubstructMatch(slot):
                        add(s, "precedent")
            for s in vocab_matches(slot, t.slot_keys[j]):
                add(s, "vocab")
        return out[: args.max_partners]

    def reagent_for(t: Template) -> str | None:
        """Реагент для шаблона, которому он обязателен: сначала ищем в пуле,
        вне --internal-only допускаем реагент из прецедента."""
        for s, m in pool.items():
            if t.can_be_reagent(s, m):
                return s
        if args.internal_only:
            return None
        # реагент из прецедента, закрывающего нужные функции
        ex = next(
            (
                ex
                for ex, fs in zip(t.examples, t.ex_functions, strict=True)
                if covers(t.required_functions, fs)
            ),
            None,
        )
        if ex is None and t.intrinsic_functions:
            return None  # превращению нужен реагент, которого нет ни в наборе, ни в прецедентах
        return canon(".".join((ex or t.examples[0])["spectators"]))

    def unit(inp, t: Template, level: int):
        """
        Все кандидаты одной пары «вещество × шаблон» в порядке перебора.
        Общего состояния не меняет (дедупликация и новые вещества — в merge),
        поэтому пары считаются в разных процессах. Возвращает (счётчики,
        события); событие — (ключ реакции, запись или None, отказ или None).
        """
        X = inp["mol"]
        name = inp["name"]
        counts = collections.Counter()
        events = []
        local_seen = set()

        def emit(combo, role, reagent=""):
            """combo: [(smiles, source)], где source='input' — само соединение;
            reagent — сопутствующий реагент (в продукт не входит)."""
            smiles_combo = [s for s, _ in combo]
            own = [made_template[s] for s in smiles_combo if made_template.get(s, ("",))[0] == t.id]
            repeat, series = 1, None
            if own:
                repeat, series = max(k for _, k, _ in own) + 1, own[0][2]
                if repeat > args.max_repeat:
                    counts["skipped_repeat"] += 1  # C12E3 -> C12E4 -> ...
                    return
            mols = [Chem.MolFromSmiles(s) for s in smiles_combo]
            for prod, provenance in run_products_traced(t.rxn, mols, keep_new_stereo=False).items():
                if prod in smiles_combo:
                    continue
                if not t.selective(provenance, mols):
                    counts["rejected_selectivity"] += 1
                    continue
                core = ".".join(sorted(smiles_combo)) + ">" + reagent + ">" + prod
                if core in local_seen:
                    continue
                local_seen.add(core)
                rxn_smiles = ".".join(smiles_combo) + ">" + reagent + ">" + prod
                if not args.no_roundtrip and not t.roundtrip(smiles_combo, prod):
                    rej = {"input_name": name, "template_id": t.id, "reaction_smiles": rxn_smiles}
                    events.append((core, None, rej))
                    continue
                # Само входное вещество тоже участник: на стадиях ≥2 это
                # промежуточный продукт, а не вещество исходного набора.
                participants = smiles_combo + ([reagent] if reagent else [])
                if all(frag_set(s) in user_frags for s in participants):
                    source = "internal"  # все участники — исходный набор
                elif all(frag_set(s) in pool_frags for s in participants):
                    source = "internal+intermediate"
                else:
                    source = (
                        ",".join(sorted({src for _, src in combo if src not in ("input",)}))
                        or "none"
                    )
                ex = t.precedent
                rec = {
                    "level": level,
                    "input_name": name,
                    "input_smiles": inp["canon"],
                    "role": role,
                    "reaction_key": core,
                    # реакция без средней части: варианты с разными
                    # сопутствующими реагентами — одна реакция (R13)
                    "reaction_core": ".".join(sorted(smiles_combo)) + ">>" + prod,
                    # серия олигомеров: C12E1..C12E3 — одна реакция для 2.1
                    "series": series or ".".join(sorted(smiles_combo)) + ">>" + prod,
                    "repeat": repeat,
                    "reagent": reagent,
                    "type_key": t.type_key,
                    "template_id": t.id,
                    "template_count": t.count,
                    "partners": ".".join(s for s, src in combo if src != "input"),
                    "partner_source": source,
                    "product": prod,
                    "reaction_smiles": rxn_smiles,
                    "needs_counterion": Chem.GetFormalCharge(Chem.MolFromSmiles(prod)) != 0,
                    "precedent_id": ex["id"],
                    "precedent_rxn": ".".join(ex["reactants"]) + ">>" + ex["product"],
                    "precedent_source": ex.get("source") or NO_SOURCE,
                    "catalyst": ".".join(t.catalysts),
                }
                events.append((core, rec, None))

        # (а) вещество — реагент: его атомы входят в продукт
        hit_slots = [i for i, s in enumerate(t.slots) if X.HasSubstructMatch(s)]
        reagent = ""
        if hit_slots and t.needs_reagent:
            reagent = reagent_for(t)
            if reagent is None:  # нужного реагента в наборе нет
                hit_slots = []
        if hit_slots:
            counts["templates_as_reactant"] += 1
        for i in hit_slots:
            lists = [
                partners_for(t, j) if j != i else [(inp["canon"], "input")]
                for j in range(len(t.slots))
            ]
            for combo in itertools.islice(itertools.product(*lists), args.max_combos):
                emit(list(combo), "reactant", reagent)
        # (б) вещество — сопутствующий реагент (основание, окислитель):
        # в прецедентах было, в продукт не вошло и выполняет нужную функцию
        if t.can_be_reagent(inp["canon"], X):
            counts["templates_as_reagent"] += 1
            lists = [partners_for(t, j) for j in range(len(t.slots))]
            for combo in itertools.islice(itertools.product(*lists), args.max_combos):
                emit(list(combo), "reagent", inp["canon"])
        return counts, events

    results, rejected = [], []
    stats = collections.defaultdict(collections.Counter)
    seen = set()
    made_by: dict[str, dict] = {}  # продукт -> первая реакция, его давшая

    def merge(inp, counts, events, new_products):
        """Последовательная сборка в порядке перебора: первая пара, давшая
        реакцию, решает её судьбу — как при обработке в одном процессе."""
        name = inp["name"]
        stats[name].update(counts)
        for core, rec, rej in events:
            if (name, core) in seen:
                continue
            seen.add((name, core))
            if rec is None:
                stats[name]["rejected_roundtrip"] += 1
                rejected.append(rej)
                continue
            results.append(rec)
            prod = rec["product"]
            if prod not in pool and prod not in new_products:
                new_products[prod] = rec
                made_by.setdefault(prod, rec)

    # Стадия 1 — исходные соединения; стадии 2..depth — только новые продукты
    # предыдущей стадии (реакции «старое + старое» уже перебраны).
    jobs = n_jobs(args.jobs)
    print(f"[apply] шаблонов: {len(templates)}; процессов: {jobs}", file=sys.stderr)
    frontier = [i for i in inputs]
    level_summary = []
    for level in range(1, args.depth + 1):
        new_products: dict[str, dict] = {}
        for inp in frontier:
            if inp["mol"] is None:
                stats[inp["name"]]["not_a_molecule"] = 1
        active = [inp for inp in frontier if inp["mol"] is not None]
        nt = len(templates)
        # Процессы создаются здесь, после обновления пула: fork видит его текущим.
        _WORKER_STATE.update(unit=unit, active=active, templates=templates, level=level)
        for k, (counts, events) in enumerate(parallel_map(_apply_worker, len(active) * nt, jobs)):
            merge(active[k // nt], counts, events, new_products)
        _WORKER_STATE.clear()
        level_summary.append(
            (
                level,
                len({r["reaction_core"] for r in results if r["level"] == level}),
                len(new_products),
            )
        )
        print(
            f"[apply] стадия {level}: реакций всего {len(results)}, новых веществ "
            f"{len(new_products)} ({time.time() - t_start:.0f} с)",
            file=sys.stderr,
        )
        for prod in new_products:
            m = Chem.MolFromSmiles(prod)
            if m is not None and m.GetNumHeavyAtoms() <= args.max_heavy_atoms:
                pool[prod], pool_level[prod] = m, level
                rec = new_products[prod]
                made_template[prod] = (rec["template_id"], rec["repeat"], rec["series"])
                pool_frags.add(frag_set(prod))
        frontier = [
            {"name": f"[стадия {level}] {p}", "canon": p, "mol": pool[p]}
            for p in new_products
            if p in pool
        ]
        if not frontier:
            break

    for inp in inputs:
        mine = [r for r in results if r["input_name"] == inp["name"]]
        stats[inp["name"]]["reactions"] = len({r["reaction_core"] for r in mine})
        stats[inp["name"]]["internal"] = len(
            {r["reaction_core"] for r in mine if r["partner_source"] == "internal"}
        )
        stats[inp["name"]]["types"] = len({r["type_key"] for r in mine})
        stats[inp["name"]]["internal_types"] = len(
            {r["type_key"] for r in mine if r["partner_source"] == "internal"}
        )

    fields = [
        "level",
        "input_name",
        "input_smiles",
        "role",
        "reaction_key",
        "reaction_core",
        "series",
        "repeat",
        "reagent",
        "type_key",
        "template_id",
        "template_count",
        "partners",
        "partner_source",
        "product",
        "reaction_smiles",
        "needs_counterion",
        "precedent_id",
        "precedent_rxn",
        "precedent_source",
        "catalyst",
    ]
    write_csv(args.out, results, fields)
    elapsed = time.time() - t_start

    # --- отчёт о покрытии ---
    uniq = {r["reaction_core"] for r in results}
    uniq_int = {r["reaction_core"] for r in results if r["partner_source"] == "internal"}
    uniq_net = {
        r["reaction_core"]
        for r in results
        if r["partner_source"] in ("internal", "internal+intermediate")
    }
    L = [
        "# Покрытие входных соединений шаблонами\n",
        f"- Шаблонов после фильтра (частота ≥ {args.min_count} или курируемый, "
        f"самопроверка ≥ {args.min_selfcheck:.0%}): **{len(templates)}**",
        "- Партнёры: "
        + (
            "**только вещества набора и промежуточные продукты** (--internal-only)"
            if args.internal_only
            else f"набор + прецеденты + словарь из **{len(vocab_freq)}** веществ"
            + (" (только покупаемые)" if allowed is not None else "")
        ),
        f"- Стадий: {args.depth}",
        "- Обратная проверка (round-trip): "
        + ("выключена" if args.no_roundtrip else f"отклонено **{len(rejected)}** кандидатов"),
        f"- Время: {elapsed:.1f} с\n",
        "## По всему набору (так сформулирован критерий ТЗ 2.1)\n",
        f"- Уникальных реакций всего: **{len(uniq)}**; типов превращений: "
        f"**{len({r['type_key'] for r in results})}**",
        f"- Вариантов с разными сопутствующими реагентами (считаются как одна реакция): "
        f"{len({r['reaction_key'] for r in results})}",
        f"- **Для критерия 2.1** (серия олигомеров одного шаблона — одна реакция, "
        f"решение 0011): **{len({r['series'] for r in results if r['partner_source'] in ('internal', 'internal+intermediate')})}**",
        f"- Все участники — исходный набор: **{len(uniq_int)}**",
        f"- Все участники — исходный набор и полученные из него промежуточные: **{len(uniq_net)}**\n",
    ]
    if args.depth > 1:
        L.append("| Стадия | Реакций | Новых веществ |")
        L.append("|---|---|---|")
        for lv, n_r, n_p in level_summary:
            L.append(f"| {lv} | {n_r} | {n_p} |")
        L.append("")
    L += [
        "## По исходным соединениям (стадия 1)\n",
        "*Реакций* — уникальные реакции; *типов* — уникальные превращения по скелетному "
        "ключу (честная мера разнообразия); *внутри набора* — все партнёры из входного "
        "набора; *как реагент / как сопутствующее* — число шаблонов, где вещество входит "
        "в продукт или только сопровождает реакцию.\n",
        "| Соединение | Шаблонов: как реагент / как сопутствующее | Реакций | Типов | "
        "Внутри набора: реакций (типов) | Отклонено round-trip |",
        "|---|---|---|---|---|---|",
    ]
    for inp in inputs:
        c = stats[inp["name"]]
        if c.get("not_a_molecule"):
            L.append(f"| {inp['name']} | ресурс, не молекула | — | — | — | — |")
            continue
        L.append(
            f"| {inp['name']} | {c['templates_as_reactant']} / {c['templates_as_reagent']} | "
            f"{c['reactions']} | {c['types']} | {c['internal']} ({c['internal_types']}) | "
            f"{c['rejected_roundtrip']} |"
        )

    if args.depth > 1:
        names = {i["canon"]: i["name"] for i in inputs if i["canon"]}

        def route(prod, depth=0):
            rec = made_by.get(prod)
            if rec is None or depth > args.depth:
                return []
            reac_part, reag_part, _ = rec["reaction_smiles"].split(">")
            steps = []
            for s in reac_part.split("."):
                if s != prod and s in made_by:
                    steps += route(s, depth + 1)
            return steps + [rec]

        def dedup_steps(steps):
            out, seen_s = [], set()
            for r in steps:
                if r["reaction_key"] not in seen_s:
                    seen_s.add(r["reaction_key"])
                    out.append(r)
            return out

        user_canon = {i["canon"] for i in inputs if i["canon"]}

        def used_inputs(steps):
            used = set()
            for r in steps:
                reac_part, reag_part, _ = r["reaction_smiles"].split(">")
                for s in reac_part.split(".") + ([reag_part] if reag_part else []):
                    if s in user_canon:
                        used.add(s)
            return used

        deep = [p for p, r in made_by.items() if r["level"] >= 2]
        routes = {p: dedup_steps(route(p)) for p in deep}
        deep.sort(key=lambda p: (-len(used_inputs(routes[p])), len(routes[p]), p))
        if deep:
            L.append("\n## Многостадийные маршруты\n")
            L.append(
                "Отсортированы по числу задействованных веществ исходного набора, "
                "затем по длине маршрута.\n"
            )
            for p in deep[: args.max_routes]:
                steps = routes[p]
                used = sorted(names.get(s, s) for s in used_inputs(steps))
                L.append(f"**{p}** — {len(steps)} стад.; из набора: {', '.join(used)}\n")
                for k, r in enumerate(steps, 1):
                    L.append(
                        f"{k}. `{r['reaction_smiles']}` — шаблон {r['template_id']}, "
                        f"прецедент {r['precedent_id']}"
                    )
                L.append("")

    if rejected:
        L.append("\n## Примеры отклонённых обратной проверкой\n")
        for r in rejected[:10]:
            L.append(f"- {r['input_name']}: `{r['reaction_smiles']}` (шаблон {r['template_id']})")

    L.append("\n## Типы превращений по исходным соединениям\n")
    for inp in inputs:
        rows = [r for r in results if r["input_name"] == inp["name"]]
        if not rows:
            continue
        L.append(f"### {inp['name']}\n")
        by_t = collections.defaultdict(list)
        for r in rows:
            by_t[(r["role"], r["type_key"])].append(r)
        for (role, tk), rs in sorted(by_t.items(), key=lambda kv: -len(kv[1])):
            ex = rs[0]
            tids = sorted({r["template_id"] for r in rs})
            tag = "" if role == "reactant" else " *(как сопутствующий реагент)*"
            L.append(
                f"- `{tk}`{tag} — {len(rs)} реакц. (шаблоны {', '.join(tids)}); "
                f"пример `{ex['reaction_smiles']}`; прецедент {ex['precedent_id']}"
            )
        L.append("")
    with open(args.report, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    print(
        f"[apply] реакций: {len(results)} (отклонено round-trip: {len(rejected)}) "
        f"-> {args.out}; отчёт -> {args.report}",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("map", help="атомный маппинг корпуса")
    m.add_argument("--corpus", required=True)
    m.add_argument("--smiles-col", default="rxn_smiles")
    m.add_argument("--id-col", default="id")
    m.add_argument(
        "--source-col", default="source", help="колонка источника (патент); переносится как есть"
    )
    m.add_argument("--out", required=True)
    m.add_argument("--batch", type=int, default=32)
    m.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda"],
        default="auto",
        help="auto — GPU, при нехватке памяти CPU",
    )
    m.set_defaults(func=cmd_map)

    e = sub.add_parser("extract", help="извлечение шаблонов")
    e.add_argument(
        "--mapped",
        required=True,
        nargs="+",
        help="один или несколько файлов после map (корпус + ручная библиотека)",
    )
    e.add_argument("--out", required=True)
    e.add_argument("--timeout", type=int, default=10, help="сек. на одну реакцию")
    e.add_argument("--max-reactants", type=int, default=3)
    e.add_argument("--max-examples", type=int, default=20)
    e.add_argument("--no-intra", action="store_true", help="отбросить внутримолекулярные")
    e.add_argument("--jobs", type=int, default=0, help="процессов; 0 — все ядра")
    e.add_argument(
        "--trusted",
        nargs="*",
        default=[],
        help="файлы из --mapped с курируемыми реакциями: их шаблоны не режутся по --min-count",
    )
    e.set_defaults(func=cmd_extract)

    a = sub.add_parser("apply", help="применение шаблонов к входным соединениям")
    a.add_argument("--templates", required=True)
    a.add_argument("--inputs", required=True, help="CSV с колонками name,smiles")
    a.add_argument("--purchasable", help="CSV с колонкой smiles: ограничить партнёров покупаемыми")
    a.add_argument(
        "--utilities",
        help="CSV name,smiles: всегда доступные ресурсы (вода, воздух), data/inputs/utilities.csv",
    )
    a.add_argument("--out", required=True)
    a.add_argument("--report", required=True)
    a.add_argument("--min-count", type=int, default=2)
    a.add_argument("--min-selfcheck", type=float, default=0.5)
    a.add_argument("--max-partners", type=int, default=25, help="кандидатов на слот")
    a.add_argument("--max-combos", type=int, default=200, help="комбинаций на шаблон и слот")
    a.add_argument("--no-roundtrip", action="store_true", help="выключить обратную проверку")
    a.add_argument("--depth", type=int, default=1, help="число стадий (продукты -> новые входы)")
    a.add_argument(
        "--internal-only",
        action="store_true",
        help="партнёры только из набора пользователя и промежуточных продуктов",
    )
    a.add_argument(
        "--max-heavy-atoms",
        type=int,
        default=60,
        help="не пускать на следующую стадию продукты крупнее этого",
    )
    a.add_argument("--max-routes", type=int, default=15, help="маршрутов в отчёте")
    a.add_argument(
        "--max-repeat",
        type=int,
        default=3,
        help="сколько раз подряд шаблон наращивает свой продукт (олигомеры)",
    )
    a.add_argument("--jobs", type=int, default=0, help="процессов; 0 — все ядра")
    a.set_defaults(func=cmd_apply)

    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
