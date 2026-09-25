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
import re
import signal
import sys
import time
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

MAP_NUM_RE = re.compile(r":\d+\]")


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


def run_products(rxn, reactant_mols, max_products=200) -> set[str]:
    """Прямой прогон шаблона; возвращает канонические SMILES продуктов."""
    out = set()
    try:
        product_sets = rxn.RunReactants(tuple(reactant_mols), max_products)
    except Exception:
        return out
    for ps in product_sets:
        for p in ps:
            try:
                Chem.SanitizeMol(p)
                smi = Chem.MolToSmiles(p)
            except Exception:
                continue
            c = canon(smi)
            if c:
                out.add(c)
    return out


# ---------------------------------------------------------------------------
# Шаг 1. Атомный маппинг
# ---------------------------------------------------------------------------


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
        from rxnmapper import RXNMapper

        mapper = RXNMapper()
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

    write_csv(args.out, out_rows, ["id", "rxn_original", "agents", "mapped_rxn", "map_confidence"])
    print(f"[map] записано {len(out_rows)} -> {args.out}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Шаг 2. Извлечение шаблонов
# ---------------------------------------------------------------------------


def prepare_for_extraction(mapped_rxn: str):
    """
    Оставляем один основной продукт (самый тяжёлый фрагмент) — генератор
    нужен для главного продукта, побочные (вода, соли) достраивает балансировщик.
    Возвращает (реагирующие_mapped, основной_продукт_mapped,
                реагирующие_без_разметки, продукт_без_разметки, сопутствующие).
    Сопутствующие — реагенты, ни один атом которых не попал в основной продукт
    (основания, кислоты, окислители). Они не входят в шаблон, но по ним
    превращение привязывается к веществу: NaOH «ведёт» нейтрализацию, хотя
    его атомов в продукте нет.
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
    spectators = [canon(r) for r in reac_frags if not (mapnums(r) & main_maps)]
    return (
        ".".join(reacting),
        main,
        [canon(r) for r in reacting],
        canon(main),
        [s for s in spectators if s],
    )


def cmd_extract(args):
    from rdchiral.template_extractor import extract_from_reaction

    trusted_paths = set(args.trusted or [])
    rows = []
    for path in args.mapped:
        for row in read_csv(path):
            row["_trusted"] = path in trusted_paths
            rows.append(row)
    signal.signal(signal.SIGALRM, _alarm)

    templates: dict[str, dict] = {}
    stats = collections.Counter()

    for n, row in enumerate(rows, 1):
        prepared = prepare_for_extraction(row["mapped_rxn"]) if row.get("mapped_rxn") else None
        if prepared is None:
            stats["skip_unparsable"] += 1
            continue
        reac_m, prod_m, reac_plain, prod_plain, spectators = prepared
        if None in reac_plain or prod_plain is None:
            stats["skip_unparsable"] += 1
            continue

        signal.alarm(args.timeout)
        try:
            t = extract_from_reaction({"reactants": reac_m, "products": prod_m, "_id": row["id"]})
        except Timeout:
            stats["skip_timeout"] += 1
            continue
        except Exception:
            stats["skip_error"] += 1
            continue
        finally:
            signal.alarm(0)

        if not t or "reactants" not in t or not t.get("products"):
            stats["skip_no_template"] += 1
            continue
        if t.get("intra_only") and args.no_intra:
            stats["skip_intra"] += 1
            continue

        forward = f"{t['reactants']}>>{t['products']}"
        rxn = AllChem.ReactionFromSmarts(forward)
        n_slots = rxn.GetNumReactantTemplates()
        if n_slots > args.max_reactants:
            stats["skip_too_many_reactants"] += 1
            continue

        # Самопроверка: прямой шаблон на исходных реагентах должен дать
        # записанный основной продукт (перебираем порядок реагентов по слотам).
        mols = [Chem.MolFromSmiles(s) for s in reac_plain]
        ok = False
        for combo in itertools.permutations(mols, n_slots) if len(mols) >= n_slots else []:
            if prod_plain in run_products(rxn, combo):
                ok = True
                break

        rec = templates.get(forward)
        if rec is None:
            rec = templates[forward] = {
                "forward": forward,
                "retro": t["reaction_smarts"],
                "n_reactants": n_slots,
                "necessary_reagent": t.get("necessary_reagent", ""),
                "count": 0,
                "selfcheck_pass": 0,
                "examples": [],
                "trusted": False,
            }
        rec["count"] += 1
        rec["selfcheck_pass"] += int(ok)
        rec["trusted"] = rec["trusted"] or row["_trusted"]
        if len(rec["examples"]) < args.max_examples:
            rec["examples"].append(
                {
                    "id": row["id"],
                    "reactants": reac_plain,
                    "product": prod_plain,
                    "spectators": spectators,
                    "agents": row.get("agents", ""),
                    "selfcheck": ok,
                }
            )
        stats["extracted"] += 1
        if n % 5000 == 0:
            print(f"[extract] {n}/{len(rows)}; шаблонов: {len(templates)}", file=sys.stderr)

    ordered = sorted(templates.values(), key=lambda r: -r["count"])
    with open(args.out, "w", encoding="utf-8") as f:
        for i, rec in enumerate(ordered, 1):
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
    if stats["extracted"]:
        print(
            f"[extract] самопроверка пройдена: {total_ok}/{stats['extracted']} "
            f"({total_ok / stats['extracted']:.0%})",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Шаг 3. Применение шаблонов к входным соединениям
# ---------------------------------------------------------------------------

# Противоионы не считаем «реагентами» при сравнении наборов веществ.
COUNTER_IONS = {"[Na+]", "[K+]", "[Li+]", "[Cl-]", "[Br-]", "[I-]", "[H+]"}


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
        # Одномолекулярный шаблон, у которого во всех прецедентах был
        # сопутствующий реагент (нейтрализация, омыление, гидролиз): без этого
        # реагента превращение не идёт, и применять его «в пустоту» нельзя.
        self.needs_reagent = (
            len(self.slots) == 1 and bool(self.spectator_sets) and all(self.spectator_sets)
        )
        return self

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
    user_frags = {frag_set(i["canon"]) for i in inputs if i["mol"] is not None}

    # Пул известных веществ: исходный набор пользователя + продукты прошлых
    # стадий. В режиме --internal-only партнёры берутся только отсюда.
    pool = {i["canon"]: i["mol"] for i in inputs if i["mol"] is not None}
    pool_level = {c: 0 for c in pool}
    pool_frags = {frag_set(c) for c in pool}

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
                add(s, "user" if pool_level[s] == 0 else "intermediate")
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

    results, rejected = [], []
    stats = collections.defaultdict(collections.Counter)
    seen = set()
    made_by: dict[str, dict] = {}  # продукт -> первая реакция, его давшая

    def emit(inp, t, combo, role, level, new_products, reagent=""):
        """combo: [(smiles, source)], где source='input' — само соединение;
        reagent — сопутствующий реагент (в продукт не входит)."""
        smiles_combo = [s for s, _ in combo]
        mols = [Chem.MolFromSmiles(s) for s in smiles_combo]
        name = inp["name"]
        for prod in run_products(t.rxn, mols):
            if prod in smiles_combo:
                continue
            core = ".".join(sorted(smiles_combo)) + ">" + reagent + ">" + prod
            if (name, core) in seen:
                continue
            seen.add((name, core))
            rxn_smiles = ".".join(smiles_combo) + ">" + reagent + ">" + prod
            if not args.no_roundtrip and not t.roundtrip(smiles_combo, prod):
                stats[name]["rejected_roundtrip"] += 1
                rejected.append(
                    {"input_name": name, "template_id": t.id, "reaction_smiles": rxn_smiles}
                )
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
                    ",".join(sorted({src for _, src in combo if src not in ("input",)})) or "none"
                )
            ex = t.examples[0]
            rec = {
                "level": level,
                "input_name": name,
                "input_smiles": inp["canon"],
                "role": role,
                "reaction_key": core,
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
            }
            results.append(rec)
            stats[name]["reactions"] += 1
            stats[name]["internal"] += int(source == "internal")
            if prod not in pool and prod not in new_products:
                new_products[prod] = rec
                made_by.setdefault(prod, rec)

    def reagent_for(t: Template) -> str | None:
        """Реагент для шаблона, которому он обязателен: сначала ищем в пуле,
        вне --internal-only допускаем реагент из прецедента."""
        for s in pool:
            fs = frag_set(s)
            if fs and any(fs <= sp for sp in t.spectator_sets):
                return s
        if args.internal_only:
            return None
        return canon(".".join(t.examples[0]["spectators"]))

    def expand(inp, level, new_products):
        X, FX = inp["mol"], frag_set(inp["canon"])
        name = inp["name"]
        for t in templates:
            # (а) вещество — реагент: его атомы входят в продукт
            hit_slots = [i for i, s in enumerate(t.slots) if X.HasSubstructMatch(s)]
            reagent = ""
            if hit_slots and t.needs_reagent:
                reagent = reagent_for(t)
                if reagent is None:  # нужного реагента в наборе нет
                    hit_slots = []
            if hit_slots:
                stats[name]["templates_as_reactant"] += 1
            for i in hit_slots:
                lists = [
                    partners_for(t, j) if j != i else [(inp["canon"], "input")]
                    for j in range(len(t.slots))
                ]
                for combo in itertools.islice(itertools.product(*lists), args.max_combos):
                    emit(inp, t, list(combo), "reactant", level, new_products, reagent)
            # (б) вещество — сопутствующий реагент (основание, окислитель):
            # в прецедентах было, но в продукт не вошло
            if FX and any(FX <= sp for sp in t.spectator_sets):
                stats[name]["templates_as_reagent"] += 1
                lists = [partners_for(t, j) for j in range(len(t.slots))]
                for combo in itertools.islice(itertools.product(*lists), args.max_combos):
                    emit(inp, t, list(combo), "reagent", level, new_products, inp["canon"])

    # Стадия 1 — исходные соединения; стадии 2..depth — только новые продукты
    # предыдущей стадии (реакции «старое + старое» уже перебраны).
    frontier = [i for i in inputs]
    level_summary = []
    for level in range(1, args.depth + 1):
        new_products: dict[str, dict] = {}
        for inp in frontier:
            if inp["mol"] is None:
                stats[inp["name"]]["not_a_molecule"] = 1
                continue
            expand(inp, level, new_products)
        level_summary.append(
            (
                level,
                len({r["reaction_key"] for r in results if r["level"] == level}),
                len(new_products),
            )
        )
        for prod in new_products:
            m = Chem.MolFromSmiles(prod)
            if m is not None and m.GetNumHeavyAtoms() <= args.max_heavy_atoms:
                pool[prod], pool_level[prod] = m, level
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
    ]
    write_csv(args.out, results, fields)
    elapsed = time.time() - t_start

    # --- отчёт о покрытии ---
    uniq = {r["reaction_key"] for r in results}
    uniq_int = {r["reaction_key"] for r in results if r["partner_source"] == "internal"}
    uniq_net = {
        r["reaction_key"]
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
    m.add_argument("--out", required=True)
    m.add_argument("--batch", type=int, default=32)
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
    a.set_defaults(func=cmd_apply)

    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
