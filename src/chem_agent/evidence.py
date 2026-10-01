"""
evidence.py — доказательная база реакции (решение 0007, roadmap R11).

  - уровень I: реакция есть в корпусе или ручной библиотеке (те же
    реагирующие вещества среди реагентов записи и тот же продукт);
  - уровень II: аналог по шаблону (прецедент того же превращения);
  - условия прецедента из meta.csv корпуса: температура, растворитель,
    время, давление, атмосфера.

Индекс строится из corpus.csv (corpus_prep), meta.csv и
data/manual/inorganic.csv и кэшируется рядом с корпусом.
"""

from __future__ import annotations

import collections
import csv
import os
import pickle
import sys

from chem_agent.template_engine import canon, frag_set

CONDITION_FIELDS = {
    "temperature_c_min": "T мин, °C",
    "temperature_c_max": "T макс, °C",
    "pressure": "давление",
    "solvent": "растворитель",
    "atmosphere": "атмосфера",
    "time_minutes_stages": "время, мин",
}


class EvidenceIndex:
    def __init__(self):
        # продукт -> [(фрагменты реагентов записи, id, источник)]
        self.by_product: dict[str, list[tuple[frozenset, str, str]]] = collections.defaultdict(list)
        self.conditions: dict[str, dict] = {}

    def _add(self, rid: str, rxn: str, source: str) -> None:
        parts = rxn.split(">")
        if len(parts) != 3:
            return
        reac = frag_set(parts[0] + ("." + parts[1] if parts[1] else ""))
        for p in parts[2].split("."):
            c = canon(p)
            if c:
                self.by_product[c].append((reac, rid, source))

    @classmethod
    def build(cls, corpus_csv: str | None, meta_csv: str | None, manual_csv: str | None):
        idx = cls()
        csv.field_size_limit(sys.maxsize)
        for path, id_col, src_col in ((corpus_csv, "id", "source"), (manual_csv, "id", "source")):
            if path and os.path.exists(path):
                with open(path, newline="", encoding="utf-8") as f:
                    for r in csv.DictReader(f):
                        idx._add(r[id_col], r["rxn_smiles"], r.get(src_col, ""))
        if meta_csv and os.path.exists(meta_csv):
            with open(meta_csv, newline="", encoding="utf-8") as f:
                for r in csv.DictReader(f):
                    cond = {CONDITION_FIELDS[k]: r[k] for k in CONDITION_FIELDS if r.get(k)}
                    if cond:
                        idx.conditions[r["id"]] = cond
        return idx

    @classmethod
    def load(cls, corpus_csv: str | None, meta_csv: str | None, manual_csv: str | None):
        """Индекс из кэша рядом с корпусом (пересобирается, если корпус новее)."""
        if not corpus_csv:
            return cls.build(None, meta_csv, manual_csv)
        cache = corpus_csv + ".evidence.pkl"
        sources = [p for p in (corpus_csv, meta_csv, manual_csv) if p and os.path.exists(p)]
        if os.path.exists(cache) and all(
            os.path.getmtime(cache) >= os.path.getmtime(p) for p in sources
        ):
            with open(cache, "rb") as f:
                return pickle.load(f)
        idx = cls.build(corpus_csv, meta_csv, manual_csv)
        with open(cache, "wb") as f:
            pickle.dump(idx, f)
        return idx

    def level(self, reactants: list[str], product: str) -> tuple[str, str, str]:
        """('I', id, источник) — реакция есть в корпусе; ('II', '', '') — аналог."""
        need = frag_set(".".join(reactants))
        for reac, rid, source in self.by_product.get(canon(product) or product, []):
            if need <= reac:
                return "I", rid, source
        return "II", "", ""
