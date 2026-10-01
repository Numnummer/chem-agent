"""
vectordb.py — векторная база реакций (ТЗ 2.2, решение 0016).

Вектор реакции — разностный отпечаток RDKit (что изменилось: «продукт минус
реагенты») и структурный отпечаток (кто участвует), L2-нормированные и
склеенные. Текстовые эмбеддинги не используются: вектор отражает структурное
изменение, а не формулировку (architecture.md).

Назначение: (1) уникальность — у каждой реакции свой вектор (критерий 2.2);
(2) поиск похожих известных процессов (корпус, ручная библиотека);
(3) дедупликация почти-дубликатов.

Хранение: <путь>.npz (векторы) + <путь>.jsonl (метаданные).
"""

from __future__ import annotations

import json

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

FP_SIZE = 2048
ELEMENTS = {
    e: i for i, e in enumerate(["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I", "Na", "*"])
}
COMP_WEIGHT = 0.1


def reaction_vector(rxn_smiles: str) -> np.ndarray | None:
    """Вектор реакции 'A.B>реагент>P' (средняя часть не учитывается)."""
    parts = rxn_smiles.split(">")
    if len(parts) != 3:
        return None
    try:
        rxn = rdChemReactions.ReactionFromSmarts(f"{parts[0]}>>{parts[2]}", useSmiles=True)
    except Exception:
        return None
    params = rdChemReactions.ReactionFingerprintParams()
    params.fpSize = FP_SIZE
    # Морган, а не пары атомов (по умолчанию): инварианты атома включают заряд,
    # иначе кислота и её анион дают один вектор (критерий 2.2 — уникальность)
    params.fpType = rdChemReactions.FingerprintType.MorganFP
    diff = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn, params)
    dv = np.zeros(FP_SIZE, dtype=np.float32)
    for k, v in diff.GetNonzeroElements().items():
        dv[k % FP_SIZE] += v
    struct = rdChemReactions.CreateStructuralFingerprintForReaction(rxn, params)
    sv = np.zeros(struct.GetNumBits(), dtype=np.float32)
    sv[list(struct.GetOnBits())] = 1.0
    out = []
    for v in (dv, sv):
        n = np.linalg.norm(v)
        out.append(v / n if n else v)
    # Состав (атомы элементов в реагентах и продукте) с малым весом: бинарные
    # отпечатки насыщаются на гомологах (C12E2 и C12E3 — один набор
    # фрагментов), а разные брутто-формулы должны давать разные векторы.
    comp = np.zeros(2 * len(ELEMENTS), dtype=np.float32)
    for side, smiles in enumerate((parts[0], parts[2])):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        for a in Chem.AddHs(mol).GetAtoms():
            i = ELEMENTS.get(a.GetSymbol(), len(ELEMENTS) - 1)
            comp[side * len(ELEMENTS) + i] += 1
    out.append(COMP_WEIGHT * comp / max(1.0, float(np.linalg.norm(comp))))
    return np.concatenate(out) / np.sqrt(2)


class VectorDB:
    def __init__(self):
        self.keys: list[str] = []
        self.meta: list[dict] = []
        self._vecs: list[np.ndarray] = []
        self._matrix: np.ndarray | None = None

    def __len__(self) -> int:
        return len(self.keys)

    def add(self, key: str, rxn_smiles: str, meta: dict | None = None) -> bool:
        if key in self._index():
            return False
        v = reaction_vector(rxn_smiles)
        if v is None:
            return False
        self.keys.append(key)
        self.meta.append({"rxn": rxn_smiles, **(meta or {})})
        self._vecs.append(v)
        self._matrix = None
        return True

    def _index(self) -> set[str]:
        return set(self.keys)

    @property
    def matrix(self) -> np.ndarray:
        if self._matrix is None:
            self._matrix = np.vstack(self._vecs) if self._vecs else np.zeros((0, 1))
        return self._matrix

    def search(self, rxn_smiles: str, k: int = 5, exclude: str | None = None):
        """[(сходство, ключ, метаданные)] по убыванию косинусной близости."""
        v = reaction_vector(rxn_smiles)
        if v is None or not len(self):
            return []
        sims = self.matrix @ v
        order = np.argsort(-sims)
        out = []
        for i in order:
            if self.keys[i] == exclude:
                continue
            out.append((float(sims[i]), self.keys[i], self.meta[i]))
            if len(out) == k:
                break
        return out

    def unique_vectors(self) -> int:
        """Число попарно различных векторов (критерий 2.2)."""
        if not len(self):
            return 0
        return len(np.unique(np.round(self.matrix, 6), axis=0))

    def save(self, path: str) -> None:
        np.savez_compressed(path + ".npz", vectors=self.matrix, keys=np.array(self.keys))
        with open(path + ".jsonl", "w", encoding="utf-8") as f:
            for key, m in zip(self.keys, self.meta, strict=True):
                f.write(json.dumps({"key": key, **m}, ensure_ascii=False) + "\n")

    @classmethod
    def load(cls, path: str) -> VectorDB:
        db = cls()
        data = np.load(path + ".npz")
        db.keys = [str(k) for k in data["keys"]]
        db._matrix = data["vectors"]
        db._vecs = list(db._matrix)
        with open(path + ".jsonl", encoding="utf-8") as f:
            db.meta = [{k: v for k, v in json.loads(line).items() if k != "key"} for line in f]
        return db
