"""Общие фикстуры: шаблоны строятся один раз на сессию из размеченных фикстур."""

import csv
from pathlib import Path

import pytest

from chem_agent.template_engine import main

FIXTURES = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).parent.parent
INPUTS = ROOT / "data" / "inputs" / "priority_set.csv"
PURCHASABLE = ROOT / "data" / "demo" / "purchasable.csv"


@pytest.fixture(scope="session")
def templates_path(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("tpl") / "templates.jsonl"
    manual = str(FIXTURES / "manual_mapped.csv")
    main(
        [
            "extract",
            "--mapped",
            str(FIXTURES / "corpus_mapped.csv"),
            manual,
            "--trusted",
            manual,
            "--out",
            str(out),
        ]
    )
    return out


@pytest.fixture(scope="session")
def run_apply(templates_path, tmp_path_factory):
    """Запускает apply и возвращает строки candidates.csv."""

    def _run(*extra: str, inputs: Path = INPUTS) -> list[dict]:
        d = tmp_path_factory.mktemp("apply")
        out, report = d / "candidates.csv", d / "coverage.md"
        main(
            [
                "apply",
                "--templates",
                str(templates_path),
                "--inputs",
                str(inputs),
                "--min-count",
                "1",
                "--out",
                str(out),
                "--report",
                str(report),
                *extra,
            ]
        )
        with open(out, encoding="utf-8") as f:
            return list(csv.DictReader(f))

    return _run


def write_inputs(path: Path, rows: list[tuple[str, str]]) -> Path:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "smiles"])
        w.writerows(rows)
    return path
