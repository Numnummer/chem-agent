"""Сквозной MVP: набор веществ -> отранжированные реакции с экономикой (решение 0016)."""

import itertools
import json

import pytest
from conftest import INPUTS

from chem_agent.pipeline import main


@pytest.fixture(scope="module")
def analysis(templates_path, tmp_path_factory):
    out = tmp_path_factory.mktemp("mvp")
    main(
        ["--inputs", str(INPUTS), "--templates", str(templates_path), "--out", str(out)]
        + ["--depth", "4", "--jobs", "1"]
    )
    return out, json.loads((out / "results.json").read_text(encoding="utf-8"))


def test_outputs_written(analysis):
    out, _ = analysis
    for name in ("results.json", "results.csv", "summary.md", "report.html", "vectordb.npz"):
        assert (out / name).exists(), name


def test_ranked_and_complete(analysis):
    _, res = analysis
    rx = res["reactions"]
    assert rx and [r["rank"] for r in rx] == list(range(1, len(rx) + 1))
    assert all(a["score"] >= b["score"] for a, b in itertools.pairwise(rx))
    for r in rx:
        assert r["evidence"]["level"] in ("I", "II")
        assert r["evidence"]["precedent_id"]


def test_economics_and_route(analysis):
    """SDS/SLES: сырьё считается по маршруту от исходного набора (сера, кислород)."""
    _, res = analysis
    priced = [
        r
        for r in res["reactions"]
        if r["economics_usd_per_t"] and r["economics_usd_per_t"]["margin"] is not None
    ]
    assert priced
    sulfated = [r for r in res["reactions"] if "OS(=O)(=O)" in r["product"] and r["route"]]
    assert sulfated, "у сульфатирования должен быть маршрут через SO3"


def test_criteria_summary(analysis):
    _, res = analysis
    s = res["summary"]
    assert s["criterion_2_2"]["unique_vectors"] == s["criterion_2_2"]["vectors"]
    assert {"criterion_2_1", "criterion_2_2", "criterion_2_3"} <= set(s)


def test_report_embeds_data(analysis):
    out, res = analysis
    page = (out / "report.html").read_text(encoding="utf-8")
    assert "<title>" in page and '"reactions"' in page
    assert res["reactions"][0]["reaction_names"].split(" ")[0] in page


def test_route_shows_what_economics_used(analysis):
    """Маршрут в отчёте — те стадии, по которым посчитана себестоимость сырья
    (в т. ч. веществ из прайса, если «сделать» дешевле «купить»)."""
    _, res = analysis
    for r in res["reactions"]:
        e = r["economics_usd_per_t"]
        if e and any(m["basis"] == "route" for m in e["raw_materials"]):
            assert r["route"], r["reaction_names"]
