"""
Регрессионные тесты химического движка.

Каждый тест фиксирует решение или найденный дефект (см. docs/decisions/).
Если тест падает после изменения кода — это химическая регрессия, а не
«устаревший тест»: разберитесь в причине, прежде чем править ожидания.
"""

import csv
import json

import pytest
from conftest import INPUTS, PURCHASABLE, ROOT, write_inputs

from chem_agent.template_engine import (
    NO_SOURCE,
    canon,
    frag_set,
    main,
    prepare_for_extraction,
    type_key,
)

with open(INPUTS, encoding="utf-8") as f:
    USER_SET = [r["smiles"] for r in csv.DictReader(f) if r["smiles"] != "resource"]


def products(rows):
    return {r["product"] for r in rows}


def load_templates(path):
    return [json.loads(line) for line in open(path, encoding="utf-8")]


# --- представление веществ -------------------------------------------------


def test_sulfur_smiles_is_not_hydrogen_sulfide():
    """SMILES 'S' — это H2S; элементарную серу пишем как S8."""
    assert canon("S") != canon("S1SSSSSSS1")


def test_counter_ions_ignored_in_fragment_sets():
    assert frag_set("[Na+].[OH-]") == frozenset({"[OH-]"})


# --- извлечение шаблонов ---------------------------------------------------


def test_all_templates_reproduce_their_precedents(templates_path):
    for t in load_templates(templates_path):
        assert t["selfcheck_rate"] == 1.0, t["id"]


def test_solvents_are_not_spectator_reagents():
    """Растворитель, не вошедший в продукт, — агент, а не сопутствующий реагент.
    Иначе «THF + вода» станет обязательным реагентом шаблона и статьёй затрат."""
    mapped = "[CH3:1][C:2](=[O:3])[OH:4].[Na+].[OH-].C1CCOC1.O>>[CH3:1][C:2](=[O:3])[O-:4]"
    *_, spectators, solvents = prepare_for_extraction(mapped)
    assert sorted(spectators) == ["[Na+]", "[OH-]"]
    assert sorted(solvents) == ["C1CCOC1", "O"]


def test_map_keeps_source_column(tmp_path):
    """Уже размеченная реакция проходит map без модели; источник не теряется."""
    corpus = tmp_path / "c.csv"
    with open(corpus, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "rxn_smiles", "source"])
        w.writerow(
            ["X1", "[CH3:1][OH:2].[CH2:3]1[CH2:4][O:5]1>>[CH3:1][O:2][CH2:3][CH2:4][OH:5]", "US1"]
        )
    out = tmp_path / "m.csv"
    main(["map", "--corpus", str(corpus), "--out", str(out)])
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    assert rows[0]["source"] == "US1"


def test_template_count_is_distinct_precedents(tmp_path):
    """R14: одна реакция, записанная дважды с разными растворителями или
    реагентами (патент и его копия в 2naoh_dataset), — один прецедент.
    Иначе фильтр --min-count 2 пропускает разовые ошибки разметки."""
    rxn = "[CH3:1][OH:2].[CH2:3]1[CH2:4][O:5]1>>[CH3:1][O:2][CH2:3][CH2:4][OH:5]"
    mapped = tmp_path / "m.csv"
    with open(mapped, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "mapped_rxn", "agents", "source"])
        w.writerow(["A", rxn, "", "US1"])
        w.writerow(["B", rxn.replace(">>", ".C1CCOC1.[Na+].[OH-]>>"), "", "без источника"])
    out = tmp_path / "t.jsonl"
    main(["extract", "--mapped", str(mapped), "--out", str(out), "--jobs", "1"])
    (t,) = load_templates(out)
    assert t["count"] == 1
    assert [e["id"] for e in t["examples"]] == ["A"]


def test_manual_templates_are_trusted(templates_path):
    trusted = [t for t in load_templates(templates_path) if t.get("trusted")]
    assert len(trusted) >= 7


def test_type_key_merges_same_transformation():
    """Спирт + эпоксид: общий и специфичный (метанол) шаблоны — один тип."""
    general = (
        "[C:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1>>"
        "[C:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]"
    )
    methyl = (
        "[C;D1;H3:1]-[OH;D1;+0:2].[C:3]1-[CH2;D2;+0:4]-[O;H0;D2;+0:5]-1>>"
        "[C;D1;H3:1]-[O;H0;D2;+0:2]-[CH2;D2;+0:4]-[C:3]-[OH;D1;+0:5]"
    )
    assert type_key(general) == type_key(methyl)


# --- применение ------------------------------------------------------------


def test_propylene_oxide_opens_at_ch2(run_apply):
    """Региохимия переносится из прецедентов: атака по CH2, вторичный спирт."""
    rows = run_apply("--purchasable", str(PURCHASABLE))
    prods = products(r for r in rows if r["input_name"] == "Оксид пропилена")
    assert canon("COCC(C)O") in prods
    assert canon("COC(C)CO") not in prods


def test_roundtrip_rejects_epoxide_deoxygenation(run_apply):
    """Слишком общий шаблон этерификации не должен «раскрывать» эпоксид с потерей O."""
    rows = run_apply("--purchasable", str(PURCHASABLE))
    bad = canon("CCCOC(C)=O")  # пропилацетат из ПО + AcOH — химическая ошибка
    good = canon("CC(=O)OCC(C)O")  # гидроксипропилацетат — правильный продукт
    prods = products(r for r in rows if r["input_name"] == "Оксид пропилена")
    assert bad not in prods
    assert good in prods


def test_false_esterification_blocked_even_without_roundtrip(run_apply):
    """Ложная «этерификация» ПО + AcOH -> пропилацетат (маппер увёл кислород
    спирта) требует ухода OH с sp3-углерода без активации — её отсекает
    решение 0012 даже без обратной проверки."""
    rows = run_apply("--purchasable", str(PURCHASABLE), "--no-roundtrip")
    assert canon("CCCOC(C)=O") not in products(rows)


def test_roundtrip_off_lets_false_positive_through(tmp_path):
    """Контроль решения 0002 на реальном шаблоне T00277 (нитрил + этилендиамин ->
    имидазолин): прямой шаблон срабатывает и на N,N-диметилэтилендиамин и
    теряет метилы; ловит это только обратная проверка."""
    templates = ROOT / "tests" / "regression" / "roundtrip_templates.jsonl"
    inputs = write_inputs(
        tmp_path / "in.csv",
        [("S8", "S1SSSSSSS1"), ("диамин", "CN(C)CCN"), ("нитрил", "N#Cc1cccc(Br)c1")],
    )
    wrong = canon("Brc1cccc(C2=NCCN2)c1")
    found = {}
    for flag in ("", "--no-roundtrip"):
        out = tmp_path / f"c{flag}.csv"
        main(
            ["apply", "--templates", str(templates), "--inputs", str(inputs), "--internal-only"]
            + [
                "--min-count",
                "1",
                "--jobs",
                "1",
                "--out",
                str(out),
                "--report",
                str(tmp_path / "r.md"),
            ]
            + ([flag] if flag else [])
        )
        found[flag] = wrong in {r["product"] for r in csv.DictReader(open(out, encoding="utf-8"))}
    assert found == {"": False, "--no-roundtrip": True}


def test_naoh_participates_as_reagent(run_apply):
    rows = run_apply("--purchasable", str(PURCHASABLE))
    naoh = [r for r in rows if r["input_name"] == "Гидроксид натрия"]
    assert any(r["role"] == "reagent" for r in naoh)


def test_neutralization_requires_base_in_set(run_apply, tmp_path):
    """Без щёлочи в наборе кислота не должна «нейтрализоваться сама»."""
    no_base = write_inputs(
        tmp_path / "no_base.csv", [("Додеканол", "CCCCCCCCCCCCO"), ("Кислород", "O=O")]
    )
    rows = run_apply("--depth", "2", "--internal-only", inputs=no_base)
    assert canon("CCCCCCCCCCCC(=O)O") in products(rows)  # окисление есть
    assert canon("CCCCCCCCCCCC(=O)[O-]") not in products(rows)  # нейтрализации нет


def test_internal_only_uses_only_the_set(run_apply):
    rows = run_apply("--internal-only")
    assert rows
    assert {r["partner_source"] for r in rows} <= {"internal", "internal+intermediate"}


def test_internal_label_means_all_participants_from_the_set(run_apply):
    """На стадиях ≥2 входное вещество — промежуточный продукт. Реакция с ним
    не может считаться «все участники из исходного набора» (критерий 2.1)."""
    rows = run_apply("--depth", "3", "--internal-only")
    user = frozenset().union(*(frag_set(s) for s in USER_SET))
    internal = [r for r in rows if r["partner_source"] == "internal"]
    assert internal
    for r in internal:
        reac, reag, _ = r["reaction_smiles"].split(">")
        assert frag_set(reac + ("." + reag if reag else "")) <= user, r["reaction_smiles"]
    # додеканол + SO3: SO3 получен из серы за две стадии, это не «внутри набора»
    sulfation = [r for r in rows if r["product"] == canon("CCCCCCCCCCCCOS(=O)(=O)O")]
    assert sulfation
    assert all(r["partner_source"] == "internal+intermediate" for r in sulfation)


def test_sles_route_found_from_priority_set(run_apply):
    """Из перечня ТЗ за 5 стадий собирается лауретсульфат натрия (SLES)."""
    rows = run_apply("--depth", "5", "--internal-only")
    assert canon("CCCCCCCCCCCCOCCOS(=O)(=O)[O-]") in products(rows)


def test_every_candidate_has_precedent(run_apply):
    rows = run_apply("--purchasable", str(PURCHASABLE))
    assert all(r["precedent_id"] for r in rows)


def test_candidate_shows_precedent_source(run_apply):
    """Источник прецедента виден в выдаче; прецедент без источника помечен явно."""
    rows = run_apply("--internal-only")
    by_rxn = {r["reaction_smiles"]: r for r in rows}
    burn = by_rxn["S1SSSSSSS1.O=O>>O=S=O"]  # ручная библиотека, MAN-S01
    assert burn["precedent_source"].startswith("Сжигание серы")
    assert all(r["precedent_source"] for r in rows)
    # демо-корпус синтетический, источников у него нет
    assert any(r["precedent_source"] == NO_SOURCE for r in rows)


def test_parallel_run_gives_identical_output(tmp_path):
    """Число процессов не влияет на результат: те же шаблоны, те же кандидаты
    в том же порядке, тот же отчёт (кроме строки времени)."""
    from conftest import FIXTURES, INPUTS

    manual = str(FIXTURES / "manual_mapped.csv")
    outs = {}
    for jobs in ("1", "4"):
        d = tmp_path / f"j{jobs}"
        d.mkdir()
        mapped = [str(FIXTURES / "corpus_mapped.csv"), manual]
        main(
            ["extract", "--mapped", *mapped, "--trusted", manual, "--out", str(d / "t.jsonl")]
            + ["--jobs", jobs]
        )
        for mode, extra in (
            ("open", ["--purchasable", str(PURCHASABLE)]),
            ("net", ["--depth", "2"]),  # открытый режим, 2 стадии: ~800 реакций
        ):
            main(
                ["apply", "--templates", str(d / "t.jsonl"), "--inputs", str(INPUTS)]
                + ["--min-count", "1", "--jobs", jobs, *extra]
                + ["--out", str(d / f"{mode}.csv"), "--report", str(d / f"{mode}.md")]
            )
        outs[jobs] = {
            name: [
                line
                for line in (d / name).read_text(encoding="utf-8").splitlines()
                if not line.startswith("- Время:")
            ]
            for name in ("t.jsonl", "open.csv", "open.md", "net.csv", "net.md")
        }
    for name in outs["1"]:
        assert outs["1"][name] == outs["4"][name], name
    assert len(outs["1"]["net.csv"]) > 100  # сравнение не на пустом выводе


def test_oligomer_series_limited_and_counted_once(run_apply, tmp_path):
    """R18, решение 0011: шаблон наращивает свой продукт не больше 3 раз
    (C12E1..C12E3, не C12E4); серия олигомеров в критерии 2.1 — одна реакция.
    Маршрут SLES (разные шаблоны на каждой стадии) сохраняется."""
    rows = run_apply("--depth", "5", "--internal-only")
    prods = products(rows)
    c12e = ["CCCCCCCCCCCCO" + "CCO" * n for n in range(1, 5)]
    assert [canon(p) in prods for p in c12e] == [True, True, True, False]
    assert canon("CCCCCCCCCCCCOCCOS(=O)(=O)[O-]") in prods
    e3 = next(r for r in rows if r["product"] == canon(c12e[2]))
    e1 = next(r for r in rows if r["product"] == canon(c12e[0]))
    assert e3["series"] == e1["series"] and e3["repeat"] == "3"
    assert len({r["series"] for r in rows}) < len({r["reaction_core"] for r in rows})


UTILITIES = ROOT / "data" / "inputs" / "utilities.csv"


def test_utilities_are_available_partners(run_apply):
    """Решение 0011: вода и воздух — всегда доступные ресурсы (ТЗ: «доступные
    ресурсы: электричество, вода, воздух»). ЭО + вода -> этиленгликоль."""
    eg = canon("OCCO")
    without = run_apply("--internal-only")
    assert not any(
        r["product"] == eg and "O" in r["reaction_smiles"].split(">")[0].split(".") for r in without
    )
    rows = run_apply("--internal-only", "--utilities", str(UTILITIES))
    hydration = [r for r in rows if r["product"] == eg and r["partners"] == "O"]
    assert hydration
    assert all(r["partner_source"] == "internal" for r in hydration)
    # ресурсы не раскрываются как входы: у воды своих строк нет
    assert not any(r["input_smiles"] == "O" for r in rows)


@pytest.mark.parametrize(
    "product, name",
    [
        ("CCCCCCCCCCCCOS(=O)(=O)[O-]", "SDS: S8 -> SO2 -> SO3 -> сульфатирование -> NaOH"),
        ("O=S(=O)([O-])CCO", "изэтионат: SO2 + NaOH -> NaHSO3 + ЭО"),
        ("O=S([O-])([O-])=S", "тиосульфат: сульфит + S8"),
        ("C=CCO", "аллиловый спирт из ПО"),
    ],
)
def test_manual_library_routes_from_priority_set(run_apply, product, name):
    """Ручная библиотека (R3, docs/manual-library-proposal.md): маршруты из
    перечня ТЗ с водой и воздухом за 4 стадии."""
    rows = run_apply("--depth", "4", "--internal-only", "--utilities", str(UTILITIES))
    assert canon(product) in products(rows), name


def test_type_key_ignores_charge_and_query_details():
    """SO3 + H2O -> H2SO4 и SO3 + OH- -> HSO4- — одно превращение (MAN-S03,
    MAN-S14): тип не должен зависеть от заряда и порядка уточнений в SMARTS.
    Иначе разнообразие для критерия 2.1 завышается."""
    water = (
        "[O;D1;H0:1]=[S;H0;D3;+0:2](=[O;D1;H0:3])=[O;H0;D1;+0:4].[OH2;D0;+0:5]"
        ">>[O;D1;H0:1]=[S;H0;D4;+0:2](=[O;D1;H0:3])(-[OH;D1;+0:5])-[OH;D1;+0:4]"
    )
    hydroxide = (
        "[O;D1;H0:1]=[S;H0;D3;+0:2](=[O;D1;H0:3])=[O;H0;D1;+0:4].[OH-;D0:5]"
        ">>[O-;H0;D1:5]-[S;H0;D4;+0:2](=[O;D1;H0:1])(=[O;D1;H0:3])-[OH;D1;+0:4]"
    )
    assert type_key(water) == type_key(hydroxide)
