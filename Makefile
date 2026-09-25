.PHONY: setup test test-all lint format demo fixtures

PY ?= python
ENGINE = $(PY) -m chem_agent.template_engine
DEMO_OUT = outputs/demo

setup:            ## установить пакет со всеми зависимостями
	$(PY) -m pip install -e ".[mapping,dev]"

test:             ## быстрые тесты (без rxnmapper), ~10 с
	$(PY) -m pytest -m "not slow"

test-all:         ## все тесты, включая маппинг
	$(PY) -m pytest

lint:             ## проверка стиля и ошибок
	ruff check src tests
	ruff format --check src tests

format:
	ruff format src tests
	ruff check --fix src tests

demo:             ## полный прогон на демо-корпусе -> outputs/demo/
	mkdir -p $(DEMO_OUT)
	$(ENGINE) map --corpus data/demo/demo_corpus.csv --out $(DEMO_OUT)/mapped.csv
	$(ENGINE) map --corpus data/manual/inorganic.csv --out $(DEMO_OUT)/manual_mapped.csv
	$(ENGINE) extract --mapped $(DEMO_OUT)/mapped.csv $(DEMO_OUT)/manual_mapped.csv \
		--trusted $(DEMO_OUT)/manual_mapped.csv --out $(DEMO_OUT)/templates.jsonl
	$(ENGINE) apply --templates $(DEMO_OUT)/templates.jsonl --inputs data/inputs/priority_set.csv \
		--purchasable data/demo/purchasable.csv --min-count 1 \
		--out $(DEMO_OUT)/candidates.csv --report $(DEMO_OUT)/coverage.md
	$(ENGINE) apply --templates $(DEMO_OUT)/templates.jsonl --inputs data/inputs/priority_set.csv \
		--min-count 1 --depth 4 --internal-only \
		--out $(DEMO_OUT)/network.csv --report $(DEMO_OUT)/network.md

fixtures:         ## пересобрать размеченные фикстуры для тестов (нужен rxnmapper)
	$(ENGINE) map --corpus data/demo/demo_corpus.csv --out tests/fixtures/corpus_mapped.csv
	$(ENGINE) map --corpus data/manual/inorganic.csv --out tests/fixtures/manual_mapped.csv
