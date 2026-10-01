.PHONY: setup test test-all lint format demo fixtures corpus eval balance templates analyze eval-economics

PY ?= python
ENGINE = $(PY) -m chem_agent.template_engine
DEMO_OUT = outputs/demo
RAW ?= data/raw/all_balanced_reactions_training.csv
LABEL ?= full-v1
CORPUS_OUT ?= outputs/$(shell date +%F)-$(LABEL)
PREP_ARGS ?=
TEMPLATES ?=
IN ?=
INPUTS ?= data/inputs/priority_set.csv
MVP_OUT ?= outputs/$(shell date +%F)-mvp
ANALYZE_ARGS ?= --depth 4

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

corpus:           ## выгрузка all_balanced_reactions -> $(CORPUS_OUT)/corpus.csv, meta.csv
	$(PY) -m chem_agent.corpus_prep --raw $(RAW) --out-dir $(CORPUS_OUT) $(PREP_ARGS)
	# без NaOH-выборок: make corpus LABEL=core-v1 PREP_ARGS="--exclude-subset NAOH 2naoh"

fixtures:         ## пересобрать размеченные фикстуры для тестов (нужен rxnmapper)
	$(ENGINE) map --corpus data/demo/demo_corpus.csv --out tests/fixtures/corpus_mapped.csv
	$(ENGINE) map --corpus data/manual/inorganic.csv --out tests/fixtures/manual_mapped.csv

eval:             ## качество на реакциях с вердиктом химика: make eval TEMPLATES=outputs/<прогон>/templates.jsonl
	$(PY) -m chem_agent.review_eval --templates $(TEMPLATES) --out-dir $(dir $(TEMPLATES))eval

balance:          ## полные уравнения: make balance IN=outputs/<прогон>/network.csv -> *_balanced.csv
	$(PY) -m chem_agent.balancer --in $(IN) --out $(basename $(IN))_balanced.csv

templates:        ## шаблоны из корпуса -> $(CORPUS_OUT)/templates.jsonl (нужен rxnmapper; ~20 мин CPU)
	$(MAKE) corpus
	$(ENGINE) map --corpus $(CORPUS_OUT)/corpus.csv --out $(CORPUS_OUT)/mapped.csv
	$(ENGINE) map --corpus data/manual/inorganic.csv --out $(CORPUS_OUT)/manual_mapped.csv
	$(ENGINE) extract --mapped $(CORPUS_OUT)/mapped.csv $(CORPUS_OUT)/manual_mapped.csv \
		--trusted $(CORPUS_OUT)/manual_mapped.csv --out $(CORPUS_OUT)/templates.jsonl

analyze:          ## MVP: набор -> отранжированные реакции: make analyze TEMPLATES=... [INPUTS=...]
	$(PY) -m chem_agent.pipeline --inputs $(INPUTS) --templates $(TEMPLATES) --out $(MVP_OUT) $(ANALYZE_ARGS)

eval-economics:   ## точность экономики на эталоне: make eval-economics TEMPLATES=...
	$(PY) -m chem_agent.economics_eval --templates $(TEMPLATES) --out outputs/$(shell date +%F)-eval-economics
