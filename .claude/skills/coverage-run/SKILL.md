---
name: coverage-run
description: Прогон шаблонного пайплайна (map → extract → apply) на реальном корпусе и сравнение покрытия с прошлым прогоном
disable-model-invocation: true
arguments: [corpus, label]
---

Прогони пайплайн покрытия на корпусе `$corpus` с меткой прогона `$label`.

1. Создай папку `outputs/<сегодняшняя дата YYYY-MM-DD>-$label/` (далее OUT).
2. Проверь формат корпуса: первые 5 строк, есть ли колонки `id` и `rxn_smiles`
   (или спроси, как они называются). Посчитай строки. Если реакции уже
   размечены (есть `:число]` в SMILES), шаг map всё равно нужен: он просто
   переложит их в нужный формат, без модели.
3. Маппинг (долго: ~90 реакций/с на CPU) — запусти в фоне и следи за логом:
   `python -m chem_agent.template_engine map --corpus $corpus --out OUT/mapped.csv`
4. Ручная библиотека:
   `python -m chem_agent.template_engine map --corpus data/manual/inorganic.csv --out OUT/manual_mapped.csv`
5. Извлечение:
   `python -m chem_agent.template_engine extract --mapped OUT/mapped.csv OUT/manual_mapped.csv --trusted OUT/manual_mapped.csv --out OUT/templates.jsonl`
   Запиши из лога: число шаблонов, число с частотой ≥ 2, число типов,
   долю самопроверки, причины пропусков.
6. Применение, два режима:
   - открытый: `apply --templates OUT/templates.jsonl --inputs data/inputs/priority_set.csv --out OUT/candidates.csv --report OUT/coverage.md`
   - внутри набора: то же с `--depth 3 --internal-only --out OUT/network.csv --report OUT/network.md`
7. Возьми выборку из 30 кандидатов (по 5 на соединение, из разных типов) и отдай
   субагенту `chemistry-reviewer` на проверку осмысленности.
8. Сравни с последним прошлым прогоном в `outputs/` (таблица по соединениям:
   реакции, типы, внутри набора, отклонено round-trip).
9. Допиши в `docs/findings.md` раздел «Прогон $label (дата)»: ключевые цифры,
   вердикт ревьюера, что изменилось, открытые вопросы. Не больше 20 строк.

Не меняй код в этом прогоне. Если нашёл дефект — опиши его в отчёте и
предложи отдельной задачей.
