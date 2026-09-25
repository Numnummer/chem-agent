#!/usr/bin/env bash
# PostToolUse: после каждого Edit/Write форматирует .py-файл и проверяет его ruff.
# Если ошибки остались — exit 2, и Claude получает их в stderr, чтобы исправить.

file_path=$(python3 -c 'import json,sys; print(json.load(sys.stdin).get("tool_input",{}).get("file_path",""))' 2>/dev/null)

case "$file_path" in
  *.py) ;;
  *) exit 0 ;;
esac

command -v ruff >/dev/null 2>&1 || exit 0   # ruff не установлен — не мешаем работе

ruff check --fix --quiet "$file_path" >/dev/null 2>&1   # сначала автоисправления
ruff format --quiet "$file_path" 2>/dev/null             # потом форматирование
if ! out=$(ruff check --quiet "$file_path" 2>&1); then   # остались ошибки -> Claude
  echo "ruff нашёл проблемы в $file_path:" >&2
  echo "$out" >&2
  exit 2
fi
exit 0
