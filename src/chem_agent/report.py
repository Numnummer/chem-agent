"""
report.py — HTML-отчёт анализа (самодостаточная страница, светлая и тёмная тема).

Сводка по критериям ТЗ, затем отранжированные реакции; у каждой —
уравнение, экономика на 1 т, маршрут от исходного набора, прецедент,
условия и похожие известные процессы. Данные встроены в страницу.
"""

from __future__ import annotations

import html
import json

STYLE = """
/* Технологический журнал: сводка сверху, ниже — список реакций-карточек,
   раскрываемых до сырья, маршрута и доказательств. */
:root {
  --paper: #f3f5f6; --sheet: #ffffff; --ink: #18222c; --muted: #5b6874;
  --rule: #d6dde2; --sulfur: #a87c00; --sulfur-soft: #f5ecc9;
  --good: #2f7a4f; --warn: #a2620b; --bad: #a83a32;
  --f-display: "IBM Plex Sans Condensed", "Arial Narrow", sans-serif;
  --f-body: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  --f-mono: "IBM Plex Mono", "SFMono-Regular", Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --paper: #11171c; --sheet: #182129; --ink: #e3e9ed; --muted: #95a3ae;
  --rule: #2b3741; --sulfur: #e2b53c; --sulfur-soft: #3a3214;
  --good: #6cc191; --warn: #e3a453; --bad: #ec8a80; color-scheme: dark } }
:root[data-theme="dark"] {
  --paper: #11171c; --sheet: #182129; --ink: #e3e9ed; --muted: #95a3ae;
  --rule: #2b3741; --sulfur: #e2b53c; --sulfur-soft: #3a3214;
  --good: #6cc191; --warn: #e3a453; --bad: #ec8a80; color-scheme: dark }
body { background: var(--paper); color: var(--ink); font: 15px/1.5 var(--f-body); }
.wrap { max-width: 1080px; margin: 0 auto; padding-inline: 16px; padding-block: 28px 64px;
  display: grid; gap: 28px; }
h1, h2 { font-family: var(--f-display); font-weight: 600; text-wrap: balance; margin: 0; }
h1 { font-size: 2rem; line-height: 1.15; }
h2 { font-size: 1.15rem; letter-spacing: .02em; }
.eyebrow { font: 600 .72rem/1 var(--f-body); letter-spacing: .12em; text-transform: uppercase;
  color: var(--muted); }
header { display: grid; gap: 10px; }
.lead { color: var(--muted); max-width: 68ch; margin: 0; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.chip { border: 1px solid var(--rule); background: var(--sheet); border-radius: 999px;
  padding: 2px 10px; font-size: .85rem; }
.chip b { font-family: var(--f-mono); font-weight: 500; color: var(--muted); }
.criteria { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); }
.crit { background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; padding: 14px;
  display: grid; gap: 4px; min-width: 0; }
.crit .num { font: 600 1.6rem/1.1 var(--f-display); font-variant-numeric: tabular-nums; }
.crit p { margin: 0; color: var(--muted); font-size: .88rem; }
.pill { justify-self: start; font-size: .75rem; font-weight: 600; border-radius: 999px; padding: 1px 8px; }
.pill.ok { color: var(--good); border: 1px solid currentColor; }
.pill.no { color: var(--bad); border: 1px solid currentColor; }
.note { font-size: .85rem; color: var(--muted); border-left: 3px solid var(--sulfur);
  padding-left: 10px; margin: 0; max-width: 80ch; }
.tools { display: flex; flex-wrap: wrap; gap: 10px; align-items: end; }
.tools label { display: grid; gap: 3px; font-size: .78rem; color: var(--muted); }
.tools input, .tools select { font: inherit; color: var(--ink); background: var(--sheet);
  border: 1px solid var(--rule); border-radius: 4px; padding: 5px 8px; }
.tools input { min-width: 0; width: min(320px, 100%); }
.count { margin-left: auto; color: var(--muted); font-size: .85rem; font-variant-numeric: tabular-nums; }
.list { display: grid; gap: 8px; }
details.rx { background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; }
details.rx[open] { border-color: var(--sulfur); }
summary { list-style: none; cursor: pointer; padding: 12px 14px; display: grid; gap: 4px 14px;
  grid-template-columns: 2.4rem minmax(0, 1fr) auto; align-items: baseline; }
summary::-webkit-details-marker { display: none; }
summary:focus-visible { outline: 2px solid var(--sulfur); outline-offset: 2px; }
.rank { font: 600 1.05rem/1 var(--f-display); color: var(--muted); font-variant-numeric: tabular-nums; }
.title { min-width: 0; }
.names { font-weight: 600; }
.eq { font-family: var(--f-mono); font-size: .84rem; color: var(--muted); overflow-wrap: anywhere; }
.side { display: grid; justify-items: end; gap: 4px; text-align: right; }
.margin { font: 600 1rem/1 var(--f-mono); font-variant-numeric: tabular-nums; }
.margin.neg { color: var(--bad); } .margin.pos { color: var(--good); }
.tags { display: flex; gap: 6px; flex-wrap: wrap; justify-content: flex-end; }
.tag { font-size: .72rem; border-radius: 4px; padding: 1px 6px; background: var(--paper);
  border: 1px solid var(--rule); white-space: nowrap; }
.tag.lvl1 { background: var(--sulfur-soft); border-color: var(--sulfur); }
.body { border-top: 1px solid var(--rule); padding: 14px; display: grid; gap: 16px;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); }
.body section { min-width: 0; display: grid; gap: 6px; align-content: start; }
.body h3 { margin: 0; font: 600 .72rem/1 var(--f-body); letter-spacing: .1em; text-transform: uppercase;
  color: var(--muted); }
.tbl { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: .85rem; font-variant-numeric: tabular-nums; }
td, th { padding: 3px 6px; border-bottom: 1px solid var(--rule); text-align: left; }
td.n, th.n { text-align: right; }
.mono { font-family: var(--f-mono); font-size: .82rem; overflow-wrap: anywhere; }
ol.route { margin: 0; padding-left: 1.2rem; display: grid; gap: 2px; }
.dim { color: var(--muted); font-size: .85rem; }
footer { color: var(--muted); font-size: .8rem; }
@media (max-width: 560px) {
  summary { grid-template-columns: 2rem minmax(0, 1fr); }
  .side { grid-column: 2; justify-items: start; text-align: left; }
  .tags { justify-content: flex-start; } .count { margin-left: 0; }
}
"""

SCRIPT = """
const data = JSON.parse(document.getElementById('data').textContent);
const fx = data.fx;
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const chem = s => esc(s).replace(/([A-Za-z\\)])(\\d+)/g, '$1<sub>$2</sub>').replace(/->/g, '→');
const money = v => v == null ? '—' : Math.round(v).toLocaleString('ru-RU');
const basis = {price: 'прайс', class: 'класс продукта', route: 'маршрут', '': 'нет цены'};
function card(r) {
  const e = r.economics_usd_per_t;
  const m = e ? e.margin : null;
  const ev = r.evidence;
  const raw = e ? e.raw_materials.map(x => `<tr><td>${esc(x.name)}</td><td class="n">${x.t_per_t.toFixed(3)}</td>
    <td class="n">${x.price == null ? '—' : money(x.price)}</td><td class="dim">${basis[x.basis] ?? ''}</td></tr>`).join('') : '';
  const cond = Object.entries(ev.conditions || {}).map(([k, v]) => `${esc(k)}: ${esc(v)}`).join('; ');
  const route = r.route.length ? `<ol class="route">${r.route.map(s => `<li class="mono">${chem(s)}</li>`).join('')}</ol>`
    : '<p class="dim">Все вещества — из исходного набора или прайса.</p>';
  const sim = r.similar_known.map(s => `<li><span class="mono">${esc(s.id)}</span> · сходство ${s.similarity}
    ${s.source ? '· ' + esc(s.source) : ''}</li>`).join('');
  return `<details class="rx" data-stage="${r.stage}" data-level="${ev.level}" data-text="${esc((r.reaction_names + ' ' + r.equation_formula).toLowerCase())}">
  <summary>
    <span class="rank">${r.rank}</span>
    <span class="title"><span class="names">${esc(r.reaction_names)}</span><br>
      <span class="eq">${chem(r.equation_formula || r.reaction)}</span></span>
    <span class="side"><span class="margin ${m == null ? '' : m < 0 ? 'neg' : 'pos'}">${m == null ? 'нет цены' : (m > 0 ? '+' : '') + money(m) + ' $/т'}</span>
      <span class="tags"><span class="tag">стадия ${r.stage}</span>
      <span class="tag ${ev.level === 'I' ? 'lvl1' : ''}">уровень ${ev.level}</span>
      ${r.catalyst ? `<span class="tag">кат. ${chem(r.catalyst)}</span>` : ''}</span></span>
  </summary>
  <div class="body">
    <section><h3>Экономика на 1 т продукта</h3>
      ${e ? `<div class="tbl"><table><thead><tr><th>Сырьё</th><th class="n">т/т</th><th class="n">$/т</th><th>цена</th></tr></thead>
      <tbody>${raw}</tbody></table></div>
      <p class="dim">Сырьё ${money(e.raw_cost)} $ · выручка ${money(e.revenue)} $ (${basis[e.product_price_basis] ?? ''})
      · кредит за побочные ${money(e.byproduct_credit)} $ · маржа ${money(e.margin)} $ ≈ ${money(e.margin_rub)} ₽</p>
      ${e.unknown_prices.length ? `<p class="dim">Нет цены: ${esc(e.unknown_prices.join(', '))}</p>` : ''}`
      : `<p class="dim">Уравнение не составлено: ${esc(r.balance_reason)}</p>`}
    </section>
    <section><h3>Уравнение и маршрут</h3>
      <p class="mono">${chem(r.equation_formula)}</p>
      ${r.balance_reason ? `<p class="dim">${esc(r.balance_reason)}</p>` : ''}
      ${route}
    </section>
    <section><h3>Доказательства</h3>
      <p>${ev.level === 'I' ? `Реакция есть в корпусе: <span class="mono">${esc(ev.in_corpus)}</span>` :
        `Аналог по шаблону <span class="mono">${esc(ev.template_id)}</span> (${ev.template_precedents} прецед.)`}</p>
      <p class="dim">Прецедент <span class="mono">${esc(ev.precedent_id)}</span> · ${esc(ev.precedent_source)}</p>
      <p class="mono">${esc(ev.precedent_reaction)}</p>
      ${cond ? `<p class="dim">Условия прецедента: ${cond}</p>` : ''}
      ${sim ? `<h3>Похожие известные процессы</h3><ul class="dim">${sim}</ul>` : ''}
    </section>
  </div></details>`;
}
const list = document.getElementById('list');
const q = document.getElementById('q'), st = document.getElementById('stage'), lv = document.getElementById('level');
function render() {
  const text = q.value.trim().toLowerCase();
  const rows = data.reactions.filter(r => (!st.value || String(r.stage) === st.value)
    && (!lv.value || r.evidence.level === lv.value)
    && (!text || (r.reaction_names + ' ' + r.equation_formula).toLowerCase().includes(text)));
  list.innerHTML = rows.map(card).join('') || '<p class="dim">Ничего не найдено — измените фильтр.</p>';
  document.getElementById('count').textContent = `${rows.length} из ${data.reactions.length}`;
}
[q, st, lv].forEach(el => el.addEventListener('input', render));
render();
"""


def render_html(summary: dict, rows: list[dict], cfg: dict, inputs: list[dict]) -> str:
    c1, c2, c3 = summary["criterion_2_1"], summary["criterion_2_2"], summary["criterion_2_3"]

    def pill(ok):
        return (
            '<span class="pill ok">выполнен</span>'
            if ok
            else '<span class="pill no">не выполнен</span>'
        )

    chips = "".join(
        f'<span class="chip">{html.escape(r["name"])} <b>{html.escape(r["smiles"])}</b></span>'
        for r in inputs
    )
    stages = sorted({r["stage"] for r in rows})
    payload = json.dumps(
        {"reactions": rows, "fx": cfg.get("fx_rub_per_usd", 1)}, ensure_ascii=False
    ).replace("</", "<\\/")
    return f"""<title>Реакции из набора сырья</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@600&family=IBM+Plex+Sans:wght@400;600&display=swap">
<style>{STYLE}</style>
<div class="wrap">
<header>
  <span class="eyebrow">Химический модуль НИОКР · прототип</span>
  <h1>Реакции из набора сырья</h1>
  <p class="lead">Система подобрала реакции между выбранными веществами (плюс вода и воздух) до
  {summary["depth"]} стадий, составила уравнения, оценила экономику на 1 т продукта и
  отранжировала результат. У каждой реакции есть прецедент из патентного корпуса или справочника.</p>
  <div class="chips">{chips}</div>
</header>
<section class="criteria" aria-label="Критерии ТЗ">
  <div class="crit"><span class="eyebrow">2.1 · валидные реакции</span>
    <span class="num">{c1["types"]} типов</span>
    <p>{c1["valid_reactions"]} уравненных реакций, {c1["series"]} серий; {c1["level_I"]} есть в корпусе (уровень I).
    Порог — 10 типов превращений.</p>{pill(c1["passed"])}</div>
  <div class="crit"><span class="eyebrow">2.2 · векторная база</span>
    <span class="num">{c2["unique_vectors"]} / {c2["vectors"]}</span>
    <p>Уникальных векторов реакций (разностный и структурный отпечатки).</p>{pill(c2["passed"])}</div>
  <div class="crit"><span class="eyebrow">2.3 · экономика</span>
    <span class="num">{c3["with_full_economics"]} реакций</span>
    <p>С полной оценкой затрат и выручки ({c3["share"]:.0%}). Точность — по методике проекта.</p>
    <span class="pill no">цены — оценка</span></div>
</section>
<p class="note">{html.escape(c3["price_note"])} Выход принят {cfg.get("default_yield", 1):.0%} на стадию;
энергозатраты в прототипе не оцениваются. Ранг — взвешенная сумма маржи, доказательности, длины
маршрута, качества уравнения и использования только исходного набора.</p>
<section style="display:grid;gap:12px">
  <h2>Отранжированные реакции</h2>
  <div class="tools">
    <label for="q">Поиск<input id="q" type="search" placeholder="сульфат, C12H26O…"></label>
    <label for="stage">Стадия<select id="stage"><option value="">все</option>
      {"".join(f'<option value="{s}">{s}</option>' for s in stages)}</select></label>
    <label for="level">Доказательность<select id="level"><option value="">все</option>
      <option value="I">I — есть в корпусе</option><option value="II">II — аналог</option></select></label>
    <span class="count" id="count"></span>
  </div>
  <div class="list" id="list"></div>
</section>
<footer>Шаблоны: {html.escape(summary["templates"])} · время анализа {summary["elapsed_s"]} с ·
цены в $/т, курс {cfg.get("fx_rub_per_usd", "—")} ₽/$.</footer>
</div>
<script type="application/json" id="data">{payload}</script>
<script>{SCRIPT}</script>
"""
