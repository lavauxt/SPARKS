# Styles and script of the unified HTML report ([sparks_report()]).
#
# They live here, as R strings, so the report is built from package code alone:
# no files under inst/, no template lookup, nothing to copy next to the HTML.
# Raw strings (R >= 4.0) keep the CSS / JavaScript verbatim; edit them in place.

.SPARKS_REPORT_CSS <- r"---(
/* SPARKS unified report -- design tokens, layout and components.
   Colours follow the validated palette: categorical hues are applied by JS in
   fixed order; chrome (surface / ink / grid) swaps between light and dark. */
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface: #fcfcfb;
  --text: #0b0b0b;
  --text-2: #52514e;
  --muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --ring: rgba(11, 11, 11, 0.10);
  --accent: #2a78d6;
  --accent-ink: #ffffff;
  --accent-soft: #e7f0fc;
  --good: #006300;
  --warn: #8a5a00;
  --radius: 10px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
    --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
}

* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0; background: var(--page); color: var(--text);
  font: 14px/1.5 var(--font); -webkit-font-smoothing: antialiased;
}
a { color: var(--accent); }
h1, h2, h3 { margin: 0; font-weight: 600; letter-spacing: -0.01em; }
h2 { font-size: 18px; } h3 { font-size: 14px; }
code, pre { font-family: var(--mono); font-size: 12px; }
pre { margin: 0; padding: 12px; overflow: auto; background: var(--page);
      border: 1px solid var(--ring); border-radius: 8px; max-height: 420px; }
.muted { color: var(--muted); } .small { font-size: 12px; }
.sr-only { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); }

/* ---- Header & tabs ---- */
.top { background: var(--surface); border-bottom: 1px solid var(--ring); }
.top-inner { max-width: 1440px; margin: 0 auto; padding: 14px 16px 0; }
.title-row { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }
.title-row h1 { font-size: 20px; }
.badge { display: inline-flex; align-items: center; gap: 6px; padding: 2px 10px;
         border-radius: 999px; background: var(--accent-soft); color: var(--text-2);
         font-size: 12px; border: 1px solid var(--ring); }
.badge.demo { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.spacer { flex: 1; }
.tabs { display: flex; gap: 2px; overflow-x: auto; scrollbar-width: thin; margin-top: 10px; }
.tab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit;
       padding: 10px 14px; cursor: pointer; border-bottom: 2px solid transparent;
       white-space: nowrap; border-radius: 6px 6px 0 0; }
.tab:hover { background: var(--accent-soft); color: var(--text); }
.tab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.tab:focus-visible, button:focus-visible, select:focus-visible, input:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
main { max-width: 1440px; margin: 0 auto; padding: 16px; }
.panel[hidden] { display: none; }

/* ---- Cards, grids, controls ---- */
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius);
        padding: 14px; min-width: 0; }
.card > header { display: flex; align-items: baseline; gap: 8px; margin-bottom: 8px; flex-wrap: wrap; }
.card > header .hint { margin-left: auto; }
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
.g3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.g-main { grid-template-columns: minmax(0, 1fr) 360px; align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text); }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }

/* ---- Legend (HTML so categories can be toggled) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }

/* ---- Tables ---- */
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 520px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 5px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl th { position: sticky; top: 0; background: var(--surface); cursor: pointer; user-select: none;
          font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); z-index: 1; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl tbody tr:hover { background: var(--accent-soft); }
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }

/* ---- Selection panel ---- */
.bars { display: flex; flex-direction: column; gap: 5px; }
.bar-row { display: grid; grid-template-columns: 110px 1fr 70px; gap: 8px; align-items: center; font-size: 12px; }
.bar-row .lab { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; color: var(--text-2); }
.bar-track { background: var(--grid); border-radius: 4px; height: 10px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px 4px 4px 4px; }
.note { font-size: 12px; color: var(--muted); }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }

/* ---- Gallery ---- */
.gallery { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.shot { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; overflow: hidden; margin: 0; }
.shot button { all: unset; display: block; cursor: zoom-in; width: 100%; }
.shot img { display: block; width: 100%; height: 190px; object-fit: contain; background: #fff; }
.shot figcaption { padding: 6px 10px; font-size: 12px; color: var(--text-2); overflow-wrap: anywhere; }
.shot .nopreview { height: 190px; display: grid; place-items: center; color: var(--muted); background: var(--page); }
dialog.lightbox { border: 0; padding: 0; background: transparent; max-width: 96vw; max-height: 96vh; }
dialog.lightbox::backdrop { background: rgba(0, 0, 0, 0.7); }
dialog.lightbox img { max-width: 96vw; max-height: 92vh; background: #fff; border-radius: 8px; display: block; }

@media (max-width: 1000px) { .g-main, .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } .bar-row { grid-template-columns: 80px 1fr 56px; } }
@media print {
  .tabs, .toolbar, .btn { display: none !important; }
  .panel[hidden] { display: block !important; }
  .card { break-inside: avoid; }
}
@media (forced-colors: active) { .bar-fill, .swatch { forced-color-adjust: none; } }

/* ---- Additions: notices, lazy images, section spacing ---- */
.notice { font-size: 12.5px; color: var(--text-2); background: var(--accent-soft);
          border: 1px solid var(--ring); border-radius: 8px; padding: 8px 12px; }
.notice.warn { border-color: var(--warn); }
.mt { margin-top: 12px; }
.shot img[data-pending] { opacity: 0; }
.shot img { transition: opacity .15s; }
)---"

.SPARKS_REPORT_JS <- r"---(
/* SPARKS unified report -- client-side app.
   Reads the small JSON index in #sparks-data. Heavy payloads (datasets, tables,
   images) live in their own inert <script> blocks and are decoded on first use.
   Every figure is drawn with the inlined plotly.js; no network access is needed. */
(function () {
  'use strict';

  // ------------------------------------------------------------------ data & utils
  var DATA = JSON.parse(document.getElementById('sparks-data').textContent);
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  var blobCache = {};
  function blob(id) {
    if (blobCache[id] !== undefined) return blobCache[id];
    var el = document.getElementById('sd-' + id), v = null;
    if (el) { try { v = JSON.parse(el.textContent); } catch (e) { if (window.console) console.error(e); } el.textContent = ''; }
    blobCache[id] = v; return v;
  }
  function imgSrc(im) {
    if (!im.embedded) return null;
    if (im._src) return im._src;
    var el = document.getElementById('si-' + im.id);
    im._src = el ? el.textContent : null;
    return im._src;
  }
  // Tables: header is in DATA.tables, columns are decoded on first use.
  function tdata(T) {
    if (!T.data) T.data = blob(T.id) || {};
    T.cols.forEach(function (c) { if (!T.data[c]) T.data[c] = []; });
    return T;
  }

  function h(tag, attrs) {
    var el = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'class') el.className = v;
      else if (k === 'text') el.textContent = v;
      else if (k.slice(0, 2) === 'on' && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else el.setAttribute(k, v === true ? '' : v);
    });
    for (var i = 2; i < arguments.length; i++) append(el, arguments[i]);
    return el;
  }
  function append(el, kid) {
    if (kid === null || kid === undefined || kid === false) return;
    if (Array.isArray(kid)) kid.forEach(function (k) { append(el, k); });
    else if (kid.nodeType) el.appendChild(kid);
    else el.appendChild(document.createTextNode(String(kid)));
  }
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
  function range(n) { var a = new Array(n); for (var i = 0; i < n; i++) a[i] = i; return a; }
  function uniq(a) { var s = new Set(), o = []; a.forEach(function (v) { if (!s.has(v)) { s.add(v); o.push(v); } }); return o; }
  function isNum(v) { return typeof v === 'number' && isFinite(v); }
  function median(a) {
    if (!a.length) return NaN;
    var s = a.slice().sort(function (p, q) { return p - q; }), m = s.length >> 1;
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  }
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function fmt(v) {
    if (v === null || v === undefined || (typeof v === 'number' && !isFinite(v))) return '';
    if (typeof v !== 'number') return String(v);
    if (Number.isInteger(v)) return String(v);
    if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  function bytes(n) {
    if (n < 1024) return n + ' B';
    if (n < 1048576) return (n / 1024).toFixed(0) + ' KB';
    return (n / 1048576).toFixed(1) + ' MB';
  }
  function download(name, text, type) {
    var blobObj = new Blob([text], { type: type || 'text/csv;charset=utf-8' });
    var a = h('a', { href: URL.createObjectURL(blobObj), download: name });
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  function csvField(v) {
    var s = v === null || v === undefined ? '' : String(v);
    return /[",\n\r\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function toCSV(cols, rowFn, n) {
    var lines = [cols.map(csvField).join(',')];
    for (var i = 0; i < n; i++) lines.push(cols.map(function (c) { return csvField(rowFn(c, i)); }).join(','));
    return lines.join('\n') + '\n';
  }
  function store(key, val) {
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) { /* storage may be unavailable */ }
    return null;
  }
  function rafThrottle(fn) {
    var queued = false;
    return function () { if (queued) return; queued = true; requestAnimationFrame(function () { queued = false; fn(); }); };
  }

  // ------------------------------------------------------------------ theme & palette
  var CHROME = {
    light: { surface: '#fcfcfb', text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', mid: '#f0efec' },
    dark: { surface: '#1a1a19', text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', mid: '#383835' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  var SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
  var DIV = ['#e34948', '#f0a3a2', '#f0efec', '#9ec5f4', '#2a78d6'];
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";
  var SEQ_SCALE = SEQ.map(function (c, i) { return [i / (SEQ.length - 1), c]; });

  var themePref = store('sparks-theme') || 'auto';
  function resolved() {
    if (themePref === 'light' || themePref === 'dark') return themePref;
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  function applyThemeAttr() {
    if (themePref === 'auto') document.documentElement.removeAttribute('data-theme');
    else document.documentElement.setAttribute('data-theme', themePref);
  }
  function tint(hex, amt, toward) {
    if (!amt) return hex;
    var n = function (s) { return [1, 3, 5].map(function (i) { return parseInt(s.slice(i, i + 2), 16); }); };
    var a = n(hex), b = n(toward);
    return '#' + a.map(function (v, i) { return ('0' + Math.round(v + (b[i] - v) * amt).toString(16)).slice(-2); }).join('');
  }
  // First eight levels use the validated palette in fixed order; further levels recycle it
  // with lighter tints (documented deviation) and are always direct-labelled.
  var colorCache = {};
  function catColors(n, theme) {
    var key = theme + n;
    if (colorCache[key]) return colorCache[key];
    var base = CAT[theme], out = [], tints = [0, 0.38, 0.62];
    for (var i = 0; i < n; i++) {
      out.push(tint(base[i % 8], tints[Math.min(Math.floor(i / 8), 2)], theme === 'light' ? '#ffffff' : '#000000'));
    }
    colorCache[key] = out; return out;
  }
  var themeRegex = null, themeMap = null;
  function reEscape(s) { return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }
  function mapTheme(obj, theme) {
    if (theme !== 'dark') return obj;
    if (!themeRegex) {
      themeMap = {};
      Object.keys(DATA.theme_map || {}).forEach(function (k) { themeMap[k.toLowerCase()] = DATA.theme_map[k]; });
      var keys = Object.keys(themeMap);
      themeRegex = keys.length ? new RegExp(keys.map(reEscape).join('|'), 'gi') : null;
    }
    if (!themeRegex) return obj;
    return JSON.parse(JSON.stringify(obj).replace(themeRegex, function (m) { return themeMap[m.toLowerCase()] || m; }));
  }

  // ------------------------------------------------------------------ plotly helpers
  var registry = [];
  var CONFIG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['autoScale2d'] };

  function baseLayout(t, o) {
    o = o || {};
    var c = CHROME[t];
    var ax = function (label, extra) {
      return Object.assign({
        title: { text: label || '', font: { color: c.secondary } },
        gridcolor: c.grid, zerolinecolor: c.axis, linecolor: c.axis,
        tickfont: { color: c.muted }, automargin: true
      }, extra || {});
    };
    return {
      title: o.title ? { text: o.title, x: 0, xanchor: 'left', font: { size: 14, color: c.text } } : undefined,
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: FONT, size: 12, color: c.secondary },
      xaxis: ax(o.xlab, o.xaxis), yaxis: ax(o.ylab, o.yaxis),
      showlegend: !!o.legend,
      legend: { orientation: 'h', y: -0.22, font: { color: c.secondary } },
      margin: o.margin || { l: 56, r: 16, t: o.title ? 44 : 16, b: 48 },
      hoverlabel: { font: { family: FONT } },
      annotations: o.annotations, shapes: o.shapes, barmode: o.barmode, dragmode: o.dragmode, bargap: o.bargap
    };
  }
  // build(theme) -> {data, layout, config}
  function plot(div, build) {
    if (div._note) { clear(div); div._note = false; }
    div._render = function () {
      var s = build(resolved());
      var cfg = Object.assign({}, CONFIG, s.config || {});
      s.layout.autosize = true;
      window.Plotly.react(div, s.data, s.layout, cfg);
      div._dirty = false;
    };
    if (registry.indexOf(div) < 0) registry.push(div);
    div._render();
    return div;
  }
  function specPlot(div, spec) {
    plot(div, function (t) {
      var s = mapTheme(JSON.parse(JSON.stringify(spec)), t);
      s.layout = s.layout || {};
      delete s.layout.width; delete s.layout.height; delete s.layout.title;   // the card header carries the title
      if (s.layout.margin) s.layout.margin.t = Math.min(s.layout.margin.t || 16, 16);
      s.config = Object.assign({}, s.config || {});
      s.data = s.data || [];
      return s;
    });
  }
  var themeListeners = [];
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
    themeListeners.forEach(function (cb) { cb(); });
  }
  // Replace a figure with a plain-text note (no stale plotly state left behind).
  function showNote(div, text) {
    if (window.Plotly) window.Plotly.purge(div);
    clear(div); div._render = null; div._note = true; div._dirty = false;
    div.appendChild(emptyNote(text));
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (d._dirty) d._render(); else if (window.Plotly && d.data) window.Plotly.Plots.resize(d);
    });
  }
  function plotBox(cls) { return h('div', { class: 'plot ' + (cls || '') }); }
  function emptyNote(text) { return h('div', { class: 'empty', text: text }); }
  function notice(text, warn) { return h('p', { class: 'notice' + (warn ? ' warn' : ''), text: text }); }

  // ------------------------------------------------------------------ tables
  function colType(T, c) {
    var arr = T.data[c];
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // Build a table object from plain columns: {col: [..]}.
  function makeTable(cols, data, extra) {
    var n = cols.length ? (data[cols[0]] || []).length : 0;
    return Object.assign({ cols: cols, data: data, n: n }, extra || {});
  }
  // T = {cols, data: {col: [..]}, n}. opts: filename, pageSize, render: {col: fn(value,i)->node}, sort: [col, 'asc'|'desc']
  function tableView(host, T, opts) {
    opts = opts || {};
    tdata(T);
    var n = T.n, pageSize = opts.pageSize || 25, page = 0, filter = '', pred = null;
    var sortCol = opts.sort ? opts.sort[0] : null, sortDir = opts.sort ? opts.sort[1] : 'asc';
    var types = {}; T.cols.forEach(function (c) { types[c] = colType(T, c); });
    var idx = range(n);
    var search = h('input', { type: 'search', placeholder: 'Filter rows…', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var tbody = h('tbody'), thead = h('thead');
    var table = h('table', { class: 'tbl' }, thead, tbody);
    var wrap = h('div', { class: 'tbl-wrap' }, table);
    var dl = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      download((opts.filename || 'table') + '.csv', toCSV(T.cols, function (c, k) { return T.data[c][idx[k]]; }, idx.length));
    } });
    host.appendChild(h('div', { class: 'tbl-bar' }, n > 10 ? search : null, count, dl, pager));
    host.appendChild(wrap);
    if (opts.note) host.appendChild(h('div', { class: 'note', text: opts.note }));

    function recompute() {
      var f = filter.toLowerCase();
      idx = [];
      for (var i = 0; i < n; i++) {
        if (pred && !pred(i)) continue;
        if (f) {
          var hit = false;
          for (var j = 0; j < T.cols.length && !hit; j++) {
            var v = T.data[T.cols[j]][i];
            if (v !== null && v !== undefined && String(v).toLowerCase().indexOf(f) >= 0) hit = true;
          }
          if (!hit) continue;
        }
        idx.push(i);
      }
      if (sortCol) {
        var arr = T.data[sortCol], dir = sortDir === 'asc' ? 1 : -1, num = types[sortCol] === 'number';
        idx.sort(function (a, b) {
          var x = arr[a], y = arr[b];
          var xm = x === null || x === undefined || (num && !isFinite(x)), ym = y === null || y === undefined || (num && !isFinite(y));
          if (xm || ym) return xm === ym ? 0 : (xm ? 1 : -1);
          return (num ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true })) * dir;
        });
      }
      page = Math.min(page, Math.max(0, Math.ceil(idx.length / pageSize) - 1));
      draw();
    }
    function draw() {
      var tr = h('tr');
      T.cols.forEach(function (c) {
        tr.appendChild(h('th', {
          scope: 'col', tabindex: '0', text: c, 'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); }
        }));
      });
      clear(thead).appendChild(tr);
      var frag = document.createDocumentFragment();
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      for (var k = start; k < end; k++) {
        var i = idx[k], row = h('tr');
        T.cols.forEach(function (c) {
          var v = T.data[c][i], r = opts.render && opts.render[c];
          row.appendChild(h('td', { class: types[c] === 'number' ? 'num' : '' }, r ? r(v, i) : fmt(v)));
        });
        frag.appendChild(row);
      }
      clear(tbody).appendChild(frag);
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows' + (T.truncated ? ' (truncated from ' + fmtInt(T.n_total) + ')' : '');
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      if (pages > 1) {
        pager.appendChild(h('button', { class: 'btn', type: 'button', text: '‹', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
        pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
        pager.appendChild(h('button', { class: 'btn', type: 'button', text: '›', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
      }
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return { setPredicate: function (fn) { pred = fn; page = 0; recompute(); }, refresh: recompute };
  }

  // ------------------------------------------------------------------ datasets
  var state = { ds: 0, listeners: [] };
  // Dataset header (always present) merged with its payload (decoded on first use).
  function dsAt(i) {
    var d = DATA.datasets[i];
    if (!d._loaded) {
      var b = blob(d.id) || { cells: [], x: [], y: [], meta: {}, genes: {}, default_genes: [] };
      Object.keys(b).forEach(function (k) { d[k] = b[k]; });
      d._loaded = true;
    }
    return d;
  }
  function ds() { return dsAt(state.ds); }
  function onDataset(cb) { state.listeners.push(cb); }
  function datasetSelect() {
    if (DATA.datasets.length < 2) return null;
    var sel = h('select', { 'aria-label': 'Dataset', 'data-role': 'dataset' }, DATA.datasets.map(function (d, i) {
      return h('option', { value: i, text: d.label + (d.kind === 'subset' ? ' (subset, ' : ' (') + fmtInt(d.n_cells) + ' cells)' });
    }));
    sel.value = state.ds;
    sel.addEventListener('change', function () {
      state.ds = +sel.value;
      $$('select[data-role="dataset"]').forEach(function (o) { o.value = state.ds; });
      state.listeners.forEach(function (cb) { cb(); });
    });
    return h('label', { class: 'field' }, 'Dataset', sel);
  }
  function catCols(d) { return Object.keys(d.meta).filter(function (k) { return d.meta[k].type === 'cat'; }); }
  function numCols(d) { return Object.keys(d.meta).filter(function (k) { return d.meta[k].type === 'num'; }); }
  function geneNames(d) { return Object.keys(d.genes); }
  function preferredGroup(d, cats) {
    cats = cats || catCols(d);
    var pref = ['cell_type', 'celltype', 'singleR_labels_main', 'SingleR_main', 'seurat_clusters', 'cluster', 'orig.ident'];
    for (var i = 0; i < pref.length; i++) if (cats.indexOf(pref[i]) >= 0) return pref[i];
    return cats[0];
  }
  var geneCache = {};
  function geneVec(d, g) {
    var key = d.id + '|' + g;
    if (geneCache[key]) return geneCache[key];
    var v = new Float32Array(d.cells.length), e = d.genes[g];
    if (e) for (var i = 0; i < e.i.length; i++) v[e.i[i]] = e.v[i];
    geneCache[key] = v; return v;
  }
  function numVec(d, key) {
    var m = d.meta[key];
    if (m && m.type === 'num') return Float64Array.from(m.values, function (v) { return v === null ? NaN : v; });
    return geneVec(d, key);
  }
  function labelsOf(d, col) {
    var m = d.meta[col];
    return m.codes.map(function (c) { return m.levels[c]; });
  }

  // ------------------------------------------------------------------ figure builders (client side)
  function volcano(div, V, o) {
    plot(div, function (t) {
      var c = CHROME[t];
      var up = { x: [], y: [], t: [] }, dn = { x: [], y: [], t: [] }, ns = { x: [], y: [], t: [] }, lab = [];
      for (var i = 0; i < V.lfc.length; i++) {
        var l = V.lfc[i], p = V.p[i];
        if (!isNum(l) || !isNum(p)) continue;
        var y = -Math.log10(Math.max(p, 1e-300)), k = ns;
        if (p < o.pCut && l >= o.lfcCut) k = up; else if (p < o.pCut && l <= -o.lfcCut) k = dn;
        k.x.push(l); k.y.push(y); k.t.push(V.genes[i]);
        if (k !== ns) lab.push({ x: l, y: y, g: V.genes[i], up: k === up });
      }
      var type = V.lfc.length > 4000 ? 'scattergl' : 'scatter';
      var mk = function (name, s, color, opac) {
        return { type: type, mode: 'markers', name: name + ' (' + s.x.length + ')', x: s.x, y: s.y, text: s.t,
          marker: { color: color, size: 6, opacity: opac }, hovertemplate: '%{text}<br>log2FC %{x:.2f}<br>-log10 p %{y:.1f}<extra></extra>' };
      };
      var anns = [];
      [true, false].forEach(function (isUp) {
        lab.filter(function (a) { return a.up === isUp; }).sort(function (a, b) { return b.y - a.y; }).slice(0, o.labelTop).forEach(function (a) {
          anns.push({ x: a.x, y: a.y, text: a.g, showarrow: false, xanchor: a.x > 0 ? 'left' : 'right', yanchor: 'bottom', font: { size: 10, color: c.secondary } });
        });
      });
      var line = { color: c.axis, width: 1, dash: 'dot' };
      var shapes = [
        { type: 'line', x0: -o.lfcCut, x1: -o.lfcCut, y0: 0, y1: 1, yref: 'paper', line: line },
        { type: 'line', x0: o.lfcCut, x1: o.lfcCut, y0: 0, y1: 1, yref: 'paper', line: line }];
      if (o.pCut > 0 && o.pCut <= 1) shapes.push({ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: -Math.log10(o.pCut), y1: -Math.log10(o.pCut), line: line });
      return {
        data: [mk('Not significant', ns, c.muted, 0.45), mk('Down', dn, DIV[4], 0.85), mk('Up', up, DIV[0], 0.85)],
        layout: baseLayout(t, { title: o.title, xlab: 'log2 fold change', ylab: '-log10 adjusted p', legend: true, annotations: anns, shapes: shapes })
      };
    });
  }

  function violin(div, d, values, groupCol, title, ylab) {
    var m = d.meta[groupCol], levels = m.levels;
    var by = levels.map(function () { return []; });
    for (var i = 0; i < values.length; i++) { var v = values[i]; if (isNum(v)) by[m.codes[i]].push(v); }
    plot(div, function (t) {
      var cols = catColors(levels.length, t);
      var data = [];
      levels.forEach(function (lv, k) {
        var v = by[k];
        if (!v.length) return;
        if (v.length > 3000) { var step = v.length / 3000; v = range(3000).map(function (j) { return v[Math.floor(j * step)]; }); }
        data.push({ type: 'violin', y: v, x: v.map(function () { return lv; }), name: lv, points: false, spanmode: 'hard', scalemode: 'width',
          width: 0.85, fillcolor: cols[k], line: { color: cols[k], width: 1 }, meanline: { visible: true }, opacity: 0.9, hoverinfo: 'y+name' });
      });
      return { data: data, layout: baseLayout(t, { title: title, xlab: groupCol, ylab: ylab || 'Expression (log-normalised)', legend: false, xaxis: { type: 'category' } }) };
    });
  }

  function dotplot(div, d, genes, groupCol) {
    var m = d.meta[groupCol], levels = m.levels, cells = d.cells.length;
    var rows = [];
    genes.forEach(function (g) {
      var v = geneVec(d, g), sum = levels.map(function () { return 0; }), pos = sum.slice(), cnt = sum.slice();
      for (var i = 0; i < cells; i++) { var k = m.codes[i]; cnt[k]++; sum[k] += v[i]; if (v[i] > 0) pos[k]++; }
      var mean = sum.map(function (s, k) { return cnt[k] ? s / cnt[k] : 0; });
      var mu = mean.reduce(function (a, b) { return a + b; }, 0) / mean.length;
      var sd = Math.sqrt(mean.reduce(function (a, b) { return a + (b - mu) * (b - mu); }, 0) / Math.max(1, mean.length - 1));
      levels.forEach(function (lv, k) { rows.push({ g: g, lv: lv, mean: mean[k], pct: cnt[k] ? 100 * pos[k] / cnt[k] : 0, z: sd > 0 ? (mean[k] - mu) / sd : 0 }); });
    });
    var lim = rows.reduce(function (a, r) { return Math.max(a, Math.abs(r.z)); }, 1);
    plot(div, function (t) {
      return {
        data: [{ type: 'scatter', mode: 'markers', x: rows.map(function (r) { return r.g; }), y: rows.map(function (r) { return r.lv; }),
          text: rows.map(function (r) { return r.g + ' / ' + r.lv + '<br>mean ' + fmt(r.mean) + '<br>' + Math.round(r.pct) + '% of cells'; }), hoverinfo: 'text',
          marker: { size: rows.map(function (r) { return 4 + 0.22 * r.pct; }), sizemode: 'diameter', color: rows.map(function (r) { return r.z; }), cmin: -lim, cmax: lim,
            colorscale: SEQ_SCALE, colorbar: { title: { text: 'Scaled mean' }, thickness: 10, len: 0.6 },
            line: { color: CHROME[t].surface, width: 1 } } }],
        layout: baseLayout(t, { ylab: groupCol, legend: false, xaxis: { type: 'category', categoryarray: genes }, yaxis: { type: 'category', autorange: 'reversed' }, margin: { l: 110, r: 16, t: 16, b: 80 } })
      };
    });
  }

  function stackedBars(div, xLabels, series, o) {
    plot(div, function (t) {
      var c = CHROME[t], cols = catColors(series.length, t);
      var data = series.map(function (s, k) {
        return { type: 'bar', name: s.name, x: xLabels, y: s.values, marker: { color: cols[k], line: { color: c.surface, width: 2 } },
          customdata: s.counts, hovertemplate: '%{x}<br>' + s.name + ': %{y:.1f}' + (o.percent ? '% (%{customdata} cells)' : '') + '<extra></extra>' };
      });
      return { data: data, layout: baseLayout(t, { title: o.title, xlab: o.xlab, ylab: o.percent ? 'Share of cells' : 'Cells', legend: true, barmode: 'stack', bargap: 0.3,
        xaxis: { type: 'category' }, yaxis: o.percent ? { range: [0, 100], ticksuffix: '%' } : {} }) };
    });
  }

  function scatterQC(div, d, xKey, yKey, groupCol) {
    var x = numVec(d, xKey), y = numVec(d, yKey), m = groupCol ? d.meta[groupCol] : null;
    plot(div, function (t) {
      var levels = m ? m.levels : ['All cells'], cols = catColors(levels.length, t);
      var tr = levels.map(function (lv) { return { x: [], y: [] }; });
      for (var i = 0; i < x.length; i++) {
        if (!isFinite(x[i]) || !isFinite(y[i])) continue;
        var k = m ? m.codes[i] : 0; tr[k].x.push(x[i]); tr[k].y.push(y[i]);
      }
      var type = x.length > 4000 ? 'scattergl' : 'scatter';
      return {
        data: levels.map(function (lv, k) {
          return { type: type, mode: 'markers', name: lv, x: tr[k].x, y: tr[k].y, marker: { color: cols[k], size: 3, opacity: 0.5 },
            hovertemplate: lv + '<br>' + xKey + ' %{x}<br>' + yKey + ' %{y}<extra></extra>' };
        }),
        layout: baseLayout(t, { xlab: xKey, ylab: yKey, legend: levels.length > 1 && levels.length <= 30 })
      };
    });
  }

  // ------------------------------------------------------------------ tabs
  var TABS = [];
  function tab(id, label, available, build) { TABS.push({ id: id, label: label, available: available, build: build, built: false }); }
  function card(title, hint, body) {
    return h('section', { class: 'card' }, h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function tablesOfKind(kinds) { return DATA.tables.filter(function (t) { return kinds.indexOf(t.kind) >= 0; }); }
  function figuresOf(tabId, group) { return DATA.figures.filter(function (f) { return f.tab === tabId && (group === undefined || (f.group || '') === group); }); }
  function figureCards(figs, cls) {
    return figs.map(function (f) {
      var box = plotBox(cls || '');
      var c = card(f.title, f.group || '', box);
      setTimeout(function () { specPlot(box, f.spec); }, 0);
      return c;
    });
  }
  function nav(id) { activate(id, true); }
  function pickCol(T, cands) { for (var i = 0; i < cands.length; i++) if (T.cols.indexOf(cands[i]) >= 0) return cands[i]; return null; }
  var RUN = (DATA.settings && DATA.settings.run) || {};
  function ranOff(flag) { return RUN[flag] === false; }
  function kv(pairs) {
    var dl = h('dl', { class: 'kv' });
    pairs.forEach(function (p) { if (p[1] === null || p[1] === undefined || p[1] === '') return; dl.appendChild(h('dt', { text: p[0] })); dl.appendChild(h('dd', { text: Array.isArray(p[1]) ? p[1].join(', ') : String(p[1]) })); });
    return dl;
  }
  function countsTable(c, label) {
    var tot = c.counts.reduce(function (a, b) { return a + b; }, 0);
    return makeTable([label, 'cells', 'percent'], {
      [label]: c.levels.slice(), cells: c.counts.slice(),
      percent: c.counts.map(function (v) { return tot ? Math.round(1000 * v / tot) / 10 : 0; }) });
  }

  // ---- Overview
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var S = DATA.summary;
    var tiles = [['cells', 'Cells'], ['samples', 'Samples'], ['conditions', 'Conditions'], ['genes', 'Genes'], ['clusters', 'Clusters'], ['tables', 'Result tables'], ['files', 'Files']];
    root.appendChild(h('div', { class: 'tiles' }, tiles.filter(function (t) { return S[t[0]] !== null && S[t[0]] !== undefined; }).map(function (t) {
      return h('div', { class: 'tile' }, h('div', { class: 'v', text: fmtInt(S[t[0]]) }), h('div', { class: 'l', text: t[1] }));
    })));
    var left = [], right = [];
    figuresOf('overview').forEach(function (f) { var b = plotBox('short'); left.push(card(f.title, f.group || '', b)); setTimeout(function () { specPlot(b, f.spec); }, 0); });
    if (DATA.datasets.length) {
      var dsRows = DATA.datasets.map(function (d) { return d; });
      left.push(card('Datasets', 'objects saved by the run', (function () {
        var host = h('div');
        tableView(host, makeTable(['dataset', 'kind', 'comparison', 'cells', 'genes', 'clusters'], {
          dataset: dsRows.map(function (d) { return d.label; }), kind: dsRows.map(function (d) { return d.kind; }),
          comparison: dsRows.map(function (d) { return d.comparison || ''; }), cells: dsRows.map(function (d) { return d.facts.cells; }),
          genes: dsRows.map(function (d) { return d.facts.genes; }), clusters: dsRows.map(function (d) { return d.facts.clusters_n === undefined ? null : d.facts.clusters_n; }) }),
          { filename: 'sparks_datasets', pageSize: 10 });
        return host;
      })()));
    }
    var marks = tablesOfKind(['markers_all', 'cluster_markers', 'deg_sig'])[0];
    if (marks) {
      tdata(marks);
      var gc = pickCol(marks, ['gene', 'Gene']), cc = pickCol(marks, ['cluster', 'group']), lc = pickCol(marks, ['avg_log2FC', 'logFC']);
      if (gc && cc && lc) {
        var byG = {};
        for (var i = 0; i < marks.n; i++) { var g = String(marks.data[cc][i]); (byG[g] = byG[g] || []).push(i); }
        var dl = h('dl', { class: 'kv' });
        Object.keys(byG).sort(function (a, b) { return a.localeCompare(b, undefined, { numeric: true }); }).forEach(function (g) {
          var ii = byG[g].filter(function (i) { return isNum(marks.data[lc][i]); }).sort(function (a, b) { return marks.data[lc][b] - marks.data[lc][a]; }).slice(0, 4);
          dl.appendChild(h('dt', { text: g })); dl.appendChild(h('dd', { text: ii.map(function (i) { return marks.data[gc][i]; }).join(', ') }));
        });
        right.push(card('Top markers', marks.label, dl));
      }
    }
    var cmp = S.comparisons || [];
    right.push(card('Run', null, kv([['Results folder', DATA.meta.results_dir_name], ['Comparisons', cmp.length ? cmp : '—'],
      ['Species', DATA.meta.species || '—'], ['Generated', DATA.meta.generated], ['SPARKS', DATA.meta.version]])));
    var jump = TABS.filter(function (t) { return t.id !== 'overview' && t.available(); }).map(function (t) {
      return h('button', { class: 'btn', type: 'button', text: t.label, onclick: function () { nav(t.id); } });
    });
    right.push(card('Jump to', null, h('div', { class: 'toolbar', style: { marginBottom: 0 } }, jump)));
    root.appendChild(h('div', { class: 'grid g2 mt' }, h('div', { class: 'stack' }, left), h('div', { class: 'stack' }, right)));
  });

  // ---- QC
  tab('qc', 'QC', function () {
    return figuresOf('qc').length > 0 || DATA.datasets.length > 0 || !!DATA.settings || !!(DATA.qc && DATA.qc.retention) ||
      DATA.images.some(function (i) { return i.category === 'QC'; });
  }, function (root) {
    var Sx = DATA.settings, top = [];
    if (Sx) {
      var q = Sx.qc || {};
      top.push(card('Filtering thresholds', 'from config_used.yaml', h('div', null, kv([
        ['Min detected genes (>)', q.min_features], ['Max detected genes (<)', q.max_features], ['Max counts (<)', q.max_counts],
        ['Max mitochondrial % (<)', q.max_mt_percent], ['Min cells per gene', q.min_cells], ['Mitochondrial pattern', Sx.mt_pattern],
        ['Gene removal pattern', Sx.gene_removal_pattern]]), h('p', { class: 'note', text: 'Doublets were removed with scDblFinder after these filters.' }))));
    }
    var R = DATA.qc || {};
    if (R.retention) {
      var hostR = h('div');
      tableView(hostR, makeTable(R.retention.cols, R.retention.data), { filename: 'sparks_cell_retention', pageSize: 15 });
      top.push(card('Cell retention by sample', 'before filtering → after QC and doublet removal', hostR));
    }
    if (top.length) root.appendChild(h('div', { class: 'grid g2' }, top));
    var figs = figuresOf('qc', '');
    if (figs.length) root.appendChild(h('div', { class: 'grid g2 mt' }, figureCards(figs)));

    if (DATA.datasets.length) {
      var host = h('div', { class: 'stack' }), sel = datasetSelect();
      root.appendChild(h('div', { class: 'mt' }, card('Per-sample QC (after filtering)', 'violins from the embedded cells; summary table from all cells', h('div', { class: 'stack' }, sel, host))));
      var draw = function () {
        var d = ds(), F = d.facts || {};
        clear(host);
        var g = d.meta['orig.ident'] ? 'orig.ident' : (d.meta.condition ? 'condition' : preferredGroup(d));
        var metrics = ['nFeature_RNA', 'nCount_RNA', 'percent.mt'].filter(function (m) { return d.meta[m] && d.meta[m].type === 'num'; });
        if (F.samples) {
          var S = F.samples, cols = ['sample', 'cells'].concat(S.condition ? ['condition'] : [], S.quant_format ? ['input format'] : []);
          var tbl = makeTable(cols, { sample: S.sample, cells: S.cells, condition: S.condition, 'input format': S.quant_format });
          var hs = h('div'); tableView(hs, tbl, { filename: 'sparks_samples', pageSize: 15 });
          var note = Sx && Sx.input && Sx.input.save_raw_matrix ? 'Raw, pre-QC count matrices were saved per sample under <comparison_group>/RawMatrix/<folder_id>/ (' + (Sx.input.raw_matrix_format || 'mtx') + ').' : null;
          host.appendChild(h('div', null, h('h3', { text: 'Samples and input data' }), hs, note ? h('p', { class: 'note', text: note }) : null));
        }
        if (F.qc) {
          var hq = h('div'); tableView(hq, makeTable(F.qc.cols, F.qc.data), { filename: 'sparks_qc_summary', pageSize: 15 });
          host.appendChild(h('div', null, h('h3', { text: 'QC summary (median, Q1, Q3)' }), hq));
        }
        if (!metrics.length || !g) { host.appendChild(emptyNote('No QC metric columns in this dataset.')); return; }
        var row = h('div', { class: 'grid g3' });
        metrics.forEach(function (m) {
          var b = plotBox('short'); row.appendChild(b);
          setTimeout(function () { violin(b, d, d.meta[m].values, g, m, m); }, 0);
        });
        host.appendChild(row);
        if (d.meta.nCount_RNA && d.meta.nCount_RNA.type === 'num') {
          var sc = h('div', { class: 'grid g2' });
          ['percent.mt', 'nFeature_RNA'].filter(function (m) { return d.meta[m] && d.meta[m].type === 'num'; }).forEach(function (m) {
            var b = plotBox('short'); sc.appendChild(b);
            setTimeout(function () { scatterQC(b, d, 'nCount_RNA', m, g); }, 0);
          });
          host.appendChild(sc);
        }
        var cc = F.counts || {}, blocks = [];
        [['Phase', 'Cell cycle phase', 'cell_cycle'], ['Assigned_Sex', 'Sex assignment', 'sex_scoring'], ['scDblFinder.class', 'Doublet class', null]].forEach(function (p) {
          if (cc[p[0]]) { var hh = h('div'); tableView(hh, countsTable(cc[p[0]], p[0]), { filename: 'sparks_' + p[0] }); blocks.push(card(p[1], 'all cells', hh)); }
          else if (p[2] && RUN[p[2]]) blocks.push(card(p[1], null, notice('Enabled in the configuration but no ' + p[0] + ' column was found in this object.', true)));
        });
        if (blocks.length) host.appendChild(h('div', { class: 'grid g3' }, blocks));
        else host.appendChild(notice('Cell cycle and sex scoring were not run for this dataset.'));
      };
      draw(); onDataset(draw);
    }
    var alevin = figuresOf('qc', 'Alevin barcode diagnostics');
    if (alevin.length) root.appendChild(h('div', { class: 'mt' }, card('Alevin barcode-level diagnostics', 'whitelisted barcodes are opaque', h('div', { class: 'grid g3' }, figureCards(alevin)))));
    var qcImgs = DATA.images.filter(function (i) { return i.category === 'QC'; });
    if (qcImgs.length) root.appendChild(h('div', { class: 'mt' }, card('Pipeline QC figures', qcImgs.length + ' images (pre / post filtering)', gallery(qcImgs, false))));
  });

  // ---- Embedding explorer
  tab('embedding', 'Embedding', function () { return DATA.datasets.length > 0; }, function (root) {
    var hidden = {}, sel = [], selSet = null, cache = {};
    var plotDiv = plotBox('tall'), legend = h('div', { class: 'legend' });
    var colorSel = h('select', { 'aria-label': 'Colour by' }), size = h('input', { type: 'range', min: 1, max: 10, step: 0.5, value: 4, 'aria-label': 'Point size' });
    var opac = h('input', { type: 'range', min: 0.1, max: 1, step: 0.05, value: 0.8, 'aria-label': 'Opacity' });
    var selBody = h('div', { class: 'stack' });
    var dsSel = datasetSelect();

    function fillColor() {
      var d = ds(), cats = catCols(d), nums = numCols(d), genes = geneNames(d);
      clear(colorSel);
      var add = function (label, items) { if (items.length) colorSel.appendChild(h('optgroup', { label: label }, items.map(function (k) { return h('option', { value: k, text: k }); }))); };
      add('Categories', cats); add('QC & numeric', nums); add('Genes', genes);
      colorSel.value = preferredGroup(d, cats) || (nums[0] || genes[0]);
      hidden = {}; sel = []; selSet = null; cache = {};
    }
    // Per (dataset, colour key): draw order, labels and label positions. Computed once.
    function prep(d, key) {
      var ck = d.id + '|' + key;
      if (cache[ck]) return cache[ck];
      var n = d.cells.length, isCat = d.meta[key] && d.meta[key].type === 'cat';
      var vals = isCat ? null : numVec(d, key);
      var ord = range(n);
      if (!isCat) ord.sort(function (a, b) { var x = vals[a], y = vals[b]; x = isFinite(x) ? x : -Infinity; y = isFinite(y) ? y : -Infinity; return x - y; });
      var pos = new Int32Array(n); ord.forEach(function (c, p) { pos[c] = p; });
      var labels = isCat ? labelsOf(d, key) : null, levels = isCat ? d.meta[key].levels : [];
      var centers = {};
      if (isCat && levels.length > 1 && levels.length <= 40) {
        var xs = levels.map(function () { return []; }), ys = levels.map(function () { return []; }), codes = d.meta[key].codes;
        for (var i = 0; i < n; i++) { xs[codes[i]].push(d.x[i]); ys[codes[i]].push(d.y[i]); }
        levels.forEach(function (lv, k) { if (xs[k].length) centers[lv] = [median(xs[k]), median(ys[k])]; });
      }
      var text = ord.map(function (i) { return d.cells[i] + '<br>' + (isCat ? labels[i] : key + ': ' + fmt(vals[i])); });
      cache[ck] = { isCat: isCat, vals: vals, ord: ord, pos: pos, labels: labels, levels: levels, centers: centers, text: text,
        x: ord.map(function (i) { return d.x[i]; }), y: ord.map(function (i) { return d.y[i]; }) };
      return cache[ck];
    }
    function draw() {
      var d = ds(), key = colorSel.value;
      if (!key) { showNote(plotDiv, 'No metadata or genes to colour by.'); return; }
      var P = prep(d, key), n = d.cells.length;
      renderLegend(P.levels);
      plot(plotDiv, function (t) {
        var c = CHROME[t], cols = P.isCat ? catColors(P.levels.length, t) : null;
        var marker = { size: +size.value };
        if (P.isCat) {
          var codes = d.meta[key].codes;
          marker.color = P.ord.map(function (i) { return cols[codes[i]]; });
          marker.opacity = P.ord.map(function (i) { return hidden[P.labels[i]] ? 0.04 : +opac.value; });
        } else {
          marker.color = P.ord.map(function (i) { var v = P.vals[i]; return isFinite(v) ? v : null; }); marker.opacity = +opac.value;
          marker.colorscale = SEQ_SCALE;
          marker.colorbar = { title: { text: key }, thickness: 10, len: 0.6 };
        }
        var trace = {
          type: n > 8000 ? 'scattergl' : 'scatter', mode: 'markers', x: P.x, y: P.y,
          customdata: P.ord, marker: marker, hoverinfo: 'text', text: P.text,
          selected: { marker: { opacity: 1 } }, unselected: { marker: { opacity: 0.12 } }
        };
        if (selSet) trace.selectedpoints = sel.map(function (i) { return P.pos[i]; });
        var anns = [];
        Object.keys(P.centers).forEach(function (lv) {
          if (hidden[lv]) return;
          anns.push({ x: P.centers[lv][0], y: P.centers[lv][1], text: lv, showarrow: false, font: { size: 11, color: c.text }, bgcolor: t === 'light' ? 'rgba(252,252,251,0.7)' : 'rgba(26,26,25,0.7)' });
        });
        var L = baseLayout(t, { xlab: d.reduction.toUpperCase() + ' 1', ylab: d.reduction.toUpperCase() + ' 2', legend: false, annotations: anns, dragmode: 'lasso',
          xaxis: { showgrid: false, zeroline: false }, yaxis: { showgrid: false, zeroline: false, scaleanchor: 'x' }, margin: { l: 48, r: 16, t: 8, b: 40 } });
        L.uirevision = d.id;
        return { data: [trace], layout: L };
      });
      if (!plotDiv._bound) {
        plotDiv._bound = true;
        plotDiv.on('plotly_selected', function (ev) {
          if (!ev || !ev.points) return;
          sel = ev.points.map(function (p) { return p.customdata; }); selSet = sel.length ? new Set(sel) : null; renderSel();
        });
        plotDiv.on('plotly_deselect', function () { sel = []; selSet = null; renderSel(); });
      }
      renderSel();
    }
    var drawSoon = rafThrottle(draw);
    var lastLevels = [];
    function renderLegend(levels) {
      lastLevels = levels;
      clear(legend);
      if (!levels.length) return;
      var cols = catColors(levels.length, resolved());
      levels.forEach(function (lv, k) {
        legend.appendChild(h('button', { type: 'button', 'aria-pressed': hidden[lv] ? 'false' : 'true', title: 'Toggle ' + lv,
          onclick: function () { hidden[lv] = !hidden[lv]; draw(); } }, h('span', { class: 'swatch', style: { background: cols[k] } }), lv));
      });
    }
    var groupSel = h('select', { 'aria-label': 'Break selection down by' });
    var groupChoice = null;
    groupSel.addEventListener('change', function () { groupChoice = groupSel.value; renderSel(); });
    function renderSel() {
      var d = ds(), n = d.cells.length;
      clear(selBody);
      if (!sel.length) {
        selBody.appendChild(h('p', { class: 'note', text: 'Use the lasso or box tool in the plot toolbar to select cells. Double-click the plot to clear. Click a legend entry to hide a category.' }));
        return;
      }
      var isSel = new Uint8Array(n); sel.forEach(function (i) { isSel[i] = 1; });
      selBody.appendChild(h('p', null, h('strong', { text: fmtInt(sel.length) }), ' cells selected (' + (100 * sel.length / n).toFixed(1) + '% of ' + fmtInt(n) + ')'));
      var cats = catCols(d);
      if (cats.length) {
        clear(groupSel); cats.forEach(function (c) { groupSel.appendChild(h('option', { value: c, text: c })); });
        groupSel.value = groupChoice && cats.indexOf(groupChoice) >= 0 ? groupChoice : preferredGroup(d, cats);
        selBody.appendChild(h('label', { class: 'field' }, 'Break down by', groupSel));
        var m = d.meta[groupSel.value], cnt = m.levels.map(function () { return 0; }), tot = cnt.slice();
        for (var i = 0; i < n; i++) { tot[m.codes[i]]++; if (isSel[i]) cnt[m.codes[i]]++; }
        var cols = catColors(m.levels.length, resolved()), bars = h('div', { class: 'bars' });
        m.levels.forEach(function (lv, k) {
          if (!cnt[k]) return;
          bars.appendChild(h('div', { class: 'bar-row', title: cnt[k] + ' of ' + tot[k] + ' ' + lv + ' cells' }, h('span', { class: 'lab', text: lv }),
            h('span', { class: 'bar-track' }, h('span', { class: 'bar-fill', style: { width: (100 * cnt[k] / sel.length) + '%', background: cols[k], display: 'block' } })),
            h('span', { class: 'small', text: (100 * cnt[k] / sel.length).toFixed(1) + '%' })));
        });
        selBody.appendChild(bars);
      }
      // Quick look: panel genes, selection vs rest.
      var gl = geneNames(d), rows = gl.map(function (g) {
        var v = geneVec(d, g), a = 0, b = 0;
        for (var i = 0; i < n; i++) { if (isSel[i]) a += v[i]; else b += v[i]; }
        var ma = a / sel.length, mb = n > sel.length ? b / (n - sel.length) : 0;
        return { g: g, a: ma, b: mb, d: ma - mb };
      }).sort(function (p, q) { return q.d - p.d; });
      if (rows.length && n > sel.length) {
        var tb = h('table', { class: 'tbl' }, h('thead', null, h('tr', null, ['Gene', 'Selected', 'Rest', 'Δ'].map(function (x) { return h('th', { text: x }); }))),
          h('tbody', null, rows.slice(0, 8).map(function (r) { return h('tr', null, h('td', { text: r.g }), h('td', { class: 'num', text: fmt(r.a) }), h('td', { class: 'num', text: fmt(r.b) }), h('td', { class: 'num', text: fmt(r.d) })); })));
        selBody.appendChild(h('div', null, h('h3', { text: 'Higher in selection' }), h('div', { class: 'tbl-wrap', style: { maxHeight: '240px' } }, tb),
          h('p', { class: 'note', text: 'Quick look over the ' + gl.length + ' genes embedded in this report (mean log-normalised expression). Run full differential expression on a selection in the explorer app.' })));
      }
      var mk = function () { return sel.map(function (i) { return d.cells[i]; }); };
      selBody.appendChild(h('div', { class: 'toolbar', style: { marginBottom: 0 } },
        h('button', { class: 'btn primary', type: 'button', text: 'Download barcodes', onclick: function () {
          var gcol = cats.length ? groupSel.value : null, g = gcol ? labelsOf(d, gcol) : null;
          var cols = gcol ? ['barcode', gcol] : ['barcode'];
          download('sparks_selection.csv', toCSV(cols, function (c, k) { return c === 'barcode' ? d.cells[sel[k]] : g[sel[k]]; }, sel.length));
        } }),
        h('button', { class: 'btn', type: 'button', text: 'Copy', onclick: function (e) {
          var btn = e.currentTarget, txt = mk().join('\n');
          var done = function () { btn.textContent = 'Copied'; setTimeout(function () { btn.textContent = 'Copy'; }, 1200); };
          if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(txt).then(done, function () { download('sparks_selection.txt', txt, 'text/plain'); });
          else download('sparks_selection.txt', txt, 'text/plain');
        } }),
        h('button', { class: 'btn', type: 'button', text: 'Clear', onclick: function () { sel = []; selSet = null; draw(); } })));
      selBody.appendChild(h('p', { class: 'note', text: 'Import the CSV in the explorer app (Selection → Import barcodes) to run DE, re-cluster or score a signature on exactly these cells.' }));
    }
    colorSel.addEventListener('change', function () { hidden = {}; draw(); });
    themeListeners.push(function () { if (document.body.contains(plotDiv)) { renderSel(); renderLegend(lastLevels); } });
    size.addEventListener('input', drawSoon); opac.addEventListener('input', drawSoon);
    onDataset(function () { fillColor(); draw(); });
    root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Colour by', colorSel),
      h('label', { class: 'field' }, 'Point size', size), h('label', { class: 'field' }, 'Opacity', opac)));
    root.appendChild(h('div', { class: 'grid g-main' }, card('Embedding', 'lasso-select cells to inspect them', h('div', null, plotDiv, legend)), h('section', { class: 'card' }, h('header', null, h('h3', { text: 'Selection' })), selBody)));
    fillColor(); draw();
  });

  // ---- Genes
  tab('genes', 'Genes', function () { return DATA.datasets.some(function (d) { return d.n_genes_embedded > 0; }); }, function (root) {
    var geneIn = h('input', { type: 'text', list: 'sparks-genes', 'aria-label': 'Gene', autocomplete: 'off' });
    var dl = h('datalist', { id: 'sparks-genes' });
    var groupSel = h('select', { 'aria-label': 'Group by' });
    var multi = h('input', { type: 'text', 'aria-label': 'Gene list', style: { minWidth: '320px' } });
    var vBox = plotBox(), dBox = plotBox();
    var dsSel = datasetSelect();
    function fill() {
      var d = ds(), g = geneNames(d);
      clear(dl); g.forEach(function (x) { dl.appendChild(h('option', { value: x })); });
      clear(groupSel); catCols(d).forEach(function (c) { groupSel.appendChild(h('option', { value: c, text: c })); });
      groupSel.value = preferredGroup(d) || '';
      geneIn.value = (d.default_genes && d.default_genes[0]) || g[0] || '';
      multi.value = (d.default_genes && d.default_genes.length ? d.default_genes : g).slice(0, 10).join(', ');
    }
    function draw() {
      var d = ds(), g = geneIn.value.trim(), grp = groupSel.value;
      if (!grp || !d.meta[grp]) { showNote(vBox, 'This dataset has no categorical metadata to group by.'); showNote(dBox, ''); return; }
      if (d.genes[g]) violin(vBox, d, geneVec(d, g), grp, g + ' by ' + grp);
      else showNote(vBox, 'Gene "' + g + '" is not embedded in this report (' + geneNames(d).length + ' genes available). Use the explorer app for any gene.');
      var list = uniq(multi.value.split(/[\s,;]+/).filter(Boolean).filter(function (x) { return d.genes[x]; }));
      if (list.length) dotplot(dBox, d, list, grp); else showNote(dBox, 'No matching genes for the dot plot.');
    }
    [geneIn, groupSel, multi].forEach(function (el) { el.addEventListener('change', draw); });
    onDataset(function () { fill(); draw(); });
    root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Gene', geneIn, dl), h('label', { class: 'field' }, 'Group by', groupSel),
      h('label', { class: 'field' }, 'Dot plot genes (comma-separated)', multi)));
    root.appendChild(h('div', { class: 'grid g2' }, card('Expression distribution', 'violin, mean marked', vBox), card('Dot plot', 'size = % of cells expressing; colour = scaled mean', dBox)));
    fill(); draw();
  });

  // ---- Markers & DE
  var DE_KINDS = ['markers_all', 'deg_sig', 'cluster_markers', 'pseudobulk_deg', 'explorer'];
  function deSets() {
    return tablesOfKind(DE_KINDS).filter(function (t) { return t.kind !== 'explorer' || t.cols.indexOf('avg_log2FC') >= 0; });
  }
  function summaryFor(T) { return DATA.tables.filter(function (s) { return s.kind === 'deg_summary' && s.summary_of === T.id; })[0]; }
  tab('de', 'Markers & DE', function () { return deSets().length > 0 || figuresOf('de').length > 0 || tablesOfKind(['deg_counts']).length > 0 || ranOff('deg'); }, function (root) {
    if (ranOff('deg')) root.appendChild(h('div', { style: { marginBottom: '12px' } }, notice('Between-condition DEG was disabled (deg$run = false). Cluster markers, if present, are shown below.')));
    var sets = deSets();
    // Prefer the full table over its significant-only subset for the same analysis.
    var RANK = { markers_all: 3, deg_sig: 2, cluster_markers: 1, pseudobulk_deg: 3, explorer: 3 };
    var keyOf = function (t) { return t.kind === 'explorer' || t.kind === 'pseudobulk_deg' ? t.id : [t.comparison, t.unit, t.grouping].join('|'); };
    var best = {};
    sets.forEach(function (t) { var k = keyOf(t); best[k] = Math.max(best[k] || 0, RANK[t.kind]); });
    var taken = {};
    sets = sets.filter(function (t) { var k = keyOf(t); if (taken[k] || RANK[t.kind] !== best[k]) return false; taken[k] = true; return true; });
    if (!sets.length) {
      if (figuresOf('de').length) root.appendChild(h('div', { class: 'grid g2' }, figureCards(figuresOf('de'))));
      else if (!ranOff('deg')) root.appendChild(notice('No cell-level DEG tables are available. This happens when comparisons had fewer than two conditions, groups had too few cells, or no test produced results.'));
      return;
    }
    var setSel = h('select', { 'aria-label': 'Result set', style: { maxWidth: 'min(520px, 100%)' } }, sets.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
    var clusterSel = h('select', { 'aria-label': 'Group' });
    var lfcIn = h('input', { type: 'number', value: 0.25, step: 0.05, min: 0, style: { minWidth: '80px', width: '90px' }, 'aria-label': 'log2FC threshold' });
    var pIn = h('input', { type: 'number', value: 0.05, step: 0.01, min: 0, max: 1, style: { minWidth: '80px', width: '90px' }, 'aria-label': 'Adjusted p threshold' });
    var topIn = h('input', { type: 'number', value: window.innerWidth < 600 ? 4 : 8, step: 1, min: 0, max: 30, style: { minWidth: '70px', width: '80px' }, 'aria-label': 'Labels' });
    var vBox = plotBox('tall'), tblHost = h('div'), sumHost = h('div'), info = h('div', { class: 'note' });
    var cur = null, view = null;
    var draw = function () {
      var T = tdata(sets[+setSel.value]), lc = pickCol(T, ['avg_log2FC', 'logFC', 'avg_logFC', 'log2FoldChange']), pc = pickCol(T, ['p_val_adj', 'padj', 'FDR', 'adj.P.Val', 'p_val', 'PValue']);
      var gc = pickCol(T, ['gene', 'Gene']), cc = pickCol(T, ['cluster', 'group', 'Group', 'contrast']);
      if (cur !== T) {
        cur = T; clear(clusterSel);
        if (cc) uniq(T.data[cc].map(String)).sort(function (a, b) { return a.localeCompare(b, undefined, { numeric: true }); }).forEach(function (g) { clusterSel.appendChild(h('option', { value: g, text: g })); });
        clusterSel.disabled = !cc;
        clear(tblHost); view = null;
        clear(sumHost);
        var S = summaryFor(T);
        if (S) tableView(sumHost, S, { filename: S.id + '_summary', pageSize: 15, note: 'Counts use the adjusted p-value column of the full table on disk with a 0.05 cutoff; no test is re-run.' });
        else sumHost.appendChild(emptyNote('No summary for this table.'));
      }
      var lfcCut = Math.max(0, +lfcIn.value || 0), pCut = Math.min(1, Math.max(0, +pIn.value || 0));
      var idx = range(T.n).filter(function (i) { return !cc || String(T.data[cc][i]) === clusterSel.value; });
      if (lc && pc) {
        var V = { lfc: idx.map(function (i) { return T.data[lc][i]; }), p: idx.map(function (i) { return T.data[pc][i]; }), genes: idx.map(function (i) { return gc ? T.data[gc][i] : String(i); }) };
        volcano(vBox, V, { lfcCut: lfcCut, pCut: pCut, labelTop: Math.max(0, +topIn.value || 0), title: (cc ? clusterSel.value + ' — ' : '') + T.label });
      } else showNote(vBox, 'This table has no fold-change / p-value columns; see the table below.');
      var nsig = lc && pc ? idx.filter(function (i) { var p = T.data[pc][i], l = T.data[lc][i]; return isNum(p) && isNum(l) && p < pCut && Math.abs(l) >= lfcCut; }).length : 0;
      info.textContent = fmtInt(idx.length) + ' genes in this group, ' + fmtInt(nsig) + ' pass the thresholds' + (T.truncated ? ' (table truncated in the report to the strongest rows; see the results folder for the full file)' : '') + '.';
      if (!view) view = tableView(tblHost, T, { filename: T.id, sort: pc ? [pc, 'asc'] : null });
      view.setPredicate(function (i) { return !cc || String(T.data[cc][i]) === clusterSel.value; });
    };
    [setSel, clusterSel, lfcIn, pIn, topIn].forEach(function (el) { el.addEventListener('change', draw); });
    root.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field' }, 'Result set', setSel), h('label', { class: 'field' }, 'Group', clusterSel),
      h('label', { class: 'field' }, '|log2FC| ≥', lfcIn), h('label', { class: 'field' }, 'adj. p <', pIn), h('label', { class: 'field' }, 'Labels / side', topIn)));
    root.appendChild(h('div', { class: 'grid g-main' }, card('Volcano', 'thresholds update live', h('div', null, vBox, info)),
      h('div', { class: 'stack' }, card('Per-group summary', 'all genes on disk', sumHost), figureCards(figuresOf('de'), 'short'))));
    root.appendChild(h('div', { class: 'mt' }, card('Results table', null, tblHost)));
    draw();
  });

  // ---- Enrichment
  tab('enrichment', 'Enrichment', function () {
    return figuresOf('enrichment').length > 0 || tablesOfKind(['gsea', 'escape']).length > 0 || !!(DATA.notes && DATA.notes.gsea_status) ||
      RUN.gsea === true || RUN.escape === true || DATA.datasets.some(function (d) { return d.facts && d.facts.escape_assays && d.facts.escape_assays.length; });
  }, function (root) {
    var notes = [];
    var st = DATA.notes && DATA.notes.gsea_status;
    if (st) {
      var hs = h('div');
      tableView(hs, makeTable(['file', 'message'], { file: st.file, message: st.text }), { filename: 'sparks_gsea_status' });
      notes.push(card('GSEA diagnostics', 'skip reasons written by the run', hs));
    }
    if (!tablesOfKind(['gsea']).length) {
      notes.push(notice(RUN.gsea === false ? 'Ranked GSEA was disabled in the configuration (gsea$run = false).' :
        'No ranked GSEA tables for MSigDB collection "' + ((DATA.settings && DATA.settings.gsea_collection) || 'H') + '". GSEA needs edgeR pseudobulk ranks, two conditions with biological replicates and gene sets; see the diagnostics and the Run tab.'));
    }
    var assays = uniq([].concat.apply([], DATA.datasets.map(function (d) { return (d.facts && d.facts.escape_assays) || []; })));
    if (assays.length) notes.push(notice('Per-cell pathway-score assays in the saved objects: ' + assays.join(', ') + '.'));
    else if (!tablesOfKind(['escape']).length) notes.push(notice(RUN.escape === false ? 'Per-cell pathway scoring was disabled (escape$run = false).' : 'No escape assay or score table is available; check the escape / msigdb packages and the Run tab.'));
    if (notes.length) root.appendChild(h('div', { class: 'stack', style: { marginBottom: '12px' } }, notes));
    root.appendChild(h('div', { class: 'grid g2' }, figureCards(figuresOf('enrichment'), 'tall')));
    var sets = tablesOfKind(['gsea', 'escape']);
    if (sets.length) {
      var host = h('div'), sel = h('select', { 'aria-label': 'Table', style: { maxWidth: 'min(520px, 100%)' } }, sets.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
      var show = function () { clear(host); var T = sets[+sel.value]; tableView(host, T, { filename: T.id, sort: T.cols.indexOf('padj') >= 0 ? ['padj', 'asc'] : null }); };
      sel.addEventListener('change', show);
      root.appendChild(h('div', { class: 'mt' }, card('Enrichment tables', null, h('div', null, h('label', { class: 'field' }, 'Table', sel), host))));
      show();
    }
    var imgs = DATA.images.filter(function (i) { return /^(Escape_Heatmap_|GSEA_Projection_)/.test(i.label); });
    if (imgs.length) root.appendChild(h('div', { class: 'mt' }, card('Enrichment figures', 'GSEA projections show group-level NES repeated across cells, not per-cell scores', gallery(imgs, false))));
  });

  // ---- Composition
  tab('composition', 'Composition', function () { return DATA.datasets.length > 0 || figuresOf('composition').length > 0 || tablesOfKind(['proportions', 'prop_test']).length > 0; }, function (root) {
    if (DATA.datasets.length) {
      var dsSel = datasetSelect(), box = plotBox('tall');
      var stackSel = h('select', { 'aria-label': 'Stack by' }), xSel = h('select', { 'aria-label': 'Compare across' });
      var modeSel = h('select', { 'aria-label': 'Scale' }, h('option', { value: 'percent', text: 'Proportions (%)' }), h('option', { value: 'count', text: 'Cell counts' }));
      var sizeSel = h('select', { 'aria-label': 'Column' }), sizeHost = h('div');
      var fill = function () {
        var d = ds(), cats = catCols(d);
        [stackSel, xSel].forEach(function (s) { clear(s); cats.forEach(function (c) { s.appendChild(h('option', { value: c, text: c })); }); });
        stackSel.value = preferredGroup(d, cats) || '';
        var alt = cats.filter(function (c) { return c !== stackSel.value; });
        xSel.value = alt.indexOf('condition') >= 0 ? 'condition' : (alt.indexOf('orig.ident') >= 0 ? 'orig.ident' : (alt[0] || cats[0] || ''));
        var cc = (d.facts && d.facts.counts) || {};
        clear(sizeSel); Object.keys(cc).forEach(function (c) { sizeSel.appendChild(h('option', { value: c, text: c })); });
        sizeSel.value = Object.keys(cc).indexOf(stackSel.value) >= 0 ? stackSel.value : (Object.keys(cc)[0] || '');
      };
      var draw = function () {
        var d = ds(), a = d.meta[stackSel.value], b = d.meta[xSel.value];
        if (!a || !b) { showNote(box, 'Needs two categorical columns.'); return; }
        var cnt = a.levels.map(function () { return b.levels.map(function () { return 0; }); });
        for (var i = 0; i < d.cells.length; i++) cnt[a.codes[i]][b.codes[i]]++;
        var tot = b.levels.map(function (_, j) { return cnt.reduce(function (s, r) { return s + r[j]; }, 0); });
        var pct = modeSel.value === 'percent';
        stackedBars(box, b.levels, a.levels.map(function (lv, k) {
          return { name: lv, values: cnt[k].map(function (v, j) { return pct ? (tot[j] ? 100 * v / tot[j] : 0) : v; }), counts: cnt[k] };
        }), { percent: pct, xlab: xSel.value, title: stackSel.value + ' composition by ' + xSel.value });
      };
      var sizes = function () {
        var d = ds(), c = d.facts && d.facts.counts && d.facts.counts[sizeSel.value];
        clear(sizeHost);
        if (!c) { sizeHost.appendChild(emptyNote('No categorical columns.')); return; }
        tableView(sizeHost, countsTable(c, sizeSel.value), { filename: 'sparks_sizes_' + sizeSel.value, sort: ['cells', 'desc'] });
      };
      [stackSel, xSel, modeSel].forEach(function (s) { s.addEventListener('change', draw); });
      sizeSel.addEventListener('change', sizes);
      onDataset(function () { fill(); draw(); sizes(); });
      root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Stack by', stackSel), h('label', { class: 'field' }, 'Compare across', xSel), h('label', { class: 'field' }, 'Scale', modeSel)));
      root.appendChild(h('div', { class: 'grid g-main' }, card('Composition', 'embedded cells (down-sampled above the report limit)', box),
        card('Group sizes', 'exact counts over all cells', h('div', null, h('label', { class: 'field' }, 'Column', sizeSel), sizeHost))));
      fill(); draw(); sizes();
    }
    var figs = figuresOf('composition');
    if (figs.length) root.appendChild(h('div', { class: 'grid g2 mt' }, figureCards(figs)));
    var pt = tablesOfKind(['prop_test']);
    if (pt.length) {
      var host = h('div'), sel = h('select', { 'aria-label': 'Proportion test' }, pt.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
      var show = function () { clear(host); var T = pt[+sel.value]; tableView(host, T, { filename: T.id }); };
      sel.addEventListener('change', show);
      root.appendChild(h('div', { class: 'mt' }, card('Proportion tests (scProportionTest)', null, h('div', null, h('label', { class: 'field' }, 'Test', sel), host))));
      show();
    }
  });

  // ---- All tables
  tab('tables', 'Tables', function () { return DATA.tables.length > 0; }, function (root) {
    var kinds = uniq(DATA.tables.map(function (t) { return t.kind; }));
    var kindSel = h('select', { 'aria-label': 'Kind' }, kinds.map(function (k) { return h('option', { value: k, text: (DATA.kind_labels[k] || k) + ' (' + DATA.tables.filter(function (t) { return t.kind === k; }).length + ')' }); }));
    var tSel = h('select', { 'aria-label': 'Table', style: { maxWidth: 'min(520px, 100%)' } }), host = h('div');
    var list = function () {
      clear(tSel);
      DATA.tables.forEach(function (t, i) { if (t.kind === kindSel.value) tSel.appendChild(h('option', { value: i, text: t.label })); });
      show();
    };
    var show = function () {
      clear(host);
      if (tSel.value === '') { host.appendChild(emptyNote('No tables.')); return; }
      var T = DATA.tables[+tSel.value];
      tableView(host, T, { filename: T.id });
    };
    kindSel.addEventListener('change', list); tSel.addEventListener('change', show);
    root.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field' }, 'Kind', kindSel), h('label', { class: 'field' }, 'Table', tSel)));
    root.appendChild(card('Result table', 'every table the pipeline wrote (large ones are truncated here; files on disk are complete)', host));
    list();
  });

  // ---- Gallery
  var lightbox = null, io = null;
  function lazyImg(img, im) {
    img.setAttribute('data-pending', '');
    var load = function () { img.src = imgSrc(im); img.removeAttribute('data-pending'); };
    if (!('IntersectionObserver' in window)) { load(); return img; }
    if (!io) io = new IntersectionObserver(function (es) { es.forEach(function (e) { if (e.isIntersecting) { io.unobserve(e.target); e.target._load(); } }); }, { rootMargin: '300px' });
    img._load = load; io.observe(img);
    return img;
  }
  function gallery(images, withFilters) {
    var host = h('div'), grid = h('div', { class: 'gallery' });
    var q = h('input', { type: 'search', placeholder: 'Filter figures…', 'aria-label': 'Filter figures' });
    var cat = h('select', { 'aria-label': 'Category' }, h('option', { value: '', text: 'All categories' }), uniq(images.map(function (i) { return i.category || ''; })).filter(Boolean).map(function (c) { return h('option', { value: c, text: c }); }));
    var grp = h('select', { 'aria-label': 'Analysis' }, h('option', { value: '', text: 'All analyses' }), uniq(images.map(function (i) { return i.group || ''; })).filter(Boolean).map(function (c) { return h('option', { value: c, text: c }); }));
    var count = h('span', { class: 'muted small' });
    var limit = 120;
    var draw = function () {
      clear(grid);
      var f = q.value.toLowerCase(), shown = images.filter(function (im) {
        return (!cat.value || im.category === cat.value) && (!grp.value || im.group === grp.value) && (!f || (im.label + ' ' + im.rel).toLowerCase().indexOf(f) >= 0);
      });
      shown.slice(0, limit).forEach(function (im) {
        var body = im.embedded ? h('button', { type: 'button', 'aria-label': 'Enlarge ' + im.label, onclick: function () { openLightbox(im); } }, lazyImg(h('img', { alt: im.label, decoding: 'async' }), im))
          : h('div', { class: 'nopreview' }, im.href ? h('a', { href: im.href, text: 'Open file (' + bytes(im.bytes) + ')', target: '_blank', rel: 'noopener' }) : h('span', { text: 'Not embedded (' + bytes(im.bytes) + ')' }));
        grid.appendChild(h('figure', { class: 'shot' }, body, h('figcaption', { text: im.label + (im.group ? ' — ' + im.group : '') })));
      });
      count.textContent = shown.length + ' figures' + (shown.length > limit ? ' (showing ' + limit + ')' : '');
      if (shown.length > limit) grid.appendChild(h('button', { class: 'btn', type: 'button', text: 'Show more', onclick: function () { limit += 120; draw(); } }));
    };
    [q, cat, grp].forEach(function (el) { el.addEventListener('input', function () { limit = 120; draw(); }); });
    if (withFilters !== false) host.appendChild(h('div', { class: 'toolbar' }, q, h('label', { class: 'field' }, 'Category', cat), h('label', { class: 'field' }, 'Analysis', grp), count));
    host.appendChild(grid); draw();
    return host;
  }
  function openLightbox(im) {
    if (!lightbox) { lightbox = h('dialog', { class: 'lightbox', 'aria-label': 'Figure preview' }); lightbox.addEventListener('click', function () { lightbox.close(); }); document.body.appendChild(lightbox); }
    clear(lightbox).appendChild(h('img', { src: imgSrc(im), alt: im.label }));
    if (lightbox.showModal) lightbox.showModal();
  }
  tab('figures', 'Figures', function () { return DATA.images.length > 0; }, function (root) {
    var embedded = DATA.images.filter(function (i) { return i.embedded; }).length;
    root.appendChild(card('Pipeline figures', embedded + ' embedded of ' + DATA.images.length + (embedded < DATA.images.length ? ' (the rest are linked)' : ''), gallery(DATA.images, true)));
  });

  // ---- Files
  tab('files', 'Files', function () { return DATA.files && DATA.files.rel && DATA.files.rel.length > 0; }, function (root) {
    var F = DATA.files;
    var T = makeTable(['file', 'kind', 'category', 'size'], { file: F.rel, kind: F.kind, category: F.category, size: F.size });
    var host = h('div');
    tableView(host, T, { filename: 'sparks_files', pageSize: 50, render: {
      file: function (v, i) { return F.href[i] ? h('a', { href: F.href[i], target: '_blank', rel: 'noopener', text: v }) : document.createTextNode(v); },
      size: function (v) { return bytes(v); } } });
    root.appendChild(card('Result files', 'links work while this report sits next to the results folder', host));
  });

  // ---- Run
  tab('run', 'Run', function () { return true; }, function (root) {
    var R = DATA.run, blocks = [];
    var D = R.diagnostics;
    if (D && D.message && D.message.length) {
      var hd = h('div');
      tableView(hd, makeTable(['type', 'message'], { type: D.type, message: D.message }), { filename: 'sparks_diagnostics', pageSize: 20 });
      blocks.push(card('Pipeline diagnostics (' + D.message.length + ')', 'warnings, errors and skipped analyses from pipeline.log', hd));
    } else blocks.push(card('Pipeline diagnostics', null, notice(R.has_log ? 'No warnings, errors or skipped analyses were logged.' : 'pipeline.log is not available.')));
    if (R.warnings && R.warnings.length) blocks.push(card('Warnings (' + R.warnings.length + ')', null, h('pre', { text: R.warnings.join('\n') })));
    if (R.config_yaml) blocks.push(card('Configuration used', null, h('pre', { text: R.config_yaml })));
    if (R.log_tail && R.log_tail.length) blocks.push(card('Pipeline log (tail)', null, h('pre', { text: R.log_tail.join('\n') })));
    if (R.session_info && R.session_info.length) blocks.push(card('Session info', null, h('pre', { text: R.session_info.join('\n') })));
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ------------------------------------------------------------------ shell
  var panels = {}, tabButtons = {}, active = null;
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && panels[x.id]; })[0] || TABS.filter(function (x) { return panels[x.id]; })[0];
    if (!t) return;
    id = t.id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) {
      t.built = true;
      try { t.build(panels[id]); } catch (e) { panels[id].appendChild(h('div', { class: 'empty', text: 'This tab could not be built: ' + e.message })); if (window.console) console.error(e); }
    } else refreshVisible(panels[id]);
    // Newly built plots in this tab need a layout pass once visible.
    setTimeout(function () { refreshVisible(panels[id]); }, 60);
    active = id;
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) { /* file:// may refuse */ } }
  }
  function init() {
    applyThemeAttr();
    document.title = DATA.meta.title;
    var tabsEl = $('#tabs'), main = $('#panels'), boot = $('#boot');
    if (boot) boot.remove();
    if (!window.Plotly) { main.appendChild(emptyNote('plotly.js failed to load; figures are unavailable.')); }
    var avail = TABS.filter(function (t) { try { return t.available(); } catch (e) { return false; } });
    avail.forEach(function (t, i) {
      var b = h('button', { class: 'tab', role: 'tab', id: 'tab-' + t.id, 'aria-controls': 'panel-' + t.id, 'aria-selected': 'false', tabindex: '-1', type: 'button', text: t.label,
        onclick: function () { activate(t.id, true); },
        onkeydown: function (e) {
          var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
          if (!k) return; e.preventDefault();
          var nx = avail[(i + k + avail.length) % avail.length]; activate(nx.id, true); tabButtons[nx.id].focus();
        } });
      tabButtons[t.id] = b; tabsEl.appendChild(b);
      panels[t.id] = h('section', { class: 'panel', role: 'tabpanel', id: 'panel-' + t.id, 'aria-labelledby': 'tab-' + t.id, hidden: true });
      main.appendChild(panels[t.id]);
    });
    var btn = $('#theme-btn');
    var label = function () { btn.textContent = 'Theme: ' + themePref; };
    label();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('sparks-theme', themePref); applyThemeAttr(); label(); rerenderAll();
    });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      var onChange = function () { if (themePref === 'auto') rerenderAll(); };
      if (mq.addEventListener) mq.addEventListener('change', onChange); else if (mq.addListener) mq.addListener(onChange);
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
)---"
