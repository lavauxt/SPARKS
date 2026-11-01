/* SPARKS unified report -- client-side app.
   Reads the JSON blob in #sparks-data, builds one lazily-rendered panel per tab and
   draws every figure with the inlined plotly.js. No network access is needed. */
(function () {
  'use strict';

  // ------------------------------------------------------------------ data & utils
  var DATA = JSON.parse(document.getElementById('sparks-data').textContent);
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

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
  function uniq(a) { return a.filter(function (v, i) { return a.indexOf(v) === i; }); }
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
    var blob = new Blob([text], { type: type || 'text/csv;charset=utf-8' });
    var a = h('a', { href: URL.createObjectURL(blob), download: name });
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
  function catColors(n, theme) {
    var base = CAT[theme], out = [], tints = [0, 0.38, 0.62];
    for (var i = 0; i < n; i++) {
      out.push(tint(base[i % 8], tints[Math.min(Math.floor(i / 8), 2)], theme === 'light' ? '#ffffff' : '#000000'));
    }
    return out;
  }
  var themeRegex = null, themeMap = null;
  function mapTheme(obj, theme) {
    if (theme !== 'dark') return obj;
    if (!themeRegex) {
      themeMap = {};
      Object.keys(DATA.theme_map || {}).forEach(function (k) { themeMap[k.toLowerCase()] = DATA.theme_map[k]; });
      themeRegex = new RegExp(Object.keys(themeMap).join('|'), 'gi');
    }
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
      delete s.layout.width; delete s.layout.height; delete s.layout.title;   // the card header carries the title
      if (s.layout.margin) s.layout.margin.t = Math.min(s.layout.margin.t || 16, 16);
      s.config = Object.assign({}, s.config || {});
      return s;
    });
  }
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
    themeListeners.forEach(function (cb) { cb(); });
  }
  var themeListeners = [];
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

  // ------------------------------------------------------------------ tables
  function colType(T, c) {
    var arr = T.data[c];
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // T = {cols, data: {col: [..]}, n}. opts: filename, pageSize, render: {col: fn(value,i)->node}, sort: [col, 'asc'|'desc']
  function tableView(host, T, opts) {
    opts = opts || {};
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
    host.appendChild(h('div', { class: 'tbl-bar' }, search, count, dl, pager));
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
          if (x === null || x === undefined) return 1;
          if (y === null || y === undefined) return -1;
          return (num ? x - y : String(x).localeCompare(String(y))) * dir;
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
      clear(tbody);
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      for (var k = start; k < end; k++) {
        var i = idx[k], row = h('tr');
        T.cols.forEach(function (c) {
          var v = T.data[c][i], r = opts.render && opts.render[c];
          row.appendChild(h('td', { class: types[c] === 'number' ? 'num' : '' }, r ? r(v, i) : fmt(v)));
        });
        tbody.appendChild(row);
      }
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows' + (T.truncated ? ' (truncated from ' + fmtInt(T.n_total) + ')' : '');
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '‹', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '›', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return { setPredicate: function (fn) { pred = fn; page = 0; recompute(); }, refresh: recompute };
  }
  function tableById(id) { return DATA.tables.filter(function (t) { return t.id === id; })[0]; }

  // ------------------------------------------------------------------ datasets
  var state = { ds: 0, listeners: [] };
  function ds() { return DATA.datasets[state.ds]; }
  function onDataset(cb) { state.listeners.push(cb); }
  function datasetSelect() {
    if (DATA.datasets.length < 2) return null;
    var sel = h('select', { 'aria-label': 'Dataset' }, DATA.datasets.map(function (d, i) { return h('option', { value: i, text: d.label + ' (' + fmtInt(d.cells.length) + ' cells)' }); }));
    sel.value = state.ds;
    sel.addEventListener('change', function () {
      state.ds = +sel.value;
      $$('select[data-role="dataset"]').forEach(function (o) { o.value = state.ds; });
      state.listeners.forEach(function (cb) { cb(); });
    });
    sel.setAttribute('data-role', 'dataset');
    return h('label', { class: 'field' }, 'Dataset', sel);
  }
  function catCols(d) { return Object.keys(d.meta).filter(function (k) { return d.meta[k].type === 'cat'; }); }
  function numCols(d) { return Object.keys(d.meta).filter(function (k) { return d.meta[k].type === 'num'; }); }
  function geneNames(d) { return Object.keys(d.genes); }
  function preferredGroup(d) {
    var cats = catCols(d);
    var pref = ['cell_type', 'celltype', 'SingleR_main', 'seurat_clusters', 'cluster', 'orig.ident'];
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
        if (l === null || p === null || !isFinite(l) || !isFinite(p)) continue;
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
      return {
        data: [mk('Not significant', ns, c.muted, 0.45), mk('Down', dn, DIV[4], 0.85), mk('Up', up, DIV[0], 0.85)],
        layout: baseLayout(t, {
          title: o.title, xlab: 'log2 fold change', ylab: '-log10 adjusted p', legend: true, annotations: anns,
          shapes: [
            { type: 'line', x0: -o.lfcCut, x1: -o.lfcCut, y0: 0, y1: 1, yref: 'paper', line: line },
            { type: 'line', x0: o.lfcCut, x1: o.lfcCut, y0: 0, y1: 1, yref: 'paper', line: line },
            { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: -Math.log10(o.pCut), y1: -Math.log10(o.pCut), line: line }]
        })
      };
    });
  }

  function violin(div, d, values, groupCol, title, ylab) {
    var labs = labelsOf(d, groupCol), levels = d.meta[groupCol].levels;
    plot(div, function (t) {
      var cols = catColors(levels.length, t);
      var by = levels.map(function () { return []; });
      for (var i = 0; i < values.length; i++) by[d.meta[groupCol].codes[i]].push(values[i]);
      var data = levels.map(function (lv, k) {
        var v = by[k];
        if (v.length > 3000) { var step = v.length / 3000; v = range(3000).map(function (j) { return v[Math.floor(j * step)]; }); }
        return { type: 'violin', y: v, x: v.map(function () { return lv; }), name: lv, points: false, spanmode: 'hard', scalemode: 'width',
          width: 0.85, fillcolor: cols[k], line: { color: cols[k], width: 1 }, meanline: { visible: true }, opacity: 0.9, hoverinfo: 'y+name' };
      });
      return { data: data, layout: baseLayout(t, { title: title, xlab: groupCol, ylab: ylab || 'Expression (log-normalised)', legend: false, xaxis: { type: 'category' } }) };
    });
    return labs;
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
    plot(div, function (t) {
      var lim = Math.max(1, Math.max.apply(null, rows.map(function (r) { return Math.abs(r.z); })));
      return {
        data: [{ type: 'scatter', mode: 'markers', x: rows.map(function (r) { return r.g; }), y: rows.map(function (r) { return r.lv; }),
          text: rows.map(function (r) { return r.g + ' / ' + r.lv + '<br>mean ' + fmt(r.mean) + '<br>' + Math.round(r.pct) + '% of cells'; }), hoverinfo: 'text',
          marker: { size: rows.map(function (r) { return 4 + 0.22 * r.pct; }), sizemode: 'diameter', color: rows.map(function (r) { return r.z; }), cmin: -lim, cmax: lim,
            colorscale: SEQ.map(function (c, i) { return [i / (SEQ.length - 1), c]; }), colorbar: { title: { text: 'Scaled mean' }, thickness: 10, len: 0.6 },
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
        yaxis: o.percent ? { range: [0, 100], ticksuffix: '%' } : {} }) };
    });
  }

  // ------------------------------------------------------------------ tabs
  var TABS = [];
  function tab(id, label, available, build) { TABS.push({ id: id, label: label, available: available, build: build, built: false }); }
  function card(title, hint, body) {
    return h('section', { class: 'card' }, h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function tablesOfKind(kinds) { return DATA.tables.filter(function (t) { return kinds.indexOf(t.kind) >= 0; }); }
  function figuresOf(tabId) { return DATA.figures.filter(function (f) { return f.tab === tabId; }); }
  function figureCards(tabId, cls) {
    return figuresOf(tabId).map(function (f) {
      var box = plotBox(cls || '');
      var c = card(f.title, f.group || '', box);
      setTimeout(function () { specPlot(box, f.spec); }, 0);
      return c;
    });
  }
  function nav(id) { activate(id, true); }

  // ---- Overview
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var S = DATA.summary;
    var tiles = [['cells', 'Cells'], ['samples', 'Samples'], ['genes', 'Genes'], ['clusters', 'Clusters'], ['tables', 'Result tables'], ['files', 'Files']];
    root.appendChild(h('div', { class: 'tiles' }, tiles.filter(function (t) { return S[t[0]] !== null && S[t[0]] !== undefined; }).map(function (t) {
      return h('div', { class: 'tile' }, h('div', { class: 'v', text: fmtInt(S[t[0]]) }), h('div', { class: 'l', text: t[1] }));
    })));
    var left = [], right = [];
    var cps = figuresOf('overview');
    cps.forEach(function (f) { var b = plotBox('short'); left.push(card(f.title, f.group || '', b)); setTimeout(function () { specPlot(b, f.spec); }, 0); });
    var marks = tablesOfKind(['markers_all', 'deg_sig', 'cluster_markers'])[0];
    if (marks) {
      var gc = pickCol(marks, ['gene', 'Gene']), cc = pickCol(marks, ['cluster', 'group']), lc = pickCol(marks, ['avg_log2FC', 'logFC']);
      if (gc && cc && lc) {
        var groups = uniq(marks.data[cc]), rows = [];
        groups.forEach(function (g) {
          var ii = range(marks.n).filter(function (i) { return marks.data[cc][i] === g; }).sort(function (a, b) { return marks.data[lc][b] - marks.data[lc][a]; }).slice(0, 4);
          rows.push(g, ii.map(function (i) { return marks.data[gc][i]; }).join(', '));
        });
        var dl = h('dl', { class: 'kv' });
        for (var i = 0; i < rows.length; i += 2) { dl.appendChild(h('dt', { text: rows[i] })); dl.appendChild(h('dd', { text: rows[i + 1] })); }
        right.push(card('Top markers', marks.label, dl));
      }
    }
    var cmp = DATA.summary.comparisons || [];
    right.push(card('Run', null, h('dl', { class: 'kv' },
      h('dt', { text: 'Results folder' }), h('dd', { text: DATA.meta.results_dir_name }),
      h('dt', { text: 'Comparisons' }), h('dd', { text: cmp.length ? cmp.join(', ') : '—' }),
      h('dt', { text: 'Species' }), h('dd', { text: DATA.meta.species || '—' }),
      h('dt', { text: 'Generated' }), h('dd', { text: DATA.meta.generated }),
      h('dt', { text: 'SPARKS' }), h('dd', { text: DATA.meta.version }))));
    var jump = TABS.filter(function (t) { return t.id !== 'overview' && t.available(); }).map(function (t) {
      return h('button', { class: 'btn', type: 'button', text: t.label, onclick: function () { nav(t.id); } });
    });
    right.push(card('Jump to', null, h('div', { class: 'toolbar', style: { marginBottom: 0 } }, jump)));
    root.appendChild(h('div', { class: 'grid g2', style: { marginTop: '12px' } }, h('div', { class: 'stack' }, left), h('div', { class: 'stack' }, right)));
  });
  function pickCol(T, cands) { for (var i = 0; i < cands.length; i++) if (T.cols.indexOf(cands[i]) >= 0) return cands[i]; return null; }

  // ---- QC
  tab('qc', 'QC', function () { return figuresOf('qc').length > 0 || DATA.datasets.length > 0 || DATA.images.some(function (i) { return i.category === 'QC'; }); }, function (root) {
    var grid = h('div', { class: 'grid g2' }, figureCards('qc'));
    root.appendChild(grid);
    if (DATA.datasets.length) {
      var box = h('div', { class: 'stack', style: { marginTop: '12px' } });
      root.appendChild(card('Per-cell QC metrics by sample', 'from the merged object', box));
      var host = h('div');
      var draw = function () {
        var d = ds(), g = d.meta['orig.ident'] ? 'orig.ident' : preferredGroup(d);
        clear(host);
        var metrics = ['nFeature_RNA', 'nCount_RNA', 'percent.mt'].filter(function (m) { return d.meta[m] && d.meta[m].type === 'num'; });
        if (!metrics.length || !g) { host.appendChild(emptyNote('No QC metrics in this dataset.')); return; }
        var row = h('div', { class: 'grid g3' });
        metrics.forEach(function (m) {
          var b = plotBox('short'); row.appendChild(b);
          setTimeout(function () { violin(b, d, d.meta[m].values, g, m, m); }, 0);
        });
        host.appendChild(row);
      };
      var sel = datasetSelect();
      if (sel) box.appendChild(sel);
      box.appendChild(host); draw(); onDataset(draw);
    }
    var qcImgs = DATA.images.filter(function (i) { return i.category === 'QC'; });
    if (qcImgs.length) root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Pipeline QC figures', qcImgs.length + ' images', gallery(qcImgs, false))));
  });

  // ---- Embedding explorer
  tab('embedding', 'Embedding', function () { return DATA.datasets.length > 0; }, function (root) {
    var hidden = {}, sel = [], selSet = null, ordCache = null;
    var plotDiv = plotBox('tall'), legend = h('div', { class: 'legend' });
    var colorSel = h('select', { 'aria-label': 'Colour by' }), size = h('input', { type: 'range', min: 1, max: 10, step: 0.5, value: 4, 'aria-label': 'Point size' });
    var opac = h('input', { type: 'range', min: 0.1, max: 1, step: 0.05, value: 0.8, 'aria-label': 'Opacity' });
    var selTitle = h('h3', { text: 'Selection' }), selBody = h('div', { class: 'stack' });
    var dsSel = datasetSelect();

    function fillColor() {
      var d = ds(), cats = catCols(d), nums = numCols(d), genes = geneNames(d);
      clear(colorSel);
      var add = function (label, items) { if (items.length) colorSel.appendChild(h('optgroup', { label: label }, items.map(function (k) { return h('option', { value: k, text: k }); }))); };
      add('Categories', cats); add('QC & numeric', nums); add('Genes', genes);
      colorSel.value = preferredGroup(d) || (nums[0] || genes[0]);
      hidden = {}; sel = []; selSet = null;
    }
    function draw() {
      var d = ds(), key = colorSel.value, n = d.cells.length, theme = resolved();
      var isCat = d.meta[key] && d.meta[key].type === 'cat';
      var vals = isCat ? null : (d.meta[key] ? Float64Array.from(d.meta[key].values) : geneVec(d, key));
      var ord = range(n);
      if (!isCat) ord.sort(function (a, b) { return vals[a] - vals[b]; });
      var pos = new Int32Array(n); ord.forEach(function (c, p) { pos[c] = p; });
      var labels = isCat ? labelsOf(d, key) : null;
      var levels = isCat ? d.meta[key].levels : [];
      renderLegend(levels, key);
      plot(plotDiv, function (t) {
        var c = CHROME[t], cols = isCat ? catColors(levels.length, t) : null;
        var marker = { size: +size.value };
        if (isCat) {
          marker.color = ord.map(function (i) { return cols[d.meta[key].codes[i]]; });
          marker.opacity = ord.map(function (i) { return hidden[labels[i]] ? 0.04 : +opac.value; });
        } else {
          marker.color = ord.map(function (i) { return vals[i]; }); marker.opacity = +opac.value;
          marker.colorscale = SEQ.map(function (col, i) { return [i / (SEQ.length - 1), col]; });
          marker.colorbar = { title: { text: key }, thickness: 10, len: 0.6 };
        }
        var trace = {
          type: n > 8000 ? 'scattergl' : 'scatter', mode: 'markers',
          x: ord.map(function (i) { return d.x[i]; }), y: ord.map(function (i) { return d.y[i]; }),
          customdata: ord, marker: marker, hoverinfo: 'text',
          text: ord.map(function (i) { return d.cells[i] + '<br>' + (isCat ? labels[i] : key + ': ' + fmt(vals[i])); }),
          selected: { marker: { opacity: 1 } }, unselected: { marker: { opacity: 0.12 } }
        };
        if (selSet) trace.selectedpoints = sel.map(function (i) { return pos[i]; });
        var anns = [];
        if (isCat && levels.length > 1 && levels.length <= 40) {
          levels.forEach(function (lv) {
            if (hidden[lv]) return;
            var xs = [], ys = [];
            for (var i = 0; i < n; i++) if (labels[i] === lv) { xs.push(d.x[i]); ys.push(d.y[i]); }
            if (xs.length) anns.push({ x: median(xs), y: median(ys), text: lv, showarrow: false, font: { size: 11, color: c.text }, bgcolor: t === 'light' ? 'rgba(252,252,251,0.7)' : 'rgba(26,26,25,0.7)' });
          });
        }
        var L = baseLayout(t, { xlab: d.reduction.toUpperCase() + ' 1', ylab: d.reduction.toUpperCase() + ' 2', legend: false, annotations: anns, dragmode: 'lasso',
          xaxis: { showgrid: false, zeroline: false }, yaxis: { showgrid: false, zeroline: false, scaleanchor: 'x' }, margin: { l: 48, r: 16, t: 8, b: 40 } });
        return { data: [trace], layout: L, config: { modeBarButtonsToAdd: [] } };
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
    var lastLegend = [[], null];
    function renderLegendNow() { renderLegend(lastLegend[0], lastLegend[1]); }
    function renderLegend(levels, key) {
      lastLegend = [levels, key];
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
      var keep = groupChoice;
      clear(groupSel); cats.forEach(function (c) { groupSel.appendChild(h('option', { value: c, text: c })); });
      groupSel.value = keep && cats.indexOf(keep) >= 0 ? keep : preferredGroup(d);
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
          var g = labelsOf(d, groupSel.value);
          download('sparks_selection.csv', toCSV(['barcode', groupSel.value], function (c, k) { return c === 'barcode' ? d.cells[sel[k]] : g[sel[k]]; }, sel.length));
        } }),
        h('button', { class: 'btn', type: 'button', text: 'Copy', onclick: function () { if (navigator.clipboard) navigator.clipboard.writeText(mk().join('\n')); } }),
        h('button', { class: 'btn', type: 'button', text: 'Clear', onclick: function () { sel = []; selSet = null; draw(); } })));
      selBody.appendChild(h('p', { class: 'note', text: 'Import the CSV in the explorer app (Selection → Import barcodes) to run DE, re-cluster or score a signature on exactly these cells.' }));
    }
    colorSel.addEventListener('change', function () { hidden = {}; draw(); });
    themeListeners.push(function () { if (document.body.contains(plotDiv)) renderSel(), renderLegendNow(); });
    size.addEventListener('input', draw); opac.addEventListener('input', draw);
    onDataset(function () { fillColor(); draw(); });
    root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Colour by', colorSel),
      h('label', { class: 'field' }, 'Point size', size), h('label', { class: 'field' }, 'Opacity', opac)));
    root.appendChild(h('div', { class: 'grid g-main' }, card('Embedding', 'lasso-select cells to inspect them', h('div', null, plotDiv, legend)), h('section', { class: 'card' }, h('header', null, selTitle), selBody)));
    fillColor(); draw();
  });

  // ---- Genes
  tab('genes', 'Genes', function () { return DATA.datasets.some(function (d) { return geneNames(d).length > 0; }); }, function (root) {
    var geneIn = h('input', { type: 'text', list: 'sparks-genes', 'aria-label': 'Gene', autocomplete: 'off' });
    var dl = h('datalist', { id: 'sparks-genes' });
    var groupSel = h('select', { 'aria-label': 'Group by' });
    var multi = h('input', { type: 'text', 'aria-label': 'Gene list', style: { minWidth: '320px' } });
    var vBox = plotBox(), dBox = plotBox();
    var dsSel = datasetSelect(), status = h('div', { class: 'note' });
    function fill() {
      var d = ds(), g = geneNames(d);
      clear(dl); g.forEach(function (x) { dl.appendChild(h('option', { value: x })); });
      clear(groupSel); catCols(d).forEach(function (c) { groupSel.appendChild(h('option', { value: c, text: c })); });
      groupSel.value = preferredGroup(d);
      geneIn.value = (d.default_genes && d.default_genes[0]) || g[0] || '';
      multi.value = (d.default_genes && d.default_genes.length ? d.default_genes : g).slice(0, 10).join(', ');
    }
    function draw() {
      var d = ds(), g = geneIn.value.trim(), grp = groupSel.value;
      status.textContent = '';
      if (d.genes[g]) { violin(vBox, d, geneVec(d, g), grp, g + ' by ' + grp); }
      else { showNote(vBox, 'Gene "' + g + '" is not embedded in this report (' + geneNames(d).length + ' genes available).'); }
      var list = multi.value.split(/[\s,;]+/).filter(Boolean).filter(function (x) { return d.genes[x]; });
      if (list.length) dotplot(dBox, d, uniq(list), grp); else showNote(dBox, 'No matching genes for the dot plot.');
    }
    [geneIn, groupSel, multi].forEach(function (el) { el.addEventListener('change', draw); });
    onDataset(function () { fill(); draw(); });
    root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Gene', geneIn, dl), h('label', { class: 'field' }, 'Group by', groupSel),
      h('label', { class: 'field' }, 'Dot plot genes (comma-separated)', multi), status));
    root.appendChild(h('div', { class: 'grid g2' }, card('Expression distribution', 'violin, mean marked', vBox), card('Dot plot', 'size = % of cells expressing; colour = scaled mean', dBox)));
    fill(); draw();
  });

  // ---- Markers & DE
  var DE_KINDS = ['markers_all', 'deg_sig', 'cluster_markers', 'pseudobulk_deg', 'explorer'];
  function deSets() {
    return tablesOfKind(DE_KINDS).filter(function (t) { return t.kind !== 'explorer' || t.cols.indexOf('avg_log2FC') >= 0; });
  }
  tab('de', 'Markers & DE', function () { return deSets().length > 0 || figuresOf('de').length > 0; }, function (root) {
    var sets = deSets();
    // Prefer the full table over its significant-only subset for the same analysis.
    var seen = {};
    sets = sets.filter(function (t) {
      var k = t.kind === 'explorer' ? t.id : [t.comparison, t.unit, t.grouping, t.kind === 'pseudobulk_deg' ? 'pb' : 'sc'].join('|');
      var rank = { markers_all: 3, deg_sig: 2, cluster_markers: 1, pseudobulk_deg: 3, explorer: 3 }[t.kind];
      if (seen[k] && seen[k].rank >= rank) return false; seen[k] = { rank: rank }; return true;
    });
    if (sets.length) {
      var setSel = h('select', { 'aria-label': 'Result set' }, sets.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
      var clusterSel = h('select', { 'aria-label': 'Group' });
      var lfcIn = h('input', { type: 'number', value: 0.25, step: 0.05, min: 0, style: { minWidth: '80px', width: '90px' }, 'aria-label': 'log2FC threshold' });
      var pIn = h('input', { type: 'number', value: 0.05, step: 0.01, min: 0, max: 1, style: { minWidth: '80px', width: '90px' }, 'aria-label': 'Adjusted p threshold' });
      var topIn = h('input', { type: 'number', value: window.innerWidth < 600 ? 4 : 8, step: 1, min: 0, max: 30, style: { minWidth: '70px', width: '80px' }, 'aria-label': 'Labels' });
      var vBox = plotBox('tall'), tblHost = h('div'), info = h('div', { class: 'note' });
      var cur = null, view = null;
      var draw = function () {
        var T = sets[+setSel.value], lc = pickCol(T, ['avg_log2FC', 'logFC', 'log2FoldChange']), pc = pickCol(T, ['p_val_adj', 'padj', 'FDR', 'adj.P.Val', 'p_val', 'PValue']);
        var gc = pickCol(T, ['gene', 'Gene']), cc = pickCol(T, ['cluster', 'group', 'contrast']);
        if (cur !== T) {
          cur = T; clear(clusterSel);
          if (cc) uniq(T.data[cc].map(String)).forEach(function (g) { clusterSel.appendChild(h('option', { value: g, text: g })); });
          clear(tblHost); view = null;
        }
        if (!lc || !pc) showNote(vBox, 'This table has no fold-change / p-value columns; see the table below.');
        var idx = range(T.n).filter(function (i) { return !cc || String(T.data[cc][i]) === clusterSel.value; });
        if (lc && pc) {
          var V = { lfc: idx.map(function (i) { return T.data[lc][i]; }), p: idx.map(function (i) { return T.data[pc][i]; }), genes: idx.map(function (i) { return gc ? T.data[gc][i] : String(i); }) };
          volcano(vBox, V, { lfcCut: +lfcIn.value, pCut: +pIn.value, labelTop: +topIn.value, title: (cc ? clusterSel.value + ' — ' : '') + T.label });
        }
        var nsig = lc && pc ? idx.filter(function (i) { return T.data[pc][i] < +pIn.value && Math.abs(T.data[lc][i]) >= +lfcIn.value; }).length : 0;
        info.textContent = fmtInt(idx.length) + ' genes in this group, ' + fmtInt(nsig) + ' pass the thresholds' + (T.truncated ? ' (table truncated in the report to the strongest rows; see the results folder for the full file)' : '') + '.';
        if (!view) {
          view = tableView(tblHost, T, { filename: T.id, sort: pc ? [pc, 'asc'] : null });
        }
        view.setPredicate(function (i) { return !cc || String(T.data[cc][i]) === clusterSel.value; });
      };
      [setSel, clusterSel, lfcIn, pIn, topIn].forEach(function (el) { el.addEventListener('change', draw); });
      root.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field' }, 'Result set', setSel), h('label', { class: 'field' }, 'Group', clusterSel),
        h('label', { class: 'field' }, '|log2FC| ≥', lfcIn), h('label', { class: 'field' }, 'adj. p <', pIn), h('label', { class: 'field' }, 'Labels / side', topIn)));
      root.appendChild(h('div', { class: 'grid g-main' }, card('Volcano', 'thresholds update live', h('div', null, vBox, info)),
        h('div', { class: 'stack' }, figureCards('de', 'short'))));
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Results table', null, tblHost)));
      draw();
    } else {
      root.appendChild(h('div', { class: 'grid g2' }, figureCards('de')));
    }
  });

  // ---- Enrichment
  tab('enrichment', 'Enrichment', function () { return figuresOf('enrichment').length > 0 || tablesOfKind(['gsea', 'escape']).length > 0; }, function (root) {
    root.appendChild(h('div', { class: 'grid g2' }, figureCards('enrichment', 'tall')));
    var sets = tablesOfKind(['gsea', 'escape']);
    if (sets.length) {
      var host = h('div'), sel = h('select', { 'aria-label': 'Table' }, sets.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
      var show = function () { clear(host); var T = sets[+sel.value]; tableView(host, T, { filename: T.id, sort: T.cols.indexOf('padj') >= 0 ? ['padj', 'asc'] : null }); };
      sel.addEventListener('change', show);
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Enrichment tables', null, h('div', null, h('label', { class: 'field' }, 'Table', sel), host))));
      show();
    }
  });

  // ---- Composition
  tab('composition', 'Composition', function () { return DATA.datasets.length > 0 || figuresOf('composition').length > 0 || tablesOfKind(['proportions', 'prop_test']).length > 0; }, function (root) {
    if (DATA.datasets.length) {
      var dsSel = datasetSelect(), box = plotBox('tall');
      var stackSel = h('select', { 'aria-label': 'Stack by' }), xSel = h('select', { 'aria-label': 'Compare across' });
      var modeSel = h('select', { 'aria-label': 'Scale' }, h('option', { value: 'percent', text: 'Proportions (%)' }), h('option', { value: 'count', text: 'Cell counts' }));
      var fill = function () {
        var d = ds(), cats = catCols(d);
        [stackSel, xSel].forEach(function (s) { clear(s); cats.forEach(function (c) { s.appendChild(h('option', { value: c, text: c })); }); });
        stackSel.value = preferredGroup(d);
        xSel.value = cats.indexOf('orig.ident') >= 0 && stackSel.value !== 'orig.ident' ? 'orig.ident' : (cats.filter(function (c) { return c !== stackSel.value; })[0] || cats[0]);
      };
      var draw = function () {
        var d = ds(), a = d.meta[stackSel.value], b = d.meta[xSel.value];
        if (!a || !b) return;
        var cnt = a.levels.map(function () { return b.levels.map(function () { return 0; }); });
        for (var i = 0; i < d.cells.length; i++) cnt[a.codes[i]][b.codes[i]]++;
        var tot = b.levels.map(function (_, j) { return cnt.reduce(function (s, r) { return s + r[j]; }, 0); });
        var pct = modeSel.value === 'percent';
        stackedBars(box, b.levels, a.levels.map(function (lv, k) {
          return { name: lv, values: cnt[k].map(function (v, j) { return pct ? (tot[j] ? 100 * v / tot[j] : 0) : v; }), counts: cnt[k] };
        }), { percent: pct, xlab: xSel.value, title: stackSel.value + ' composition by ' + xSel.value });
      };
      [stackSel, xSel, modeSel].forEach(function (s) { s.addEventListener('change', draw); });
      onDataset(function () { fill(); draw(); });
      root.appendChild(h('div', { class: 'toolbar' }, dsSel, h('label', { class: 'field' }, 'Stack by', stackSel), h('label', { class: 'field' }, 'Compare across', xSel), h('label', { class: 'field' }, 'Scale', modeSel)));
      root.appendChild(card('Composition', 'computed from the merged object in this report', box));
      fill(); draw();
    }
    var figs = figuresOf('composition');
    if (figs.length) root.appendChild(h('div', { class: 'grid g2', style: { marginTop: '12px' } }, figureCards('composition')));
    var pt = tablesOfKind(['prop_test']);
    if (pt.length) {
      var host = h('div'), sel = h('select', { 'aria-label': 'Proportion test' }, pt.map(function (t, i) { return h('option', { value: i, text: t.label }); }));
      var show = function () { clear(host); var T = pt[+sel.value]; tableView(host, T, { filename: T.id }); };
      sel.addEventListener('change', show);
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Proportion tests (scProportionTest)', null, h('div', null, h('label', { class: 'field' }, 'Test', sel), host))));
      show();
    }
  });

  // ---- All tables
  tab('tables', 'Tables', function () { return DATA.tables.length > 0; }, function (root) {
    var kinds = uniq(DATA.tables.map(function (t) { return t.kind; }));
    var kindSel = h('select', { 'aria-label': 'Kind' }, kinds.map(function (k) { return h('option', { value: k, text: (DATA.kind_labels[k] || k) }); }));
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
    root.appendChild(card('Result table', 'every table the pipeline wrote (large ones are truncated here)', host));
    list();
  });

  // ---- Gallery
  var lightbox = null;
  function gallery(images, withFilters) {
    var host = h('div'), grid = h('div', { class: 'gallery' });
    var q = h('input', { type: 'search', placeholder: 'Filter figures…', 'aria-label': 'Filter figures' });
    var cat = h('select', { 'aria-label': 'Category' }, h('option', { value: '', text: 'All categories' }), uniq(images.map(function (i) { return i.category || ''; })).filter(Boolean).map(function (c) { return h('option', { value: c, text: c }); }));
    var grp = h('select', { 'aria-label': 'Analysis' }, h('option', { value: '', text: 'All analyses' }), uniq(images.map(function (i) { return i.group || ''; })).filter(Boolean).map(function (c) { return h('option', { value: c, text: c }); }));
    var count = h('span', { class: 'muted small' });
    var draw = function () {
      clear(grid);
      var f = q.value.toLowerCase(), shown = images.filter(function (im) {
        return (!cat.value || im.category === cat.value) && (!grp.value || im.group === grp.value) && (!f || (im.label + ' ' + im.rel).toLowerCase().indexOf(f) >= 0);
      });
      shown.slice(0, 120).forEach(function (im) {
        var body = im.src ? h('button', { type: 'button', 'aria-label': 'Enlarge ' + im.label, onclick: function () { openLightbox(im); } }, h('img', { loading: 'lazy', src: im.src, alt: im.label }))
          : h('div', { class: 'nopreview' }, h('a', { href: im.rel, text: 'Open file (' + bytes(im.bytes) + ')', target: '_blank', rel: 'noopener' }));
        grid.appendChild(h('figure', { class: 'shot' }, body, h('figcaption', { text: im.label })));
      });
      count.textContent = shown.length + ' figures' + (shown.length > 120 ? ' (showing 120)' : '');
    };
    [q, cat, grp].forEach(function (el) { el.addEventListener('input', draw); });
    if (withFilters !== false) host.appendChild(h('div', { class: 'toolbar' }, q, h('label', { class: 'field' }, 'Category', cat), h('label', { class: 'field' }, 'Analysis', grp), count));
    host.appendChild(grid); draw();
    return host;
  }
  function openLightbox(im) {
    if (!lightbox) { lightbox = h('dialog', { class: 'lightbox', 'aria-label': 'Figure preview' }); lightbox.addEventListener('click', function () { lightbox.close(); }); document.body.appendChild(lightbox); }
    clear(lightbox).appendChild(h('img', { src: im.src, alt: im.label }));
    if (lightbox.showModal) lightbox.showModal();
  }
  tab('figures', 'Figures', function () { return DATA.images.length > 0; }, function (root) {
    var embedded = DATA.images.filter(function (i) { return i.src; }).length;
    root.appendChild(card('Pipeline figures', embedded + ' embedded of ' + DATA.images.length + (embedded < DATA.images.length ? ' (the rest are linked)' : ''), gallery(DATA.images, true)));
  });

  // ---- Files
  tab('files', 'Files', function () { return DATA.files.length > 0; }, function (root) {
    var cols = ['file', 'kind', 'category', 'size'];
    var T = { cols: cols, n: DATA.files.length, data: {
      file: DATA.files.map(function (f) { return f.rel; }), kind: DATA.files.map(function (f) { return f.kind; }),
      category: DATA.files.map(function (f) { return f.category || ''; }), size: DATA.files.map(function (f) { return f.size; }) } };
    var reports = DATA.files.filter(function (f) { return f.kind === 'report'; });
    if (reports.length) {
      root.appendChild(h('div', { style: { marginBottom: '12px' } }, card('Pipeline reports', 'generated by the run; open alongside this report',
        h('div', { class: 'toolbar', style: { marginBottom: 0 } }, reports.map(function (f) { return h('a', { class: 'btn', href: f.rel, target: '_blank', rel: 'noopener', text: f.rel.split('/').pop() }); })))));
    }
    var host = h('div');
    tableView(host, T, { filename: 'sparks_files', pageSize: 50, render: {
      file: function (v) { return h('a', { href: v, target: '_blank', rel: 'noopener', text: v }); },
      size: function (v) { return bytes(v); } } });
    root.appendChild(card('Result files', 'links work when this report sits in the results folder', host));
  });

  // ---- Run
  tab('run', 'Run', function () { return true; }, function (root) {
    var R = DATA.run, blocks = [];
    if (R.warnings && R.warnings.length) blocks.push(card('Warnings (' + R.warnings.length + ')', null, h('pre', { text: R.warnings.join('\n') })));
    if (R.config_yaml) blocks.push(card('Configuration used', null, h('pre', { text: R.config_yaml })));
    if (R.log_tail && R.log_tail.length) blocks.push(card('Pipeline log (tail)', null, h('pre', { text: R.log_tail.join('\n') })));
    if (R.session_info && R.session_info.length) blocks.push(card('Session info', null, h('pre', { text: R.session_info.join('\n') })));
    if (!blocks.length) blocks.push(emptyNote('No run information was found in the results folder.'));
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ------------------------------------------------------------------ shell
  var panels = {}, tabButtons = {}, active = null;
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && x.available(); })[0] || TABS.filter(function (x) { return x.available(); })[0];
    if (!t) return;
    id = t.id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) { t.built = true; try { t.build(panels[id]); } catch (e) { panels[id].appendChild(h('div', { class: 'empty', text: 'This tab could not be built: ' + e.message })); if (window.console) console.error(e); } }
    else refreshVisible(panels[id]);
    // Newly built plots in this tab need a layout pass once visible.
    setTimeout(function () { refreshVisible(panels[id]); }, 60);
    active = id;
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) { /* file:// may refuse */ } }
  }
  function init() {
    applyThemeAttr();
    document.title = DATA.meta.title;
    var tabsEl = $('#tabs'), main = $('#panels');
    var avail = TABS.filter(function (t) { return t.available(); });
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
      if (mq.addEventListener) mq.addEventListener('change', onChange);
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
