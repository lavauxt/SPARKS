// Browser smoke test for SPARKS_all_results.html (needs `playwright` and a Chromium).
//
//   node tests/browser/report-smoke.js path/to/SPARKS_all_results.html [screenshot_prefix]
//
// Opens the report from disk (file://, no server), visits every tab, makes a lasso
// selection on the embedding, colours by a gene, toggles the theme and the legend,
// changes the volcano thresholds and checks there is no console error. Exits 1 on failure.
const { chromium } = require('playwright');

const file = process.argv[2];
const shot = process.argv[3] || null;
if (!file) { console.error('usage: node report-smoke.js report.html [screenshot_prefix]'); process.exit(2); }

const failures = [];
const check = (ok, msg) => { console.log((ok ? 'ok   ' : 'FAIL ') + msg); if (!ok) failures.push(msg); };

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1360, height: 900 } });
  const errors = [];
  page.on('console', m => { if (['error', 'warning'].includes(m.type())) errors.push(m.type() + ': ' + m.text()); });
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  const external = [];
  page.on('request', r => { if (!/^(file|data|blob|about):/.test(r.url())) external.push(r.url()); });

  await page.goto('file://' + require('path').resolve(file));
  await page.waitForSelector('html[data-ready="1"]', { timeout: 30000 });

  // every tab renders something, and plots draw in the tabs that have them
  const tabs = await page.$$eval('.tab', els => els.map(e => e.id.replace('tab-', '')));
  check(tabs.length >= 5, 'tabs found: ' + tabs.join(', '));
  for (const t of tabs) {
    await page.click('#tab-' + t);
    await page.waitForTimeout(700);
    const info = await page.evaluate(() => {
      const p = document.querySelector('.panel:not([hidden])');
      return { text: p ? p.innerText.length : 0, plots: p ? p.querySelectorAll('.js-plotly-plot .main-svg').length : 0 };
    });
    check(info.text > 20, `tab "${t}" has content (${info.text} chars, ${info.plots} plot layers)`);
    if (shot) await page.screenshot({ path: `${shot}_${t}.png` });
  }
  for (const t of ['embedding', 'qc', 'de']) {
    if (!tabs.includes(t)) continue;
    await page.click('#tab-' + t); await page.waitForTimeout(700);
    const n = await page.locator('.panel:not([hidden]) .js-plotly-plot').count();
    check(n > 0, `tab "${t}" draws at least one plotly figure (${n})`);
  }

  // lasso selection on the embedding
  if (tabs.includes('embedding')) {
    await page.click('#tab-embedding'); await page.waitForTimeout(700);
    const box = await page.locator('#panel-embedding .js-plotly-plot .nsewdrag').first().boundingBox();
    const x0 = box.x + box.width * 0.05, y0 = box.y + box.height * 0.05;
    await page.mouse.move(x0, y0); await page.mouse.down();
    for (const [dx, dy] of [[0, 0], [0.9, 0], [0.9, 0.9], [0, 0.9], [0, 0]]) {
      await page.mouse.move(x0 + box.width * dx, y0 + box.height * dy, { steps: 6 });
    }
    await page.mouse.up(); await page.waitForTimeout(800);
    const text = await page.locator('#panel-embedding').innerText();
    check(/selected/i.test(text), 'lasso selection reports the selected cells');

    // colour by a gene, back to a category, legend toggle
    const colour = page.locator('#panel-embedding select[aria-label="Colour by"]');
    const genes = await colour.locator('optgroup[label*="ene"] option').evaluateAll(o => o.map(x => x.value));
    if (genes.length) {
      await colour.selectOption(genes[0]); await page.waitForTimeout(600);
      check(await page.locator('#panel-embedding .js-plotly-plot').count() > 0, 'colouring by gene ' + genes[0]);
    }
    const cats = await colour.locator('option').evaluateAll(o => o.map(x => x.value));
    await colour.selectOption(cats[0]); await page.waitForTimeout(500);
    const btn = page.locator('#panel-embedding .legend button').first();
    if (await btn.count()) {
      await btn.click(); await page.waitForTimeout(300);
      check((await btn.getAttribute('aria-pressed')) === 'false', 'legend click hides a group');
    }
  }

  // theme: auto -> light -> dark
  await page.click('#theme-btn'); await page.click('#theme-btn'); await page.waitForTimeout(400);
  check((await page.evaluate(() => document.documentElement.getAttribute('data-theme'))) === 'dark', 'theme toggle reaches dark');
  if (shot) await page.screenshot({ path: `${shot}_dark.png` });

  // volcano thresholds change the counts
  if (tabs.includes('de')) {
    await page.click('#tab-de'); await page.waitForTimeout(800);
    const input = page.locator('#panel-de input[aria-label="log2FC threshold"]');
    if (await input.count()) {
      const before = await page.locator('#panel-de .note').first().innerText();
      await input.fill('3'); await input.dispatchEvent('change'); await page.waitForTimeout(500);
      const after = await page.locator('#panel-de .note').first().innerText();
      check(before !== after, 'volcano thresholds update the summary');
    }
  }

  check(external.length === 0, 'no network requests (' + external.slice(0, 3).join(', ') + ')');
  check(errors.length === 0, 'no console errors or warnings' + (errors.length ? ': ' + errors.slice(0, 3).join(' | ') : ''));
  await browser.close();
  if (failures.length) { console.error(`\n${failures.length} check(s) failed`); process.exit(1); }
  console.log('\nall checks passed');
})().catch(e => { console.error(e); process.exit(1); });
