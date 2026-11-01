// Browser smoke test for the Shiny explorer (needs `playwright` and a Chromium).
//
//   Rscript -e 'd <- tempfile(); SPARKS::sparks_demo_results(d); SPARKS::sparks_explorer(d, port = 7733, launch.browser = FALSE)'
//   node tests/browser/explorer-smoke.js http://127.0.0.1:7733 [screenshot_prefix]
//
// Lasso -> save as category -> differential expression -> same analysis again (must be a
// cache hit: no new history row) -> colour by gene -> build the unified report.
const { chromium } = require('playwright');

const url = process.argv[2] || 'http://127.0.0.1:7733';
const shot = process.argv[3] || null;
const failures = [];
const check = (ok, msg) => { console.log((ok ? 'ok   ' : 'FAIL ') + msg); if (!ok) failures.push(msg); };

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 950 } });
  const errors = [];
  const notFound = [];
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('response', r => { if (r.status() === 404) notFound.push(r.url()); });

  await page.goto(url, { waitUntil: 'load' });
  await page.waitForSelector('#explorer-embedding .main-svg', { timeout: 90000 });
  await page.waitForTimeout(1500);
  check(true, 'embedding drawn');

  const box = await page.locator('#explorer-embedding .nsewdrag').first().boundingBox();
  const x0 = box.x + box.width * 0.55, y0 = box.y + box.height * 0.55;
  await page.mouse.move(x0, y0); await page.mouse.down();
  for (const [dx, dy] of [[0, 0], [0.25, 0], [0.25, 0.4], [0, 0.4], [0, 0]]) {
    await page.mouse.move(x0 + box.width * dx, y0 + box.height * dy, { steps: 6 });
  }
  await page.mouse.up(); await page.waitForTimeout(1500);
  const summary = await page.locator('#explorer-sel_summary').innerText();
  check(/\d+/.test(summary) && !/^\s*0\b/.test(summary), 'lasso selected cells: ' + summary.replace(/\s+/g, ' ').slice(0, 80));

  await page.fill('#explorer-cat_name', 'Region A'); await page.click('#explorer-save_cat'); await page.waitForTimeout(800);
  check(/Region A/.test(await page.locator('#explorer-cats_ui').innerText()), 'selection saved as a category');

  await page.click('#explorer-run_de');
  await page.waitForSelector('#explorer-res_volcano .main-svg', { timeout: 180000 });
  check(true, 'differential expression ran and the volcano is drawn');
  if (shot) await page.screenshot({ path: `${shot}_de.png` });

  const rows0 = await page.locator('#explorer-history tbody tr').count();
  await page.click('a[data-value="Embedding"]'); await page.waitForTimeout(500);
  await page.click('#explorer-run_de'); await page.waitForTimeout(5000);
  await page.click('a[data-value="Analyses"]'); await page.waitForTimeout(800);
  const rows1 = await page.locator('#explorer-history tbody tr').count();
  check(rows0 >= 1 && rows1 === rows0, `repeating the analysis is a cache hit (history rows ${rows0} -> ${rows1})`);

  await page.click('a[data-value="Embedding"]');
  const genes = await page.evaluate(() => {
    const s = document.querySelector('#explorer-color_by').selectize;
    return Object.keys(s.options).filter(k => /^[A-Z0-9]+$/.test(k));
  });
  if (genes.length) {
    await page.evaluate(g => document.querySelector('#explorer-color_by').selectize.setValue(g), genes[0]);
    await page.waitForTimeout(2500);
    check(await page.locator('#explorer-embedding .main-svg').count() > 0, 'coloured by gene ' + genes[0]);
  }

  await page.click('a[data-value="Analyses"]'); await page.waitForTimeout(500);
  await page.click('#explorer-build_report');
  await page.waitForSelector('#explorer-report_link a', { timeout: 180000 });
  check(true, 'unified report built from the app: ' + (await page.locator('#explorer-report_link a').innerText()));

  // Some distro-packaged shiny / DT builds (e.g. Debian's r-cran-shiny, r-cran-dt) ship without
  // two optional static assets that the page asks for: shiny's `strftime-min.js` (date
  // formatting; the explorer has no date input) and DT's bootstrap-5 stylesheet. They are
  // reported as warnings, not failures. Any other 404 or console error still fails the test.
  const KNOWN_MISSING = /strftime|dt-core-bootstrap5/;
  const skew = notFound.length > 0 && notFound.every(u => KNOWN_MISSING.test(u));
  const real = errors.filter(e => !(skew && (/Failed to load resource.*404/.test(e) || /GLOBAL is not defined/.test(e))));
  if (notFound.length) console.log('404s: ' + notFound.join(', '));
  if (skew) console.log('warn  ignored missing optional static assets of this shiny/DT install');
  check(real.length === 0, 'no browser console errors' + (real.length ? ': ' + real.slice(0, 3).join(' | ') : ''));
  await browser.close();
  if (failures.length) { console.error(`\n${failures.length} check(s) failed`); process.exit(1); }
  console.log('\nall checks passed');
})().catch(e => { console.error(e); process.exit(1); });
