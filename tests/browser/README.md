# Browser smoke tests

These are not run by `devtools::test()`. They drive a real Chromium with
[Playwright](https://playwright.dev/) and check what unit tests cannot: that the
report's JavaScript runs, plots draw, selection and theming work, and that the
Shiny explorer wires the lasso, the analyses and the targets cache end to end.

```sh
npm install playwright && npx playwright install chromium   # once
```

**Report** (offline, from disk):

```sh
Rscript -e 'd <- tempfile(); SPARKS::sparks_demo_results(d); cat(SPARKS::sparks_report(d, quiet = TRUE))'
node tests/browser/report-smoke.js /path/printed/above/SPARKS_all_results.html
```

**Explorer** (start the app, then point the test at it):

```sh
Rscript -e 'd <- tempfile(); SPARKS::sparks_demo_results(d); SPARKS::sparks_explorer(d, port = 7733, launch.browser = FALSE)' &
node tests/browser/explorer-smoke.js http://127.0.0.1:7733
```

Both exit with status 1 if a check fails. An optional last argument is a prefix
for screenshots. SPARKS must be installed (not just `load_all()`-ed) for the
explorer to run analyses in a background process; otherwise they run in-process,
which is fine for the test.
