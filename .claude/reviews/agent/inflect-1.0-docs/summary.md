# Linus review - agent/inflect-1.0-docs

## Scope

Reviewed the current working tree against `HEAD` (`9b781e5`). The change
replaces three overlapping vignettes with one task-based guide, simplifies the
pkgdown navigation, rewrites the README and primary reference prose, updates
documentation contract tests, and regenerates the site. Runtime R function
bodies are unchanged.

Validation completed before review:

- source vignette and pkgdown site rendered successfully;
- 582 package tests passed with no failures, warnings, or skips;
- `R CMD check` completed with 0 errors, 0 warnings, and 0 notes.

The initial review was invoked without `--fix`, so no findings were auto-fixed
at that stage. The follow-up execution described below clears all four findings.

## Findings

### [FIXED] Generated-documentation gates fail open

Files: `tests/testthat/test-package-rename.R:32`,
`tests/testthat/test-release-version.R:44`

The tests discard missing rendered files before checking them, while the
stale-slug loop scans only surviving authored sources. A build can omit the new
article or retain retired pages in `docs/`, `search.json`, or `sitemap.xml` and
still pass.

Why this is bad taste: the repository tracks generated documentation, so its
existence and route closure are part of the release contract. Optional checks
turn a known stale-site failure mode into a silent success.

The fix: assert every required rendered file exists, then assert every retired
slug is absent from generated article files, search records, and the sitemap.

### [FIXED] Candidate selection depends on row order

Files: `vignettes/fastINFLECT.Rmd:113`, `README.md:73`

Both examples select `subset(result$selection, partition_available)[1, ]`
immediately after explaining that candidate rows are not rankings. This makes
table order an undocumented preference. When no candidate partition is
available, one-row indexing produces an all-`NA` pseudo-row instead of a clear
failure.

Why this is bad taste: the example hides a selection policy inside indexing and
handles the empty case accidentally.

The fix: require an explicit method, verify that exactly one matching row is
available, and stop with an instruction to add the fitted `k` to `set.i` when
the partition was not tested.

### [FIXED] Pass-rate wording overstates the result

File: `README.md:89`

The README says a higher `qc_pass_rate` means fewer pairs were flagged. The
denominator changes with `k` and unresolved entries remain in that denominator.
A higher fraction can coexist with a larger absolute number of failed or
unresolved pairs.

Why this is bad taste: the sentence turns a relative rate into an unsupported
absolute claim in the part of the documentation meant to teach interpretation.

The fix: say only that a larger fraction passed the selected criterion, then
direct readers to failed and unresolved counts for the candidate partition.

### [FIXED] Partition extraction stops before interpretation

File: `vignettes/fastINFLECT.Rmd:117`

The guide assigns `partition` but never displays or interprets it. The result
map also omits `criterion_pass`, the matrix used to calculate the displayed
aggregate for the selected `uniform.test` mode.

Why this is bad taste: the consolidated guide reaches the central result and
then leaves the reader with an unused object and an incomplete map between the
aggregate rate and its decisions.

The fix: show a compact partition summary, explain that positions are SOM nodes
and values are metacluster labels, and include `criterion_pass` beside the
component matrices.

## Fix verification

All four findings are fixed in the current working tree:

- authored-source and checkout-only generated-site tests now require every
  contract file, the exact canonical route, nonempty search paths, and absence
  of all four retired slugs; the source package receives one explicit
  generated-site skip because `_pkgdown.yml` and `docs/` are excluded;
- README and vignette select `candidate_method <- "inflection"`, require exactly
  one row, require an available partition, and require the character `k` key;
- both guides display the SOM-node map and label counts, `criterion_pass`, its
  component matrices, `failure_reason`, and the selected criterion-summary row;
- both guides define `qc_pass_rate` as `100 * passed / total`, retain unresolved
  decisions in `total`, and limit interpretation to the fraction that passed.

The rendered example executes the directly tested inflection candidate at
`k = 9`. It displays 375 SOM nodes, node-label counts
`48, 6, 93, 23, 96, 26, 36, 39, 8`, and the selected combined-criterion summary:
205 passed, 83 failed, 0 unresolved, 288 total, and 71.18056% passed.

Fresh validation completed with the R 4.5.1 environment:

1. `Rscript -e 'devtools::document()'` completed.
2. The standalone vignette rendered after `pkgload::load_all(".")`.
3. `Rscript -e 'pkgdown::build_site()'` completed with the locally cached
   Pandoc 3.8.3, followed by `Rscript data-raw/clean-pkgdown-search.R`.
4. `git diff --check` passed.
5. Focused `package-rename|release-version` tests passed: 46 passed, 0 failed,
   0 warnings, 0 skips.
6. The full checkout suite passed: 600 passed, 0 failed, 0 warnings, 0 skips.
7. `R CMD build .` produced `fastINFLECT_1.0.0.tar.gz`; the verified tarball
   contains zero `.claude` paths.
8. `R CMD check --no-manual fastINFLECT_1.0.0.tar.gz` ended `Status: OK`: 0
   errors, 0 warnings, and 0 notes. Its installed-package suite reported 573
   passed and the one expected generated-site skip. Repository index lookups
   were unavailable in the execution environment, so the existing packet-local
   `multimode` library was exposed to the check.

Three narrow, authorized release-hygiene amendments support these gates:
`docs/**/*.html whitespace=-blank-at-eof` suppresses only pkgdown's consistent
second terminal newline in generated HTML;
`docs/articles/*.md whitespace=-blank-at-eol` suppresses canonical R-console
padding in generated article Markdown; and `^[.]claude$` keeps the execution
packet out of source packages. Pkgdown's own Markdown conversion reproduced the
article byte-for-byte. Isolated Git probes confirmed that HTML trailing
whitespace, trailing whitespace in non-generated Markdown, generated-article
blank lines at EOF, and non-generated blank lines at EOF remain errors.
After this follow-up amendment, plain `git diff --check` passed, the actual
untracked article audit emitted no whitespace finding, the rebuilt tarball
again contained zero `.claude` paths, and `R CMD check --no-manual` again ended
`Status: OK` with 573 passes and the one expected generated-site skip.

## What looks good

### [GOOD TASTE] Direct navigation

File: `_pkgdown.yml:12`

The navbar exposes one `Get started` route and one task-grouped reference index.
It removes article-index hopping and keeps benchmarks outside the user path.

### [GOOD TASTE] Task-shaped guide

File: `vignettes/fastINFLECT.Rmd:26`

The guide follows four concrete jobs: run, read, inspect, and report. The
sequence is local, short, and easier to navigate than the three documents it
replaces.

### [GOOD TASTE] Build discipline

The source vignette, generated site, full test suite, and installed-package
check all pass. The change does not hide a documentation rewrite behind broken
examples or stale generated output.

## Calibrated non-findings

- README overlap with the vignette is normal for a package quick start. The
  concrete inconsistencies above matter; duplication by itself does not.
- Pkgdown search-text negation loss appears generator-owned and was not shown to
  originate in this change.
- Retired URL redirects may be useful if those routes were deployed, but this
  branch is not merged into `origin/main`; deployment history must be confirmed
  before treating redirects as a release blocker.
- Raw YAML substring assertions and the one-use `schedule` variable are minor
  maintenance issues, not user-facing defects.

## Overall taste rating

**GOOD**

The navigation remains simple, the inspection policy is explicit, and the
release gates now fail closed without weakening source-package checks. The
rendered example connects the candidate, inspected partition, selected decision
matrix, component evidence, and bounded aggregate interpretation.

Fix status: 4 of 4 findings fixed. The Linus documentation review is cleared.
