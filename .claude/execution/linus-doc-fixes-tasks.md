# Task ledger: Linus documentation fixes

Plan: `.claude/execution/linus-doc-fixes-plan.md`

- [x] T1 — Strengthen authored-source and generated-site documentation
  contracts, including required-file, canonical-route, and retired-route
  assertions.
- [x] T2 — Replace row-order selection in README and vignette with the explicit
  checked `inflection` inspection example.
- [x] T3 — Show and interpret the node-to-metacluster map in README and
  vignette.
- [x] T4 — Add and explain `criterion_pass`, component matrices, and
  `failure_reason`.
- [x] T5 — Show selected-criterion passed/failed/unresolved/total counts and
  correct the `qc_pass_rate` interpretation in both guides.
- [x] T6 — Regenerate roxygen, vignette HTML, pkgdown, and cleaned search
  output.
- [x] T7 — Pass `git diff --check`, focused documentation tests, and the full
  test suite.
- [x] T8 — Build the source tarball and pass `R CMD check --no-manual` with
  0 errors, 0 warnings, and 0 notes.
- [x] T9 — Update the existing Linus review summary with truthful fix status
  and fresh validation evidence.

The executor must mark each item complete only after retaining corresponding
command output or direct file evidence in its final report.

## Authorized generator amendment

During T7, canonical pkgdown 2.2.0 generation was found to add a second final
newline to refreshed `docs/reference/*.html` files. The existing generated HTML
uses this form consistently, but Git reported the newly refreshed formerly
stale page as `blank-at-eof`. The coordinator authorized the narrow
`.gitattributes` rule `docs/**/*.html whitespace=-blank-at-eof`; it leaves all
other whitespace checks active. An isolated Git probe confirmed the HTML
attribute value, continued to flag trailing whitespace in HTML, and continued
to flag `blank-at-eof` outside generated HTML.

The first T8 source check also showed that the execution packet directory
`.claude` entered the tarball and produced a package-structure note. The
coordinator authorized the escape-safe `^[.]claude$` entry in `.Rbuildignore`
without broadening the package ignore policy. A first literal-backslash form did
not match and was replaced after an explicit `grepl()` check. T8 therefore
requires a fresh source build, tar listing check, and package check.

A parent audit after T9 found canonical R-console padding in the untracked
`docs/articles/fastINFLECT.md`. Re-running pkgdown's own `convert_md()` with the
pinned Pandoc 3.8.3 reproduced the article byte-for-byte, including those line
endings. The authorized `docs/articles/*.md whitespace=-blank-at-eol` rule is
limited to generated article Markdown. An isolated Git probe confirmed that it
suppresses the canonical article padding while still flagging trailing
whitespace in non-generated Markdown and `blank-at-eof` in generated article
Markdown. The follow-up plain `git diff --check` passed, the actual no-index
article audit emitted no whitespace finding, the rebuilt tarball contained zero
`.claude` paths, and the full installed-package check again ended `Status: OK`
with 573 passes and the one expected generated-site skip.
