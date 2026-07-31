# Execution handoff: clear Linus documentation findings

Repository: `/exports/para-lipg-hpc/mdmanurung/INFLECT`

Branch: `agent/inflect-1.0-docs`

Frozen plan: `.claude/execution/linus-doc-fixes-plan.md`

Task ledger: `.claude/execution/linus-doc-fixes-tasks.md`

Source review: `.claude/reviews/agent/inflect-1.0-docs/summary.md`

## Executor contract

1. Read the frozen plan, task ledger, and source review completely before
   editing.
2. Inspect current files and preserve every existing working-tree change.
3. Execute every task in order, updating the task ledger as evidence is
   obtained.
4. Use `apply_patch` for authored source and test edits. Generated outputs may
   be rewritten only by their canonical R build commands.
5. Use only the R 4.5.1 binaries under
   `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin`.
6. Do not change runtime function bodies, package APIs, benchmarks, release
   version, or retired-vignette policy. Do not reset, clean, stage, commit, or
   push.
7. If a gate fails, diagnose and fix only within the frozen scope, regenerate
   dependent outputs, and rerun every invalidated gate.
8. Return a concise report containing files changed, each validation command
   and outcome, exact test counts when available, the final `R CMD check`
   status, and any residual concern. Do not declare completion without all
   completion criteria in the frozen plan.
