# MECHANISM — HARDENED CONVENTIONS (authoritative, read first)

> **This file is the single source of truth for how the mechanism system works.** It was hardened
> 2026-07-19 after a session that built the fold pipeline, the CS-format sync, and a 158-node graph.
> A future session — including a weaker model — MUST NOT deviate from this via a wrong grep, a stale
> file, a hand-rolled fold, or an isolation leak.
>
> **FIRST ACTION, every session:** run `python3 scripts/mechanism_hardening_check.py`. If it FAILS,
> fix the violation at its source — **never work around a failing check** (a worked-around guard is
> worse than none). `mechanism_preflight.py` runs it for you.
>
> If anything below conflicts with a doc you found by grep, **this file wins** — the other is stale
> (see §7). Do not trust a docstring over what the code does + what the hardening check asserts.

## 1. ISOLATION INVARIANTS (non-negotiable — the reason this fork exists)
- **Never write CS's graph.** `data/ANCHOR_GRAPH.json` is CS's (coordinator reads it). Biology lives ONLY
  in `data/MECHANISM_ANCHOR_GRAPH.json`. The hardening check asserts the two graphs share **zero** node
  ids — a collision is an isolation breach.
- **Never `git push`.** `upstream` = cad-to-simulation, FETCH-ONLY (pull via `merge -X ours`). Push stays
  mechanically disabled.
- **Knowledge lives in-repo.** `bt_memory/`, this doc, the graph, `WAVE_PLAN.md`. Never load or write
  `~/.coordinator*/…/memory` (it bleeds to other sessions/accounts). `context_loader.py --memdir bt_memory`.
- **Neutral dialect off-repo.** Anything another process/session reads (`/tmp` tags, a question relayed
  to CS, commit trailers) carries NO biology terms and NO mechanism identity. Protect identity: a
  relayed question is stripped to the neutral, inherited-machinery core (see the bt_memory note).
- **Commit trailer:** `Co-Authored-By: mechanism (coordinator coordinator 4.8, 1M ctx) <noreply@anthropic.com>` —
  never the inherited CS coordinator trailer.

## 2. THE CANONICAL FOLD PATH (the only way biology enters the graph)
```
research subagent (return-only, MT-JSON-BEGIN … MT-JSON-END sentinel)
  -> mechanism_fold.extract_cell_json(raw)          # ONE shared extractor; refuses, never fabricates
  -> mechanism_fold.fold_one(...)                   # builds HYPOTHESIS submission
  -> fold_gate_v2(submission)                      # CS's blessed gate; HYPOTHESIS short-circuits to ALLOW
  -> upsert a CANONICAL node into data/MECHANISM_ANCHOR_GRAPH.json
```
- **Use `fold_gate_v2`, not v1, not a hand-rolled gate.** v2 = v1 receipt + ATOMS coverage/triviality;
  a HYPOTHESIS folds as hypothesis unchanged. The check asserts `mechanism_fold` imports `fold_gate_v2`.
- **Do NOT re-implement extraction.** `extract_cell_json` is the one shared reader (CS itself had 4
  fragile copies — we keep one). It prefers the `MT-JSON-BEGIN/END` sentinel, then a ```json fence,
  then the largest balanced-brace object, and **REFUSES** (routes to `data/body_twin/fold_queue.jsonl`)
  rather than folding a thin/garbage extraction.

## 3. THE NODE SCHEMA (anchor_graph_tools-compatible → CS-combinable)
Every node: `{id, claim, type, status, evidence, depends_on, load, risk, regime_note, actionable_by,
cluster, cert_design}`.
- `status` ∈ `{PROVEN, REFUTED, DEFERRED, OPEN, ASSUMED}` (CS's `VALID_STATUS`). A hypothesis =
  **`OPEN`** ("designed, not measured"). The rich mechanism grade (`MEASURED-B` etc.) is a **derived
  scorecard label** (`mechanism_scorecard.py`), stored in `cert_design.mechanism_grade`, **never** in
  `status` — else the scorecard and CS validate both break.
- `type` ∈ `{THEORY, EMPIRICAL, ENGINEERING, DATA-ANCHOR, GOAL}`; `risk` ∈ `{LOW, MED, HIGH, SOURCE-COND}`.
- The graph wrapper is `{_meta, nodes}`. Authoritative builder: `scripts/mechanism_build_anchor_graph.py`.

## 4. WAVE DISCIPLINE (firing research subagents)
- Subagents are **return-only**: prompt each "no Write/Edit/memory; return ONE line
  `MT-JSON-BEGIN {…minified…} MT-JSON-END`; cite DOI/PMID per number; honest gap = PASS; decorrelated
  anchor required." QC **every** output — it is a HYPOTHESIS, not truth.
- **CS-dialect sanitizer:** agents sometimes echo CS reduced-rep dialect (`V_d`, `SELTO`, `manifold`,
  `reduced-rep`, `sigma_min`, `gradient-projection`) into a proposed_cell. `mechanism_fold._sanitize`
  scrubs it from the folded claim; the periodic overview re-scrubs. The check FAILS on residual dialect.
- **Couplings → edges:** capture each cell's `couples_to`; resolve prose→node-ids into
  `cert_design.couples_to_ids` (the "everything-hangs-together" densification).
- **Video compute:** `nice -n 19 ffmpeg -hwaccel cuda -c:v av1_cuvid`; NEVER `cv2.VideoCapture`
  (froze the box once). See the bt_memory note.

## 4b. PIPELINE STAGES + THE CROSS-CUTTING DISCIPLINE
> Supersedes an earlier draft that listed four co-equal "agent roles" — that OVER-SPLIT it: the
> literature-scout is really *acquisition* (a paper is just a source-type), and verification is a
> cross-cutting *function*, not its own stage. The clean model is **three stages + one discipline.**

The coordinator drives THREE stages with subagents, and verification wraps all of them:

1. **ACQUIRE** (subagents, weak, massively parallel) — bring the world in machine-readable. Two source
   kinds, same stage: (a) the **cartographer** lands DATA sources (datasets/video/files) → `SOURCE_REGISTRY`,
   routing by type; (b) the **literature-scout** lands KNOWLEDGE — return-only research agents each emitting
   one `MT-JSON-BEGIN … END` block per hole. Both are acquisition. (The waves this session were the scout.)
2. **DESIGN** (the COORDINATOR owns it — strong model) — decide the cert: which occluded hidden state to
   corner, which ≥2 **decorrelated** legs (different physics/instruments) meet on it, and which **held-out**
   external anchor checks it (never a tautology). A scout's `proposed_cell` is a *suggestion the coordinator
   owns + QCs*, not a standalone role. A designed-not-measured cell folds via `mechanism_fold → fold_gate_v2`
   as `status:OPEN`. ⚠ At 100× the coordinator can't be the sole designer — design decomposes into
   (weak propose) + (MACHINE admission gate: anchor-independence / `coupling_admission_precheck`) + (strong
   adjudicate the residual), and strong-model leverage migrates to *building the gates*. See
   `docs/MECHANISM_SCALING_DOCTRINE.md`.
3. **MEASURE / CERTIFY** (subagents + compute) — **this is where the decorrelation is actually MEASURED**:
   build the legs' forward-models, run them against acquired data, check they agree on the hidden state vs
   the held-out anchor. Objective — the verdict is the number (`fold_gate_v2` §5 RESULT form), not a model's
   judgement. Advances `OPEN → ASSUMED/PROVEN/REFUTED`. (Currently un-done for all cells — the graph is a
   *design frontier*, not yet certified.)

**VERIFICATION — cross-cutting discipline, not a stage.** Wraps every stage: QC every subagent output (it is
a hypothesis), an adversarial refute-leg on a DIFFERENT context, and the receiver's OWN disk re-run
(agent-GREEN never counts, R46). CS spawns it as a verify-subagent at intake; the coordinator applies it
throughout.

Sequence: **acquire → design → measure/certify**, verification throughout. Model class per stage:
acquire + measure = weak/parallel; design = strong (→ increasingly machine-gated at scale); verify = machine
gates + a medium adversary.

**Naming convention** — label every spawned subagent `<role>: <topic>` (e.g. `hypothesizer: iron-metabolism`,
`executor: knee-cell`, `verifier: <claim>`, `cartographer: <source>`) so the FleetView distribution of roles is
legible at a glance. NB: `watertight-researcher` is the agent-TYPE (the method baked into its system prompt),
NOT the role — one type serves several roles. (Today all waves are `hypothesizer:` because only the ACQUIRE
stage runs; `executor:` / `verifier:` appear once execution starts.)

## 5. RESULT-FORM CONTRACT (when you EXECUTE a cell, not just fold a hypothesis)
A measured RESULT must be gate-runnable on the first `fold_gate_v2(claim, atoms_report=report)` call:
- Claim = cge fields (`decisive_number`, `decisive_number_path`, `artifacts` full-path, `legs`,
  `claimed_scope{n,substrate}`) + v1 receipt (`verdict∈{positive,negative,partial}`,
  `rerun.{command, artifact, checks[{key,expected[,tol]}]}`). Mirror `data/claims/I-VLM-DATASET.json`.
- ATOMS block lives in the REPORT; `artifact`-exists key is `artifact` (not `path`); decisive atom =
  tight inequality; mirror `reports/probes/compressed_domain_probe.json`.
- Evidence = a `reports/probes/<name>.json` on disk (or a fold-ledger id / `docs/*.md` doc-pointer),
  **not** a raw `agent_outputs/*.json` path.
- Verification = an independent leg + **your own disk re-run** (agent-GREEN never counts, R46).

**Absorbed from CS's 2026-07-19 fold-path audit (booked doctrine — apply when building the execute path;
detail in [[cs-fold-execute-path-hardening-absorbed-dig-fix-klobber-own-atoms-hyp-failopen]]):**
- **The cell writes its OWN ATOMS block** into `reports/probes/<name>.json` during report generation — never
  rely on a post-hoc inject. A reproducer rerun regenerates the report and clobbers hub-applied atoms (the
  "vanishing-ATOMS loop"); self-written atoms survive regeneration.
- **An atom never mutates what it measures** — round-trip against a *pinned* copy. (CS's VLM determinism atom
  deleted its own artifacts as step 1 → every QC run destroyed the thing it measured.)
- **KLOBBER-GUARD:** any grade daemon calls the gate with `execute=False` — it must never auto-run a claim's
  `rerun.command`. (We have not pulled CS's `grade_watcher.py`; when we build our execute path, guard it from
  the start.)
- **`_HYP` fail-open:** `fold_gate.py`'s `_HYP` regex (`HYPOTHES*`, `awaiting-qc`, `speculativ`) is OR'd into the
  hypothesis flag, so a measured RESULT whose text merely *contains* those words is silently downgraded to an
  un-gated hypothesis. Keep those words out of the fields `_HYP` scans on a RESULT fold (`fold_gate.py:49,138,140`).
- ✅ Absorbed as CODE 2026-07-20: the `_dig` ambiguity-fix in `scripts/fold_gate.py` (literal-dot keys like
  `PF_Boron_fs0.5` no longer false-negative). `fold_gate_v2.py` was already at parity with upstream.

## 6. THE MACHINE-CHECKABLE GUARD
`scripts/mechanism_hardening_check.py` asserts §1–§3 mechanically and exits 1 with the exact violation on
any drift. It is the self-enforcing backstop. Extend it — never bypass it — when a new invariant is added.

## 7. AUTHORITATIVE vs STALE FILES (avoid the wrong-file trap)
Trust, in order: **this doc → the hardening check → the code it checks (`mechanism_fold.py`,
`mechanism_build_anchor_graph.py`, `mechanism_scorecard.py`) → the bt_memory notes below.** A doc that
predates 2026-07-19 and disagrees with this file is stale — do not follow it. (A full CS-style
stale-sweep of all `docs/` is the planned second layer, to be mirrored from CS's own review.)

## 8. bt_memory pointers (the durable reasons behind the above)
- `cs-has-no-live-robust-json-extractor-…` — why extraction is ours, refuse-don't-fabricate.
- `mechanism-sync-contract-fold-gate-v2-and-result-form-cs-verified` — the gate + RESULT form.
- `neutral-dialect-applies-to-any-question-relayed-to-another-session-…` — identity/dialect off-repo.
- Recall: `context_loader.py "<task>" --memdir bt_memory`, then `bt_memory/MEMORY.md`.
