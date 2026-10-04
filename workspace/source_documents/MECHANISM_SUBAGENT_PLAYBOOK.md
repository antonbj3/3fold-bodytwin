# MECHANISM SUBAGENT PLAYBOOK — how the coordinator fans out to USE the archive

You (mechanism) are one session; the sibling CS project runs ~10 lanes. You compensate by **spinning subagents
constantly**. This is the concrete how. The work is **using the 449-source archive to build and CLOSE cert-cells**
— NOT hunting new data (discovery is done; 26 waves). New-data hunting is Pattern D only, for named FENCED holes.

## HARD RULES — apply to EVERY subagent, no exceptions
1. **Return-only.** Prompt every subagent: *"Do NOT use the Write tool. Do NOT write any file or memory. Return your
   result as your FINAL MESSAGE (raw JSON). Your output is data to the coordinator."* Subagents don't get memory
   auto-injected and have no memory mandate — this keeps their file tools idle so nothing lands in the wrong place.
   YOU (coordinator) persist the vetted result into `bt_memory/` + the graph. Never `~/.coordinator`.
2. **Output is a HYPOTHESIS — QC before trust.** Subagents are both false-positive AND false-negative prone. Run the
   QC checklist (below) on every return before folding anything. A subagent claiming "verified" is not verification.
3. **Isolation in the prompt:** this repo only (tool_find/recall → here); claims → the mechanism graph; never write
   coordinator memory; neutral terms in any /tmp state. (This is `docs/MECHANISM_PREAMBLE.md` — prepend it or its essence.)
4. **Measured ≠ constructed.** A leg value must be MEASURED from the data, never a plausible number. A model's
   prediction is a FILLED prior, never a decorrelated leg against its own training data (common-mode).
5. **Pace:** one wave per cron tick (a few concurrent, cap ~5); heavy data jobs → `resource_gate.py check` + flock
   `/tmp/gpu.lock` first. Commit each vetted unit; never push.
6. **⚠ `/mnt/data_root` silently corrupts bulk small-file writes** (loop35 mount — found 2026-07-18: `unzip -d`
   of 11k PNGs silently zero-filled ~0-40% of them, no error, source zip clean). Every execution subagent that
   touches archives MUST **read files in-memory** (`zipfile.ZipFile(...).read(name)` + `io.BytesIO`) — never
   bulk-extract many small files onto /mnt/data_root. Large SINGLE files (h5ad/csv/tar/gctx) stream fine; the
   hazard is bulk many-small-file extraction. MD5/shape-verify anything you do extract. And avoid MONOLITHIC
   multi-GB no-range downloads (L1000 21GB, ToxCast 7.5GB) inside a tick — they need a dedicated long download,
   not a tick-bound subagent. [[mnt-data-root-silently-corrupts-bulk-small-file-writes]].
7. **⚠ No background-and-defer.** Measured 2× (both permutation-heavy microbiome agents): a subagent that
   backgrounds a long script (big `n_perm`, nested CV) hits the turn boundary and returns a **placeholder**
   ("waiting for the monitor to notify me…") instead of the result JSON. That placeholder is a NON-RESULT —
   never fold it. Two defenses: (a) PROMPT every compute-heavy subagent *"do NOT background a long job and
   defer to a monitor; run it synchronously — reduce n_perm (e.g. to 50-200) if needed — and return the
   ACTUAL numbers in your FINAL message"*; (b) when a placeholder still comes back, `SendMessage` the agent to
   load any checkpointed partial results or re-run the decisive statistic once at small n_perm and return the
   JSON (the nudge has recovered the real result both times). [[subagent-backgrounds-long-job-returns-placeholder-nudge-for-actual-result]].

## THE QC CHECKLIST (coordinator runs this on every subagent return)
- **Source-ids exist:** every cited id is a key in `data/MECHANISM_SOURCE_REGISTRY.json`
  (`scripts/mechanism_qc_cell_designs.py`, or grep). A missing id = fabrication/typo → reject or fix to nearest verified.
- **Decorrelation is real, not common-mode:** the ≥2 legs use DIFFERENT physics/data. Two legs sharing a model, a
  training set, or a backbone = n_eff≈1 fake decorrelation → flag, don't credit.
- **Anchor is external + held-out:** the anchor is an INDEPENDENT measured GT not used in the fit. An anchor derived
  from the same model/data it checks = a tautology gate → certifies nothing → reject.
- **Regime stated:** the claim names the regime it holds in (a fixed threshold on a noisy statistic is regime-blind).
- **Honest-negative = PASS:** a booked negative with a mechanism (no anchor, no 2nd leg, regime-gate closed) is a
  correct result, not a failure — keep it FENCED, don't force a cert.
- **Numbers reproduce:** for execution results, re-run the leg statistic yourself on a sample; confirm held-out split.

## PATTERN A — cert-DESIGN (map archive → a cell's decorrelated legs + anchor)   [wave-1 used this]
When: a build-goal cell has no cert_design yet, or a STUB needs filling.
Spawn (return-only): give the cell's hidden state, the ANCHOR node schema, the cert model (≥2 decorrelated legs +
held-out external anchor + occluded-truth→consequence rule), and point at `MECHANISM_SOURCE_REGISTRY.json` + `WAVE_PLAN.md`.
Ask for: `{node{...}, cert_design{hidden_state, occluded, legs[{leg,measures,data_source_ids,decorrelation}],
anchor{measured_gt,source_id,held_out}, composition_mode, couples_to_inherited}, honest_gaps}` — cite exact source-ids.
QC: full checklist. Fold into `cell_designs_wave1.json` → `mechanism_build_anchor_graph.py` → re-QC → commit.

## PATTERN B — cert-EXECUTION (measure a cell's legs, fold vs the anchor)   [how a cell CLOSES]
When: a cell is designed (OPEN) and its data is reachable. This is the actual measurement that advances status.
Spawn (return-only): one subagent per leg (or per cell if light) — *load the cited source, compute the leg's
statistic in the stated regime, return the number + how it was computed + the held-out split used.* Then a fold step:
compare the legs on the shared hidden state, check against the held-out anchor, emit a measurement claim with a
`target{file,object_id,field_path,value}` for `measurement_claim_sync` to fold → status advances OPEN→(PROVEN|REFUTED).
QC: reproduce one leg's number; confirm the anchor was genuinely held out; confirm agreement isn't self-selected.
Cheapest genuine first-close: **MOL-SAMESAMPLE-PATCHSEQ-CALIBRATION** (LOW risk, open DANDI 000020/000023).

## PATTERN C — adversarial VERIFY (try to REFUTE a cell before trusting a pass)
When: a cell reports a pass, or a design looks too clean (swept-clean = a hypothesis until an independent red-team
fails to break it). Spawn N skeptics (return-only), each with a distinct lens, each prompted to REFUTE (default to
refuted=true if uncertain): (i) common-mode — do two legs secretly share a model/dataset? (ii) tautology-anchor — is
the anchor derived from what it checks? (iii) regime-blind — does the threshold survive a regime shift? (iv) held-out
— was the anchor actually excluded from the fit? Kill the pass if ≥majority refute. Book survivors.

## PATTERN D — targeted gap-fill (FENCED holes ONLY — the only new-data hunt)
When: a cell has a NAMED missing leg/anchor (an honest_gap in the inventory, or a FENCED hole in `MECHANISM_HOLES_SEED.jsonl`).
Spawn `watertight-researcher` (inherits the method): find an OPEN dataset for THAT specific gap. Discipline: band
BT-HOLD, honest-negative = PASS (a confirmed "no open source" is a correct FENCED result), NEVER fabricate a URL,
flag licenses. Do NOT broaden into general discovery — we already have the archive; this is surgical.

## LEDGER
Book every closed cell, every honest-negative, and every refuted design in `bt_memory/` (+ `LEDGER.jsonl`) so
mechanism1 (cl) and mechanism2 (cl2) both inherit it. In-repo only.
