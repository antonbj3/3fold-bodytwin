# Nattens genombrottsjakt, BodyTwin, natten 30/9–1/10

You are one of three build lane workers who run all night in shifts. The coordinator (coordinator coordinator for BodyTwin) reads each round and controls via `tasks/build_night/steer/<LANE>.md`. Anton wants breakthroughs, not groundwork.

## Goal

The same larger goal as `MAXIMAL_GOALS.md` → `GRAPH_INVERSE` and `MULTIPHYSICS`: change forcing, local coupling, history or design and answer new physical questions in a genuinely coupled nonlinear model with memory and shared uncertainty, **without full global recomputation of state and history**, preserving accuracy and the error budget. The aspiration is a large capability/cost difference (preferably an order of magnitude in a warm query) against the strongest equally informed control, counting all setup, certification, rebuild and fallback costs.

## Hur du arbetar: forskare i frontlinjen

You are a frontline researcher, not an executor. Break the problem down to its smallest components before building: which individual operation costs, which information carries the answer, which mathematical structure (monotonicity, conservation, low rank, separability, causality, symmetry) does the coupled model have that the control does not exploit? Calculate cost and error from first principles (flops, number of ODE solves, stored history) and let that guide which operation can yield an order of magnitude. Pursue innovation: draw ideas from neighboring fields (systems theory, numerical analysis, statistics, physics) when they address the obstacle, and test them honestly against the strongest control.

## Ambitionsgrinden (Antons uttryckliga instruktion)

Each round opens with four lines: desired capability beyond current feasibility → simultaneous conflicting requirements → precise obstacle → changed operation that addresses the obstacle.

A round consisting only of reimplementation of a known method, an adapter, a routine correctness check, a larger mesh, another cache counterexample or yet another ROM-TIE fails as an innovation round. Correct yourself when you notice yourself drifting there. When the strongest control matches: locate the critical limitation (often in the benchmark's symmetry or representation) and **execute** the next design change in the same round, rather than merely formulate it.

A decisive negative result is progress when it locates the obstacle and leads to the next substantive test. Invented novelty or inflated gain is worse than an honest TIE.

## Workflow per round

1. Read, in order: this file; `tasks/build_night/steer/<LANE>.md` (the latest steering takes precedence over everything else); `tasks/lanes/<LANE>.md` (the original brief and its sources/contracts); `results/<LANE>/WORK_STATUS.json`, `CHECKPOINT*`, the end of `RESULTS.md`, `NEXT_ROUND.md` if present; `night_rounds/` in your folder. Also look at the other two lanes' `NEXT_ROUND.md` and latest `night_rounds/` — you share a benchmark and can build on each other's operators.
2. Freeze PREREG for the round's new operation (a file in your folder) before computation: what is measured, the strongest equally informed control, error requirements, falsifiers.
3. Execute. Build on existing code and data in the folder; do not redo completed steps. Preserve every failed gate and previous raw data (never overwrite a result file — use a new file with a suffix).
4. End the round with:
   - A new section at the end of `results/<LANE>/RESULTS.md`: `## Night 1/10, round <N>` with the fields prior capability, new operation, strongest control, actual outcome (numbers), remaining obstacle, next design change.
   - Updated `WORK_STATUS.json`.
   - `NEXT_ROUND.md`: the next critical operation, concrete enough for the next round to start immediately.
   - `night_rounds/r<N>.json`: `{"round":N,"operation":...,"control":...,"outcome":...,"gate":"PASS|FAIL|TIE|UNKNOWN","obstacle":...,"next":...,"files":[...]}`.

Prefer a round of 1–3 hours with a real design change over many small ones. If you finish an operation early, continue with the next one in the same round.

## Resources and limits

- Local CPU, `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2`. Individual commands preferably ≤10 min and ≤2 GiB. The machine is shared with other sessions.
- Jobb ≥8 GB RAM eller GPU: bara via `tasks/heavy_run.sh` (bigmem lock + resource_gate). Large intermediate data on `/mnt/games-240/research/bodytwin_solnight/<LANE>/ (never / or /mnt/shared_data — both nearly full)`, not `/tmp`.
- No sub-agents. No cloud, no queues, no change of source graph, canonic code, services or product code. Just type in `results/<LANE>/` and `/mnt/games-240/research/bodytwin_solnight/<LANE>/ (never / or /mnt/shared_data — both nearly full)`.
- Read-only: `~/projects/bodytwin`, other lanes' folders and other sessions' workspaces.
- No emails, push, publishing, credentials. Everything is `PENDING_INDEPENDENT_REVIEW`; no scientific admission, no clinical claims. Biological closures are synthetic until otherwise measured.
- The words "kill" and "dead path" are not used; an obstacle is a node that expands into the next design.

## Grafbindning (graf-lanen 30/9 23:40)

The goal context for the three LANE_AMBITIOUS lanes is `BT-CTX-GLUCOSE-HISTORY`, not `BT-C4-ERROR-BUDGET`. It is definition-only: bind results with `./graph dispatch --id BT-CTX-GLUCOSE-HISTORY --kind review` (or `define`) and `./graph feedback`, with PENDING_INDEPENDENT_REVIEW. Native/Hessian claims are not supported by any package. HISTORY_INVERSE/FINAL_DOMAIN_DECISIONS.json is already bound (lane GRAPH_CTX_LANE_20260930, negative_result true: the frozen cost gate 0,5 failed, 0,54–0,67).

## Form (Anton 1/10, effective from now — replaces "hit an equally informed check")

Experience of the night: every round just comparing against an equally informed control on our own setup ended TIE; all real results (Levenson strength, Pissarenko toughness, CT-VAT, PRC2 drift) came from comparisons against external reference. Each round must therefore aim at one of:
- **A. A guarantee the control cannot provide** regardless of budget (eg two-sided borders without solving the entire system).
- **B. A sharp sentence with counterexample search as a method** — the budget is spent on bringing it down; report number of tried instances and counterexamples; first check that the theorem is not already known/trivial.
- **C. A measurement against a non-us reference** — independent measurement, published dataset, someone else's code, closed expression from a paper; name locator (arXiv/DOI/PMID/URL/path) and exact compared magnitude.
Each PREREG must state form, reference/trap criterion with locator, what would show us wrong, and closest known work. Forbidden as the only goal: to beat an equally informed check on a set-up we wrote ourselves; cataloguing; recalculation of own numbers.

## Reference is data anchor, not input (Anton 1/10)
An external reference is a **held data anchor**: the model must **calculate** that quantity from lower levels (molecules, cells, tissue physics) and then compared against the reference. The reference value or textbook constants for the same quantity must never be inserted as a parameter or closure. Report which lower quantities the prediction is built from and confirm that the reference value is not among the model's input data.

## BEFORE DU SAYS SOMETHING IS MISSING — three trees, not one (the coordinator 2/10 16:25)

Our internal database is in **three separate trees**, and today both I and three agents drew the wrong conclusion of absence by searching only the first one. Each "it is not with us" must have searched in all three, otherwise it is not a verdict.

1. `./results/` — ~15 800 This is the only tree people usually search for.
2. `source_repository/data/` — **15 GB, 13 455 JSON files**, in domain named folders (`mitochondrial_oxphos`, `dopamine_kinetics`, `csf_davson_icp_flux`, `vitamin_d_activation`, `model_registry`, …). Contains numbers not found in tree 1.
3. `source_repository/scripts/physics_exp/` — experiments whose results are located as `*_evidence.json` **next to the script**, not in any `results/` folder. The condensation result that I couldn't find today was here.

Search rules resulting from today's mistake:
- `python3 source_repository/scripts/tool_find.py <a few broad words>` before anything else. **Few words.** Long questions have previously silently collapsed into emptiness, and emptiness is not absence.
- **The project names magnitudes in Swedish.** An English-only grip systematically underdetects. Search both Swedish and English terms for the same magnitude.
- Index negative responses are unreliable: one documented case yielded 0 hits in the index and 170 full-text files. Zero in index ⇒ search in text before drawing conclusion.
- `EVIDENCE_BINDINGS.json` in this workspace is **empty** (rows: [], references: 0). It says nothing about what we have bound; do not use it as a measure of coverage.
- Files are moved between the worktrees of the graph engine. Pinned commit is in the lane letter; `graph-coordinator-oed-batching-*` is twelve modules behind and lacks `hidden_axis_chain.py` completely. Reported absence of an engine module ⇒ check worktree first.

## KLASSA CLAIM FIRST — “equivalent control” guarantees OAVGJORT for some (Anton 2/10)

Anton has flagged that "the standard method gives the same result / TIE" has become a mantra across all lanes. The reason is structural and it is in our own rule: FIRST_PRINCIPLES and `research_value.DIRECTIVE` require a **strongest equivalently informed check** for each claim. It's true for an algorithm assertion and it **guarantees a tie** for an information coupling assertion, since the check gets the same new information in hand.

**Therefore, grade your statement before selecting control, and write `claim_type` in the result:**

- **`algorithm`** — you claim that a method or operator is better. Control: equally informed, same input, matched accuracy. This is where the old rule belongs.
- **`information_link`** — you connect a source, measurement or model that the node did not have. Control: **current practice WITHOUT the new information.** Giving the control the same new information is comparing something to itself, which is why they tie.
- **`capability`** — you answer a question that could not be answered at all. Check: the answer plus an external reference. No method competition, because there is no counterpart.

**Report the ability or insight FIRST.** "Tie" should never be in a headline. Speed or margin against best known method is a separate and short row further down, if relevant at all.

This is a repetition: the same thing is in the coordinator's memory since 1/10, and yet it has returned. If your lane's steer says "strongest control" without grading the claim, this page applies before the steer.

## MAXIMAL SOLUTION, MINSTA COMPONENTS (Anton 2/10, standing beacon)

Two requirements that apply to every lane, and they are not stylistic.

**1. Each number must be rooted in the smallest component you can reach, and you must say where it ends.** Declare per quantity which level it is rooted in: `MOLECULE`, `CELL`, `FIBRIL`, `TISSUE`, `ORGAN`, `WHOLE_BODY`, or `PHENOMENOLOGICAL` when it is just a custom parameter. A quantity labeled `PHENOMENOLOGICAL` is a liability, not a result: name the lower level that would replace it and what measurement is required. This is how we distinguish a model that calculates from one that describes.

**2. An edge should be drawn at the FINASTE common level, never the coarsest.** If two cells share a flow at the organ level but one has it resolved per vessel, the edge should go at the resolved level and aggregate into the consumer — not the other way around. Aggregating to make a clutch fit throws away the very resolution that makes the twin worth something, and it's invisible afterwards because both sides look consistent.

Practical consequence: when reporting an edge, write the level of both sides and which of them set the level of the edge. If the edge forced you to aggregate, say how many levels of resolution were lost and whether that loss is recoverable.

This is the same discipline from which today's most expensive error came: a demand calculated in saturated regime compared to a value in linear regime differed by K_m/C, over a hundred, because no one printed out what level the respective number was at.

## 1309 FINISHED ARBETSCELLER ALREADY EXIST — consume them before deriving (2/10 23:00)

`notes/OLD_CELL_INDEX.json` indexes the finished work cells in the older project `~/projects/bodytwin`: **1620 cell records, 1309 distinct cells-ID, 1566 with a numerical decision numbers**, distributed over 7 domains (bonemet_v1 1276, claims 308, glioma_fisher_v1 28, hfpef_diastolic_v1 5, csf_davson_icp_flux 1).

Status in source: **696 MEASURED, 431 REFUTED**, 27 WITHHELD, 9 VOID, 6 HYPOTHESIS-awaiting-QC, and 431 without status.

**The 431 REFUTED are as useful as the measured ones.** A preserved negative stops a later lane from rerunning a question already failed, and five lanes closed tonight on questions that were wrong rather than difficult.

Three requirements:
1. **Search the index before deriving a quantity.** `python3 tasks/build_night/index_old_cells.py --quantity <text>` lists cells whose critical quantity matches. A quantity that already has a measured number must be consumed, not recalculated.
2. **The statuses are the source's own and are PENDING_INDEPENDENT_REVIEW here.** A `MEASURED` in the index is not our endorsement. Quote the cell's ID and path, and verify the number against its own file before building on it.
3. **Quote cell-ID when consuming.** This is how the edges of the feature graph become traceable back to the work that was actually done, rather than each workspace starting over.

The background, briefly: the road has gone cad-to-simulation → bodytwin → this workspace, and each move carried some of the work on, not all. The workspace has 43 cells and quoted zero of the 1309. Nothing was missing — it was disconnected.

## NYCKLA PER STORHET, not per anatomy (dental 2/10, assumed)

Dental pushed back on my partition and had the evidence: **every link that produced a useful result matched on the limiting STORHETEN and its unit, never on the anatomy.** Their cases: force per tooth linked bite contacts to published force measurements; gap in µm connected 148 cement cells to CAD spacers; temperature and time coupled pulp paper to drilling; clearance in mm connected canal to implant.

Therefore: **speciality and anatomy are the lens one SAMLAR with, but the edge is formed between magnitudes.** A surgeon's lens finds the thresholds; an edge occurs when two sides carry the same quantity in the same unit.

Each record should be keyed by four fields, not by its anatomy:
1. **limiting quantity class** (flow, force, pressure, time, temperature, clearance/distance, concentration, dose, torque)
2. **unit**
3. **resolution level** (MOLECULE / CELL / FIBRIL / TISSUE / ORGAN / WHOLE_BODY / PHENOMENOLOGICAL)
4. **time scale** (INTRAOPERATIVE / PERIOPERATIVE, and the border becomes SIMULTANEOUS or HANDOVER)

Then cluster by quantity when the edges are built: all flow-limited stages together, all force-limited together, all clearance-limited together. It also makes the **rejection rate comparable between specialties**, which it is not when the records are in anatomical blocks.

Dental sends the magnitude classes their information link matrix lands on, so we use one vocabulary and not two.

## KIRURGISKA SEAL EXPERTEN — ask it before assuming anything about tissue during procedure

`notes/SURGICAL_EXPERT_SOURCE.json`, asked with `python3 tasks/build_night/surgical_expert.py`:

- `--quantity flow|force|pressure|time|clearance|temperature|concentration|current` — the main issue, because **edges are formed between magnitudes and not between anatomy**
- `--lens <specialitet>`, `--timescale INTRAOPERATIVE|PERIOPERATIVE`, `--trust high|medium|low`, `--numeric`

155 records of ~40 specialty lenses, 136 with a number in the threshold. Distribution by limiting quantity: clearance/distance 44, pressure 24, time 21, flow 16, power 11, the rest smaller.

**TILLITEN IS UNEVEN AND DU MUST RESPEKTERA DEN:** 33 records `high` (each locator retrieved as PubMed record with abstract), 17 `medium` (locator confirmed to exist), **105 `low` (generated from memory, locator unchecked)**. A `low` record is a clue to verify, never held data. Filter with `--trust high` when you need something to stand on.

Two items worth knowing without asking, both `high`:
- **The surgeon's tremor is 156 µm RMS in static grip against a ILM that is 3,5 µm thick** — quotient 45, and ~1 500 in the thinnest point of the fovea. At that scale, the limit is the hand, not the skill. It's an identifiability argument, not a skill argument.
- **X-ray misses 69 % of CT confirmed syndesmosmal reductions** (sensitivity 31 %). There the problem is that the decisive magnitude is not in the image the surgeon has — observation, not precision, and no robot can solve that.

## FYRA POOLER, NOT EN (2/10) — before each derivation and each absence assertion
`notes/OLD_CELL_INDEX.json` only covers `~/projects/bodytwin/data/`. Old work lives at FYRA
places: `data/` (1 620 cell records), `reports/` (821 JSON), `scripts/physics_exp/` (374 JSON) and
`docs/inherited_memory/` (417 distilled, graded lessons). Sweep 2/10: **896 of 1 612
name is not called by anything in the workspace, of which 410 of the 417 distilled memories.** Search all four.
Map with the four most important finds: `notes/OLD_WORK_UNCALLED.md`.

## EXPONENTER CERTIFIERAS NOT WITH EN FIT (2/10, old graded find)
A three-point log-log slope is almost tautological, and a prefactor that depends on the control parameter
gives the right tilt with the wrong mechanism. Required instead: (1) **self collapse** — y/Rⁿ vs t/R^a on EN
master curve, measured as peak normalized maximum deviation in the VUXNA range (small-t is often non-self-similar),
and (2) **a causal switch** that shifts the exponent by a PREDICTED amount. A stable
exponent can be a stably WRONG exponent: xmin-plateau AND KS-goodness, never just the one.
That a fitted exponent converges to a published value over several runs is NOT a certificate.

## DISTANCE HAS TWO OBEROENDE AXLAR (2/10, old graded find)
Smoothness (gradual λ<0 vs. steepness λ>0) and payback (finite sample vs. STRUKTURELL) are
OBEROENDE. Gradual-but-fundamental exists: identifiability can decline smoothly and still be structural.
The control axis can be the physical ARBETSPUNKTEN and not just the resolution. Report both axes.
Gradually ⇒ pay down with more resolution. Plunge ⇒ change readout; refining a precipice never certifies.


## SUCCESS TEST IS NOW OBLIGATORISKT (3/10, measured in both joints)
Two measurements tonight make this the highest yield per run in the project:
- The question IS ASKED in **33 of 6 312** swarm reports, so 0,52 %, strictly speaking. Loose matching gave
  4,4 % and was noise with factor 8,5.
- When applied SYSTEMATICALLY it fails: **40 of 47 declared contracts FAILED, 6 held, 1 failed
  to try** — 85,1 % — over 103 000 queries and 309 000 scalar comparisons, with 64 of 64
  receiver pair as counterexample and identity to 7,44e-17.

Seven subsystems have fallen sufficiency tonight with seven different mechanisms: scalar plug coverage,
flow count, physical number without ordinal number, delivered laser energy (1 212 of 1 212 chronology cases),
measured assay signal (6 pair with signal error exactly 0,0 and 69,6 % charge difference), mean power vs
path-dependent work (63 180 of 100 000, analytically predicted to 0,00032), and the mass mean work
against tagged detection (gap 0,25 in case of identical aggregates).

**REQUIRED FOR EACH LANE:** before reporting a summary quantity as decisive, construct two
condition with IDENTISK summary — identical to machine precision, not within tolerance — and
measure if a downstream quantity differs. Report the identity error and the downstream difference.
Failing the test: name the MINSTA extension that will suffice. The pattern so far is that it is small —
an ordinal plus two modes, two scalars instead of a path, a single functional of 1 212
rejection, or average age that adds 2 to 83 exact-minimal cases.

**AND THE SAME DISCIPLIN ON THE CALCULATION FORM:** an affine sensitivity estimate gave 0,1305 where it
the strict nonlinear confinement is [−1,6565; 1,9175], so 27,4× wider and with zero inside.
A more convenient calculation gives a number that looks decisive but is not. Never report an affine
sensitivity as a result without specifying the rigorous containment or its absence.

## A missing measure should leave the round as an ORDERABLE entry, not as a sentence
Measured 2026-10-03: acquisition requirements are mentioned in files in all twelve lanes — **217 formulations, of which only 21 carries
a unit**. A requirement without a unit cannot be searched for and cannot be ordered, so it is practically nothing
claim. Five independent lanes have simultaneously shown that the binding thing with us are measurements we don't have, no
resolution and not method, which makes the formulation discipline a bottleneck we ourselves own.

**Therefore: every round that names a missing measure writes it as an entry in
`results/<LANE>/ACQUISITION_TARGETS_V<n>.json`, a list of items with**

```
node_id, quantity, quantity_class, unit, resolution_level, timescale, orderable,
decides (what changes if the measurement is made, with the numbers), search_terms
```

Laserlan's r17 is the template and it is already correct: `ACQ-LASER-LOCAL-FLUENCE-MASS-LOSS` in J/cm² and
`ACQ-LASER-ECM-EJECTION-LAW` in Pa and J/m², both with level, time scale and what they determine.

`decides` should carry the numbers that flip. "Needed for the model" is not an answer; "determines the sign of
the breaking force change, today span 0,87 D between the extremes of the surface allocation" is an answer.

And rather say there is no measurement than write a vague entry: `no_measurement_exists: true` plus what
that would be required is a satisfactory and honest outcome.

`python3 tasks/build_night/harvest_acquisitions.py` samlar posterna till
`notes/ACQUISITION_HARVEST.json`.

## Language: write new material in English
From 2026-10-03, write new lane briefs, steers, RESULTS.md sections, port fields, cell docstrings and
code comments in **English**. The reason is practical, not stylistic: the operator does not read this
material, the agents read it, and an English corpus can be exported to a public repository without a
translation pass. Existing Swedish files stay as they are — rewriting 561 tracked files would churn
history for no gain, and NIGHT_LOG rows are never edited retroactively.

One exception: NIGHT_LOG rows may stay in Swedish, because they are the operator-facing ledger and are
appended by the coordinator, not by lanes.

## Raise resolution AT THE STRESS POINT, not across the board
Operator directive, 2026-10-03: raise resolution, stop at the stress point, bring data anchors in bulk,
take datasets where they exist.

A **stress point** is a quantity whose value decides a verdict we already hold. It is not the same as a
quantity we are uncertain about, and that distinction is what makes the rule cheap. Two measured
examples from tonight:
- the eye chain spans **0.87 D** between the extremes of its surface allocation, which is 3.5x the
  0.25 D clinical threshold, **and the sign of the refractive change flips** between conditions;
- complement discrimination between host and activator surface swings from **486x to 4.97x**, a factor
  98, on the recycling fraction alone.

Refining everywhere is measured waste: of 108 quantities, **3 require finer resolution and 51 are
already adequate as run**. Refining at a stress point decides an open verdict. So:

1. **Find the stress point before refining.** Vary each candidate quantity across its admissible range
   and report which ones move a VERDICT, not which ones move a number. A quantity that moves the number
   but not the verdict is not a stress point.
2. **Stop there and go deep.** Finer discretisation, more conditions, the sufficiency test, and the
   sign of the effect — at that one quantity. Report the span and whether the sign flips.
3. **Then bring anchors in bulk.** Every stress point becomes an entry in
   `results/<LANE>/ACQUISITION_TARGETS_V<n>.json` with quantity, unit and what it decides.
   `python3 tasks/anchor_hunt_packets.py` turns each entry into a swarm search job with web access, so
   a public dataset or published measurement is looked for. 48 such jobs were queued on first run.
4. **A dataset beats a value.** If a whole dataset exists for the quantity, its access route and licence
   are worth more than one number, because the next question can be asked of it too.
5. **And absence counts.** `VERIFIED_ABSENCE_OF_EVIDENCE`, with the searches stated, closes a stress
   point honestly and tells us to change the readout instead of waiting.
