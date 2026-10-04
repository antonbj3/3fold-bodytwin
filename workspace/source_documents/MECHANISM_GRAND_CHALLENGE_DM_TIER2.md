# MECHANISM — Grand Challenge DM Tier-2 unblock attempt: DIAGNOSED-GAP (honest-negative)

**Goal:** get the full muscle+ligament (Tier-2) knee-contact cert on Grand Challenge subject DM, to complete the
flagship same-subject in-vivo validation beyond the muscle-free Tier-1 lower bound (`MECHANISM_GRAND_CHALLENGE_DM.md`).

**Outcome: the unblock did NOT succeed. Tier-2 on DM remains blocked.** The attempted fix (prescribe the
poorly-marker-observable secondary knee/PF DOF, so SO resolves only the primary flexion) still yields a
non-converging static optimization on every frame. Evidence, parsed by the coordinator directly from the run's
own SO log (`data/msk_smoketest/DM_ngait_og1/so/run_so.log`, 4899 lines):
- **0 frames converged**; **247 "[warning] The model appears too weak for static optimization"** + 498
  "OPTIMIZATION FAILED / Ipopt: Maximum iterations exceeded (status -1)".
- Muscles pinned at activation→1 (upper bound) across the set; final-frame Performance ~3.1e8, constraint
  violation ~3.2e6 — a wholly infeasible solve, not a marginal one.

**Root cause (confirms the flagship's diagnosis, now with the prescribe-DOF variant ruled out):** DM's TKA
knee+patellofemoral carry a 149-bundle `Blankevoort1991Ligament` passive apparatus (aggregate scalar tension up
to ~18,700 N at a pose) that the muscles+reserves must balance every instant; prescribing the secondary DOF
removes the marker-observability problem but not this passive-force load. The model is genuinely too stiff/weak
for SO as shipped.

**What Tier-2 on DM would actually need (the concrete path, for a future run):** model surgery beyond a
prescribe-DOF fix — e.g. (a) reduce/soften or selectively disable the ligament apparatus during the SO solve
(then add its contribution back analytically at the cut), (b) weld rather than prescribe the secondary DOF to
remove them from the optimization entirely, and/or (c) subject-specific reserve re-tuning. Each changes the
model materially and must be validated separately; this is a project, not a parameter tweak.

**Process note:** the agent running this stalled mid-debugging its own Tier-2 ligament-tension lookup and was
stopped; this write-up is authored by the coordinator from the on-disk SO log, so the honest-negative is
recorded rather than lost.

> ~~**Net:** the flagship **Tier-1 same-subject in-vivo confirmation stands unchanged** (twin over-predicts DM's
> measured knee contact force by 1.5–1.8× at the muscle-free floor; true value ≥ that since muscle/ligament forces
> only add).~~ **SUPERSEDED 2026-07-28 — DIRECTION REVERSED.** `MECHANISM_GRAND_CHALLENGE_DM.md` retracted this
> headline on 2026-07-22 (units bug: the raw eTibia CSV `Fx,Fy,Fz` columns are **pounds-force**, not Newtons; the
> ×4.4482 lbf→N conversion was never applied — independently confirmed by a neutral adjudication agent in
> `MECHANISM_DM_UNITS_ADJUDICATION.md` via a mass-and-literature-independent JCF/GRF ratio check, 2.27–2.83× in
> all 7 trials under the lbf reading vs a mechanically-impossible 0.51–0.64× under the Newton reading). Corrected,
> DM's true peak is **≈296 %BW** (not 58–68 %BW), so the twin's muscle-free Tier-1 floor (104 %BW) **UNDER**-predicts
> DM by roughly **2.7×**, not over-predicts by 1.5–1.8×. This Tier-2 unblock attempt itself (the SO-infeasibility
> diagnosis) is unaffected by the units bug — only the "Net" framing line above is corrected. The separate
> cross-cohort over-prediction (LaiArnold static-optimization on subject2 vs OrthoLoad, `MECHANISM_TRUST_LEDGER.md`)
> is a different model/subject/quantity and stands on its own, unretracted.

The full Tier-2 number on DM is a **diagnosed-gap** — achievable only with the model surgery above.
**Confidence tier:** diagnosed-gap (method-blocked); the SO-infeasibility evidence is directly measured from the log.
