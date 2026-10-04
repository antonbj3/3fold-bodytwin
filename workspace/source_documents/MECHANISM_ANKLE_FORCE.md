# MECHANISM ANKLE FORCE — talocrural contact + Achilles tendon force vs published anchor (2026-07-21)

## ⚠ CORRECTION (2026-07-21) — READ THIS BEFORE THE NUMBERS BELOW

`validate_ankle_force.py` calls `static_opt_knee.knee_crossing_muscles_and_forces` directly,
unmodified, for BOTH the ankle-r cut (this doc's own headline) AND the knee-r reproducibility cut
(§2.6) — the SAME shared function `docs/MECHANISM_STATIC_OPT.md`'s correction banner describes. The
self-computed ankle contact-force number was the EXACT NEGATIVE of the correct muscle-force
direction, masked by a vector-norm step. Full diagnosis/fix: `docs/MECHANISM_SIGN_BUG_AUDIT.md`,
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`. Fixed at its source (one-line branch swap in
`static_opt_knee.py`). Re-ran `validate_ankle_force.py` after the fix to confirm directly:

| | published (bug artifact) | **corrected** (re-measured, direct re-run) | ratio vs published anchor |
|---|---:|---:|---:|
| self-computed ankle contact (R − Σmuscle-crossing) | 286.68 %BW | **484.49 %BW @ t=0.58s** | 0.601 → **1.016** |

**The corrected number now matches this doc's OWN already-published, bug-immune official
`opensim.JointReaction` cross-check (484.48 %BW) to 0.0014%** — the tightest agreement of any joint
in this family, and further confirmation the correction (not the original) is right.

**This joint's honest headline is DIFFERENT IN CHARACTER from the knee/hip correction**: the ankle's
anchor (§0) is a weaker-tier published cadaveric/model-literature band (390-540 %BW, not a direct
in-vivo instrumented measurement), and the corrected number **lands almost exactly on it** (ratio
1.016, both self-computed and JointReaction) rather than over-predicting it — there is no
"over-prediction vs in-vivo" story to tell here the way there is at the knee/hip, precisely because
this anchor was never an in-vivo one. The Achilles/triceps-surae numbers (320.26 %BW scalar sum,
318.83 %BW vector-magnitude cross-check, §5.2) are **UNCHANGED, verified unaffected**: a uniform sign
flip applied identically to all 3 triceps-surae muscles negates their summed vector, but
`np.linalg.norm` of a negated vector equals the norm of the original — this specific downstream
quantity is structurally blind to this bug, confirmed by re-running (320.2649... %BW before and
after, bit-identical).

**§2.6/§7's "FREE reproducibility gate" (comparing this script's own knee-r-cut re-derivation
against the already-published knee cert numbers) now correctly reports a 67.7% MISMATCH on the
self-computed candidate** when re-run — EXPECTED, not a new defect: `PRIOR_KNEE_SELF_COMPUTED_PCT_BW
= 233.20` (line ~156) is a hardcoded historical constant describing the OLD buggy knee number; this
script's own re-derivation now correctly returns 391.10 (the corrected number), so it no longer
matches that stale constant. This hardcoded constant should be updated to 391.10 in a follow-up so
this script's own internal gate exits clean again (not done here — out of this remediation's
explicit scope, flagged in `docs/MECHANISM_SIGN_BUG_REMEDIATION.md`).

**§5.1's "striking ratio similarity... corroborates the knee cert's diagnosis" (differentiation
scheme) is SUPERSEDED** — same reattribution as the knee/hip docs: the shared ~1.68-1.69× ratio
recurring at knee AND ankle corroborates a shared CODE bug (this exact function, called identically
at both joints), not shared differentiation scheme.

Everything below this banner is preserved as the ORIGINAL, pre-correction analysis (historical
record) — read the headline table, §2.6, §5, and §7 with the correction above in mind.

---

Replicates the VALIDATED knee method (`docs/MECHANISM_JOINT_FORCE_VALIDATION.md` +
`docs/MECHANISM_STATIC_OPT.md`: Static Optimization + a Newton/virtual-work method immune to the
OpenSim 4.6 ID-tool bug) for the **ankle (talocrural)**. Every number below is machine-measured this
session (`scripts/msk/validate_ankle_force.py`, exit 0), not recalled. Isolation respected:
`.venv-msk` only, LabValidation data + a sibling agent's already-computed Static-Optimization/
JointReaction output read in place (never written to), no git commit/push.

## Headline result [SUPERSEDED — see correction banner above; preserved for the historical record]

| quantity | method | peak (%BW) | vs published |
|---|---|---:|---:|
| Ankle (talocrural) contact force | self-computed (R − Σmuscle-crossing) | ~~286.68~~ **→ 484.49 (corrected)** | ~~0.601~~ **→ 1.016** |
| Ankle (talocrural) contact force | official `opensim.JointReaction` | **484.48** | ratio 1.016 |
| Ankle (talocrural) contact force | PURE kinematics+GRF reaction (no muscles) | 107.50 | ratio 0.226 |
| **Achilles/triceps-surae tendon force** [unaffected, verified] | SO tension, scalar sum (gasmed+gaslat+soleus) | **320.26** | ratio 0.821 |
| PUBLISHED ankle/hindfoot anchor (3 independent studies, walking) | — | 390–540, mid 476.7 | 1.0 |
| PUBLISHED Achilles anchor (walking) | Giddings et al. 2000 | 390.0 | 1.0 |

**[SUPERSEDED framing, kept verbatim for the record]** Both muscle-inclusive ankle-contact candidates
and the Achilles estimate land within the pre-registered plausible band (ratio 0.30–3.00× the
published point) — PASS, not a knife-edge. The official-JointReaction candidate (484.48 %BW) lands
within 2% of the published midpoint (476.7 %BW); the self-computed candidate (286.68 %BW) and the
Achilles estimate (320.26 %BW) undershoot their anchors by a similar ~18–20%, in the same direction
and rough magnitude as the already-published knee cert's own self-computed-vs-JointReaction gap (§5
below) — a consistent, not cherry-picked, pattern across two independent joints. **[Corrected: the
self-computed candidate now matches JointReaction to 0.0014%, ratio 1.016 — it does not undershoot;
only the Achilles estimate's 0.821 ratio was ever a genuine (unaffected) undershoot.]**

## 0. HONEST TIER — read this before the numbers (stated up front, per the task's own framing)

**This is a WEAKER anchor tier than the knee/hip OrthoLoad certs.** No in-vivo instrumented-implant
program has ever existed for the human ankle — verified absent again this session (OrthoLoad's
`data/external/orthoload/` ships only `hip/knee/shoulder/spine_vbr/spine_fixator/cams_knee`, never
ankle), consistent with the pre-existing honest gap already logged in
`docs/MECHANISM_MSK_BUILD_PLAN.md` (lines 222–223, 292–294): *"Ankle (talocrural+subtalar): NONE
PUBLIC — verified, nothing found in OrthoLoad or elsewhere. No in-vivo anchor exists."* The
comparison below is therefore **model-vs-model + cadaveric/quasi-static/finite-element literature**
(1977–2000 studies), not a direct in-vivo measurement — a plausibility/ballpark check, not a
validation with the OrthoLoad knee cert's epistemic weight.

## 1. Published anchor — verified LIVE this session, not recalled

`WebSearch` was unavailable (session-wide quota exhausted — a pre-existing constraint already
flagged in `docs/MECHANISM_MUSCLE_AUDIT.md`). Verification instead used the publishers' own
structured metadata/abstract APIs directly (`curl` to NCBI eutils `esearch`/`esummary`/`efetch` and
the Crossref REST API) — machine-fetched real abstract text with DOI+PMID cross-confirmed, not a
narrated search snippet and not training-data recall.

**Ankle/hindfoot joint contact force, level walking — three independent studies, different
eras/methods/authors, all real, all quoting a specific number:**

| source | method | quantity | peak (%BW) |
|---|---|---|---:|
| Stauffer, Chao & Brewster 1977. *Clin Orthop Relat Res* 127:189-96. DOI [10.1097/00003086-197709000-00027](https://doi.org/10.1097/00003086-197709000-00027), PMID [912978](https://pubmed.ncbi.nlm.nih.gov/912978/) | quasi-static, high-speed cine film + force plate, normal subjects | ankle joint compressive force | **500** |
| Procter & Paul 1982. *J Biomech* 15(9):627-34. DOI [10.1016/0021-9290(82)90017-3](https://doi.org/10.1016/0021-9290(82)90017-3), PMID [7174695](https://pubmed.ncbi.nlm.nih.gov/7174695/) | 3-D force-equilibrium model incl. post. tibial + peroneal muscles ("Mark II"), 7 adult male subjects | talocrural (Tc.) joint resultant peak force | **390** |
| Giddings, Beaupré, Whalen & Carter 2000. *Med Sci Sports Exerc* 32(3):627-34. DOI [10.1097/00005768-200003000-00012](https://doi.org/10.1097/00005768-200003000-00012), PMID [10731005](https://pubmed.ncbi.nlm.nih.gov/10731005/) | contact-coupled finite-element foot model + cineradiography + force plate | peak **talocalcaneal** (subtalar, not talocrural proper — see caveat below) joint load | **540** |

Verbatim quotes (machine-fetched PubMed abstracts, not paraphrased): Stauffer — *"Compressive force
across the ankle joint rose to about 5 times body weight during the latter part of stance phase."*
Procter & Paul — *"Tc. joint force = 3.9 [times body weight]."* Giddings — *"The model predicted peak
talocalcaneal and calcaneocuboid joint loads of 5.4 and 4.2 body weights (BW) during walking."*

**Precision disclosure** (forced adversary on my own anchor construction, not blurred): Giddings'
most-quoted number (5.4×BW) is technically the **talocalcaneal (subtalar)** joint, anatomically
adjacent to and closely coupled with the talocrural joint but not the identical joint this script's
`ankle_r` cut computes. Restricting to the two studies that measure the talocrural joint itself
(Stauffer 500, Procter & Paul 390) gives a narrower mean of **445 %BW**; including Giddings'
closely-related hindfoot number (as the task's own stated "~4.5-5.5×BW" ballpark implies) gives the
**390–540 %BW band, mean 476.7 %BW** used as this script's comparison point throughout. Both framings
are reported; neither is hidden.

**Achilles/triceps-surae tendon force, walking** — the one verified walking-specific number found
this session: Giddings et al. 2000 (same paper), *"The maximum predicted Achilles tendon forces were
3.9 and 7.7 BW for walking and running"* → **390 %BW walking anchor**.

**Komi PV et al.'s in-vivo buckle-transducer work** (*J Biomech* 1990, DOI
[10.1016/0021-9290(90)90038-5](https://doi.org/10.1016/0021-9290(90)90038-5), PMID
[2081741](https://pubmed.ncbi.nlm.nih.gov/2081741/); *Clin Sports Med* 1992, PMID
[1638639](https://pubmed.ncbi.nlm.nih.gov/1638639/)) is real, directly relevant, and methodologically
a **stronger** anchor class than any model estimate (a literal surgically-implanted tendon-force
transducer) — but its abstract's only specific number ("as high as 9 kN, corresponding to 12.5 times
body weight") is **NOT walking-specific**; the companion 1992 paper explicitly attributes the highest
AT loads to "certain activities (e.g., hopping)." **Deliberately NOT used as a walking point
estimate** — misattributing an extreme-activity peak to this trial's activity would be a real error,
not caution (the same discipline this repo already applied when excluding OrthoLoad's gait-aid-
assisted trials from the knee cert's walking comparison). Cited here as qualitative support only:
direct in-vivo evidence that triceps-surae forces of several ×BW are physiologically real and
measurable, not a model artifact.

**Running cross-check (directional context only — `walking1` is a walking trial, NOT compared
directly)**: Giddings 2000 running Achilles 7.7×BW / ankle-adjacent joint loads 7.9–11.1×BW; Scott &
Winter 1990 (*Med Sci Sports Exerc* 22(3):357-69, DOI
[10.1249/00005768-199006000-00013](https://doi.org/10.1249/00005768-199006000-00013), PMID
[2381304](https://pubmed.ncbi.nlm.nih.gov/2381304/)) — *"Achilles tendon force: 6.1-8.2 BW... ankle
bone-on-bone compressive force: 10.3-14.1 BW."* Two independent running estimates agree with each
other and are higher than the walking numbers in the expected direction — the literature is
internally consistent, supporting (not proving) the walking numbers used above.

## 2. Method (identical architecture to the validated knee cert, replicated not reinvented)

1. **PURE kinematics+GRF Newton's-law reaction force** at the ankle-r free-body cut — bodies distal
   to `ankle_r`, found live by BFS on the model's own joint tree (`talus_r, calcn_r, toes_r` —
   confirmed, not assumed). Bug-immune by construction: never calls
   `InverseDynamicsTool`/`InverseDynamicsSolver`/`JointReaction`'s generalized-force routine (the
   ~1000-1800× knee/hip defect documented in `docs/MECHANISM_MSK_ELASTIC_BAND.md` §4).
2. **Static Optimization muscle tensions — reused READ-ONLY**, not re-run. The sibling knee cert's
   SO/JointReaction output (`data/msk_smoketest/subject2_walking1/static_optimization/{so,jr}/`,
   same model/trial/OpenCap templates) already contains everything this script needs: the converged
   SO solution does not depend on which joint is analyzed downstream, so re-running the ~158-frame
   Ipopt optimization would reproduce identical numbers at real compute cost for zero new
   information. This script instead **independently re-verifies** the sibling's convergence/
   activation gates itself (never trusting the prior doc's prose) before building anything on them,
   and logs file provenance (mtime, size) for the record — see §4.
3. **Self-computed ankle contact force** = R_vec(ankle-r cut) − Σ(ankle-crossing muscle force
   vectors, from LIVE path geometry × SO tension) — same subtraction logic as the knee cert. All 11
   anatomically-expected ankle-crossing muscles (tibant/edl/ehl/tibpost/fdl/fhl/perlong/perbrev/
   gasmed/gaslat/soleus) were detected by the geometry-based crossing detector with **zero
   assumptions** — and, forced further (§4), detected in **100% of all 158 frames**, not merely
   "at least once."
4. **Achilles/triceps-surae tendon force** = gasmed_r+gaslat_r+soleus_r SO tension: SCALAR sum
   (primary — the merged tendon carries one axial tension, the literature-comparable convention) +
   VECTOR-magnitude-of-summed-force-vectors (cross-check on how parallel the 3 muscles actually pull,
   reusing the same live-geometry crossing detector, no extra machinery).
5. **Official `opensim.JointReaction` cross-check**, read from the same reused `.sto` — its template
   already computed reactions for `joint_names=ALL`, so `ankle_r` columns were sitting in the file
   unread until this script (the knee cert only ever extracted `walker_knee_r`).
6. **FREE reproducibility gate**: the SAME reused SO/JR files, run through the SAME code pointed at
   the KNEE-r cut instead, must reproduce the ALREADY-PUBLISHED knee cert numbers — an
   over-determination check on this script's own correctness, decorrelated from whether the
   ankle-specific literature anchor is right.

## 3. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| GRF sanity (peak vertical R/L GRF) | 85–140 %BW | 108.94 / 107.64 %BW | PASS |
| Whole-body Newton residual (RMS, x/y/z) | < 8 %BW | 2.75 / 4.10 / 1.17 %BW | PASS |
| SO convergence + activation sanity (re-verified live, not trusted from prior doc) | 0 NaN, bounds, covers t=1.50s | 0 NaN, [0.010, 0.627], covered | PASS |
| Ankle-r crossing-muscle anatomical anchor | 11/11 expected muscles found | 11/11 found | PASS |
| **Ankle-r crossing-muscle per-frame coverage** (forced adversary, not just "found at least once") | triceps-surae detected in 100% of 158 frames | 1.0 / 1.0 / 1.0 | PASS |
| Sensitivity to smoothing window (70–130 ms) | peak range < 10 pct-pts | **0.02 pct-pts** | PASS |
| Free reproducibility (knee-r cut vs already-published knee cert, 3 candidates) | < 1.5% relative diff | 0.000% / 0.000% / 0.001% | PASS **[re-run post-fix: self-computed candidate now MISMATCHES this row's stale 233.20 target at 67.7% diff — EXPECTED, target constant is now stale, see top banner]** |
| Ankle contact regime sanity (ratio vs published midpoint) | in (0.30, 3.00) | 0.601 (self-computed), 1.016 (JointReaction) | PASS |
| Achilles regime sanity (ratio vs published point) | in (0.30, 3.00) | 0.821 | PASS |

Overall pipeline validity: **PASS**. Anatomical crossing-detector anchor: **PASS**.

## 4. Sibling SO/JR output reuse — provenance (not silently assumed fresh)

| file | size | mtime (UTC) |
|---|---:|---|
| `walking1_StaticOptimization_force.sto` | 351,133 B | 2026-07-21T11:06:40 |
| `walking1_StaticOptimization_activation.sto` | 302,140 B | 2026-07-21T11:06:40 |
| `walking1_JointReaction_ReactionLoads.sto` | 541,880 B | 2026-07-21T11:06:43 |

Logged, not just asserted fresh. The convergence/activation gate (§3) was re-run against these exact
files this session (`activation_max=0.627`, `0 NaN`, `158/158 frames`, matching
`docs/MECHANISM_STATIC_OPT.md`'s own published numbers exactly) — an independent re-check, not a trust
of the prior doc's prose. One inherited, already-disclosed limitation carries over unchanged since
it is the same SO run: **pelvis residual force peak = 172.96 N (22.6% BW), above the pre-registered
"good" band** (no RRA was run first) — flagged in `docs/MECHANISM_STATIC_OPT.md` §4/§8.2 already, not
a new finding, and joint-level (hip/knee/ankle/subtalar) reserve usage stays low (≤12.5% of a small
optimal_force) throughout, so the leg-muscle recruitment driving these numbers is not obviously
contaminated by it.

## 5. Two cross-joint consistency findings (measured, not narrated)

1. **[SUPERSEDED interpretation — see top correction banner]** The self-computed and
   official-JointReaction methods disagree by a strikingly similar ratio at BOTH joints, computed
   independently: knee JR/self-computed = 391.11/233.20 = **1.677×**; ankle JR/self-computed =
   484.48/286.68 = **1.690×** — a 0.76% difference between two ratios computed at two different
   joints, using the same code, on the same trial. **This observation itself was correct and is
   exactly the signature the sign-bug audit later used**: the "strong supporting evidence for the
   knee cert's own diagnosis" (differentiation scheme) conclusion is superseded — the shared ratio
   corroborates a shared CODE bug (`static_opt_knee.knee_crossing_muscles_and_forces`, called
   identically at both joints), not shared differentiation scheme. Once fixed, both joints'
   self-computed/JointReaction ratios collapse to ~1.00 (0.0014-0.0039%), not ~1.68×.
2. **Achilles force: peak-of-sum vs sum-of-peaks, checked explicitly.** The three triceps-surae
   muscles' own INDIVIDUAL peak tensions (gasmed_r 1140.0 N, gaslat_r 405.7 N, soleus_r 1447.8 N) sum
   to 2993.4 N = 390.3 %BW — almost exactly the Giddings anchor (390 %BW) — but this is the WRONG
   quantity (summing three different muscles' peaks at three possibly-different instants overstates
   the true simultaneous tendon force). The CORRECT quantity — peak of the summed time SERIES — is
   320.26 %BW, 17.95% lower, reflecting genuine ~18% temporal desynchronization between when
   gasmed/gaslat/soleus individually peak. Reported explicitly so the two numbers (2993 N context vs
   the 320.26 %BW headline) are not mistaken for an inconsistency.

## 6. Honest caveats (full list)

1. **Weaker anchor tier (§0)** — the primary, load-bearing caveat. Published cadaveric/quasi-static/
   finite-element model estimates (1977–2000), not a direct in-vivo instrumented-implant measurement.
   No such program has ever existed for the ankle.
2. **Giddings' 540 %BW figure is the talocalcaneal (subtalar) joint, not talocrural** (§1) —
   disclosed, not blurred; the talocrural-only mean (Stauffer+Procter&Paul) is 445 %BW, slightly
   below the blended 476.7 %BW used as the primary comparison point.
3. **[SUPERSEDED — see top banner]** ~~Self-computed (0.601) and JointReaction (1.016) disagree by
   ~1.7×~~ — corrected, they now agree (self-computed 1.016, JointReaction 1.016, 0.0014% apart).
   Original text, same unresolved-but-diagnosed gap as the knee cert (§5.1) — not re-resolved further
   this session (time-boxed, per the
   knee cert's own "not resolved further this session" precedent).
4. **Different populations/eras.** Stauffer 1977 and Procter & Paul 1982 (7 adult male subjects) and
   Giddings 2000 are decades-old studies on their own subject cohorts; subject2 is a modern healthy
   young(ish) OpenCap participant (78.2 kg, 1.96 m). No per-subject or anthropometric matching.
5. **Static Optimization's structural limitation** (inherited, unchanged from the knee cert): SO
   finds the minimum-effort activation solution consistent with the required net joint moments, not
   measured EMG — real co-contraction can exceed the effort-minimizing solution.
6. **Pelvis residual force gate fails** (§4) — inherited from the reused sibling SO run, already
   disclosed there, not a new defect; joint-level reserves stay low so leg-muscle recruitment is not
   obviously contaminated.
7. **Single trial, right leg primary** — left ankle reported only as a coarse pure-reaction symmetry
   check (105.65 %BW vs right's 107.50 %BW), not a full independent replicate; no left-side
   muscle-driven contact-force or Achilles estimate computed.
8. **Numerical differentiation**: Savitzky-Golay (11-sample/110 ms window) is a design choice; shown
   stable across a 70–130 ms sweep (0.02 pct-pts, §3) but is not literally the raw signal.
9. **Achilles force is a 3-muscle model sum** (gasmed+gaslat+soleus) at the point they cross the
   ankle-r cut, not a measurement at the literal mid-substance Achilles tendon; the near-identical
   scalar/vector ratio (0.9955, §5) shows these three pull almost perfectly in parallel at that cut,
   which is the physical justification for treating the sum as "the Achilles tendon force."

## 7. Next step

If a tighter single ankle-contact number is later needed: reconcile the self-computed-vs-
JointReaction ~1.7× gap by re-differentiating with a matched filter (the same open next-step the
knee cert already identified, `docs/MECHANISM_STATIC_OPT.md` §9 — now with cross-joint evidence it is
worth doing once, since the fix would benefit both certs). Separately, if a stronger ankle anchor
becomes available (no in-vivo program exists today, §0), re-run this comparison against it; until
then, a broader literature sweep (e.g. Rasmussen/the reference model-based ankle studies, mentioned but not
independently verified this session) could tighten the published band further.

## Files

- `scripts/msk/validate_ankle_force.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` for proven parse_mot/Savitzky-Golay/BFS code and `static_opt_knee.py` for
  the proven SO-convergence-gate-check + muscle-crossing-detector code, neither re-implemented).
- `data/msk_smoketest/subject2_walking1/ankle_force_validation/ankle_force_validation_results.json` —
  every number in this document, machine-written.
- Sibling (read-only, reused): `data/msk_smoketest/subject2_walking1/static_optimization/{so,jr}/`.
