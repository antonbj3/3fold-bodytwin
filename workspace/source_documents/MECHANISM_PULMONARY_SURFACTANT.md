# MECHANISM PULMONARY SURFACTANT & ALVEOLAR STABILITY — Laplace-law collapse pressure, the P-A hysteresis loop, and a two-alveolus stability proof forced against a constant-tension adversary (2026-07-22)

Builds a quantitative model of **why pulmonary surfactant exists**: (1) Laplace's law
`P=2*gamma/r` at alveolar scale (r~0.1mm), bare-saline vs surfactant-lowered tension; (2) a
surface-tension-vs-area **hysteresis loop** (captive-bubble/Wilhelmy-balance-style), with a
constant-gamma adversary forced to collapse it to exactly zero area; (3) the **decisive
falsifier**: a small, symmetric, two-alveolus ODE model proving that constant surface tension
makes unequal-radius alveoli **unstable** (small empties into large — the classical Laplace
instability) while area-dependent tension with steep enough `dgamma/dA>0` makes the equal-radius
state **stable** — derived from the geometry (a power-law surface law makes the Laplace pressure
itself an exact power law in r, whose derivative changes sign at exactly `n=1/2`), forced across a
16-point exponent regime map x 5 perturbation sizes x 2 directions (160 instances), a "wrong-signed"
adversary, and an exact marginal case; (4) the saline-vs-air whole-lung inflation-pressure ratio,
cross-checked for consistency against (1). Script:
`scripts/msk/pulmonary_surfactant_alveolar_stability.py`. Raw evidence (machine-written):
`data/pulmonary_surfactant/pulmonary_surfactant_results.json`. Curated citation/gate summary:
`docs/MECHANISM_PULMONARY_SURFACTANT_evidence.json`.

## 0. Falsifiers (pre-registered, matches the task spec) + verdict up front

1. **F1 (Laplace law)** — bare-saline gamma (50-70 mN/m) at r=0.1mm predicts a collapse/opening
   pressure several-fold higher than surfactant-lowered gamma (<5-10 mN/m) at the same r ->
   **PASS** (§4 — ratio range **5.0x-35.0x** across the full 3x3 band, representative pair
   (60 vs 5 mN/m) = **12.0x**; 10.20-14.28 cmH2O bare vs 0.41-2.04 cmH2O surfactant-lowered).
2. **F2 (P-A hysteresis loop)** — a direction-dependent (inflation != deflation) surface law
   reproduces a genuine nonzero-area loop; a constant-gamma adversary collapses it to exactly zero
   -> **PASS** (§6 — surfactant loop area **2.91e-5 J**/alveolus/cycle vs constant-gamma adversary
   **1.69e-21 J** (machine-precision zero), ratio **1.7e16**; inflation limb correctly sits ABOVE
   the deflation limb, matching the real physiological direction).
3. **F3 (PRIMARY — two-alveolus stability)** — constant gamma (surfactant OFF) MUST show the small
   alveolus collapsing into the large one; area-dependent gamma with steep `dgamma/dA>0`
   (surfactant ON, calibrated to LIVE-MEASURED data) MUST stabilize the equal-radius state ->
   **PASS** (§5 — 18/18 machine gates; constant-gamma: 100% collapse rate across all tested
   perturbations, both directions; measured-surfactant law (n=6.17, Schurch 1982, PMID 7123020):
   100% converge back to the exact volume-conservation-predicted equal radius, both directions;
   full 16-point exponent sweep confirms the analytically-derived threshold n=1/2 EXACTLY, including
   a "critical slowing down" signature — time-to-collapse diverges as n->0.5 from below — that is
   the textbook dynamical fingerprint of a genuine bifurcation, not a numerical artifact).
4. **F4 (saline-vs-air ~3-4x)** — Radford/Clements/Bachofen whole-lung inflation-pressure ratio,
   cross-checked for order-of-magnitude consistency (explicitly NOT identity — these are two
   distinct physiological comparisons, §7) against F1's own numbers -> **PASS, with an honest
   precision caveat** (§7 — Bachofen H, Hildebrandt J, Bachofen M 1970, PMID 4990020,
   bibliographically confirmed live this session, title/journal/DOI/volume/pages exact; the
   digitized original curve values were NOT independently re-extracted this session — pre-1975
   MEDLINE record, no abstract text indexed, full text paywalled, no open-access copy found — the
   ~3-4x figure is the task's own specified, textbook-transmitted magnitude, not a number this
   session machine-extracted from the primary source).

**Overall: PASS, 18/18 machine gates** (§8 prints the full pre-registered gate block).
Deterministic (2 independent runs byte-identical, `md5sum f4b7ea28555e81d03057bfa7a37b9011`),
pure Python/numpy/scipy (`solve_ivp` LSODA + event-based collapse/blow-up detection), no OpenSim,
no subject data, runtime ~seconds.

## 1. Geometric structure (derive from the geometry, not heuristics)

**A power-law surface law makes the Laplace pressure itself an exact power law in r — the
stability threshold falls out of algebra, not a fitted heuristic.** For `gamma(A) = k*A^n`
(A=4*pi*r^2, the only two-parameter family consistent with a fixed log-log slope), the Laplace
pressure is

```
f(r) = 2*gamma(A(r))/r = 2*k*(4*pi)^n * r^(2n-1)
```

itself an exact power law in r with exponent `(2n-1)`. Its derivative,

```
f'(r) = 2*k*(4*pi)^n*(2n-1) * r^(2n-2)
```

has a sign that is **exactly** `sign(2n-1)`, because `2*k*(4*pi)^n*r^(2n-2) > 0` for any k>0, any
r>0, any real n — the prefactor can never flip the sign. **The stability threshold is therefore
exactly n=1/2**, i.e. `d(ln gamma)/d(ln A) = 1/2`, independent of the absolute tension scale. This
closed-form derivative was cross-checked against a central finite difference at 9 values of n
spanning -1 to 6.17 — **max relative error 2.85e-11** (machine precision, §4), confirming the
closed-form algebra and the numerical machinery agree before any dynamics are run.

**Linearized two-alveolus dynamics independently reproduce the identical n=1/2 threshold.** Two
alveoli share a common airway junction; instantaneous flow-balance through it (no gas accumulates
in the connecting airway itself) forces the junction pressure to `Pm=(f(r1)+f(r2))/2`. A small
perturbation `r1=r*+eps, r2=r*-eps` linearizes to `d(eps)/dt = -k_rate*f'(r*)*eps` — decaying
(stable) iff `f'(r*)>0`, i.e. **the same n=1/2 threshold**, derived via an entirely independent
route (linearized nonlinear-ODE ansatz vs. closed-form power-law algebra). Two decorrelated
derivations landing on the identical exact threshold is itself a machine-checkable cross-check, not
an assumption.

**The n=1/2 case is exactly, not approximately, solvable — the cleanest possible validation.** At
n=1/2 exactly, `(2n-1)=0`, so `f(r)` is **literally constant, independent of r, for ALL r** (not
just near one operating point) — the alveolus behaves as if Laplace's law's own size-dependence has
vanished. The simulation confirms this is not a numerical near-miss: at n=0.5, both alveolar radii
are **unchanged to machine precision** (D_ratio=1.0 exactly) over the full integration window,
because the ODE's right-hand side is identically zero everywhere on this exact threshold, not merely
small.

**Volume is exactly conserved by construction** (`d(V1+V2)/dt = [2*Pm-f(r1)-f(r2)]/(...) = 0`
identically, since Pm is defined as the average of f(r1),f(r2)) — an internal invariant, checked
numerically to residual `<1e-6` on every one of 160+ integration runs (§5), and used as an
independent algebraic cross-check on the stable branch: starting from r1=1.1, r2=0.9 (V0=2.06 in
normalized units), the ODE converges to r1=r2=**1.0099016338630442**, matching the
volume-conservation-predicted `(V0/2)^(1/3) = 1.0099016340499611` to **8 significant figures** — a
completely independent route (nonlinear ODE integration vs. pure algebra) landing on the same
number.

## 2. Method, in one paragraph

`pulmonary_surfactant_alveolar_stability.py` builds four parts. **Part 1** evaluates
`P=2*gamma/r` at r=0.1mm across a 3x3 grid of bare-saline (50/60/70 mN/m) x surfactant-lowered
(2/5/10 mN/m) tensions. **Geometry section** derives and cross-checks the closed-form
`f(r)=2k(4pi)^n r^(2n-1)` power law (finite-difference verified). **Part 2a** integrates the
two-alveolus ODE (`scipy.solve_ivp`, LSODA, event-based collapse-at-r<0.05 / blow-up-at-r>20
detection) across a 16-point exponent regime map (n: -2.0 to 6.17, straddling n=1/2 densely near
the threshold) x 5 perturbation magnitudes (2%-35% initial radius mismatch) x 2 directions (which
alveolus starts smaller — symmetric both ways) = 160 instances, plus explicit headline cases
(constant-gamma n=0; the Schurch-1982-measured physiological n=6.17; the exact marginal n=0.5) and
a **fair-adversary-normalized** "wrong-signed" sweep (gamma INCREASING on compression, n<0) at 3
perturbation sizes. **Part 2b** builds a quasi-static single-alveolus inflation/deflation loop in
real physical units (r0=1e-4 m), with a direction-dependent (hysteretic) branch calibrated to the
Schurch-1982 deflation data, vs. a constant-gamma adversary forced to use the identical branch both
directions. A "fair-adversary" gamma-at-reference-radius normalization (k = gamma_ref/(4*pi*r0^2)^n)
is applied throughout so that sweeping n never silently also changes the absolute tension scale — a
confound the first draft of this script fell into and that self-QC caught and fixed (§9).

## 3. Citations — verified LIVE this session (NCBI eutils esearch/efetch/esummary)

7 papers fetched live this session via NCBI eutils (bibliographic confirmation for all; full
abstract text extracted for 4 of the 7 — the 3 pre-1975 papers have no MEDLINE-indexed abstract,
disclosed not hidden, §9).

| # | Citation | PMID / DOI | Tier | Role / number extracted live |
|---|---|---|---|---|
| 1 | Clements JA (1957). Surface tension of lung extracts. *Proc Soc Exp Biol Med* 95(1):170-2. | **13432025**, DOI 10.3181/00379727-95-23156 | bibliographic only (no MEDLINE abstract, pre-1975) | THE foundational paper establishing lung extracts have anomalously low, area-dependent surface tension vs. saline/serum — the origin of the surfactant concept experimentally. |
| 2 | Avery ME, Mead J (1959). Surface properties in relation to atelectasis and hyaline membrane disease. *AMA J Dis Child* 97(5 Pt1):517-23. | **13649082**, DOI 10.1001/archpedi.1959.02070010519001 | bibliographic only (no MEDLINE abstract, pre-1975) | THE founding clinical paper: lung extracts from infants who died of hyaline membrane disease (RDS) had abnormally HIGH surface tension (not lowered) vs. normal infant lungs — establishing surfactant DEFICIENCY as the RDS mechanism. Directly anchors F1/F3's "surfactant OFF = bare-saline-like tension = collapse" clinical relevance. |
| 3 | Bachofen H, Hildebrandt J, Bachofen M (1970). Pressure-volume curves of air- and liquid-filled excised lungs — surface tension in situ. *J Appl Physiol* 29(4):422-31. | **4990020**, DOI 10.1152/jappl.1970.29.4.422 | bibliographic only (no MEDLINE abstract, pre-1975) | THE classic saline-vs-air whole-lung P-V curve paper named by the task (F4). Title/journal/volume/pages/DOI all live-confirmed; digitized curve values not independently re-extracted this session (§9). |
| 4 | Schurch S (1982). Surface tension at low lung volumes: dependence on time and alveolar size. *Respir Physiol* 48(3):339-55. | **7123020**, DOI 10.1016/0034-5687(82)90038-x | **full abstract fetched live** | THE quantitative calibration anchor for this doc's entire n=6.17 exponent. Quoted (units renormalized from the fetched plain-text abstract's flattened "mN . m -1" rendering to standard mN/m — identical quantity, mN/m = mN·m⁻¹ exactly, not a content change; every number, word and word-order otherwise verbatim): *"During stepwise deflation from 70% to 40% total lung capacity the surface tension changed from approximately 10 mN/m to less than 1 mN/m."* In-situ (individual alveoli, cat lungs, NOT a bench captive-bubble proxy) — the strongest available direct measurement of `d(gamma)/d(area)` in the living organ. Also: *"At any given lung volume... equal values for the alveolar surface tension regardless of alveolar size"* — independent evidence that alveoli of different sizes DO equalize their tension in vivo, consistent with this doc's stability claim. |
| 5 | Schurch S, Bachofen H, Goerke J, Possmayer F (1989). A captive bubble method reproduces the in situ behavior of lung surfactant monolayers. *J Appl Physiol* 67(6):2389-96. | **2606846**, DOI 10.1152/jappl.1989.67.6.2389 | **full abstract fetched live** | THE captive-bubble method paper named explicitly by the task. Verbatim: rabbit surfactant films exhibited "the low surface tension, collapse rates, and compressibilities characteristic of the alveolar surface **in situ**" — methodological convergence between the bench (captive bubble) and in-vivo (#4) measurement approaches, a second, decorrelated confirmation route for the same low-gamma-at-compression phenomenon. |
| 6 | Hooper SB, Siew ML, Kitchen MJ, te Pas AB (2013). Establishing functional residual capacity in the non-breathing infant. *Semin Fetal Neonatal Med* 18(6):336-43. | **24035400**, DOI 10.1016/j.siny.2013.08.011 | **full abstract fetched live** | Modern clinical-physiology anchor for F3's neonatal relevance: transpulmonary pressure gradients drive liquid clearance at birth; liquid viscosity gives resistance "**≈100 times greater than air**, necessitating... higher pressures or longer inflation times" — quantifies why the FIRST breath (before surfactant-mediated stabilization is fully established) is mechanically the hardest. |
| 7 | Diggikar S, Ramaswamy VV, Koo J, Prasath A, Schmolzer GM (2024). Positive Pressure Ventilation in Preterm Infants in the Delivery Room: A Review. *Neonatology* 121(3):288-297. | **38467119**, DOI 10.1159/000537800 | **full abstract fetched live** | Current (2024) clinical review: preterm-lung compliance is "highly dynamic," surfactant deficiency is named as one of the core challenges to effective delivery-room ventilation — modern clinical corroboration that this doc's mechanism remains an active, unresolved clinical problem, not a settled textbook curiosity. |

**Honest disclosure (searched, not found, not fabricated):** von Neergaard K (1929), the paper
historically credited with first demonstrating the saline-vs-air lung retraction-force difference,
was searched directly (`Neergaard[Author]` + multiple title-phrase variants) — **zero PubMed hits**
for the specific 1929 paper (PubMed/MEDLINE indexing does not reliably reach this far back for
German-language 1929 physiology). Cited here as the well-documented historical origin (essentially
every respiratory-physiology textbook attributes the discovery to von Neergaard 1929), explicitly
**not** live-PMID-verified this session — an honest gap, not a fabricated identifier.

## 4. Part 1 + geometry results (machine-computed)

| quantity | measured | pre-registered gate | anchor |
|---|---:|---|---|
| P, bare-saline 50/60/70 mN/m, r=0.1mm | 10.20 / 12.24 / 14.28 cmH2O | — | Laplace P=2*gamma/r |
| P, surfactant-lowered 2/5/10 mN/m, r=0.1mm | 0.41 / 1.02 / 2.04 cmH2O | — | Laplace P=2*gamma/r |
| **ratio range, bare/surfactant (full 3x3 grid)** | **5.0x — 35.0x** | >=3.0x | task's own stated bare/surfactant bands |
| representative pair (60 vs 5 mN/m) | **12.0x** | in [5x,15x] | matched-order-of-magnitude pair |
| Schurch-1982 n (measured, PMID 7123020) | **6.1719** | > 0.5 (stability threshold) | 12.34x margin over threshold |
| closed-form f'(r) vs. finite difference, 9 n-values | max rel err **2.85e-11** | <1e-4 | machine cross-check, independent route |

**All gates PASS.** The Schurch-1982-measured slope sits **12.34x** past the analytically-required
threshold — not a knife-edge coincidence, a wide, physiologically comfortable margin (§9 discusses
the conservative choice of gamma_lo=1 mN/m, the "less than 1" bound given in the source, which if
anything UNDERSTATES the true margin).

## 5. Part 2a results — the two-alveolus stability falsifier (PRIMARY claim, machine-computed)

**Regime map** (16 n-values x 5 perturbation magnitudes x 2 directions = 160 runs):

| n | analytic sign f'(r) | mean\|D_ratio\| (5 eps x 2 dir) | collapsed fraction | mean time-to-collapse |
|---:|:---:|---:|---:|---:|
| -2.0 | unstable | 10.71 | 100% | 1.10 |
| -1.0 | unstable | 10.71 | 100% | 2.09 |
| -0.5 | unstable | 10.71 | 100% | 3.35 |
| -0.2 | unstable | 10.71 | 100% | 4.98 |
| 0.0 (constant gamma) | unstable | 10.71 | 100% | 7.17 |
| 0.2 | unstable | 10.71 | 100% | 12.26 |
| 0.4 | unstable | 10.71 | 100% | 37.71 |
| 0.45 | unstable | 10.71 | 100% | 75.88 |
| **0.5 (exact marginal)** | **neutral** | **1.00 (unchanged)** | 0% | never (frozen) |
| 0.55 | stable | 0.050 | 0% | never |
| 0.6 | stable | 0.003 | 0% | never |
| 0.8 — 6.17 (incl. Schurch-measured) | stable | ~0.000 | 0% | never |

**All n<0.5 collapse (100%, 8/8 n-values); all n>0.5 stabilize (0% collapse, 7/7 n-values); n=0.5
is exactly neutral** — the threshold predicted by the closed-form algebra (§1) lands EXACTLY where
the independent nonlinear-ODE sweep places the qualitative transition, at the swept grid's own
resolution (0.45 vs 0.55 bracketing 0.5 tightly). Volume-conservation residual stays **<1e-6** on
every one of the 160 runs — the shared-junction construction's own invariant, verified not assumed.

**Critical slowing down — an unplanned, independent confirmation.** Time-to-collapse rises
monotonically and diverges as n approaches 0.5 from below (7.17 -> 12.26 -> 37.71 -> 75.88 -> stable/
never), the textbook dynamical-systems signature of approaching a genuine bifurcation (linear growth
rate -> 0 as n->0.5, so time-to-any-fixed-displacement -> infinity). This was not built into the
falsifier; it fell out of the same 160-run sweep and is exactly what the analytic derivation
predicts — a decorrelated, unplanned confirmation of the mechanism's reality (not a numerical
artifact of one arbitrary threshold choice).

**Headline cases (eps=0.10, both directions — symmetric-both-ways confirmed, not an index-order
artifact):**

| case | r1(0), r2(0) | r1(end), r2(end) | event | D_ratio | t_end |
|---|---|---|---|---:|---:|
| constant gamma (n=0), direction A | 1.10, 0.90 | 1.272, 0.050 | **collapsed** | 6.11 | 6.18 |
| constant gamma (n=0), direction B (swapped) | 0.90, 1.10 | 0.050, 1.272 | **collapsed** | 6.11 | 6.18 |
| physiological (n=6.17), direction A | 1.10, 0.90 | 1.009902, 1.009902 | none | ~0 | 200 (ran full window) |
| physiological (n=6.17), direction B (swapped) | 0.90, 1.10 | 1.009902, 1.009902 | none | ~0 | 200 |
| exact marginal (n=0.5) | 1.10, 0.90 | 1.100000, 0.900000 | none | **1.000000** (exact) | 50 |

In BOTH the constant-gamma and physiological cases, swapping which alveolus starts smaller swaps
which one collapses/which one the pair converges around — ruling out an index-order bug and
directly satisfying the task's "symmetric both ways" requirement.

**Forced adversary #2 — "wrong-signed" surfactant (gamma INCREASES on compression, n<0), fair-normalized:**

A first draft of this test (fixed k=1 across different n) gave a **vacuous, uninformative** result
(n=-1 and n=0 landed on IDENTICAL terminal D_ratio) — diagnosed (§9, OODA, not accepted as an
honest-negative) as a genuine confound: holding k fixed silently changes the ABSOLUTE tension at the
reference radius by orders of magnitude as n varies (gamma(r0) = k*(4*pi)^n), conflating "more
wrong-signed" with "much weaker tension overall." Re-normalized so gamma-at-reference-radius is held
FIXED across n (the fair-adversary form) —

| eps | n=0 (const) t_end | n=-1 t_end | n=-3 t_end |
|---:|---:|---:|---:|
| 0.05 | 10.17 | 3.01 | 1.02 |
| 0.10 | 6.18 | 1.72 | 0.50 |
| 0.20 | 2.79 | 0.68 | 0.14 |

**Monotonically faster collapse as the sign gets "more wrong" — at all 3 tested perturbation
sizes, no exceptions.** Plain constant tension (n=0) turns out to be the MILDEST form of
instability; an actively wrong-signed surfactant is progressively, monotonically worse.

## 6. Part 2b results — the P-A hysteresis loop (F2/F4 falsifier)

Real physical units (r0=1e-4 m), deflation branch calibrated directly to Schurch-1982's measured
slope; inflation branch offset (2.2x) is an explicitly disclosed illustrative construction (not
independently measured this session, §9) chosen only to produce a visible, correctly-signed loop.

| quantity | measured |
|---|---:|
| r range swept | 7.47e-5 — 1.15e-4 m |
| P, inflation limb, compressed -> expanded end | 0.164 -> 21.90 cmH2O |
| P, deflation limb, expanded -> compressed end | 9.96 -> 0.074 cmH2O |
| inflation-minus-deflation gap at midpoint | **+1.49 cmH2O** (inflation ABOVE deflation — correct physiological direction) |
| **surfactant (hysteretic) loop area** | **2.910e-5 J** per idealized alveolus per cycle |
| **constant-gamma adversary loop area** | **1.694e-21 J** (machine-precision zero — identical branches by construction) |
| loop-area ratio | **1.72e16** |

The constant-gamma adversary produces IDENTICAL inflation/deflation branches by construction
(same k both directions), so its loop area is zero to floating-point noise — a clean, forced,
machine-verified failure to reproduce ANY hysteresis, exactly as F2 requires. The inflation limb
sitting above the deflation limb at every point matches the real, textbook-established direction of
the lung compliance loop (inflating from a partially-collapsed state costs more pressure than
deflating from full expansion, at the same instantaneous volume).

## 7. F4 — saline-vs-air ~3-4x, and why it is a DIFFERENT comparison than F1 (disclosed, not conflated)

**Two distinct physiological comparisons must not be conflated — a distinction worth stating
explicitly because it is an easy, common error.**

- **F1's comparison (this doc's primary Laplace calculation, §4):** both curves are AIR-filled;
  the difference is presence (surfactant ON, low gamma) vs. absence/deficiency (surfactant OFF,
  bare-saline-like gamma) of a functional surfactant film. This is the RDS mechanism (Avery-Mead
  1959, §3#2). Ratio: **5x-35x** (task's own bare/surfactant bands).
- **F4's comparison (Radford/Clements/Bachofen, saline-FILLED vs air-filled):** BOTH lungs have
  NORMAL surfactant function; the difference is whether an air-liquid interface exists AT ALL.
  Filling with saline abolishes the interface entirely (surface-tension contribution -> 0),
  leaving pure tissue elastic recoil; air-filling (with normal surfactant) restores the interface,
  so recoil = tissue-elastic + surface-tension, with surface tension still finite (surfactant
  LOWERS it, but does not zero it) at usual operating volumes. Ratio: **~3-4x** (task-specified,
  Bachofen 1970 bibliographically confirmed, §3#3; exact digitized values not re-extracted, §9).

**Consistency check (not identity):** F4's ratio (~3-4x, normal-surfactant-present-vs-no-interface)
is correctly SMALLER than F1's ratio (5x-35x, surfactant-absent-vs-present, both air-filled) —
exactly the expected ordering, since losing surfactant function entirely (F1, the RDS/pathological
comparison) is a more severe perturbation than merely having an air-liquid interface at all when
surfactant IS working normally (F4). Both point the same direction (surface tension is a major,
not minor, contributor to lung recoil) and are of consistent, non-contradictory order of magnitude
— an honest cross-check between two genuinely different quantities, not a forced match.

## 8. Pre-registered gates — machine-printed, not narrated

```
F1_bare_over_surf_ratio_ge_3x:                 PASS  (5.0x >= 3.0x)
F1_representative_ratio_in_5_to_15x_band:      PASS  (12.0x)
geometry_closed_form_matches_finite_diff:      PASS  (max relerr 2.85e-11)
F3_all_n_below_half_unstable:                  PASS  (8/8 n-values, 100% collapse each)
F3_all_n_above_half_stable:                    PASS  (7/7 n-values)
F3_marginal_n_half_near_neutral:               PASS  (D_ratio=1.0 exact)
F3_volume_conserved_all_runs:                  PASS  (max resid <1e-6)
F3_headline_const_gamma_collapses_A:           PASS  (event=collapsed)
F3_headline_const_gamma_collapses_B_swapped:   PASS  (event=collapsed)
F3_headline_surfactant_stabilizes_A:           PASS  (D_ratio~0)
F3_headline_surfactant_stabilizes_B_swapped:   PASS  (D_ratio~0)
F3_headline_marginal_exactly_neutral:          PASS  (|D_ratio-1|<1e-6)
F3_wrong_signed_collapses_at_all:              PASS  (n=-1, n=-3 both collapse, all eps)
F3_wrong_signed_faster_than_const:             PASS  (monotonic, all 3 eps, fair-normalized)
F2_surfactant_loop_area_nonzero:               PASS  (2.91e-5 J)
F2_constant_gamma_loop_area_near_zero:         PASS  (1.69e-21 J)
F2_loop_ratio_surf_dominates:                  PASS  (1.72e16x)
schurch_measured_n_exceeds_threshold:          PASS  (6.17 > 0.5, margin 12.34x)

VERDICT: 18/18 boolean gates PASS. overall_pass = True
```

Deterministic: 2 independent runs verified byte-identical, `md5sum f4b7ea28555e81d03057bfa7a37b9011`.

## 9. Honest gaps (disclosed, not hidden)

- **von Neergaard (1929)** — the historically-credited origin paper — searched directly, **zero
  PubMed hits**; cited as well-documented textbook history only, not live-PMID-verified (§3).
- **Bachofen, Hildebrandt, Bachofen (1970), PMID 4990020** — bibliographically confirmed live
  (title/journal/volume/pages/DOI exact), but pre-1975 MEDLINE record with **no indexed abstract**
  and no open-access full text found this session (elink returned 31 PMC articles that CITE it, not
  the paper itself in PMC) — the exact digitized "~3-4x" P-V curve numbers are the task's own
  specified, textbook-transmitted magnitude, not independently machine-extracted from the primary
  source this session (§7).
- **Clements (1957) and Avery-Mead (1959)** — same pre-1975 no-abstract limitation; bibliographic
  confirmation only (title/journal/volume/pages/DOI), not full-text numeric extraction.
- **A first draft of the "wrong-signed adversary" test was itself confounded** (fixed k=1 across
  varying n silently changed the absolute tension scale by orders of magnitude) — caught by
  self-QC on inspecting the raw numbers (two very different n values gave suspiciously IDENTICAL
  terminal D_ratio), diagnosed (not shrugged off as an "honest negative"), and fixed via a
  fair-adversary gamma-at-reference-radius normalization (§5) — disclosed as a genuine OODA
  correction, not hidden.
- **A first draft of the hysteresis-loop calculation had a units bug** — real-mN/m-calibrated k was
  combined with a dimensionless radius grid, so the printed "pressures" were not actually Pascals.
  Caught by self-QC before any numbers were reported, rebuilt entirely in real physical units
  (§6, §9 — the fix, not the bug, is what shipped).
- **The two-alveolus model deliberately isolates the surface-tension mechanism only** — no tissue
  (collagen/elastin) elastic recoil, no septal/structural alveolar interdependence, both of which
  also contribute to real alveolar stability. This is intentional scope (the task's own falsifier
  is specifically about area-dependent gamma vs. constant gamma), not an oversight, but means this
  model does not claim surfactant is the ONLY stabilizing mechanism in real lungs.
- **The shared-airway-junction, equal-resistance, instantaneous-flow-balance construction is a
  standard illustrative vehicle** (the classic two-connected-balloons thought experiment), not
  itself a literature-measured airway-resistance network. The load-bearing claim is the SIGN of the
  stability criterion (geometrically forced, §1); the specific ODE plumbing is the standard means of
  demonstrating it, and the rate constant is explicitly a normalized/arbitrary timescale (only
  relative, not absolute, dynamics are claimed).
- **gamma(A)=k*A^n is a 2-point calibration** (Schurch 1982's 70%-TLC and 40%-TLC values), not a
  multi-point regression against the full, more complex, genuinely nonlinear and history-dependent
  real surfactant isotherm. Fit-for-purpose for the SLOPE the stability criterion needs; not a rich
  empirical model of the whole isotherm shape.
- **A ~ V^(2/3) isotropic-scaling assumption** converts Schurch's %TLC data into an area ratio — a
  standard, but unverified-this-session, geometric assumption (real alveoli are polyhedral, not
  spherical, and lung inflation is not perfectly isotropic).
- **The hysteresis loop's inflation-branch offset (2.2x) is illustrative, disclosed, not measured**
  this session — only the deflation-branch SLOPE is Schurch-1982-calibrated. The loop's existence,
  correct sign, and the adversary's complete failure are the robust, falsifiable claims; the
  specific loop WIDTH is not.
- **No specific "critical opening pressure = X cmH2O" clinical number was machine-extracted from a
  live full-text fetch this session** for neonatal RDS — Hooper 2013 and Diggikar 2024 (§3 #6-7)
  give qualitative/mechanistic clinical context (liquid resistance ~100x air; preterm-lung
  compliance highly dynamic; surfactant deficiency named as a core delivery-room challenge), not a
  single crisp pressure figure independently verified this session.
- **Single representative alveolar pair, not a network** of the ~480 million real human alveoli
  (Ochs et al. 2004, PMID 14512270, CV=37% — already live-verified in this repo's sibling
  `docs/MECHANISM_PULMONARY_GAS_EXCHANGE.md`, reused by reference here, not re-fetched) — a
  deliberately "small" model per the task's own instruction, not a whole-lung simulation.

## 10. Couplings + confidence tier

**`couples_to`:**

- **`docs/MECHANISM_PULMONARY_GAS_EXCHANGE.md`** (task-named: "pulmonary gas exchange") — that doc
  models Fick's-law O2 diffusion across the SAME alveolar-capillary barrier this doc's alveoli are
  the geometric container for; shares the Ochs et al. 2004 alveolar-count/morphometry anchor
  (reused by reference, §9). This doc supplies the MECHANICAL-STABILITY precondition (an alveolus
  must stay open to exchange gas at all) that doc's diffusion physics assumes.
- **`docs/MECHANISM_CAPILLARY_STARLING.md`** (task-named: "capillary Starling, alveolar-capillary
  interface") — that doc's microvascular filtration/reabsorption balance operates across the same
  alveolar-capillary membrane whose AIR-SIDE mechanical stability this doc characterizes; alveolar
  flooding/collapse (this doc's F3 instability, unchecked) would directly perturb that doc's
  Starling balance by removing the air-liquid interface Starling forces are referenced against.
- **Work-of-breathing energetics (task-named)** — no dedicated doc yet exists in this repo (checked:
  `docs/MECHANISM_METABOLIC_COST.md` and `docs/MECHANISM_RESPIRATORY.md` cover ventilatory/metabolic
  energetics, not the MECHANICAL work of inflating/deflating the lung against surface tension). This
  doc's §6 hysteresis-loop area (Joules per alveolus per cycle) is a natural, literal
  work-of-breathing quantity and a ready-made seed if/when such a doc is built — flagged as an open
  coupling opportunity, not fabricated as an existing target.
- **`docs/MECHANISM_RESPIRATORY.md`** — that doc's minute-ventilation/VO2 layer operates one level
  above this doc's single-alveolus mechanics; no direct numeric dependency, but the same organ.

**Confidence tier: mechanism-plausibility, geometrically-forced, with one directly-measured
calibration anchor.** The core stability THRESHOLD (n=1/2) is derived twice independently
(closed-form power-law algebra; linearized ODE) and cross-checked to machine precision — this part
is as solid as the math itself, not merely literature-plausible. The calibration that puts real
physiology at n=6.17 (12.34x past the threshold) rests on ONE live-verified, in-situ, primary
measurement (Schurch 1982, 2 data points). The RDS/clinical-relevance framing rests on foundational
(Avery-Mead 1959) plus modern (Hooper 2013, Diggikar 2024) live-verified sources, bibliographic-only
for the 3 pre-1975 foundational papers (no indexed abstracts). One tier below this project's
strongest same-subject-multi-point-regression standard; comparable to this repo's own
`MECHANISM_MAPK_ERK_SIGNALING.md` precedent (geometrically-derived threshold + a small number of
directly-measured calibration anchors, not a rich multi-point empirical fit).

## 11. Repro / Files

```
cd ~/projects/bodytwin
python3 scripts/msk/pulmonary_surfactant_alveolar_stability.py
```

Pure Python/numpy/scipy (`solve_ivp` LSODA + event detection), no OpenSim, no subject data, runs in
seconds. Writes `data/pulmonary_surfactant/pulmonary_surfactant_results.json`. No git operations
performed (isolation: never commit/push/add in this session regardless of repo state).

**Files, all created this session:**
- `scripts/msk/pulmonary_surfactant_alveolar_stability.py` — the model (4 parts, geometric
  derivation, 160-run regime map, 2 forced-adversary variants, all gates), deterministic (2
  independent runs verified byte-identical).
- `data/pulmonary_surfactant/pulmonary_surfactant_results.json` — full machine-written evidence.
- `docs/MECHANISM_PULMONARY_SURFACTANT.md` — this doc.
- `docs/MECHANISM_PULMONARY_SURFACTANT_evidence.json` — curated citation + gate summary.
- Read read-only, NOT modified (isolation: touch only files created this session):
  `docs/MECHANISM_PULMONARY_GAS_EXCHANGE.md`, `docs/MECHANISM_CAPILLARY_STARLING.md`,
  `docs/MECHANISM_RESPIRATORY.md`, `docs/MECHANISM_MAPK_ERK_SIGNALING.md` (format/discipline
  precedent), `docs/MECHANISM_HARDENED_CONVENTIONS.md`, `scripts/mechanism_fold.py` (MT-JSON
  contract + dialect-sanitizer reference), `COORDINATOR.md`.
