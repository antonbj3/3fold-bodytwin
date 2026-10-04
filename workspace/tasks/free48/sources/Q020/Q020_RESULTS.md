BT-HX-Q020
# Can vessels–interstitium–lymph explain both swelling and substance retention?

**First executable mechanistic model, 4 compartments (plasma / interstitium / cells / lymph),
1 macromolecule (albumin).** Files: `PREREG.md` + `PREREG.sha256`, `model.py`, `test_model.py`,
`results.json`. Runtime ~1,5 min per scenario sweep, 0 threads, < 200 MB.

## 1. Verdict

**Partly — and the answer is no for the quantitative part.** The same parameter set
can *reproduce both quantities*, but **not with the amount and sign the question
suggests within the model’s validity domain**, and the two are **not the same event**:

| Quantity | Model’s response at ΔP_c = +5 mmHg | Reference |
|---|---|---|
| **Swelling** ΔV_ISC | **+3.42 L** (11.07 → 14.41 L, +31 %), P_ISC −2.11 → **+1.06 mmHg** | P_ISC > 0 in edema: Guyton 1965, `10.1161/01.res.16.5.452` (VERIFIED summary) |
| **Substance retention** ΔM_alb | **−2.19 g** (143.3 → 162.7 g stored, but ΔM against start −2.2 g) | no primary source looked up — **UNVERIFIED** |
| Retention ratio R_ret | **−0.64 mg/ml** (negative) | prereg P1 required ≈ +11.3 mg/ml |

Swelling is thus robust and reproduces the primary source’s qualitative requirement (P_ISC becomes
positive). **Retention is not robust**: in this parameter window, interstitial
albumin amount is nearly constant, and the net change is *negative* (dilution) at
ΔP_c ≤ ~5 mmHg. The retention branch is reached only at loads that simultaneously drain plasma
(see §5) — i.e. the model predicts that retention and swelling are **linked to
the same threshold, but not the same load**.

## 2. What was built on

No internal nodes: `~/projects/bodytwin`, `notes/`, `results/MAP2/`, `MECHANISM_ANCHORGRAPH.json`
do **not** exist in the sandbox. Q020’s inputs **K01, K04, K09 are unavailable** — this is
the job’s largest single limitation. Replacement: published literature (web allowed).
Full parameter and source table: `PREREG.md` §3 and `model.py` `PARAMS`.

## 3. The mechanism (why a joint answer is possible)

The core is three first-principles relationships, not an empirical swelling formula:

1. **Starling with glycocalyx coupling** (`capillary()`, `model.py:200`):
   `J_net = K_f[(P_c−P_ISC) − σ_g(E)·(π_p−π_i)]`, where `E` (the glycocalyx’s remaining fraction)
   controls both `σ_g`, `σ_s` **and** `K_f` — and `E` in turn depends on `J_net`.
   This is the feedback making *retention* and *swelling* the same mechanism.
2. **Directing interstitium** (`P_isc`, `model.py:160`): `P_ISC = a·ln(V_ISC/V_k)`.
   Guyton 1965 (`10.1161/01.res.16.5.452`, verified summary) shows that
   compliance is low at negative P_ISC and high just above atmospheric pressure.
3. **Saturated lymph with a knee** (`lymph()`, `model.py:225`):
   `Q_L = K_lymph·P_pump·exp(β·ΔP_L/(1+R_mob))`, and `R_mob` rises 9.96·10⁴ times
   above the knee — basis: Guyton, Scheel & Murphree 1966, Circ Res 19:412-419,
   `10.1161/01.res.19.2.412`: mobility in subcutaneous tissue falls **"usually decreased
   more than 100,000-fold"** when P crosses atmospheric pressure (VERIFIED summary).

The core argument for Q020: **water inflow is large and protein inflow is small**
(`≈(1−σ_s)·c_p·J_net`, `σ_s = 0.90`), while **protein exits via lymph proportionally
to `c_i`**. Because lymph flow is saturated but protein import is linear in `J_net`,
a *final* equilibrium protein amount arises without active storage. This is the shared
mechanism. Its failure to emerge in practice is because `J_clear` and dilution
catch up throughout the loads the model can handle.

## 4. Preregistered criteria — outcomes (all in `results.json:prereg_verdict`)

| Criterion | Bound | Value | Outcome |
|---|---|---|---|
| P1 R_ret = c_ISC and invariant | spread < 15 %, all > 0 | R_ret = −12.32 / −0.64 mg/ml, spread 180 %, **negative** | **FAILS** |
| P2 ΔV spread over ΔP_c | ≥ 3.0 | 2.88 (only 2 of 3 loads physical) | **FAILS** (marginally) |
| P3 P_ISC,ss > 0 at +5 mmHg | > 0 | +1.06 mmHg | **HOLDS** |
| P4 mobility ratio | ≥ 100 | 9.96·10⁴ | **HOLDS** |
| P5 threshold in the window | 4–20 mmHg | ΔP_c,krit = **9 mmHg** (no physical solution ≥ 9) | **HOLDS** |
| P6 baseline flows | 0.2–4 mL/min | J_net 0.444, Q_L 0.450 | **HOLDS** |

Counter-tests: **N1 without lymph** → no physical solution (plasma −4.83 L): lymph is
load-bearing, not cosmetic. **N2 σ_s=0** → ΔM_alb = −59.5 g (protein washed out),
showing that `σ_s` sets the retention branch. **N3 damaged glycocalyx** and
**N0 linear P–V** both give less spread than the main model, so P2/P4 carry weight.

## 5. The model’s validity limit (most important finding)

With a **closed water budget without infusion**, a persistent ΔP_c drains plasma:
the limit lies between **+8 and +9 mmHg** (V_ISC 17.6 L, V_plasma 0.03 L at +8).
ΔP_c = +10 mmHg gives **no** physical solution. This is not a numerical error — it is
the model’s response to Q020’s question presupposing a fluid supply (kidney/intake) that
the frozen scenario did not contain. Preregistration is thus **incomplete at the
only point deciding retention’s sign**. I deliberately have **not** inserted a
renal compensator after the outcome — that would be result steering.

Sensitivity ±50 % (`results.json:sensitivity`), ΔV_ISC at ΔP_c = +5:
- **K_f0 (filtration coefficient): controlling** — both signs give no physical solution
  (the limit is sharp, not soft).
- **a_stiff (interstitial stiffness): ΔV 3.03 / 4.13 L** (−11 % / +21 %).
- **β (lymph pump gain): ΔV 3.41 / 3.42 L** (±0.05 %) — **negligible**;
  lymph coverage is determined by P_pump and the knee, not β.

## 6. Tests and unit checks

`python3 test_model.py` → **29 PASS, 0 FAIL**. Analytical limiting cases:
`V(t) = V_0 + K_f·ΔP·t` against RK4 (rel. error 1.6e-16, T5a); P–V inverse
`V(P) = V_k·exp(P/a)` round-trip (2.7e-15 mmHg, T2); exact water balance
`dV_p+dV_i+dV_c ≡ 0` (3.5e-18 mL/min, T3); equilibrium balance `J_net = Q_L` (rel 0.004)
and `J_p = retention+klarering` (rel 0.025, convergence-limited).
Unit check: 16 constitutive relationships are checked by machine (`unit_check`).
**Two real errors were found by the checks, not by me:** (i) dimensional error in the lymph equation
(`P_pump [mmHg]` was used as flow) and (ii) cell flow was subtracted from both plasma
and interstitium — violation of mass conservation. Both logged in `PREREG.md` §8.

## 7. What the next resolution step requires

1. **The scenario must be completed with fluid supply** (renal excretion + intake,
   τ ≈ 6–12 h) and run as *euvolemic conjunction*. Without it, retention’s
   sign cannot be tested across the load interval that is physiologically interesting.
   This is blocking and must be done **before** any retention number is reported.
2. **Look up** the three values that are **UNVERIFIED** and govern the answer:
   human whole-body lymph flow (L/day), lymph’s P/Q capacity curve and lymph’s
   protein permeability (σ_L). K01/K04/K09 must be made available — without them,
   the question is not tested against BodyTwin’s own basis.
3. **Rerun N2/N3/N0 over the completed scenario** and require them to fail.
   With the current load interval, none gave an unambiguous no, which is
   evidence of weakness for P2.
4. **Sensitivity sweep over the entire parameter grid** instead of 6 points:
   because K_f0’s ±50 % gives a sharp limit, the boundary’s shape must be mapped.
5. **Two macromolecules** (albumin + fibrinogen) are needed to distinguish
   selective retention from nonselective dilution.

## 8. What is NOT claimed

That the model describes any patient’s edema. It is 4 compartments, 1 macromolecule,
zero measurement data. `V_ISC0 = 11 L`, `c_ISC0 = 15 g/L`, `K_f0 = 0.05 mL/min/mmHg`,
`P_pump = 0.5 mmHg`, `β = 0.05 mmHg⁻¹` are **assumptions** (see `PREREG.md` §3, marked
ASSUMPTION/UNVERIFIED), not measurements. No verdict words have been used in the answer above.
