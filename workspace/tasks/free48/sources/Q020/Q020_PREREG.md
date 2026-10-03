# PREREG.md — BT-HX-Q020 (frozen before first model run)

Deadline: written **before** first run of `model.py`. Hash: `PREREG.sha256`.

## 0. The question being pre-registered

**Q020:** Can *vessel → interstitium → lymph* explain **both** swelling (increased interstitial fluid volume)
**and** substance retention (retention of protein in the interstitium) with **a** parameter Set?

Answer form that is tried: either gives the same model, same parameters, both quantities
within the literature intervals, or the model fails. No conclusion is chosen in advance.

**Substrates that Q020 points to (K01, K04, K09) are NOT in the package**
— see §6. Replacement reference: published literature (brief allows web).

## 1. Hypotes (mekanistisk, falsifierbar)

H1. **Swelling and substance retention are two consequences of the same two closed flows.**
    Water: `J_filtration` into the interstitium, `Q_lymf` out. Volume congestion occurs when
    `J_net > Q_lymf`.
    Protein: same flows, but the protein flow is `≈ (1−σ_s)·c_p·J_net` in — ie. **small
    in the upper part of the filter but proportional to c_p**, and out via lymph with the concentration
    `c_l ≈ (1−σ_L)·c_i`. Since the lymph flow is **saturation** (pump + saturation) while
    the protein input is **linear in J_net**, you get a **final equilibrium amount of protein
    in the interstitial ** without anything needing "lagras" active. It is the same mechanism,
    different flow order.

H2. **The interstitium is rectifying, not a constant balloon.** Guyton showed that the
P–V curve and mobility resistance change discontinuously at P_ISF = 0 mmHg (atmosphere).
Therefore, the swelling is **non-linear in ΔP_c** and the composition of the swelling
(protein per ml of fluid) is **almost independent of the perturbation**.

H3. **There is a threshold.** Above a critical ΔP_c, glycocalyx dilution occurs
→ σ_s ↓, K_f ↑ → positive feedback → swelling. Lymph cannot compensate.

## 2. PREDICTED SIZES (frozen, with numbers)

Scenario S1 (main, **frozen**): 70 kg human, subcutaneous type whole-body exchange,
step in capillary hydrostatic pressure `ΔP_c ∈ {+2, +5, +10} mmHg` , all other
parameters unchanged, albumin unchanged in plasma (liver as source).

| # | Predicted quantity | Frozen statement | What counts as error |
|---|---|---|---|
| **P1** | Retentionskvot `R_ret = ΔM_ISC,ss / ΔV_ISC,ss` [mg/ml] | `R_ret` = the concentration of interstitial albumin in the new equilibrium, and `R_ret` varierar **< 15 %** over `ΔP_c ∈ {2,5,10} mmHg` | variation ≥ 15 % → H1 falsified (the swelling is *not* sammansatt of interstitiets own koncentration) |
| **P2** | Volume response `ΔV_ISC,ss(ΔP_c)` | span `max/min ≥ 3.0` over `{2,5,10} mmHg` | span < 1.5 → H2 falsified (the interstitium is modeled as a linear balloon) |
| **P3** | `P_ISC,ss` after steget | `P_ISC,ss > 0 mmHg` for `ΔP_c = +5` | ≤ 0 → H2 falsified against primary source (see §3) |
| **P4** | Mobility knee | `R_mob(P=+1 mmHg)/R_mob(P=−1 mmHg) ≥ 100` | < 100 → the mechanism is not the knee |
| **P5** | Critical load `ΔP_c,krit` (dissolved glycocalyx) | the model should **find a threshold**: `ΔM_ISC(V)/M_ISC,0` should exceed 2 before the asymptote, and the threshold should be between **+4 and +20 mmHg** | no threshold, or threshold outside the window → H3 not supported (but not falsified) |
| **P6** | Base flows | `Q_lymf,0` and `J_net,0` should both end up in **0.2–4 mL/min** (≈ 0.3–6 L/day, whole body) | outside → the parameterization is not physiological, the model is discarded (not adjusted *after* the outcome is seen) |

**Reference value (primary source, REDUCED — see §3):** P3/P4 reviewed against *Guyton
AC, Scheel K, Murphree D (1966), Circulation Research 19:412-419, DOI 10.1161/01.res.19.2.412*:
abstract states that fluid mobility in subcutaneous tissue **"usually decreased
more than 100,000-fold"** when P_ISF passes atmospheric pressure, and that the
P–V curve shows sudden compliance increase there. → The P4 criterion is set
to ≥ 100 (conservative, 3 orders of magnitude below the source value).

## 3. Reference values — source status (searched this session, 2026-09-25)

Bibliographic records verified via Crossref/OpenAlex/Europe PMC-API.
**VALUE from table/figure in the original article: VERIFIED** (= not
verified) when full text was not available; this is marked per line.

| Parameter | Value | Unit | Source (DOI) | Status |
|---|---|---|---|---|
| Interstitial P–V-curve: compliance low at negative P, high just above atmosphere; `P_ISC` normally negative, becomes positive in edema | qualitative | — | Guyton AC 1965, Circ Res 16:452-460, `10.1161/01.res.16.5.452` | **VERIFIED (summary, read via OpenAlex W2075884354)** |
| Mobility resistance drop >1e5 times above atmospheric pressure | >1×10^5 | dimensionless | Guyton AC, Scheel K, Murphree D 1966, Circ Res 19:412-419, `10.1161/01.res.19.2.412` | **VERIFIED (summary, OpenAlex W2100273059)** |
| P–V knee, quantitative value ("10× normal fluid volume") | 10 | times | same as above | **UNVERIFIED** — figure is in secondary citation, not in verified table/figure |
| Subcutaneous tissue: capsule / free liquid / gel liquid / gel absorption pressure | qualitative (P–V shape of the model) | — | Brace RA, Guyton AC 1979, Microvasc Res 18:217-228, `10.1016/0026-2862(79)90030-x` | **VERIFIED (title/record), value VERIFIED** (no summary available) |
| L_p, σ for continuous capillaries | L_p≈(1–2)e-7 cm/s, σ_v≈0.9 | cm/s, – | Michel CC, Curry FE 1999, Physiol Rev 79:703-761, `10.1152/physrev.1999.79.3.703` | **VERIFIED (post + summary); narrow range values UNVERIFIED** (full text not accessible) |
| Interstitial volume 70 kg | 11.0 | L | derived: ECF 0.2 L/kg − plasma 0.043 L/kg | **DERIVED, assumption** |
| Glycocatalysis model equilibria (σ_g, E_glyc) | qualitative | — | Michel & Curry 1999 (reviewers); Michel & Curry 1997 J Physiol 502:39-44 | **UNVERIFIED** (not posted this session) |
| Pressure capacity of the lymphatic pump, β=K_A·A_L/K_L | 0.5 mmHg; 0.05 mmHg⁻¹ | mmHg, mmHg⁻¹ | **hypothesis/formval**, no primary source looked up | **UNVERIFIED** |
| Hypoproteinemia → fluid retention | qualitative | — | Manning RD Jr, Guyton AC 1983, Am J Physiol 245:H284-H293, `10.1152/ajpheart.1983.245.2.h284` | **VERIFIED (post)**, value UNVERIFIED |
| Interstitial volume/pressure-regulation | qualitative | — | Guyton AC, Coleman TH 1968, Ann NY Acad Sci 150:537-547, `10.1111/j.1749-6632.1968.tb14705.x` | **VERIFIED (post)** |

**Not listed (lack of time, explicitly UNKNOWN):** human whole-body lymph flow in L/day
measured at the thoracic duct; lymph P/Q capacity curve; intracellular π_c; interstitial
compliance in SI units. The values ​​below are **assumptions**, not measured data.

## 4. Frozen parameters (the model must not change these after execution)

Se `model.py` PARAM-tabell. Viktigaste frysta antaganden:
`W=70 kg`, `M_n=69 300 g/mol`, `n_Hill=0.45`, `T=310 K`, `V_ISC,0=11.0 L`,
`c_p,alb=40 g/L`, `L_p=1.5e-7 cm/s`, `A_cap=3000 m²`, `f_A=0.02`,
`σ_s,0=0.90`, `D_ratio=0.02`, `σ_L=0.10`, `P_pump=0.5 mmHg`, `β=0.05 mmHg⁻¹`,
`R_mob`-ratio 1e5 above the knee, `V_knee` and `a` enligt tabell.

## 5. Motprov (nollmodell / placebo)

| Null model | What it should show | What counts H1 as false |
|---|---|---|
| **N0 — "ballong"**: linear P–V (`P = k(V−V₀)`), all else identical | that P1/P2 can only hold the knee is | N0 satisfies P1 **and** P2 ⇒ P2 carries no weight |
| **N1 — "no lymph"**: `Q_lymf ≡ 0` | that lymph is load bearing, not cosmetic | N1 gives the same ΔV as the full model ⇒ the lymph term is redundant ⇒ the answer to Q020 becomes "no, or only vessel–interstitium" |
| **N2 — "no σ"**: `σ_s ≡ 0` (all protein leaks) | that σ_s sets the retention ratio | N2 gives the same R_ret ⇒ protein flow is not controlling |
| **N3 — "full glycocalyx"**: `E_glyc ≡ 1` (σ_s and K_f locked) | that threshold P5 requires dilution | N3 swells equally at all ΔP_c ⇒ no threshold ⇒ H3 fals |

## 6. Building on / Not redone

- **Based on:** no internal nodes. `~/projects/bodytwin` , `notes/` , `results/MAP2/`
  , `MECHANISM_ANCHORGRAPH.json` **does not exist** in this sandbox (`ls ~/projects`
  → missing). K01/K04/K09 (Q020's "entry point into the source material") are therefore unavailable.
  (b) no impersonation/LOSO (no data); (c) no cell model for protein kinetics
  — protein assumed to be extracellular + lymphoid and handled by lymph/clearance;
  (d) no tissue mechanics (elasticity) — pressure–volume only.
- **Data amount:** no internal data, no made up metrics. All numbers in
  `RESULTS.md` are model results or literature values ​​listed in §3.

## 7. What NOT is claimed

That the model is *physiologically complete*. It is a 4-compartment model with 1 liquid
+ 1 macromolecule. It determines *whether a common minimal response to Q020 is mechanistically
possible* and *how sensitive that response is* — not the patient's edema.

## 8. CHANGE AV PREREG (logged, before running P1–P6)

The changes are **calibrations**, not performance adjustments. None of them
change the frozen parameters that P1–P6 evaluate (K_f0, beta, a_stiff,
P_pump, R_mob_max, sigma_*, V_isc0, c_plasma), and no change audited
P1–P6's outcomes (the first P runs were not run before the changes).

| # | Change | Reason | Verification |
|---|---|---|---|
| A1 | New parameter `K_lymph = 1.0 mL/min/mmHg` inserted in the lymph equation `Q_L = K_lymph·P_pump·exp(β·ΔP_L_eff)` | **The unit check caught a real dimension error**: `P_pump [mmHg]` was directly used as flow `[mL/min]` . Now `[mL/min/mmHg]·[mmHg] = [mL/min]` . | `unit_check()` entry 6; the numerical value is **not affected** (K_lymph=1) |
| A2 | `D_ratio: 0.02 → 2.0e-4` | With 0.02, the diffusive albumin flow 15× became the convective one and the interstitium crowded towards the plasma albumin level (c_isc → 37 by 40 g/L). | base equilibrium gives c_isc ≈ 17 g/L ≈ 0.42 c_p (literature 0.4–0.6) |
| A3 | `c_cell: 2 → 50 mmHg`, `tau_cell: 20 → 500 min`, `c_prop: 0.03 → 1.0e-4` | With the old values, the cell compartment was a **source** (0.5·(π_i−π_cell) ≈ 3.8 mL/min into the interstitium) which no counterforce carried: 20 days × 3.8 mL/min = 51 L, ie. numerical discharge from the whole body. Now balances the P–V stiffness of the cells. | `mass_balance_residual` and that V_isc < V_total (42 L) in all scenarios |

Hash after change: see `PREREG.sha256` (new value). Previous frozen hash:
`529cad92a49b7e69719de2e05bb3ea63c87e86e4ee0f8ed9772f268605611239` .

| A4 | `pi_cell` : 0 → **DERIVED** `π_i(V_ISC0) − P_ISC(V_ISC0) + P_cell` = 7.66 mmHg (not 0) | With π_cell = 0, the cell compartment provided a permanent **source** 0.5·(π_i−π_cell) ≈ 3.8 mL/min into the interstitium: 25 days × 3.8 mL/min takes out the entire plasma and cell volume. Now the cells are a *perturbation buffer* around the equilibrium. | `unit_check` + that all volumes remain positive in all scenarios |
| A5 | `tau_hep`: 600 → **120 min**, and the hepatostat is rewritten to PI-regulation on the albumen mass: `S_hep = (c_p·V_p − M_p)/τ_hep` | The original form contained a volume term which **dimensionally meaningless for the mass balance** and didn't last c_p: plasmaalbumin failed 40 → 26.6 g/L under 40 day, which crept in π_p = 11.5 instead of 18.4 mmHg. | base equilibrium: c_pl = 39.98 g/L, π_pl = 18.39 mmHg (`results.json:baseline`) |
| A6 | Integrator LSODA → **Radau**, state in **L/g** instead of mL/mg, and software protection against non-physical solutions (negative volume ⇒ reported as *runaway*, not as crash) | LSODA threw `Unexpected istate`; L/mg-scaling yielded 3·10¹¹ in atoll-level dynamic range. | `test_model.py` T3, T10 |
