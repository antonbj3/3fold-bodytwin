# Contract for combined and sequential interventions (v0)

Agent O4, 2026-09-22. This is a research contract for model experiments and no clinical recommendations. No treatment or patient is intended.

## 1. Minimal representation

An intervention is a record with the following fields. All fields are required. An unknown value is written as `UNKNOWN`, never as a default value.

| Field | Content | Example (Lemaire cell) |
|---|---|---|
| `id` | stable identity | `INT-PTH-CONT` |
| `target` | model quantity and cell, with exported function | `bone_rankl_opg_lemaire2004.rhs`, input `I_P` |
| `change_kind` | `parameter` / `input_signal` / `boundary_condition` / `mechanism_swap` / `geometry` | `input_signal` |
| `change` | value and unit, or a function of time | 1000 pM/dag |
| `schedule` | start, stop, shape (step/pulse/ramp), repetition | step day 20–80 |
| `provenance` | `measured` / `calibrated` / `literal_cited` / `assumption` (same vocabulary as `thickness_v1`) | `literal_cited` (Lemaire 2004 Fig. 2 protocol) |
| `link_evidence` | evidence that the change acts through this particular target. Distinguish **interventional** evidence from **associative** | model postulate; no external calibration in BodyTwin |
| `validity` | model validity range (dose, time, population) | the cell's own parameters; no bone mass |
| `observables` | which outputs are assessed, with unit | osteoclast C (pM), osteoblast B (pM), AUC excess (pM·dag) |

A **combination** is a list of interventions plus the fields:

- `order`: simultaneous or a specified sequence, with gaps in days.
- `interaction_hypothesis`: expected sign of nonadditivity, or `UNKNOWN`.
- `feedbacks`: which feedbacks in the model can give interaction. In the Lemaire cell, both PTH and OPG act through the RANKL denominator (Φ_L), so interaction is structurally expected.

## 2. Comparisons that are always run

Six arms on the same solver, with the same tolerances and the same initial state:

1. **base**: no interventions. Everything is compared against the base arm and not against the initial state, because the model can drift.
2. **A** alone.
3. **B** alone.
4. **A+B** simultaneously.
5. **A→B**.
6. **B→A**.

Rapporteras per observabel:

- Effect: Δ = arm − base.
- Nonadditivity: Δ(A+B) − Δ(A) − Δ(B).
- Order effect: Δ(A→B) − Δ(B→A).
- Noise floor: base drift and solver tolerance. An effect below 10× the floor is reported as zero.
- Sensitivity to dose in each arm. The **uncalibrated** arm strength must be swept before sign or magnitude is interpreted.

Associative relationships from observed data must not be entered as `link_evidence` of interventional type. A graph distance or a data connection is no evidence of a treatment effect.

## 3. Minimum runnable example (RUN 2026-09-22)

Cell: `src/bodytwin/cells/musculoskeletal/bone_rankl_opg_lemaire2004.py` (VERIFIED-FRESH). Its `rhs(t, y, I_P, I_L, I_O)` takes three intervention inputs. Note: the module runs its self-tests on import, so stdout is captured.

```python
import sys, io, contextlib, numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, "<current_bodytwin>/src/bodytwin/cells/musculoskeletal")
with contextlib.redirect_stdout(io.StringIO()):
    import bone_rankl_opg_lemaire2004 as L
ARMS = {"bas": [], "A_PTH": [("I_P",1000,20,80)], "B_OPG": [("I_O",1000,20,80)],
        "A+B": [("I_P",1000,20,80),("I_O",1000,20,80)],
        "A->B": [("I_P",1000,20,80),("I_O",1000,80,140)],
        "B->A": [("I_O",1000,20,80),("I_P",1000,80,140)]}
te = np.linspace(0, 300, 3001)
for arm, sched in ARMS.items():
    def f(t, y, s=sched):
        kw = {"I_P":0.0,"I_L":0.0,"I_O":0.0}
        for k,a,t0,t1 in s: kw[k] += a if t0 <= t <= t1 else 0.0
        return L.rhs(t, y, **kw)
    sol = solve_ivp(f, (0,300), L.y0, method="LSODA", rtol=1e-10, atol=1e-14, t_eval=te, max_step=0.5)
    print(arm, np.trapezoid(sol.y[2]-L.y0[2], te))
```

Result, osteoclast C, AUC excess in pM·dag:

| Arm | C-AUC | C-topp/C0 | B-AUC |
|---|---|---|---|
| bas | 1,5e-5 (drift) | 1,00006 | 4,1e-5 |
| A (PTH) | 0,1659 | 4,04 (cellens eget +304 % reproduceras) | 0,0588 |
| B (OPG, 1000 pM/dag, **ANTAGEN dos**) | −0,0028 | min 0,947 | −0,0014 |
| A+B | 0,1503 | 3,82 | 0,0548 |
| A→B | 0,1637 | 4,04 | 0,0578 |
| B→A | 0,1629 | 4,04 | 0,0574 |

- Nonadditivity (C-AUC): −0,0127. This is about 8 % of the A effect and about 850× the base drift.
- Order effect: +0,0008. This is about 50× drift, but small compared with the effects.

**Interpretation within the model.** Simultaneous OPG dampens the PTH peak more than OPG's own effect. This agrees with both acting in the same RANKL denominator. The result says no more than that. It is a structural model result. The OPG dose is not calibrated, and the model lacks a bone mass state, pharmacokinetics and mechanical load.

## 4. Known gaps, which cells can fill them, and countertests

| Lucka | Befintlig byggsten | Motprov |
|---|---|---|
| Bone mass/BMD state | `bone_remodeling_transient_arithmetic` (`transient_from_acf_step`, `persistent_rate_from_imbalance`) | null case: no interventions gives unchanged mass |
| Mechanical load as intervention | `bone_wolff_law_mechanostat.mechanostat_rate`, `bone_remodeling` (strain) | without coupling to RANKL/OPG, the load arm should give exactly zero C effect |
| Calibration of dose and effect | missing. Requires published dose–response data for PTH/OPG analogs | dose sweep; the sign of nonadditivity should be stable |
| Synthetic flag through the chain | `result_envelope_v1` (today only connected to `metabolic_cost`) | a synthetic input must give `synthetic_demo` in output. This is not the case today for `bone_remodeling`/`bmu_turnover_kinetics` |
| Rehabilitation or motion as intervention | `motor_unit_recruitment`, `force_velocity_power_rfd`, `muscle_memory_substrate_function_dissociation` | requires time course for training dose, which is missing |
