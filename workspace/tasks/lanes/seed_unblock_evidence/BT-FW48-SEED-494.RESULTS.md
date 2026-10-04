BT-FW48-SEED-494

# RESULTS — M494 coevolution suite with sealed frozen anchor suite

Status: `NOT_EXECUTED_DESIGN_REVIEW` (definition/review) + one executed synthetic matched
mechanism test. Reports remain **PENDING_INDEPENDENT_REVIEW**. No automatic scientific admission.
No EXTERNALLY_MEASURED leaf created.

## 1. Recovered observable, units, input domain, old failures

- **Observable:** `newfamily_error` — held-out coverage residual `1 - max_{m in S} <a_m, z>^2` on two
  unseen task families. Dimensionless, [0,1], lower better. Design stand-in for the brief's
  `newfamilyfel`.
- **Cost:** `benchmarkcost` — counted utility evaluations (suite tasks x planner subset evaluations
  + generator candidate evaluations). Planning units, not measured wall time.
- **Input domain:** latent factor space R^6; 9 operators (6 axes + 3 pair combinations), all cost 1;
  mechanism budget B; two unseen families built from held-out combination directions `c14, c25, c36, c16`.
- **Regime:** next program generation, frozen deterministic replay (cold suite build + warm
  re-selection). Physical tau, sampling/follow-up window, protocol: **UNKNOWN**.
- **Old failure preserved (Q005, inherited, not re-measured):** no confirmed molecular winner;
  `f_u` score `1.93 mL/min/cost`, RMSE reduction `9.397 %` just under the frozen 10 % threshold;
  `K_OCT2` negative in the frozen ranking; `P_MATE` absent at steady state. Q005 supplied the
  measurement discipline but **no task-family generator distribution** — that absence is a missing
  input here.

## 2. What was executed (frozen, PREREG.md before measurement)

Deterministic synthetic operator-coverage test, `m494_coevo_test.py`, seed 49420260930.
Budget B=3, suite N=40, matched strong static control; decision tolerance `delta_dec = 0.02`,
accuracy resolution `epsilon_y = 0.005`. Four arms: `static_strong`, `cosmetic`,
`coevo_sealed` (candidate), `coevo_unsealed` (new-operation-disabled ablation).

## 3. First actual result — candidate rejected at the frozen budget

| arm | selected operators | newfamily_error | benchmarkcost |
|---|---|---|---|
| static_strong (N=40) | e1,e2,e5 | 0.6350 | 5200 |
| cosmetic | e1,e6,c25 | 0.4684 | 5200 |
| **coevo_sealed (candidate)** | e1,e3,c36 | **0.5661** | 22380 |
| coevo_unsealed (ablation) | e1,e2,e5 | 0.6350 | 22380 |
| static_strong, **matched total cost** (n=173) | e1,e2,e5 | **0.5535** | 22490 |

- **M1 (matched suite size N=40):** candidate gain vs static_strong `+0.0689` (> delta_dec), but the
  candidate does **not** separate from the `cosmetic` nuisance arm (cosmetic error is lower by 0.0977)
  and the ablation `coevo_unsealed` is a tie with static at 0.6350.
- **M2 (matched total benchmarkcost):** strong static control with n=173 tasks reaches `0.5535`,
  **beating** the candidate by `0.0126`. The candidate is not admitted.
- **Decision: `REJECT_AT_FROZEN_B=3`.** No separation, a stronger matched-cost baseline, and a
  nuisance arm that (by sampling luck) covers an unseen direction better.

## 4. Amendment — changed construction (post-hoc, not part of the frozen decisive test)

The diagnosed obstruction is **coverage capacity**: with B=3 mechanism slots and 4 unseen
directions, suite construction cannot change the achievable error floor, so the operation has no
room to matter. Raising the budget (same data, same matched-cost rule):

| budget | coevo_sealed err | matched-cost static err | matched-cost gain |
|---|---|---|---|
| 3 (frozen) | 0.5661 | 0.5535 | **-0.0126** (reject) |
| 4 | 0.5365 | 0.5179 | -0.0186 (reject) |
| 5 | 0.3395 | 0.4807 | **+0.1412** (separates) |

So the operation's value is **conditional**, appearing only when mechanism slots exceed the number
of unseen directions. This is a conjecture from a synthetic replay, not physical validation.

## 5. Validity, unit definitions, failed criteria (preserved)

- The result is a design-level mechanism test on an **assumed** unseen-family distribution. It is
  not physical validation, not a universal certificate, and no confidence intervals are built from
  these correlated replays (planned replay levels, not independent samples).
- Failed criteria kept: candidate failed M2; `cosmetic` nuisance arm was not dominated;
  `coevo_unsealed` tied `static_strong` at B=3. `epsilon_y`, `delta_dec` unchanged.
- Unit definitions: `newfamily_error` dimensionless residual; `benchmarkcost` counted evaluations.
- **UNKNOWN** sigma/observation for `newfamilyfel` and `benchmarkcost` (L3); transfer error to the
  recursive planner (L4); physical tau/protocol. UNKNOWN stays unknown.

## 6. What remains incomplete

- A calibrated unseen-task-family generator (the decisive missing input) and the consumer's declared
  risk/tolerance.
- Transfer/reduction bound to the recursive operator/method planner; native graph nodes
  `BT-C4-ERROR-BUDGET`, `BT-IM1-OBSERVABILITY` remain open prerequisites.
- The B=4/B=5 separation is an amendment on synthetic data; it needs independent execution and an
  independent reviewer before any admission.

## 7. Consumer and next test

Consumer: recursive operator/method planner (Gate G2). Next test per `results.json:next_test`:
re-run with a declared unseen-family generator and B > number of unseen directions, and require
`err(static_matched_cost) - err(coevo_sealed) > delta_dec` with valid coverage and no lost primitive
coverage. Same job/source identity `BT-FW48-SEED-494`; do not relabel this synthetic replay as an
executed v6 physical test.
