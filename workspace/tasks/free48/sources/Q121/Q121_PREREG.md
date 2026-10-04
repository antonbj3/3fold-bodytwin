# BT-HX-Q121 — Preregistration

## Scope and frozen question

This is a first runnable lumped, mechanistic model of colonic gas and contents. It asks whether production, absorption, dissolution, microbial consumption, retention and rectal evacuation can be separated from the same observed gas amount, pressure and volume. It is not a clinical prediction and it does not fit individual measurements.

Primary scenario: 15 g lactulose delivered to the proximal colon over 4 h in a healthy-adult reference case; follow-up for 24 h. The model predicts cumulative H2-equivalent excretion in rectal gas plus breath, in mL at standard temperature and pressure, and separately predicts gas volume, pressure, dissolved gas and the flow partition.

## Reference and criterion

The primary reference is:

- Stefan U. Christl, Peter R. Murgatroyd, Glenn R. Gibson and the collaborator H. Cummings (1992), “Production, metabolism, and excretion of hydrogen in the large intestine”, *Gastroenterology* 102(4 Pt 1):1269–1277, DOI `10.1016/0016-5085(92)90765-q`.
- Web lookup used: the Europe PMC abstract page, https://europepmc.org/article/med/1551534, and the publisher DOI page, https://doi.org/10.1016/0016-5085(92)90765-q. The abstract reports total H2 excretion of `227.0 +/- 60.7 mL/24 h` after `15 g` lactulose in a whole-body calorimeter study of 10 healthy subjects. The number is a study mean, not a target for calibration.
- The reference is `227.0 mL H2/24 h`; the secondary mechanistic reference is Hammer (1993), *Gut* 34:818–822, DOI `10.1136/gut.34.6.818`, whose abstract reports a `76 mL/6 h` hydrogen accumulation threshold and absorption efficiency changing from about 90% at low gas load to 20% at high load.

Frozen primary acceptance criterion: the deterministic nominal prediction of cumulative rectal-plus-breath H2-equivalent excretion at 24 h must lie between `113.5` and `454.0 mL` (within a factor of two of 227.0 mL). This criterion is evaluated once after the first run and is not changed afterward. A result outside the interval is reported as failed, not repaired by post-hoc tuning. The secondary Hammer values are reported as qualitative checks and are not silently substituted for the primary endpoint.

## Hypotheses and signatures

1. Production: a +/-50% change in the fermentation rate must change cumulative gas production and total excretion, while holding other parameters fixed.
2. Absorption: a +/-50% change in the blood-wall transfer coefficient must change the breath/rectal partition and leave total gas production unchanged.
3. Retention: a +/-50% change in motility must change the residence-time distribution, gas volume and pressure, and therefore rectal flow; it must not be equivalent to changing production.
4. Dissolution: gas-to-liquid transfer must be separately visible as dissolved amount and must respond to liquid volume and gas partial pressure.

A candidate that cannot distinguish these interventions fails the mechanistic test even if its single aggregate output happens to match the reference.

## Null models and falsification

- Production-only placebo: set absorption, dissolution exchange, microbial consumption, proximal return and rectal evacuation to zero. It is a mass-balance control, not the proposed model.
- Absorption-only placebo: set fermentation to zero and inject a fixed upstream gas flux. It tests whether an output can be explained without production.
- The full model is rejected for this first run if the gas-state balance residual exceeds `1e-6` relative to cumulative production, any state becomes materially negative, the analytic limit test fails, or the primary criterion fails without an explicit uncertainty statement.

## Builds on

- Node `BT-HX-Q121`; inputs: `inputs/QUESTION.md`, `inputs/NIGHT_PREAMBLE.md`, `BRIEF.md`.
- The interrupted-session log (`agent.log`) contains only inspection steps and no reusable model output.
- No earlier `PREREG.md`, `model.py`, `test_model.py`, `RESULTS.md` or `results.json` was present at the start of this session.
- Public source checked for the reference and mechanistic priors: Christl et al. (1992), DOI above; Hammer (1993), DOI `10.1136/gut.34.6.818`; Gibson et al. (1990), “Alternative pathways for hydrogen disposal during fermentation in the human colon”, DOI `10.1136/gut.31.6.679`.

## Not redone

No patient-specific 3D geometry, MRI/CT, acoustic data, tissue perfusion, microbial community measurements, methane-speciation fit, or internal BodyTwin data are used. The model is intentionally a one-dimensional three-compartment reduction: a future geometry-resolved version can replace the compartment laws without changing the conserved quantities.
