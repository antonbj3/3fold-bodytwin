# Executable surgical research chain

`tool + path + tissue layers → supplied injury/gap → radius-resolved bleeding/seal → O2/FV history → U/I/M direction/bridge inventories → conditional strength`. RESPONSE R4's D/H/A reservoir runs in parallel with explicitly missing chemistry→rupture law.

Requires Python3.10+, NumPy and SciPy in the existing environment. No external services or networks are used. The package contains frozen Q036/FV, coagulation and platelet definitions; their origin, hash and changes are in SOURCE_MANIFEST and CLOSURE_MANIFEST. A MEASURED label refers only to the source’s stated species/cohort/assay. SYNTHETIC is a chosen closure/transfer. UNKNOWN is always null.

From the workspace root:

```bash
results/LANE_SURGICAL_SYNTHESIS/run_chain.sh --out results/LANE_SURGICAL_SYNTHESIS/my_evidence_run
results/LANE_SURGICAL_SYNTHESIS/run_chain.sh --mode scenario --out results/LANE_SURGICAL_SYNTHESIS/my_scenario_run
results/LANE_SURGICAL_SYNTHESIS/run_chain.sh --resume results/LANE_SURGICAL_SYNTHESIS/r1/nominal_complete_v4/checkpoint21.npz --out results/LANE_SURGICAL_SYNTHESIS/my_restart
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 python results/LANE_SURGICAL_SYNTHESIS/tests_r1.py --out results/LANE_SURGICAL_SYNTHESIS/my_tests
```

Each output destination must be new. Default mode `evidence` leaves biological outcomes null and lists all missing/nonempirical ports. `scenario` is an explicit opt-in; an unknown consumed port stops that mode too. Native rupture stress in Pa remains null. A custom config is supplied with `--config <json>`; create it from `scenario_config()`/`evidence_config()` and change value, unit, status, source, scope and joint ID together. See `r1/nominal_complete_v4/summary.json` for full config and closure manifest.

Python-API:

```python
from surgical_chain.scenario import scenario_config
from surgical_chain.chain import run
cfg = scenario_config()
cfg['parameters']['bridge_fraction']['value'] = 0.5
conditional = run(cfg, mode='scenario')
```

Each Port has `value, unit, status, uncertainty, source, scope, joint_id`. Unit errors are rejected; no implicit mm→m, Γ→strength, strain→ischemia or tracer→mass conversion occurs. Layer Γ, gap and biological/perfusion width are separate inputs. Edge radius/speed are documented but lack a calibrated prediction law to these ports. `incision.work` is cleavage-only; total toolwork is unknown because contact/release/pre-tip-work is not known.

The blood model uses per-bin Poiseuille flow and state-dependent shear and the frozen15-species/2-platelet operators. The plug’s persistent seal-state is paid once; erosion is zero in this declared scenario. Seal reduces assumed atmospheric conductance and never restores severed vessels. Vessel/flow maps and the shear-rate law are biologically uncalibrated.

O2 uses FV permeability and stored oxygen inventory; cardiacQ036 cells are a synthetic skin transfer. The earlier FV C-state feeds back to O2; U/I/M tensors are downstream observers of F/O2 and do not feed back to this C. D/H/A is an isolated fixedC10-reference branch and does not determine native strength. Main strength uses the explicit affine75% reference. These different closures must not be presented as an empirically identified joint collagen law.

Checkpoint21 contains grid, full lateFVstate, early12state/hazard, Q_U/Q_I/Q_M, born/removed/initial-reference, D/H/A and absolute time; JSON adds laws/rates/boundaries/config. Replay of the prefix is unnecessary. The spatial suffix is solved again; no certifiedlocalquery or10×gain is demonstrated. The same conventional operator can restart equally well. Full strength/restart maxdiff7,11e−15pp and chemistry1,39e−17 in the frozen run.

Separate observable curves are in `controls.SeparateCurve`: they interpolate their own calibration domain; new interventions and extrapolations return unknown until an extra closure is specified. Rich control with the same physical state/operator has TIE. Measured species/cohort marginals have no joint posterior; three shared scenario directions are sensitivity, not CI.

Existing main run: `r1/nominal_complete_v4/`. Regressions: `r1/tests_final_v2/verification.json`, `TESTS_R1_final_v2.log`. Numerical refinement: `REFINEMENT_R1.json`. The first seal counterexample and all initial outputs are preserved. Graph receipts are created with the lane's own `../graph`; they are local review receipts awaiting coordinator import. All artifacts PENDING_INDEPENDENT_REVIEW.
