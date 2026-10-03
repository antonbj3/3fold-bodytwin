BT-HX-Q033

**Based on:** `BRIEF.md`, `inputs/QUESTION.md`, `inputs/NIGHT_PREAMBLE.md` and the two verified primary sources `10.1055/a-2150-0587` (Lackmann, Table 1/Figure 7) as well as `10.3389/fmed.2024.1481645` (Barthe, Figure 2/4). The aborted session had no model files left. `PREREG.sha256` verifies the frozen pre-registration.

**Model:** `model.py` describes an elliptical orifice with `g(y)=2b*sqrt(1-(y/a)^2)`, exponential orifice stress, integrated apposition force, `A=pi*a*b`, Darcy flow and evaporative boundary flux. The default case is 46 × 13 mm. The source mechanical reference is 4,2–6,0 N in fresh carcass data.

**Derived Results** (`results.json`): `F_close=3,5718 N` (1,8400 N bias + 1,7318 N surge force), `epsilon_max=0,650`, `A_exposed=4,6967×10^-4 m²`, `Q_cut=2,4030×10^-10 m³/s` against `Q_intact=4,7097×10^-11 m³/s` (`5,1022×`). The evaporation factor is exactly `5,9×`; it is a verified barrier damage analog, not a measurement of surgical wound fluid. Tip stress becomes 1,463 MPa with the assumed tip radius.

**Sensitivity:** `modulus_pa` ±50 % gives 2,7059–4,4377 N; `prestress_pa` ±50 % gives 2,6518–4,4918 N. `barrier_damage_multiplier` ±50 % gives evaporation ratio 2,95–8,85; the total flux only changes ±0,16 % because the assumed Darcy path dominates. Area and opening are unchanged by these parameter changes.

**Verification:** all frozen criteria in `PREREG.md` pass, including placebo `b=0`; five tests in `test_model.py` pass, one of which is an analytical limit case. `py_compile` passes; ruff/mypy/pyright is not present in the environment. No measurement data has been found.

**Next resolution step:** measure before/after geometry, force–deformation, pressure and flow on the same tissue sample, with uncontrolled incision and intact placebo. Then adjust `E`, `sigma0`, `beta`, road and barrier factor; only then can the sensitivity be interpreted specimen-specifically.
