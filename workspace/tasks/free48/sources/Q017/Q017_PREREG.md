# BT-HX-Q017 — preregistration

## Status and directory check

This is a first executable mechanistic model of a simultaneous heart–kidney transient. The directory check before the work found only `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, `inputs/QUESTION.md` and empty `ALLOW_WEB`; no `PREREG.md`, `PREREG.sha256`, `model.py`, `test_model.py`, `RESULTS.md` or `results.json` remained from an interrupted session. The interrupted lane log shows that predecessors had only read the directory and had no executable results.

## Locked hypothesis

A measurable transient requires at least five interconnected, but not fully_segmented, feedbacks:

1. **Fluid → preload → cardiac output:** extracellular volume changes venous filling and therefore cardiac output.
2. **Cardiac output → arterial pressure → renal flow:** a compatible arterial pressure gives renal blood flow through afferent and efferent resistance.
3. **Renal flow → GFR and urine:** the glomerular pressure difference gives filtration; the reabsorption fraction gives urine and feeds back to fluid balance.
4. **Pressure and GFR → TGF/myogenic afferent tone:** delayed feedback limits renal flow variation.
5. **Pressure/volume → sympathetic activity, Ang II and natriuretic activity:** fast neural feedback and slower hormonal feedback change vascular resistance, renal resistance and reabsorption.

The hypothesis is that a pressure–flow–volume–organ response transient cannot be determined by `MAP = CO × TPR` alone or by a GFR curve alone. In the first run it quantifies only the reduced core; individual RAAS cascades, ADH, AQ2 and nephron segments are not reused as new segment models.

## Builds on

The following are read-only from `source_repository/` and built on without copying another agent's code:

- `data/MECHANISM_ANCHOR_GRAPH.json`: `INT-CARDIORENAL-PRESSURE-AXIS`, `ORG-KIDNEY-NEPHRON`, `MODEL-RAAS-BLOOD-PRESSURE`, `MODEL-ARTERIAL-WINDKESSEL`, `SYS-BAROREFLEX-AUTONOMIC`, `ORG-RENAL-FLUID-ELECTROLYTE`, `ORG-CARDIAC-PUMP-MECHANICS`.
- `scripts/msk/renal_filtration.py` and `data/renal_filtration/renal_filtration_results.json`: Davies–Shock GFR/ERBF/ERPF anchors; the nephron's existing glomerular resistance and GFR mathematics.
- `scripts/msk/nephron_transport_cell.py`: existing 80–180-mmHg sweep and separate nephron segments, which are not connected here.
- `scripts/msk/raas.py`, `docs/MECHANISM_RAAS.md` and the RAAS result: pressure–natriuresis, Ang II/aldosterone and the slow feedback direction.
- `scripts/msk/baroreflex.py`: the baroreflex's fast timescale and direction; no existing implementation is copied.
- `scripts/msk/venous_return.py` and `scripts/msk/fluid_compartments.py`: compliance/volume limits and the existing open boundary between filling and pressure.
- `scripts/msk/cardiac_output_geometric.py`: geometric pump model supporting preload/flow, not rerun here.
- `FUNCTIONALITY_SEARCHSPACE.md` in `results/HYPERREALISM_ROADMAP_20260924`: K01 = conservative transport, K09 = dynamic joint execution, K14 = individual/population; Q017 defines that pressure, flow, fluid amount and organ response must be read simultaneously.

## Not redone

- No OpenSim/the reference model run, no internal human data profile and no invented transient.
- No full RAAS cascade, no ACE/Ang II kinetics, no adrenergic feedback with separate organs, and no full natriuretic BNP/ANP cascade are reused as if they were verified here.
- No segmented PCT/TAL/DCT/CD model, no concentration gradient and no ADH–AQP2 dynamics are rebuilt. They are explicitly the next resolution step.
- No pulse-wave model or spatial renal architecture replaces the Windkessel; this is a mean pressure/flow model.
- The volume load used here is a controlled synthetic protocol, not a measurement.

## Predicted quantity and primary reference

**Primary predicted quantity:** the relative change in cardiac output from before to the end of a controlled 3 L/0,9 % NaCl infusion over 3 h:

\[
 r_{CO} = \frac{CO_{end}-CO_0}{CO_0}.
\]

**Reference value:** `r_CO_ref = 0.147`, that is, 14,7 percent. The primary source is Kumar et al., *Critical Care* 8:R128–R136 (2004), DOI `10.1186/cc2844`, PMID 15153240, Table 1, group 1 with 24 healthy men. The group received 3 L normal saline over 3 h; cardiac index increased from 2,90 to 3,32 L/min/m², i.e. +14,7 ± 2,4 percent. Table 1 also gives MAP 80,7 → 79,6 mmHg and TPR 599 → 491 dyn·s·cm⁻⁵. The open full text and table are verified via PMC/Europe PMC. The source is a primary human study; the comparison uses the relative change, not the study participants' absolute cardiac index.

This is a compatibility and direction test, not a claim that the model reproduces exactly the same study participants or their autonomic contractility.

### Frozen primary criterion

\[
 \left|r_{CO,pred}-0.147\right|\leq 0.050
\]

with `r_CO,pred > 0`. The interval is 9,7–19,7 percent. The limit is frozen before the first run and must not be changed after the outcome is known.

### Secondary, fixed controls

1. GFR and renal blood flow at rest must match the Davies–Shock table values used for initialization within 1 percent. This is an initialization and unit check, not an independent transient validation.
2. The quasi-static renal pressure sweep 80, 85, ..., 180 mmHg must keep GFR within ±20 percent of GFR₀ with full TGF/myogenic feedback. The limit is frozen here and is not a new measurement; it is a mechanistic limiting-flow check of the model.
3. The cumulative fluid balance must have relative residual <10⁻⁹ at all saved times.
4.inga_states may become NaN or negative where the amount is defined as non-negative.
5. The full model must give at least three channels that are not identical: `P_art`, `renal_flow`, `GFR/urine`, `ECFV` and hormone activity must have separate time series.

## The model's core mechanism

All time series use seconds. `V` is extracellular fluid amount in mL, `P_a` and `P_v` are arterial and venous pressure in mmHg, `R_a` is afferent renal resistance in mmHg·s/mL, and `Q`, `Q_r` and `G` are flows in mL/s. `S`, `A` and `N` are relative, dimensionless sympathetic, Ang II and natriuretic activities around rest.

### Heart and pressure

\[
 P_v=P_{v0}+f_s\frac{V-V_0}{C_{sys}}
\]

\[
 R_{sys}=R_{sys0}(1+k_{A}A+k_SS-k_NN)
\]

\[
 Q=clip(Q_0+k_{pre}(P_v-P_{v0})+k_{con}S+q_{ext},0,Q_{max})
\]

\[
 C_a\frac{dP_a}{dt}=Q-\frac{P_a}{R_{sys}}.
\]

This is a mean-pressure Windkessel/Frank–Starling reduction: compliance stores pressure, preload changes venous filling and cardiac output, and vascular resistance changes pressure. `f_s=0.30` is the fraction of a volume change treated as stressed, sensitive circulatory volume in this reduction.

### Renal flow and filtration

\[
 P_r=P_a-P_v
\]

\[
 Q_r=clip\left(\frac{P_r}{R_a+m_RR_e},0,Q_{r,max}\right)
\]

\[
 P_{gc}=P_a-Q_rR_a
\]

\[
 G=K_f\,[P_{gc}-P_{bow}-\pi_{gc}]_+
\]

Urine and reabsorption:

\[
 f_{reabs}=clip(f_{reabs,0}+k_{reabs,A}A+k_{reabs,S}S-k_{reabs,N}N,0.97,0.999)
\]

\[
 U=G(1-f_{reabs}).
\]

`Q_r` is a blood flow and `G` is a filtrate flow; they must not be equated. `P_bow` and `pi_gc` are held as local glomerular pressures/oncotic pressures in the reduction. Full nephron transport is not claimed here.

### Feedbacks

\[
 R_{a,target}=R_{a0}\exp[k_{myo}(P_r-P_{r0})+k_{TGF}(G/G_0-1)+k_{aff,S}S-k_{aff,N}N]
\]

\[
 \tau_a\frac{dR_a}{dt}=R_{a,target}-R_a.
\]

This is local myogenic/TGF-like negative feedback. The negative GFR part makes increased distal delivery raise afferent resistance; the pressure response makes resistance increase at high pressure and decrease at low pressure.

\[
 \tau_S\frac{dS}{dt}=S^*-S,\quad
 \tau_A\frac{dA}{dt}=A^*-A,\quad
 \tau_N\frac{dN}{dt}=N^*-N
\]

med

\[
 S^*=clip[-k_{S,P}(P_a-P_{a0})+k_{S,V}(V_0-V)/V_{scale},-1.5,1.5]
\]

\[
 A^*=clip[k_{A,P}(P_{a0}-P_a)+k_{A,V}(V_0-V)/V_{scale}+k_{A,S}S,-0.5,2.0]
\]

\[
 N^*=clip[k_{N,V}(V-V_0)/V_{scale}+k_{N,P}(P_a-P_{a0}),-0.5,2.0].
\]

`S` and `N` are fast feedbacks, `A` a slow feedback. They are combined feedback modes, not separate measured hormone concentrations.

### Fluid balance and cumulative counters

\[
 \frac{dV}{dt}=I_{in}(t)-U-L_{other}
\]

\[
 \frac{dI_{cum}}{dt}=I_{in}(t),\quad
 \frac{dU_{cum}}{dt}=U,\quad
 \frac{dL_{cum}}{dt}=L_{other}.
\]

The supplied volume and urine outflow are separate account volumes; fluid cannot disappear through hidden feedback.

## Fryst protokol

All arms start from the same initial state and have a 300 s baseline before the perturbation. The main arm is `saline_load`: 3 000 mL is infused evenly over 10 800 s from t=300 s to t=11 100 s. The simulation ends at 14 400 s. This matches the Kumar study's dose and timescale but is not the same study participants.

Kontrollarmar:

- `no_renal_autoregulation`: TGF/myogenic coefficients are set to zero; everything else unchanged.
- `no_neural_hormonal`: S, A and N are held at zero and their effects on TPR, afferent resistance and reabsorption are switched off.
- `no_volume_pressure`: `f_s=0`, so fluid changes but venous filling/pressure does not.
- `afterload_step`: `R_sys` is multiplied by 1,20 from t=300 s.
- `renal_resistance_step`: `R_e` is multiplied by 1,25 from t=300 s.

No control arm is used to move the primary criterion after the run.

## Frysta parametrar

| parameter | value | unit | source or assumption |
|---|---:|---|---|
| `P_a0` | 100 | mmHg | Bikia 2024, Table 2, in-vivo MBP 100 ± 12 |
| `Q0` | 4.9 | L/min | Bikia 2024, Table 2, in-vivo CO 4.9 ± 1.2 |
| `C_a` | 1.4 | mL/mmHg | Bikia 2024, Table 2, in-vivo arterial compliance |
| `R_sys0` | 100/(4.9×1000/60) | mmHg·s/mL | derived from the same MAP/CO anchor; not freely parameterized |
| `P_v0` | 5 | mmHg | representing operating-point; not a measurement for Q017 |
| `C_sys` | 68.0 | mL/mmHg | Maas et al. 2012, 0.97 mL/mmHg/kg × 70 kg rounded; ICU population, assumption |
| `V0` | 17690 | mL | Bhave & Neilson 2011, Table 3, 73-kg man ECF; review-tier, not transient data |
| `f_stressed` | 0.30 | fraction | Magder reference in existing `venous_return.py`; reused as an assumption |
| `G0` | 122.8/60 | mL/s | Davies & Shock 1950, Table II, 20–29-year group |
| `Q_r0` | 1076.8/60 | mL/s | Davies & Shock 1950, Table II, effective renal blood flow |
| `P_gc0` | 60 | mmHg | local glomerular working point; assumption from existing nephron model |
| `P_bow` | 18 | mmHg | local Bowman's pressure; assumption |
| `pi_gc` | 32 | mmHg | glomerular oncotic pressure; assumption |
| `R_e0` | `(60-5)/Q_r0` | mmHg·s/mL | derived from `P_gc0=P_v+Q_rR_e` |
| `R_a0` | `(P_a0-P_v0)/Q_r0-R_e0` | mmHg·s/mL | derived from baseline flow and `P_gc0` |
| `K_f` | `G0/(P_gc0-P_bow-pi_gc)` | mL/(s·mmHg) | derived from the Starling equation |
| `f_reabs0` | 0.992 | fraction | controlled assumption; gives about 1 mL/min baseline urine |
| `L_other` | 0.005 | mL/s | controlled non-renal loss; not measured |
| `k_pre` | 8.0 | (mL/s)/mmHg | assumed preload–CO coefficient; not fitted to Kumar after the run |
| `k_con` | 6.0 | mL/s per aktivitetsenhet | assumed small sympathetic pump term |
| `k_A,k_S,k_N` on TPR | 0.8, 0.3, 0.4 | fraction per activity unit | reduced RAAS/NP effect; assumptions |
| `k_myo` | 0.015 | mmHg⁻¹ | assumed myogenic sensitivity |
| `k_TGF` | 0.8 | fraction per GFR deviation | assumed TGF sensitivity; direction from Ito & Abe 1996 |
| `k_aff,S/N` | 0.08, 0.05 | fraction per activity unit | assumed afferent tone effects |
| `tau_a` | 20 | s | assumed local feedback time |
| `tau_S,tau_N,tau_A` | 5, 30, 300 | s | assumed fast/medium/slow feedback times |
| `k_reabs,A/S/N` | 0.015, 0.008, 0.010 | fraction per activity unit | assumed hormonal reabsorption effect |
| `V_scale` | 500 | mL | controlled normalization of volume deviation |
| `bolus` | 3000 over 10800 from 300 | mL, s | Kumar 2004 study dose; synthetic protocol |
| `Qmax,Qrmax` | 150, 40 | mL/s | numerical/physiological safety limits; assumptions |
| `rtol, atol, max_step` | 1e-8, 1e-9, 5 | 1, 1, s | numerical setting |

No parameter may change after `PREREG.sha256` is written. The run's parameter list is written to `results.json` with the same values.

## Sensitivity plan

After the main run, three governing parameters are varied one at a time by ±50 percent, without simultaneous variation:

1. `k_pre` (volym/preload till cardiac output),
2. `k_myo` (renal autoregulation),
3. `tau_A` (slow RAAS lag).

For each parameter, `r_CO`, peak and final numbers for `P_a`, `Q_r`, `G`, `U`, `V` and the relative change against the main case are recorded. No ranking is chosen after the outcomes are known; the largest absolute percentage change is reported.

## Analytical limiting case and countertests

- With autofeedback switched off and fixed `R_a,R_e`, the numerical derivative of `Q_r(P_r)` must match `1/(R_a+R_e)` within 10⁻⁸.
- The baseline case without a bolus must keep `V`, `P_a`, `G` and `U` within 1 percent after 600 s.
- `K_f=0` must give `G=U=0` without creating a negative amount.
- An incorrect pressure zero with `R_a=0` must not pass as a numerical kidney case; this is used only as a diagnostic placeholder and is not part of the main criterion.

## Executability and delivery

`model.py` must implement the equations, parameter/unit checks, mass balance, main protocol, control arms and sensitivity. `test_model.py` must at least test the analytical pressure–flow limiting case, baseline balance, non-negative flows and unit checks. `results.json` must contain all final and peak counts with parameter and source status. `RESULTS.md` must distinguish source, derivation and assumption and must not claim that any internal transient has been measured.

## Boundary

This is a mechanistic core model, not a clinical risk model. It does not represent renal venous congestion, disease state, age, sex, training,RAS mutation details or a full pressure–natriuresis curve. Actual validation requires synchronous measurements of MAP/CVP or PCWP, cardiac output, renal blood flow, GFR, urine/volume and RAAS/NP signals under the same perturbation.
