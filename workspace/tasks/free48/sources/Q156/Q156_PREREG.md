# BT-HX-Q156 — preregistration

Frozen before the first model run: 2026-09-25. This is a first runnable mechanistic test, not a fit to an unreported tracer curve.

## Hypothesis

A one-dimensional porous slab with a CSF reservoir, a perivascular-space (PVS) channel, and an interstitial-fluid (ISF) domain can reproduce the approximately 600 s diffusion-only half-filling time reported for a 3-kD paravascular tracer when the same 250 micrometre length scale and `D = 1e-10 m^2/s` are used. Adding PVS-only advection or reversible binding is not required to meet this first criterion; either mechanism is accepted as necessary only if it changes the prediction in a way that can later be distinguished by held-out spatial and temporal data.

## Frozen prediction and reference

Primary predicted quantity: `t50_ISF`, the first time at which the distal ISF sensor (the last tissue cell) reaches 50% of the imposed CSF step concentration. The volume-mean ISF concentration is reported separately, together with PVS concentration, CSF concentration, exchange flux, blood/lymph clearance flux, total mass balance, and the effective spreading time.

Published reference:

- Byung-Ju Jin, Alex J. Smith, and Alan S. Verkman, 2016, *Journal of General Physiology* 148(6):489–501, DOI `10.1085/jgp.201611684`.
- Source page: `https://pmc.ncbi.nlm.nih.gov/articles/PMC5129742/`.
- Cited location: Results, “How large a hydrostatic pressure gradient is needed to drive advective flow in the ECS?”, Fig. 2C discussion.
- Value: diffusion-alone half-filling/half-emptying time approximately **600 s** for a **3-kD tracer** with `D = 1e-10 m^2/s`; the paper describes this as similar to the experimental paravascular-to-parenchymal time course. The geometric reference is the 250 micrometre arteriole–venule spacing used in that model.
- Verification: **verified from the full-text PMC page by webfetch on 2026-09-25**. This is a published comparison value, not an internal or newly collected measurement.

## Frozen acceptance criterion

The primary diffusion-only prediction passes if all of the following hold:

1. `t50_ISF / 600 s` is between **0.5 and 2.0** (factor two).
2. Concentrations remain finite and non-negative, and the maximum relative mass-balance residual over the run is below `1e-6`.
3. The diffusion-only run is reported separately from the PVS-advection and binding alternatives; no mechanism is declared required from an unmeasured fit.

A failed primary criterion is reported as failed, not repaired by parameter tuning. UNKNOWN is used for a requested transport contribution that the present run does not measure.

## Mechanistic scope

The model is conservative on a 1-D axial domain. CSF is a well-mixed reservoir at the proximal end; the PVS is a longitudinal fluid channel; ISF is a parallel porous domain; and explicit exchange and sink terms represent the vascular/lymphatic boundary. Diffusion is present in both PVS and ISF. Advection is allowed only in the PVS and is pressure-linked through a local Darcy relation; ISF bulk advection is zero. Reversible binding is represented as a free-to-bound exchange in ISF and is disabled in the primary run.

The finite-volume equations are:

```text
d(eps_i V c_i)/dt = div(-A eps_i D_i grad(c_i)) + J_PVS->ISF - J_blood,i - J_bind
d(eps_p V c_p)/dt = div[-A eps_p D_p grad(c_p) + A u_p c_p] - J_CSF->PVS - J_PVS->ISF - J_blood,p - J_out,p
d(V b)/dt = J_bind
d(V_CSF c_CSF)/dt = -J_CSF->PVS - J_out,CSF
J_bind = k_on V_i c_i - k_off V b
u_p = max(0, k_h (P_a - P_v) / (mu L))
```

Here concentrations are per fluid volume, `b` is bound tracer per total volume, face fluxes are in mol/s, and all source/sink terms are included in the mass balance. The PVS–ISF exchange uses the harmonic interfacial volume `2 V_p V_i/(V_p+V_i)` times `k_exchange (c_p-c_i)`, so equal concentration is the exchange equilibrium even when fluid fractions differ. The numerical scheme is explicit Euler with a conservative face-flux update.

## Parameters and assumptions

| parameter | symbol | baseline | unit | source or assumption |
|---|---:|---:|---|---|
| axial length | `L` | `250e-6` | m | Jin et al. 2016 model geometry |
| number of cells | `N` | `40` | 1 | numerical resolution; 6.25 micrometre cells at baseline |
| cross-sectional area | `A` | `1e-9` | m2 | unit-area normalization; scales flux and reservoir choice |
| ISF volume fraction | `eps_i` | `0.20` | 1 | Jin et al. Fig. 3; measured-model comparison value |
| PVS volume fraction | `eps_p` | `0.10` | 1 | explicit modelling assumption |
| ISF effective diffusivity | `D_i` | `1e-10` | m2/s | Jin et al. parameter table for 3-kD tracer |
| PVS effective diffusivity | `D_p` | `1e-9` | m2/s | fluid-channel assumption |
| CSF volume | `V_CSF` | `1e-9` | m3 | well-mixed reservoir assumption |
| CSF exchange rate | `k_CSF` | `1.0` | 1/s | prescribed-step boundary conductance; CSF reservoir feeds the proximal PVS interface |
| PVS–ISF exchange rate | `k_ex` | `1e-2` | 1/s | explicit membrane/interface assumption |
| ISF blood clearance | `k_b,i` | `0` | 1/s | primary diffusion comparison isolates filling |
| PVS blood/lymph clearance | `k_b,p` | `0` | 1/s | primary diffusion comparison isolates filling |
| PVS pressure difference | `DeltaP` | `0` | mmHg | null advection model |
| PVS hydraulic conductivity | `k_h` | `1.9e-15` | m2 | calibrated assumption for a `1 mmHg` alternative |
| dynamic viscosity | `mu` | `1e-3` | Pa s | water-like fluid assumption |
| binding association rate | `k_on` | `0` | 1/s | disabled in primary run |
| binding dissociation rate | `k_off` | `0` | 1/s | disabled in primary run |
| CSF outflow rate | `k_out,CSF` | `0` | 1/s | closed source reservoir in primary run |

Dimensions are checked in `model.py`: concentration times volume is mol; `eps*D*grad(c)*A` and `k*V*c` are mol/s; `k_h*DeltaP/(mu*L)` is m/s; reaction terms are mol/s.

## Controls and falsification

The pure-diffusion run is the null model. The PVS-advection run changes only `DeltaP`; it does not add a generic ISF bulk flow. The binding run changes only `k_on` and `k_off`. A proposed added pathway is rejected if it cannot be identified by held-out spatial/temporal data or if its net flux violates the mass balance. The present run cannot establish sleep dependence, pressure dependence in vivo, or a universal glymphatic mechanism.

## Builds on

- Q156, `inputs/QUESTION.md`, especially its stated `eps*c` conservation equation and the warning that the existing Peclet fragment is not general evidence for CSF/glymphatic transport.
- `src/bodytwin/cells/nervous/glymphatic_peclet_bound.py`, as identified by the brief; no code was imported and the external repository was not modified.
- The conservative 1-D advection–diffusion comparison in Jin, Smith & Verkman (2016), DOI `10.1085/jgp.201611684`.

## Not redone

No internal data, no claimed tracer measurement, no three-dimensional vascular geometry, no measured pressure field, no fitted binding parameters, and no universal sleep/glymphatic mechanism. The next resolution step is simultaneous 3-D geometry plus time-resolved tracer concentrations and exchange-flux constraints.
