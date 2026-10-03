BT-HX-Q156

## Result
Source verified: Jin, Smith & Verkman (2016), *Journal of General Physiology*, DOI `10.1085/jgp.201611684`, Fig. 2C Text, `https://pmc.ncbi.nlm.nih.gov/articles/PMC5129742/`. It states approximately `600 s` diffusion-only half-filling for 3-kD tracer (`D=1e-10 m^2/s`); source is published comparison, not internal measurement.

**Diffusion-only (hypothesis):** `t50_ISF=261.17337075283876 s` distal ISF; reference ratio `0.43528895125473127`. The prediction is thus faster than the reference. Integrated flows: CSF→PVS `7.499437541955433e-14 mol`, PVS→ISF `4.999625027939717e-14 mol`, blood/lymph `0 mol`. Mass balance residual `6.377562522783843e-13`; all concentrations are finite and non-negative.

**Approval:** NO. The criterion requires `0.5 <= t50/600 <= 2.0`; `0.43528895125473127` falls below the limit despite good mass balance. This is a robust failure, not an UNKNOWN.

**Options:** PVS advection (`1 mmHg`) gives no `t50` within `7200 s`, exports `2.7551151492305117e-12 mol` through PVS outlet and has residual `3.881650037568785e-15`. Binding (`k_on=1e-3`, `k_off=2e-4`) gives `t50=252.38482683547525 s` and residual `2.1920286232655233e-14`. No alternative transport is identified here as necessary.

## Sensitivity
Largest `t50` span at ±50 %: `k_exchange` (`367.4172995447236` → `235.60378519593647 s`, relative `-0.5046973738893494`), `L` (`231.1883012470413` → `331.41164072183716 s`, `0.38374256604300655`), `D_isf` (`265.50682383364114` → `259.13038166591457 s`, `-0.02441459536761467`). This makes yield and geometry larger control parameters than diffusivity in this particular model.

## Next resolution step
Simultaneously measure time-resolved CSF-, PVS- and distal ISF concentrations as well as build-up fluxes at known tracer size; simultaneously measure wall geometry and PVS/CSF pressure. A pressure or sleep perturbation must distinguish PVS advection from diffusion. Current data do not allow to select unique pathway or quantify binding.
