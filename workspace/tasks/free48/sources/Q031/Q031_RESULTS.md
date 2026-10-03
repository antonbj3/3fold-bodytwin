BT-HX-Q031

# Resultat

## Kort svar

The frozen limit is that the membrane potential **dynamics** must be coupled to the volume if the voltage clamp changes the volume more than 2 % of the initial volume, `chi_V = max_t |V_full(t) - V_clamp(t)| / V0 > 0,02`. In the base case, `chi_V = 0,008394857` became; thus, full dynamic voltage coupling does not meet this requirement for this particular protocol. The voltage must still exist in the ion flows through the Nernst law. The conclusion is not universal: the pumping density +50 % gave `chi_V = 0,018841727`, close to the limit.

## Vad som byggdes

`model.py` is a coupled ODE-model with membrane potential `Vm`, cell volume `V`, three ion quantities, water quantity, cumulative ion fluxes and a small countercharge reserve. Na/K-pump moves 3 Na out and 2 K in per cycle. The ion fluxes are given by Nernst potentials and linear conduction; the water flow is given by RT-baserad osmotic pressure balance. The surface scales with the volume, so the same geometry changes capacitance, conduction, pump flow and water flow.

`PREREG.md` was frozen before the run. Its current SHA-256 is in `PREREG.sha256` and the same hash is in `results.json` under `prereg_sha256`. The local directory had no previous model, test, or results files to continue from.

## Source, derivation and assumptions

**Source.** Farinas and Verkman, Biophysical Journal 71:3511–3522 (1996), DOI `10.1016/S0006-3495(96)79546-2` , Figure 6 and results text pp. 3516–3518. There, for the MDCK cell layer approximately double cell volume after the 300 to 150 mOsm step, `t50 cirka 6 s` , `Pf = 6,1 ± 2 × 10^-4 cm/s` , initial volume `1,7 pL` and surface/volume `6750 cm^-1` are reported. Jakobsson 1980, DOI `10.1152/ajpcell.1980.238.5.c196` , used only as structural support for transport plus bulk electroneutrality and isotonicity.

**Hderivation.** Pf was converted to a pressure coefficient by `Lhyd = 1,65 × 10^-11 m s^-1 Pa^-1`. Placebo gives `V∞/V0 = 2,000000000`, exactly from 300/150 mOsm. Other values ​​in `parameter_table` are assumptions or assumed constants; they are not claimed measurements.

**Important uncertainty.** Pf and surface/volume come from a cell layer, not a free sphere. The hydraulic surface efficiency `0,025` is an explicit bearing geometry assumption. Therefore, the reference gate is a reproducibility check, not external validation.

## Basresultat

- Final volume of the full model: `V30/V0 = 1,920899334` .
- HHalf-life to half the final volume change: `t50 = 6,469760674 s`.
- Idealisk placebo: `V∞/V0 = 2,000000000`.
- Initial membrane potential: `-80,376985 mV`; terminal potential: `-74,745154 mV`; largest transient deviation: `5,631831 mV`.
- Terminal quantities: `NNa = 1,617280222 × 10^-14 mol`, `NK = 2,378644335 × 10^-13 mol`, `NCl = 1,603629320 × 10^-14 mol`.
- Terminala koncentrationer: `[Na,K,Cl] = [4,952583,72,841014,4,910780] mM`; total osmolaritet `155,587193 mM`.
- Voltage clamp final volume: `1,929294191`; largest volume difference `chi_V = 0,008394857`.

The reference port `1,8 ≤ V∞/V0 ≤ 2,2` and `3 s ≤ t50 ≤ 10 s` were passed. The frozen 2 % limit for mandatory dynamic stress coupling was not passed in the base case.

## Balanser

All frozen balance requirements on `10^-6` were passed:

- highest relative ion abundance error: `1,113685984 × 10^-14`;
- highest relative water balance error: `8,792443335 × 10^-15`;
- highest relative charge balance error: `1,954750880 × 10^-15`;
- highest relative bulk electroneutrality error: `3,783673394 × 10^-10`;
- Greatest countercharge reserve: `1,848063561 × 10^-6` of initial osmolarity amount.

The charge control used is end-to-end tracked ion charge plus integrated net ion flux. The capacitance `Cm A dVm/dt = -(I_p + I_ion)` is reported as voltage law, not as an additional net exchange of ion charge.

## Sensitivity, accurate ±50 %

| Parameter | -50 % | +50 % | Reading |
|---|---:|---:|---|
| Hydraulisk permeabilitet | `V∞/V0=1,737448`, `t50=9,523521 s`, `chi_V=0,001989833` | `V∞/V0=1,974359`, `t50=4,682843 s`, `chi_V=0,013928960` | Controls time scale; -50 % misses reference port. |
| K lead | `V∞/V0=1,929967`, `t50=6,547059 s`, `chi_V=0,011859088` | `V∞/V0=1,916176`, `t50=6,429798 s`, `chi_V=0,003905153` | Small effect on final volume and t50. |
| Pump density | `V∞/V0=1,923359`, `t50=6,489584 s`, `chi_V=0,000620641` | `V∞/V0=1,919018`, `t50=6,455674 s`, `chi_V=0,018841727` | Small effect on volume, but largest effect on the neighborhood of the coupling boundary. |

No sensitivity run violated the balance requirement. The results are a temporal model, not a simultaneous measurement finding, as no simultaneous measurement data were included.

## What failed or cannot be determined

- The negative base coupling conclusion applies only to the generic model and the chosen protocol; it is not a universal cell biological limit.
- The pump density +50 % reaches close to the 2 % limit, so the conclusion is sensitive to an assumed pump value.
- `Lhyd` and the hydraulic surface efficiency are not independently measured for this particular model.
- The channels are constantly open; channel activation, ATP-kinetik, ion activities, hydrostatic pressure, membrane mechanics and the cell type's other ions are missing.
- No new measurement data was created. The reference match must not be described as validation.

## Next resolution step

1. Replace the generic leads and pump value with the measured or directly identified parameters of the cell type.
2. Measure `Vm`, at least two ion amounts/concentrations and volume in the same cell during the same osmotic step; use these as a common failure scenario.
3. Identify hydraulic permeability and exposed surface directly instead of combining cell layer Pf with `0,025`.
4. Map `chi_V` over a parameter field and combined perturbations to find the biological limit, not just the one-point limit.
5. Add channel gating and pressure/elasticity first if the data show that their time constant is close to the water or stress relaxation.
