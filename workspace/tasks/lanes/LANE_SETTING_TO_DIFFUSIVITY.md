# LANE_SETTING_TO_DIFFUSIVITY — make the generator setting decisive, not gating

Results directory `results/LANE_SETTING_TO_DIFFUSIVITY/`. Surgery, decisions outside the eye.

## The saturated state
`tasks/assembly/incision_setting_decision.py` is built and recorded as edge `T-E26`. It reads
published thermal damage depth **0,56–1,70 mm at 1–3 s** (PMID 12642261,
doi `10.1177/03635465030310021601`) as `δ = √(D_eff·t)` and obtains

- `D_eff = 3,1360e−07 m²/s` from 0,56 mm at 1 s
- `D_eff = 9,6333e−07 m²/s` from 1,70 mm at 3 s
- thus **2,24× to 6,88×** soft tissue’s pure conduction 1,4e−07 m²/s

The energy balance is verified by the coordinator: evaporation 2260 J/m, fracture 4 J/m (**0,177 %**),
heating 37→100 °C exactly **+10,0 %**, out of a 1 mm × 1 mm kerf.

## The obstacle, precisely
Two things make the decision weaker than it looks:

1. **The setting does not decide, it gates.** The dwell limit is identical at settings 30 and
   60, because nothing published that we have links the setting to `D_eff`. The setting enters
   only as refusal above 60, where the damage correlation drops from R = 0,95/0,98/0,92 to R = 0,73.
2. **The reference constrains `D_eff` only within a factor 3,0719**, making the dwell envelope **2,0719×**
   the conservative dwell itself. A decision whose envelope is twice its value hardly chooses anything.

## The operation
Build the calibration `setting → D_eff` from published data, or show that it cannot be
identified from depth-versus-time alone. Candidates already named in the swarm’s material:

| source | what it carries |
|---|---|
| PMID 12642261, doi `10.1177/03635465030310021601` | damage depth 0,56–1,70 mm at 1–3 s human menisci |
| PMID 39161807, doi `10.1016/j.helvyon.2024.e35266` | thermal spread, maximum temperature, time to baseline temperature in liver/muscle/skin |
| PMID 40442741, doi `10.1186/s12893-025-02969-8` | same quantities, independent measurement |
| PMID 30294056, doi `10.1007/s00466-017-1529-6` | EOS for vapour mass fraction against T and p, validated on ex vivo pig liver at **multiple power levels** |

The last one is key: it is validated at multiple power levels, which is exactly the
setting-dependent dimension that is missing.

## Strongest control
A single `D_eff` for all settings, i.e. the current decision. If the setting-dependent
calibration does not narrow the dwell envelope below 1,0× the conservative dwell, it has not bought
anything, and then it does not enter.

## Falsifier
If `D_eff` cannot be identified from depth-versus-time at several power levels — for example, if power and
dwell time only appear as a product in each published record — then the setting is in principle
gating rather than deciding, and that is the result. Write out the quantity that would break
the product.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no biological validation claim, no excluded_category material, nothing that
names an individual. Each source with PMID or DOI, volume and pages.
