# LANE_SETTING_TO_DIFFUSIVITY — make the generator setting decisive, not gating

Results catalog `results/LANE_SETTING_TO_DIFFUSIVITY/`. Surgery, decisions outside the eye.

## The saturated mode
`tasks/assembly/incision_setting_decision.py` is built and posted as edge `T-E26`
. It reads published thermal damage depth **0,56–1,70 mm at 1–3 s** (PMID 12642261,
doi `10.1177/03635465030310021601` ) as `δ = √(D_eff·t)` and gets

- `D_eff = 3,1360e−07 m²/s` out of 0,56 mm at 1 s
- `D_eff = 9,6333e−07 m²/s` out of 1,70 mm at 3 s
- so **2,24× to 6,88×** soft tissue pure conduction 1,4e−07 m²/s

The energy balance is verified by the coordinator: evaporation 2260 J/m, fracture 4
J/m (**0,177%**), heating 37→100 °C exactly **+10,0%**, out of a 1 mm × 1 mm notch.

## Obstacle, exactly
Two things make the decision weaker than it looks:

1. **The setting does not decide, it gates.** The dwell limit is identical when setting 30 and
   60, because nothing published we have links the setting to `D_eff`. The setting enters
   only as the refusal above 60, where the damage correlation drops from R = 0,95/0,98/0,92 to R = 0,73.
2. **Reference pins `D_eff` only within a factor 3,0719**, making the stay case **2,0719×** themselves
   the conservative break. A decision whose envelope is twice its value hardly chooses anything.

## Operationen
Bygg kalibreringen `setting → D_eff` from published data, or show that it cannot be
identify out of depth-versus-time alone. Candidates already named in the swarm's material:

| source | what it carries |
|---|---|
| PMID 12642261, doi `10.1177/03635465030310021601` | damage depth 0,56–1,70 mm at 1–3 s human menisci |
| PMID 39161807, doi `10.1016/j.helvyon.2024.e35266` | thermal dispersion, maximum temperature, time to baseline temperature i liver/muscle/skin |
| PMID 40442741, doi `10.1186/s12893-025-02969-8` | same quantities, independent measurement |
| PMID 30294056, doi `10.1007/s00466-017-1529-6` | EOS for vapor mass fraction against T and p, validated on ex vivo pig liver at **multiple effects** |

The last one is key: it is validated at multiple power levels, which
is precisely the setting-dependent dimension that is missing.

## Strongest control
A single `D_eff` for all settings, i.e. current decision. About the setting dependent
the calibration does not narrow the stop envelope below 1,0× the conservative stop it has not bought
something, and then it won't go in.

## Falsifier
If `D_eff` cannot be identified from depth-versus-time in the case of several effects — for example, if effect and
dwell time only appears as a product in each published record — that is the principle of the setting
gating and not deciding, and that's the result. Print the quantity that would break
produkten.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no claim of biological validation, no excluded_category material, nothing that
names an individual person. Every source with PMID or DOI, volume and pages.
