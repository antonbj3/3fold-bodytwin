BT-HX-Q090
# Result

## Status

First runnable mechanistic model is complete. `PREREG.md` was frozen before execution and the checksum is `a84ba74c54f80b296bccb00c26806236509770a3c24057c4e9eefe4f196c8ca4`. No measured load distribution was found in the work catalogue; the results below are therefore model predictions, not claimed measurement data.

## Building on

`BRIEF.md`, `inputs/QUESTION.md` and `inputs/NIGHT_PREAMBLE.md`. There was no previous model or results file to reuse. The model couples a rigid jaw plate with two molar clearances and two bilateral TMJ paths. PDL fibers are homogenized into an axial modulus; PDL and alveolar bone are in series. TMJ disc and two CART bearings are in series. Darcy damping provides a reversible kinetic term. Momentum balance, energy balance and dimensional checks are included.

## Prediction against reference

The control case is a single `100 N` tooth path. The model gives `39,218977 µm` tooth displacement and `500 kPa` PDL mean pressure. The reference anchor is `0,08 mm` (`80 µm`) at nominal `100 N`, according to **Bien & Topp (1970), “Changes in the periodicity of the tooth mobility pattern during mastication”, Journal of Periodology — VERIFIED, UR MINNET**. The ratio is `0,490237`. The frozen arithmetic criterion is quotient within `[0,5; 2,0]`; it becomes `PROXY_FAIL`. Since the source and load protocol are not verified, external validity is `UNKNOWN`, not an alleged defeat.

## Base load distribution

Symmetric scenario: `F_close = 100 N`, `M_external = 0`.

| Magnitude | Left | Right | Total / Share |
|---|---:|---:|---:|
| Tooth contact force | 2,469193 N | 2,469193 N | 4,938387 N / 4,938387 % |
| TMJ force | 47,530807 N | 47,530807 N | 95,061613 N / 95,061613 % |
| Tooth displacement | 2,795887 µm | 2,795887 µm | — |
| PDL pressure | 12,345967 kPa | 12,345967 kPa | — |
| TMJ angle proxy | 0,008010° | 0,008010° | — |

The translation of the jaw is `2,795887 µm`; the symmetric pitch is `0°`. This is a model load distribution, not a measurement. The work is `0,139760 mJ` and stored elastic energy `0,139760 mJ`.

## Sensitivity

Each parameter was changed separately with `±50 %`.

1. `t_PDL`: largest impact, maximum absolute relative change `94,856507 %`. The tooth ratio becomes `9,622768 %` at `0,5×` and `3,330935 %` at `1,5×`.
2. `E_PDL`: second largest, `48,698176 %` maximum absolute relative change.
3. `E_TMJ_disc`: `0,958975 %`; the map layers dominate the assumed TMJ serial path.

The mechanical placebo case — all compliances scaled equally — preserves load shares with a maximum difference of `1,387779×10^-17`.

## Verification

The equilibrium, dimensionality, analytical limit case, transient and placebo criteria are `PASS`. The transient control reaches static end position with post-ramp inertia ratio `0,001444945`. The model run and tests are reproduced with:

```text
python3 model.py
python3 -m unittest -v
```

All tables, figures, sources and statuses are in `results.json`; code and tests are in `model.py` and `test_model.py`.

## Next resolution step

Measure the individual's `CT/microCT` geometry, PDL fiber angles and ply thicknesses as well as simultaneous force–motion filming. Use the values in a coupled 3D contact/tissue FE model with open/closed TMJ status, and calibrate against the same loading protocol as the reference. Until then, the load distribution must be reported as `UNKNOWN` externally and used only as a hypothesis-generating model.
