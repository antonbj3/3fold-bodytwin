# PROOF_LANE_TWO_TIMESCALE — can two timescales be distinguished from one with a slope correction?

## The measured state, calculated by the coordinator today
I built `tasks/assembly/incision_setting_decision.py` and edge `T-E26`. The energy balance is
verified from a 1 mm × 1 mm kerf, water’s heat of vaporization 2,26 MJ/kg and soft tissue’s specific
heat 3,6 kJ/kg/K: vaporization **2260 J per meter of incision**, mechanical fracture at G = 1 kJ/m² **4 J/m**
(**0,177 %**), heating 37→100 °C exactly **+10,0 %**. The kerf is self-consistent: `2(w+d) = 4e−3
m²/m` and `w·d = 1 mm²` force `w = d = 1 mm`.

Published thermal injury depth 0,56–1,70 mm at 1–3 s (PMID 12642261,
doi `10.1177/03635465030310021601`), read as `δ = √(D_eff·t)`, gives `D_eff` 3,1360e−07 to
9,6333e−07 m²/s, i.e. **2,24× to 6,88×** soft tissue’s pure conduction 1,4e−07 m²/s.

Injury volume against generator setting in skeletal muscle: **R = 0,73 (p = 0,008)** over the entire 10–120, but
**R = 0,95 / 0,98 / 0,92 (p ≤ 0,001)** for depth / radius / volume restricted to 10–60.

## The question
The decision’s falsifier is a signature: if the chain has **two** timescales — diffusion building
the weakness front, evaporative vaporization maintaining it — injury size against log(setting) should have
slope ≈ 0,5 in the conduction-limited region and **saturate** (slope → 0) when the vaporization front
takes over. A pure diffusion model gives 0,5 across the entire interval.

But the source itself admits that **a small, sliding slope correction can reproduce a flattening
slope**. Thus: determine exactly under what conditions the two-timescale model is **identifiable** against a
one-timescale model with a free slope correction, and under which it is not.

I want the conditions as theorems with a measurable quantity per condition, and especially:

- which functional property of the break **no** sliding slope correction can mimic
  (curvature sign, second derivative, a finite moment, an order property — state which);
- how many measurement points and what spread are required to distinguish them on the published
  intervals 10–60 and 10–120;
- whether the difference between R = 0,95 and R = 0,73 is itself a sufficient witness, or only compatible
  with both models.

## The strongest control
The one-timescale model with a free slope correction, fitted on the same published points. If
the two-timescale model cannot be distinguished from it, the falsifier in my decision is empty, and then I should
know it.

## Falsifier
Construct a one-timescale model yourself with a slope correction that reproduces every property you
identify as impossible to mimic. If you succeed, your own condition is wrong and should be rejected, not patched.

## Warning
A width-based rounding rule fails — your own `PROOF_LANE_PRECISION_PROP` let one fail on
all 1 200 constructed boundary cases, and the right rule is that both envelope endpoints must round
the same way. Use the endpoint rule when stating precision.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no statement of biological validation, nothing that
names an individual person. Every published source with PMID or DOI, volume and pages.
