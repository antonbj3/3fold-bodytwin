# Steering LANE_SURGICAL_RESPONSE — round 4 (coordinator, 1/10 03:15)

Review of round 3: **the first mechanistic result tonight that beats curve fitting on a held-out world curve** — U/I/M collagen inventories with selective proteolysis, amount fitted only to day 28, chemistry only to HP day 0–10, **no strength fit**: Levenson held-out days RMSE 8,31 pp against Hill 14,91 and previous native 25,25 (verified r3/crosslink_port/summary.json). Strictly: an equally informed per-observable control still gives TIE, and the nominal strength gate FAIL. But the chemistry misses: HP day 42 is predicted as 0,020 against measured 0,043 mol/mol (≈52 % too low) — constant LOX-driven maturation cannot carry the late increase.

## Obstacles and changed operation

The late crosslink increase requires a process that continues after LOX activity (and new synthesis) has declined. Candidates from first principles:
1. **Immature → mature crosslink conversion** (divalent reducible → trivalent HP/LP) with its own slow kinetics, independent of LOX after initiation.
2. **Mechanosensitive LOX/maturation**: loading of the scar (increasing as strength grows) drives continued maturation — feedback strength → crosslink → strength.
3. **Selective turnover**: immature collagen is degraded and replaced, raising HP per collagen without new LOX.
Derive what each candidate predicts for HP over time and for strength at day 61/90 **before** seeing HP42; freeze; judge against HP42 and Levenson's late days. Falsifier: the candidate hitting HP42 but ruining strength RMSE is wrong.

## Koppling

SKIN_TOUGHNESS_GAP shows that mm-scale fiber bridging is the only mechanism reaching skin tear toughness. The scar's fracture strength over time may partly be bridging that rebuilds; use their PORTS if available.

Strongest control: per-observable curves with the same data. Gain = one parameter set hitting the HP curve including day 42 and the strength curve without strength fitting.
