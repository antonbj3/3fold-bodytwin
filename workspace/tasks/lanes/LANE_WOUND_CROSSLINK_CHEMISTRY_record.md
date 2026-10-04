# LANE_WOUND_CROSSLINK_CHEMISTRY

Den kirurgiska kedjans kvarvarande lucka (1/10). Resultatmapp `results/LANE_WOUND_CROSSLINK_CHEMISTRY/`.

## The location

- `results/LANE_SURGICAL_RESPONSE/` R3: U/I/M collagen inventories yield Levenson strength RMSE 8,31 pp **when measured HP is input**; R4: HP inventory rises ×2,22 dag 10→42 while collagen mass stands still (×0,93) — pure selective degradation excluded; no pre-registered candidate (immature→mature conversion, mechanosensitive maturation, selective turnover) passes early chemistry + HP42 + potency; The LH3 glycosylation variant hits but only illustrativeely.
- `results/LANE_SURGICAL_SYNTHESIS/` R3: with common collagen inventory, the integrated chain provides 21,8 pp autonomously; 8,3 requires external HP. Metric ranking in MEASUREMENT_SPEC_R3.md.
- `results/LANE_SURGICAL_BINDINGS/` PORTS: crosslink density ↔ fiber level strength (synthetic).

## Ability

A chemical mechanism for cross-link evolution in healing skin wounds that predicts HP (and preferably DHLNL/HLNL, LP, pyrrole) over time **from independent data**, so that the chain reaches the component level without external HP entry.

## Do like this

1. **Data Hunt First**: published time series of crosslinking chemistry in healing skin or incisional wounds (reducible divalent DHLNL/HLNL, mature HP/LP, glycosylation, LH2/LH3/LOX expression/activity), from the same or comparable models as Levenson (rat, incision). Freeze numbers with table/figure and art.
2. Break chemistry down to its smallest constituents: lysyl hydroxylation (LH1/LH2) → LOX oxidation → divalent → trivalent (spontaneous, depends on nearby telopeptide/helix positions and time) → glycosylation (LH3/GLT25D) that controls pathway selection. What rate steps can support a late HP increase without de novo synthesis?
3. Preregister candidates against data not used in alignment (eg align to divalents, predict HP; or align early, predict late). Judge against RESPONSE's HP42 = 0,043 mol/mol and Levenson's strength via the SYNTHESIS chain (`results/LANE_SURGICAL_SYNTHESIS/surgical_chain/`, run it with your chemistry port).
4. Strongest control: per-observable curves (Hill of HP) with the same data. Gain = a chemical mechanism that hits the HP course and via the chain the strength curve without strength fit.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
