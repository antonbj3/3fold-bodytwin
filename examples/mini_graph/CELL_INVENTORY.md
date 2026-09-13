# BODYTWIN CELL INVENTORY — the bio build-goal map

Synthesized from the cell-design catalog into the dependency graph. **2 cells, adversarially audited (none MEASURED yet — verdicts are on the DESIGN).**

## Certifiability (design-audit verdicts)

| bucket | count | meaning |
|---|---|---|
| REAL | 0 | design survives adversarial verify + external fact-check → ready to EXECUTE |
| WEAKENED | 0 | fixable defect, residual noted per cell |
| RE-SCOPED | 0 | honest weaker claim (archive lacks the ideal instrument) |
| REFUTED | 0 | anchor/leg defect not yet resolved |
| FENCED / HONEST-NEG | 0 | no valid anchor in archive (=PASS, non-promotable) |
| un-audited | 2 | (honest-neg no-anchor cells / node-0) |

**REAL (execute these first):** 

Each cell = an occluded hidden state cornered by ≥2 decorrelated legs vs a held-out external anchor. `occluded✓` = needs a consequence leg.


## foundational  (1 cells)

### MODEL-EXAMPLE-MEMBRANE-POTENTIAL  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **claim:** A Goldman-Hodgkin-Katz model with measured ion concentrations and permeability ratios reproduces a measured resting membrane potential within 10 mV.
- **hidden state:** intracellular ion activities
- **legs:** {'name': 'GHK from measured concentrations', 'source_ids': ['SRC-CONC-1']} · {'name': 'electrogenic pump contribution', 'source_ids': ['SRC-PUMP-1']}
- **anchor (held-out):** independently measured resting potential in the same preparation  ·  **mode:** agreement
- **regime:** adult mammalian, 37 C
