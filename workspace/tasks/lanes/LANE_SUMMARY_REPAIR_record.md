# LANE_SUMMARY_REPAIR — repair the cells whose missing quantity they already carry

Result directory `results/LANE_SUMMARY_REPAIR/`.

## Why this lane exists, and why it is cheap
Thirty-two cells fail a sufficiency test: the summary their readout uses does not determine the readout.
The analysis that established this also named, per cell, the minimal extension — and for a handful the
named quantity **already exists as a variable inside the cell**. Those are repairable today with no
measurement at all. I verified three by hand rather than trusting a substring match, which is how I
found that a fourth was a false hit in prose rather than a variable:

| cell | mismatch type | the quantity the readout needs | verified present as |
|---|---|---|---|
| Q009 | linear split plus time | bound concentration | `cp_bound`, `auc_plasma_bound`, `d_bound`, `initial_tissue_bound_mg_l` |
| Q022 | co-factor | strain | `strain`, `strain_array`, `initial_strain` |
| Q036 | time | injury block | `injury_block` |
| Q005 | co-factor | `P_MATE` | named in `Q005_model.py` — verify as a variable before using |
| Q168 | time | portal TMA | named in `Q168_model.py` — verify as a variable before using |
| Q121 | locator | proximal gas | **NOT a variable**; the match was prose. Do not attempt this one. |

The characterisation behind it is about locality, not nonlinearity: a summary suffices if and only if
what the readout actually asks for sits at the same place, the same instant and the same weighting, with
no second independent coordinate. Failures are of exactly three kinds — wrong locator, wrong time, missing
co-factor — and each independent mismatch costs one scalar. Note that the rule TIES against an equally
informed control (the run's own `summary_form` label separates all 44 cases), so do not report the rule
as new; report the repairs.

## Do this
1. **Repair one cell at a time, smallest first.** Extend the summary the readout consumes with the
   quantity already present, change nothing else, and re-run the cell's own sufficiency test.
2. **Report the verdict both ways per cell**: did it fail before, does it pass after, and the identity
   error and downstream difference in both cases. A repair that does not flip the verdict is a negative
   result and must be reported as one — it would mean the named extension was the wrong one.
3. **And say which direction the old gate was wrong in.** If the unrepaired cell passed its scalar gate
   while failing the distributional one, it was FALSELY SAFE, and that is the hazardous direction. Count
   falsely safe and falsely cautious separately; a symmetric count hides the one that costs a claim.
4. **Do not touch the five cells whose extension is unavailable** (`pore diameter`, `contact area`,
   `graph distance`, `peak density`, `angles`). Those are acquisition items; write them in the standard
   schema with quantity, unit and what they decide, and leave them.
5. **Verify the quantity is a variable, not a word**, before each repair. One of six candidates failed
   exactly that check.

## Control and falsifier
- **Control:** the same cell with its current summary, on the same inputs. Equally informed by
  construction — same model, same data, only the summary the readout consumes differs.
- **Falsifier:** if extending the summary with the named quantity does not change the verdict in at
  least half the repaired cells, the per-cell naming is unreliable and the whole repair list needs
  re-deriving rather than applying. Report that as the headline if it happens.
- **Forbidden:** changing a cell's physics or parameters — only what the readout consumes; reporting a
  repair without the before-and-after identity error; reporting the locality rule as a new result when
  it ties against the existing label.

## Delivery
`PORT.json` with one row per repaired cell carrying both verdicts, both identity errors, the downstream
difference, and the falsely-safe or falsely-cautious direction; plus the acquisition items for the five
cells that cannot be repaired.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
