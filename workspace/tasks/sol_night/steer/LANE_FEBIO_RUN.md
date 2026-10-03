# LANE_FEBIO_RUN — run the external model that has never been run

Resultatmapp `results/LANE_FEBIO_RUN/`.

## Why, and what is already settled
Seed 1: external biological models must be connectable so that more questions become askable. The blocker has
been recorded as untested interoperability — a solver connection was rejected at 0 % valid frames, and
FEBio has never been run.

**But connection is no longer the hard part, and that was measured tonight.** Another external model was run
against a held-out measurement: 39 converged state frames, 100 matched force points, and a peak error of
**40,89 %** (1722,0983 N versus 1222,2613 N). Geometry was ruled out as an error source (CAD proxy shifts the error
6e-5), sensor direction too (5° cone takes 40,9 % to 39,5 %), and timestamps quantitatively (0,881 N of
499,84 N, hence 0,18 %).

**And the selection criterion is new:** an external model reproducing a quantity we already compute removes no
rank deficiency and therefore makes no new question askable. Its usefulness lies in whether it carries an observable we
lack. Measured: 13,3 % of 5 984 swarm reports make an explicit claim about structural
unidentifiability, versus 0,8 % complaining of parameter uncertainty.

## Do this
1. **Determine FIRST which observable FEBio would carry that we do not have.** Read
   `notes/RANK_LIFT_OBSERVABLES.json` — 70 jobs name the observable that would remove their
   rank deficiency. Match against what FEBio actually computes. If there is no match, the run is not
   worth anything beyond showing that the connection works, and **then you must say so before you run**.
2. Run the smallest run that gives a quantity with a unit, against a reference answer that is not ours. Report
   convergence, number of valid frames, and error against the reference answer — in that order.
3. State the provenance class for each number: INDEPENDENT_MEASUREMENT / PEER_REVIEWED / STANDARD / VENDOR /
   UNSOURCED. A vendor number is not a measurement.

## Control and falsifier
- **Control:** our own cell for the same quantity, with the same inputs. The gain is an observable we did not
  have, not a prettier number.
- **Falsifier:** if FEBio only reproduces quantities we already compute, it is irrelevant to seed 1
  however well it connects — report it as a selection result, it is useful.
- **Forbidden:** counting a successful run as the finding; calibrating against our own number.

## Delivery
`PORT.json`: the observable FEBio carries, its unit and level, the rank deficiency it removes if any, and
the error against an external reference answer.

# New round (coordinator, 3/10 00:25)
You refuted a universal source claim with 136 counterexamples out of 144, which is a real finding — but seed 1’s blocker remains unchanged: no external model ran. Changed operation, and do the selection step FIRST: read notes/RANK_LIFT_OBSERVABLES.json and say which observable FEBio carries that lifts a named rank deficiency. If none exists, say so and close the lane — it is a selection result and it is useful. If one exists, run the smallest run giving that quantity against external ground truth, and report convergence, valid frames and error against ground truth in that order. Compare with how the other external model went: 39 converged frames, 100 matched force points, 40,89 % peak error. Prohibited: counting a successful run as the finding.
