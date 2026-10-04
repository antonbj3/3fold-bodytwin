# LANE_ANCHOR_PREDICTION

Anton 1/10: "ground truth is only data anchors, our cell simulations should be able to compute this without constants in a textbook." Results directory `results/LANE_ANCHOR_PREDICTION/`.

## Anchors (verified or to verify)

- **ATP flux in C2C12 myoblasts, basal: 41,7 pmol ATP/min/µg protein** (Mookerjee, Gerencser, Nicholls & Brand, J Biol Chem 292:7189, 2017, PMC5409486 — the coordinator has opened the page and seen the quote verbatim). Relevant families: MITOSTRESS (`tasks/free48/sources/MITOSTRESS/`), Q077 (cristae/ATP), Q088 (oxygen/ATP-PCr).
- More from `results/BT-GRAF-FACIT1-*/FACIT_BINDING.json` and `results/LANE_EXTERNAL_FACIT_HUNT/FACIT_INDEX.json` (when it exists) — **only those where the system matches** (same species/tissue/regime; FACIT1-03 steel back-stress is an example of the wrong system).

## Task per anchor

1. Break the quantity down to lower levels: which molecular/cellular mechanisms determine it (e.g. for ATP flux: oxygen consumption × P/O ratio, glycolytic lactate flux × ATP/lactate, proton leak, mitochondrial content per µg protein, enzyme capacities)? Which of these exist in our source models as a mechanism, and which are there as a constant?
2. Build or connect a computation that **predicts** the anchor from the lower levels, where each input is a more fundamental quantity (with source and uncertainty) — **the anchor value and textbook constants for the same quantity must not be among the inputs**. Preregister the prediction before comparing.
3. Compare against the anchor with propagated uncertainty. Report: predicted value ± uncertainty, anchor, ratio, and which lower quantity carries the most uncertainty (= next measurement or next mechanism to build).
4. Where the source model has the quantity as a constant: write a minimal proposal for a mechanism that replaces the constant (in the results directory; the sources are changed by the coordinator after review).

Form C with first-principles requirements. No internal data.
