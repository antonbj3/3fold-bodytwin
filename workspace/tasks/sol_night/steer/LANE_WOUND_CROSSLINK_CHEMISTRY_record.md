# Direction LANE_WOUND_CROSSLINK_CHEMISTRY — round 2 (the coordinator, 1/10 07:30)

Review of round 1: good data hunt — two primary workbooks (55 sheets, 5 genotypes) frozen. Decisive diagnosis: with measured available site capacity the ceiling becomes **0,0242 HP per collagen** (verified r1/CHAIN_SCORE_v1.json), while RESPONSE's target HP42 = 0,043 — no faster rates can reach the target. Either the chemistry is wrong, or capacity and target are measured against **different references**.

## Changed operation: reconcile the references first (first principles)

The night's source review showed that unit/reference errors are common. HP is reported in the literature as mol HP per mol collagen (triple helix), per α-chain, per 1000 amino acids or per hydroxyproline — factors up to ~3 (chains per helix) and more. Do this:
1. Trace the exact reference and method (HPLC/fluorescence, hydrolysis, standard) for RESPONSE's HP42 = 0,043 (results/LANE_SURGICAL_RESPONSE/WORLD_TARGETS_R3*.json and the source paper) and for every number in your workbooks. Make a table reference → conversion factor → value in a shared unit.
2. Recalculate the capacity ceiling and target in the same unit. If the gap 0,0242 against 0,043 disappears (or reverses) it is a reference error — report it clearly; it then also affects RESPONSE and SYNTHESIS.
3. Only if the gap remains: test the chemistry (available pairs, glyco pathway, product fate) against data not used in fitting.
