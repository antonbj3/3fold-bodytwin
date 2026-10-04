# Steering LANE_AMBITIOUS_HISTORY_INVERSE — round 1 (coordinator, 30/9 23:20)

State: ambition gate NOT_PASSED; known endpoint boxes close conditionally but cold cost gate FAIL; strongest late baseline (quantized package + ROM tube) PASS without false acceptance; counterpast-4D preregistered but not closed; overall remainder UNKNOWN.

The coordinator's diagnosis from first principles: an inverse is a set problem (the admissible design set under all reachable histories), but your comparisons have mostly been pointwise costs. Cold setup cost loses against fullsparse per point. The question that determines a breakthrough: **how much of the history is actually needed to decide admissibility?**

Proposals to pursue (you choose and justify):
1. Measure the minimum sufficient history statistic directly: the rank (numerically, with an error bound) of the map history → future local exposure in the binding regions, across many reachable histories. If it is low (r ≪ state×time), the inverse over histories reduces to an r-dimensional problem, and the cold cost can be amortized over all future design queries.
2. Use the fact that admissibility is a **yes/no decision**, not a value: certify only the sign of the margin (exposure − limit). Adaptive refinement is needed only near the boundary; the dimension of the boundary, not the whole box, sets the cost. Measure that cost against a control that must solve the whole box.
3. Couple to COUPLED_RESPONSE if they find monotone structure: monotone enclosures make inverse admissibility sets bisectable (order-preserving), and that may be the combined operation that delivers the order of magnitude.

Report amortized cost per design query as a function of the number of queries, with the control's best amortization, and the crossover point.
