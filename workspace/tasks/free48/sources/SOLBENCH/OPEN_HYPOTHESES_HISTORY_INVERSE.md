# Direction LANE_AMBITIOUS_HISTORY_INVERSE — round 1 (the coordinator, 30/9 23:20)

The state: ambition gate NOT_PASSED; known endpoint boxes close conditionally but cold cost gate FAIL; strongest late baseline (quantized package + ROM tube) PASS without false acceptance; counterpast-4D preregistered but not closed; overall remainder UNKNOWN.

The coordinator's diagnosis from first principles: an inverse is a set problem (the admissible design set under all reachable histories), but your comparisons have mostly been pointwise costs. Cold setup cost loses against fullsparse per point. The question that decides breakthrough: **how much of the history is actually needed to decide admissibility?**

Proposals to attack (you choose and justify):
1. Measure the smallest sufficient history statistic directly: the rank (numerically, with an error bound) of the mapping history → future local exposure in the binding regions, over many reachable histories. If it is low (r ≪ state×time) the inverse over histories is reduced to an r-dimensional problem, and the cold cost can be amortized over all future design questions.
2. Use the fact that admissibility is a **yes/no decision**, not a value: certify only the sign of the margin (exposure − limit). Adaptive refinement is needed only near the boundary; the dimension of the boundary, not the entire box, sets the cost. Measure that cost against a control that must solve the entire box.
3. Connect to COUPLED_RESPONSE if they find monotonic structure: monotonic enclosures make inverse admissibility sets bisectable (order-preserving), and it may be the combined operation that gives the order of magnitude.

Report amortized cost per design question as a function of number of questions, with the control's best amortization, and the crossover point.
