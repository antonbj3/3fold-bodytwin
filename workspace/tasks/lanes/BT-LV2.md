# BT-LV2 — certified case selection for muscle paths (wrapping transitions)
First read tasks/lanes/_PREAMBLE.md (binding). Authorized to run without asking. Light CPU. Run after BT-LV1 if possible (reuse its interval tools).
ANTON'S DIRECTION: a proven geometric filter that determines (a) when the fast computation safely selects the right case, (b) when higher precision is needed, (c) when the material still leaves several possibilities. An arithmetic proof alone does not resolve wrong side conditions — anatomically motivated side conditions are also needed.
MATERIAL: results/X1b (knee with moving axis + wrapping; moment arm jumps when the path changed sides; forced side as post hoc), results/W1 (geodesic wrapping on cylinder, fixed side per muscle), results/K1 (straight path vs wrapping 2,2 %).
ASSIGNMENT:
1. PREREG (hashed). Formalize case selection: straight line / wrapping right / wrapping left as a function of joint angle and geometry; derive the case boundaries (tangent conditions) and an interval filter that classifies each (angle, geometry uncertainty) as SAFE case X / NEEDS PRECISION / AMBIGUOUS (several cases consistent with the material — e.g. landmark/atlas errors from H2b/AT1).
2. Side conditions: formulate anatomically motivated conditions (which side the muscle can pass on, source) and show that they remove X1b's jumps; distinguish this from numerical correctness.
3. Test on X1b's muscles over angle sweeps: fraction safe/precision/ambiguous; moment arm continuity; comparison against X1b raw and forced side; cost (how often expensive computation is needed). Countertest: the filter should flag ambiguous when geometry uncertainty is increased to AT1's origin error (12–19 mm).
Output only results/BT-LV2/. Finish with RESULTS.md.
