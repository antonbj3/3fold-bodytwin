### Q168 — Measurement choices that distinguish microbial, barrier and host causes

**Question:** Which combination of local and systemic measurements can distinguish microbial production, barrier passage and host metabolism when the same final curve can be explained in several ways?

**Resolution target:** Measurement operator at niche, barrier and organ port level with relevant temporal resolution.

**Physical basis:** Causal mass-balance chains and explicit observation operators; identifiability is tested before a mechanistic conclusion.

**Mechanism:** Formulate competing causal models and measurable interventions/negative controls; this node designs future comparisons and reports no new biological effect.

**Mathematical starting point:** `Δ_A=E[Y|do(A=a_1)]−E[Y|do(A=a_0)] for a predefined perturbation A and outcome Y; the estimand has Y's unit. Identification requires specified confounders, measurement errors and transport between populations.`

**Detailed reference to develop:** Synthetic observations from a coupled full-order model and independent longitudinal, paired samples; actual interventions are later research decisions.

**Laws/parameters to determine:** Measurement sensitivity, confounders, biological variation and correlated errors.

**Outputs:** Rank in local sensitivity matrix; Prediction error in outcome units; Information gain and measurement cost

**Reuse:** No direct code anchor mapped in this limited review.

**Next work/data:** Predefine at least two competing mechanisms, measurement sites and a practically measurable perturbation.

**Comparison:** Previous graph context reported that predicted microbial metabolites did not improve host outcomes beyond raw taxonomy; original numerics not re-reviewed. Use the same samples/cohorts for taxonomy, metabolite and barrier measurement and held-out cohorts.

**Test that can reject the proposal:** Selected observations do not distinguish models on held-out cohorts, do not improve host outcomes beyond the taxonomy baseline, or the causal estimand lacks identification conditions.

**Cross-links to test:** Q060, Q148, Q152, Q165, Q167. Mechanistic development proposal without a specific external source in this review.

