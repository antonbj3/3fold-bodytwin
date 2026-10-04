# Breakthrough hunting 2026-10-03 — characterization of sufficiency

Status: **PENDING_INDEPENDENT_REVIEW**. Everything below is calculated from files
on disk (`results/LANE_SUFFICIENCY_SWEEP/PORT_R1_V1.json`, `PORT_R9_V1.json` , `RESULTS.md`
). No new physics claim; the statement is about the SHAPE OF THE GATE.

## 0. What was already in the files

`PORT_R1_V1.json` carries per contract its own field `summary_form` as
well as, for each fold, `extension.factorization` (the exact source equation)
and for each hold `hold_argument` . So I haven't had to code the shape
myself to get an independent label — the run already has one.

That label separates the output **exactly**:

| `summary_form` starts with `local_` | FAILED | HELD |
|---|---:|---:|
| yes (`local_constitutive_scalar`, `local_geometric_scalar`) | 0 | 6 |
| no (38 other form names) | 38 | 0 |

44 testable contracts (BIORESP = NOT_TESTABLE). Zero misclassifications.
This is lane driving's own label, not my post-coding.

## 1. Table: case, form of reading, expansion that was enough, failed/held

`nl` = the readout is non-linear in the state/distribution (coded out of
`extension.factorization` , i.e. out of the actual source equation). `koloc`
= the readout's representative is collocated with the summary's (same
point, same moment, same weight, no other independent coordinate).

| Cell | Contract | Form (of execution) | nl | koloc | Failure mechanism | Extension | Outcome |
|---|---|---|---:|---:|---|---|---|
| Q021 | o2_content_to_partial_pressure | local_constitutive_scalar | 1 | 1 | — | none | HELD |
| Q026 | local_ligand_to_occupancy | local_constitutive_scalar | 1 | 1 | — | none | HELD |
| Q052 | porosity_to_effective_diffusion | local_constitutive_scalar | 1 | 1 | — | none | KEPT |
| Q080 | receptor_total_to_total_equilibrium_binding | local_constitutive_scalar | 1 | 1 | — | none | KEPT |
| Q146 | plasma_amount_to_elimination_rate | local_constitutive_scalar | 0 | 1 | — | none | KEPT |
| SURG_INCISION | indentation_to_contact_area | local_geometric_scalar | 1 | 1 | — | none | KEPT |
| Q005 | renal_clearance_to_cell_retention | flux_or_clearance_to_inventory | 1 | 0 | companion factor | 1 (P_MATE) | FAILED |
| Q009 | total_tissue_to_future_free | total_pool_to_partition | **0** | 0 | linear division +  time | 1 (bound) | FAILED |
| Q012 | steady_concentration_to_half_life | steady_value_to_transient | 0 | 0 | time | 1 (CEO) | FAILED |
| Q013 | exposure_integral_to_endpoint_effect | time_integral_to_local_response | 1 | 0 | time weight | 1 (endpoint conc.) | FAILED |
| Q014 | forcing_integral_to_memory | time_integral_to_memory | **0** | 0 | time weight flat→exp | 1 (exponentially weighted drive) | FAILED |
| Q017 | volume_to_filtration | inventory_to_local_flux | 1 | 0 | locator + companion factor | 2 | FAILED |
| Q019 | total_na_to_terminal_na_flux | total_pool_to_local_flux | **0** | 0 | locator | 1 | FAILED |
| Q020 | interstitial_volume_to_osmotic_flux | inventory_to_composition | 1 | 0 | companion factor ×2 | 2 | FAILED |
| Q022 | pressure_to_muscle_perfusion | pressure_to_flow_with_geometry | 1 | 0 | companion factor | 1 (strain) | FAILED |
| Q031 | ion_amount_to_concentration | inventory_to_concentration | 1 | 0 | companion factor | 1 (V) | FAILED |
| Q036 | flow_to_reperfusion_derivative | current_state_to_memory_derivative | 1 | 0 | time | 1 (injury_block) | FAILED |
| Q043 | pooled_spikes_to_force | total_events_to_weighted_response | **0** | 0 | weight (pool→unit) ⟧|⟧ no small found | FAILED |
| Q044 | dipole_strength_to_electrode_potential | magnitude_to_directional_response | **0** | 0 | non-injective reduction | 2 (angles) | FAILED |
| Q049 | fluid_content_to_boundary_drainage | total_inventory_to_boundary_flux | **0** | 0 | locator bulk→stripe | 1 (surface print) | FAILED |
| Q052 | porosity_to_permeability | volume_fraction_to_geometry_response | 1 | 0 | factor | 1 (pore diameter) | FAILED |
| Q054 | prescribed_load_to_contact_peak | total_force_to_local_pressure | 1 | 0 | companion factor | 1 (contact area) | FAILED |
| Q058 | capacity_to_capacity_derivative | current_state_to_memory_derivative | 1 | 0 | time | 1 (fast fatigue) | FAILED |
| Q077 | area_volume_to_transport_time | bulk_geometry_to_topology | 1 | 0 | locator (topology) | 1 (graph distance) | FAILED |
| Q080 | receptor_total_to_peak_activation | total_population_to_spatial_peak | 1 | 0 | weight mean→sup | 1 (peak density) | FAILED |
| Q084 | mean_diffusivity_to_passage_time | mean_transport_coefficient_to_boundary_response | 1 | 0 | weight arithm.→resistance | 1 functional | FAILED |
| Q088 | atp_to_functional_score | current_energy_to_other_inventory | **0** | 0 | companion factor | 1 (PCr/ion) | FAILED |
| Q090 | tooth_displacement_to_force | current_displacement_to_viscous_force | 1 | 0 | time (speed) | 1 (speed) | FAILED |
| Q100 | glottal_area_to_wall_pressure | geometry_to_lagged_flow_pressure | 1 | 0 | companion factor (lag) | 1 (flow) | FAILED |
| Q107 | gastric_volume_to_pylorus_flow | total_inventory_to_local_flux | 1 | 0 | locator | 1 (atrial volume) | FAILED |
| Q115 | cell_count_to_barrier_conductance | total_population_to_age_weighted_response | 1 | 0 | weight (exp age) | 2 (M1, M2) | FAILED |
| Q121 | total_gas_to_proximal_pressure | total_inventory_to_local_pressure | 1 | 0 | locator | 1 (proximal gas) | FAILED |
| Q127 | myonuclei_count_to_density | inventory_to_concentration | 1 | 0 | companion factor | 1 (fiber area) | FAILED |
| Q140 | total_drug_to_mate_secretion | total_pool_to_local_flux | 1 | 0 | locator (4 segment) | no small found | FAILED |
| Q154 | total_substrate_to_brain_influx | total_pool_to_partition | 1 | 0 | locator + companion factor | 2 | FAILED |
| Q156 | isf_mean_to_local_diffusive_flux | mean_field_to_local_gradient | **0** | 0 | weighted mean→gradient | 1 (local diff.) | FAILED |
| Q160 | porosity_to_biofilm_permeability | volume_fraction_to_geometry_response | 1 | 0 | companion factor | 1 (EPS fraction) | FAILED |
| Q168 | plasma_tmao_to_next_derivative | current_state_to_memory_derivative | 1 | 0 | time | 1 (portal TMA) | FAILED |
| IMMUNITY | C3b_to_complement_derivative | current_state_to_regulatory_response | 1 | 0 | time | 1 (FactorH) | FAILED |
| MITOSTRESS | area_volume_to_transport_time | bulk_geometry_to_topology | 1 | 0 | locator | 1 (graph distance) | FAILED |
| SOLBENCH | regional_means_to_mean_disposal | mean_fields_to_nonlinear_sink | 1 | 0 | weight (mean of nonlinear sink) | no small found | FAILED |
| SURG_COLLAGEN | denaturation_fraction_to_cutting_ratio | fraction_to_hidden_material_response | 1 | 0 | companion factor | 1 (cross link) | FAILED |
| SURG_HEALING | flow_to_reperfusion_derivative | current_state_to_memory_derivative | 1 | 0 | time | 1 (injury_block) | FAILED |
| SURG_HEMOSTASIS | raw_luminal_flow_to_wound_flow | total_flow_to_hidden_network_gain | 1 | 0 | locator | 1 (upstream scale) | FAILED |
| BIORESP | — | UNCONSTRUCTED | — | — | — | — | NOT TESTABLE |

Expansion size among the 38 precipitates: **1 scalar in 29, 2 scalars
in 5, no small found in 4**. No case required more than two.

## 2. The hypothesis fails — plainly

The hypothesis was: *sufficiency fails exactly when the reading is a non-linear functional
of a distribution whose relevant moments are more than the rank of the summary.*

It falls on **both** ranks, and the counterexamples are in the table above:

- **Nonlinearity is not sufficient for precipitation.** Five out of six directions
  are readouts that are strictly non-linear in state: Q026 occupancy `C/(K_D+C)`
  , Q021 oxygen curve inversion, Q052 `D_eff = D0·phi/[1+b(1−phi)]²`
  , SURG_INCISION `A = pi·a0²/(1−delta/h0)` . They still hold.
- **Nonlinearity is not necessary for precipitation.** Nine precipitations
  have a readout that is **linear** in the state: most clearly
  Q156 `J = −A·eps·D·(c1−c0)/dx` (linear in the concentration
  field, falls on 192 pairs), Q009 `free = total − bound` , Q049
  `Jtop = −Gtop·P0` , Q019, Q043, Q044, Q014, Q012, Q088.
- **The moment language describes only one case.** Out of 38 precipitates,
  only Q115 is actually repaired by moments of a distribution (M1, M2 of
  `exp(−age/tau)` ). In 21 cases, the missing coordinate is a **parameter
  or other simultaneous quantity** (P_MATE, volume, contact area, fiber
  area, cross-link density, orientation angles), not a moment. In Q080
  the erring quantity is a **supremum** — not a moment of any order.

As a classifier: non-linearity gives 30/44 correct (14 wrong: 5 non-linear holds,
9 linear folds). It is not only blurry, it is almost independent of the outcome.

## 3. The sharper characterization

What differs is not the curvature of the functional
but **where its representative sits**:

> **The collocation criterion.** Write the readout `Y` as a functional of
> the state and indicate its representative: which test function (location, moment,
> weight) it is integrated against, and which quantities multiply it.
> The summary `S` is **sufficient** if and only if `Y = g(S)` for some
> measurable `g` on the entire allowed fiber — which occurs just when `Y`:s
> representative lies in `S`:s span: **same place, same moment, same weight,
> and no other independent coordinate in the product**. `g` may be as nonlinear
> as you like.
>
> Precipitation occurs at **representative mismatch**, and it has exactly three types, which
> all are visible directly in the form of the readout without any counterexample being constructed:
> 1. **Locator Mismatch** — `S` integrates over a volume/population, `Y` is evaluated
> in a point, on a line, in a segment, on a graph distance or in a
>    supremum. (12 fall)
> 2. **Time mismatch** — `S` is now, `Y` is the derivative, the transient, the
> exponentially weighted memory or future. (9 case)
> 3. **Companion factor mismatch** — `Y` is a product/ratio where `S` is a factor and the
> the second factor is an independent coordinate at the same place and time
> (volume, area, permeability, velocity, flow, angle). (17 case)

**Why the expansions are small is now a consequence, not a coincidence.** A mismatch
is a *direction*: a test function that is missing from the span. Therefore, repairing
it requires one (1) new functional per mismatch type. The five two-scalar cases are
exactly those with two simultaneous mismatches (Q017 locator+cofactor, Q154 ditto,
Q020 two cofactors, Q044 two angles, Q115 two moment orders). It predicts the **size**
of the expansion, not just its existence: *number of scalars = number of independent
representative mismatches*. It is true in 34/34 cases where minimality is proven;
in 4 cases (Q043, Q140, Q084, SOLBENCH) the mismatch is an entire weight functional
and no finite scalar expansion was found — which the criterion also predicts,
since a weight over a continuum is not a direction but a space.

**Applied to a NEW case without a counterexample:** write down the
equation of the readout; mark its place, moment and weight function;
highlight the summary; count the mismatch. Zero mismatch ⇒ sufficient.
k mismatch ⇒ folds, and the minimum expansion is k scalars (or
a functional if the mismatch is a weight over a continuum).

**Warning against superstition:** my `koloc` coding was done after the states
were read, so 44/44 is alignment in the sample. The equally informed check
is the run's OWN `summary_form` label, set before my analysis, which separates
identically (0/38 vs 6/6). The next episode is the only honest test.

## 4. Three predictions written BEFORE the outcome

The sample set is the random selection in `PORT_R9_V1.json` (n=30 out of frame 260),
which is a different set than the 45 adversarial contracts and was not used to formulate
the criterion. I read `cell` , `quantity_key` , `question` and `locators` but **not**
`status` or `R9_adjudication` . I selected the three lines whose `plan` text does
not leak the value (many other lines say outright "HOLD ..." or "Repair ...").

**Prediction A — draw 20, Q052 `outlet_flux_normalized` → `concentration_mid`
.** The summary is a fringe flux at the outlet; the readout is the concentration
in an inner center point. Locator mismatch (stripe → inner dot), a mismatch.
*Predicted: FAIL, minimum expansion 1 scalar.* Risk: the line says that
the fiber should be constructed "if constructible", so an undecided `UNKNOWN_SEARCH`
is a possible outcome without the criterion being tried.

**Prediction B — draw 27, Q020 `K_f_ss` → `J_net` .** `J_net = K_f · (net driving pressure)`
. The summary is one factor, the other is an independent
simultaneous coordinate at the same place and time. Cofactor
mismatch, a mismatch. *Predicted: FAIL, minimum expansion
1 scalar (the effective driving pressure).*

**Prediction C — draw 28, Q021 `middle_vq_mean` → global oxygen readout.**
The summary is a regional average over the middle region; the readout is
a global aggregation across all regions through a curved content curve. Weight
mismatch (regional weight → global weight). *Predicted: FAIL.* But this
is the case where I'm least sure: Q021 held in the adversarial set, the
full regional simulation couples the flows, and if the allowed fiber is
limited to source-generated stationary solutions, the intermediate mean
can lock the rest — then the outcome is FAIL or undecided, and the criterion
needs the "the fiber must be free in the mismatch direction" condition.

### The outcome (read after the above was written)

| # | Case | Predicted | Actual | Expansion Predicted | Expansion Actual | Verdict |
|---|---|---|---|---|---|---|
| A | draw 20, Q052 outlet flow → center concentration | FAILED | **FAIL** (19/19 counterex., identity error 0,0, downstream 14,617 %) | 1 | **2** (`process_state`, `permeability_override_m2`) | binary right, **size wrong** |
| B | draw 27, Q020 `K_f_ss` → `J_net` | FAILED, 1 scalar = driving pressure | **FAIL** (10/10, identity error 0,0, downstream 93,75 %) | 1 | **1** (`effective_pressure_drive_mmHg`) | right, incl. name of greatness |
| C | draw 28, Q021 regional V/Q-mean → global oxygen reading | FAILED (least secure) | **FAIL** (15/15, identity error 0,0, downstream 27,75 %) | not specified | 1 (`o2_uptake_ml_min`) | binary right, **mechanism type wrong** |

**3/3 binary right. 1/2 size right. 2/3 mechanism type correct.** Two genuine flaws:

- **A:** I counted one mismatch (locator: outlet → inner point)
  but the fiber carried two, because the summation *itself* is a
  compound flux (conductance × driving force). Procedure count
  rule correction: when the summary is a flow or a quota, it always
  carries a cofactor mismatch in addition to its locator mismatch.
  It was not written in the criteria before the test.
- **C:** I predicted precipitation via *weight mismatch* (regional → global weight).
  The reason for the line itself is different: the normalization removes the
  absolute scale, so what is missing is a cofactor (`o2_uptake_ml_min`), not
  a weight. The right judgment for the wrong reason — the criterion's type classification
  is thus not self-evident to apply to normalized quantities.

### Hold in the omitted quantity — the symmetric test

The four sides of the random selection (read after the predictions) all have `koloc = 1`:

| Case | Readout | Why collocated |
|---|---|---|
| draw 4 + 29 | Q146 AUC → `F_abs_total` | `F = AUC·CL/dos`, fixed conversion: pointwise function of the summary itself, same window, same compartment |
| draw 5 | Q058 `max_error` → `error_pass` | the predicate is `1[S ≤ f_max]` — the threshold acts **on the summary**, not on the field |
| draw 8 | Q168 `precursor_um` → the precursor derivative | the right term is the end in the precursor at fixed kinetics |

These two lines are more important than the precipitations, because they kill two tempting
coarser rules that the original series of measurements could have given rise to:

- **"Threshold always fails" is false.** draw 5 is a threshold and it holds.
  The CEM43 case fell not because it was a threshold but because the threshold
  is evaluated on the history of the underlying field, not on the summary.
- **"A derivative always fails" is false.** draw 8 is a derivative and it holds, while
  Q168 `plasma_tmao_to_next_derivative` fell. The difference is exactly the span
  condition: is the right-hand side of the ODE closed in the summation or not.

The same sharpening applies to the eleven main cases: `delivered laser energy → removed
mass` (the ablation threshold acts on the local fluence distribution, not on the
total energy), `surface temperature history → temperature at the target depth` (locator),
`physical count without ordinal numbers → 192/192` (weight: disordered count → ordered), and
`surface allocation → 0,8724 D` (co-factor). Zero of the eleven have a co-located
representative — the criterion predicts all eleven the refutations, but it is
post notice, not a test.

## 5. What is actually new, and against what control

The lack of summaries is not new — it was already measured in 45 + 20 cases. What is new
is that **the precipitation can be calculated from the form of the readout before the
test**: the three blind cases were decided correctly without any counterexample being
constructed, and in one case the name of the erring quantity was also predicted.

The equally informed check: the run's own `summary_form` field, set by the lane before
my analysis, separates the 44 adversarial cases identically (`local_*` ⇔ HELD, 0 errors).
**So the criterion doesn't beat the lane's own label on that quantity.** What's new
is three things that label doesn't provide: (i) a rule that can be applied to a
row that *does* have no shape label — which is exactly what random selection in
R9 is, (ii) a prediction of the **size** of the expansion from the number of mismatches,
which is right in 34/34 adversarial cases with proven minimality but **wrong
in 1 out of 2** blind tests, and (iii) an explanation of why the expansions
are small that is not statistical: a mismatch is a direction, not a space.

The criterion says nothing about biological validity, and nothing about
prevalence: 0,882 is adversarial sampling, and random sampling is
still not certified below half (95 % interval [0,485; 0,958]).

## 6. Measured vs. estimated — the generator and the number

Calculated from `extension.coordinates` in the collapsed contracts, ie.
the quantities the sufficiency test names as decision-making and which the
pipeline today delivers as a frozen value instead of a measurement:

- **39** named quantities in the adversarial selection refutations
- **22** in the random selection refutation
- **61** distinct (cell, quantity) pairs in the Union, **59**
  distinct quantity names, spread across **37** of the 45 cells
- hence **14** already has an orderable measurement spec (`orderable: true` in
  `ACQUISITION_TARGETS_V9.json`, 20 nodes total)
- **4** precipitates have no scalar candidate at all — the gap is an
  entire weight functional (Q043, Q140, Q084, SOLBENCH) and is therefore
  not a measured-versus-estimated decision but a model fault

**The answer to the question: 61 candidates**, of which 14 already orderable.

The limit of the claim: that *clinical practice* uses an estimate right there is
verified for **one** of the 61 (posterior cornea, 0,282 vs. 0,626 D in 69 patients).
For the other 60 it is verified that **our** pipeline uses a frozen value and that
the quantity is decision-making; if the practice estimates or measures it is
UNKNOWN per row. The list is thus a generator of *questions for practice*, and
it is computable — it falls out of the sufficiency test without extra work.

## 7. Falsifier

The criterion is passed if any of the following occurs:
1. a read with zero representative mismatch falls on a pair with identity errors
   on machine precision;
2. a readout with k ≥ 1 mismatch holds on a free fiber (the fiber must be free
   in the mismatch direction — the condition C showed is needed);
3. a minimum extension that is strictly greater than the number of independent mismatches after
   that the calculation rule in A has been applied.

Point 3 is already close: test A required a calculation rule that was
not written before the test. The size prediction should be treated
as **untested**, the binary prediction as **3/3 on blind amount**.
