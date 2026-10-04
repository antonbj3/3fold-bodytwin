# Action queue

The night log is a journal: a row there proves that something was found, not that something was done. That led
to findings being reported to the operator instead of worked on, which he pointed out on 2026-10-03.
This file is the queue. A row moves to DONE with the number that shows it is done, never on a
claim that it is done.

## OPEN

| # | what | evidence needed to close | owner |
|---|---|---|---|
| A3 | 162 consumable observables (1 consumed: the colon's gas regime, see D9) from 55 papers, 118 outside the eye, in `LANE_READ_CREATIVE/CONSUMABLE.json` | through the unit and range gate, then admitted edges with numbers in the text | me |
| A4 | 13 orderable measurement targets and 7 dataset routes from the same round | those that can be retrieved without ordering anything: retrieved and consumed | the swarm |
| A8 | CLOSED 2026-10-04: 30 STALE were 5 distinct (same alert per hourly run) and all five were LOCATOR drift, not poor evidence; 24 unchecked carried a third evidence form (`fil :: /json/pekare = text`) that the check did not read | the net: OK 149, UNRESOLVED_DECLARED 6, STALE 0, UNCHECKED 0; `relocate_stale_evidence.py` moved 4 line quotes + 1 pointer + 15 to the prose half of the same document, 2 downgraded with reasons; negative control: false number in E0078 still flagged STALE | done |
| A6 | the detail layer: harvest DONE — 312 inventoried, 99 edge proposals, 3 082 external references, 244 excluded records removed | edge proposals through the unit and range gate, then admitted; reference list crossed against our quantities | me |
| A7 | PARTLY CLOSED 2026-10-04: `np.trapz` disappeared in numpy 2.4.1 and hit EIGHT active source cells, not just native_glucose | six direct calls replaced with `np.trapezoid` in Q019, Q043 (2), Q140, Q160, native_glucose; Q013 line 159 bound the name eagerly; Q009 line 158 had `getattr(np,"trapezoid", np.trapz)` where the DEFAULT is evaluated immediately and thus threw exactly the error getattr was meant to avoid. All seven now import from their own directory and native_glucose runs the whole Dalla Man model. opensim is not on PyPI (`No matching distribution`), so myofascial_transmission cannot run here — but the DATA exists in the readable repo and needs no opensim, see A14 | partly done |

## KLART

| # | what | the number that closed it |
|---|---|---|
| D1 | `from_vector` rotated every axis 90° | round trip exact across 36 axes; toric 0,2817 → 0,2808 |
| D2 | the net's tracked copy 34 edges behind | 56 of 56 in both copies, sync in the hourly cycle |
| D3 | the queue 92 % completed work | 13 919 → 1 047 rows, pruning in the cycle |
| D4 | Sol lanes without web search | `tools.web_search=true` in all three drivers |
| D10 | 99 edge proposals admitted | the net 56 → 155 edges, external references 15 → 27, gate 99 of 99, one verified line by line |
| D9 | first harvested observable consumed into a decision | gas-regime decision outside the eye: threshold 76 mL/6h separates two regimes, single-coefficient control 77,8 % high in every row |
| D8 | three edges reported as UNCHECKED | all three already carried `evidence_unresolved` with a note; the error was in my check that did not read the field. Now **UNRESOLVED_DECLARED** (searched for, does not exist) is distinguished from **UNCHECKED** (unchecked), and UNCHECKED is zero |
| D7 | the laser refutation was mine, and it was wrong | quantities, denominators and observation operators differ; and I wrote the unit **µm⁻¹** where the source table says **µM**, thus read a concentration as an inverse length. The spread 2033× is real and does not affect the witness. Edge back to UNKNOWN, original annotation preserved verbatim |
| D6 | Q005 stored result 1000× wrong in both concentrations | the code's own run gives Met 0,20734 and Cr 0,31685 against stored 207,342 and 316,849; CL_Met identical 29,5852, thus only the concentration conversion. Old preserved as `.stale_20261003_factor1000`. Assay constant untouched. **The cell's own `unit_check` still passes** — the gate does not check what was wrong |
| D5 | memory reservation double the consumption | 1500 → 1100 for our rows, 1500 for dental's after their p95 1319 |
| A9 | Q052's transport scale: velocity 88,9–629,8 mm/s comes from `∇p=2,0e4 Pa/m` which the file itself tags "assumption", and `k_geometry` from a calibration the file itself says is not independent validation | scale rederived from the same lowest level as permeability, or a line saying the model is at scaffold rather than tissue scale; retuning ∇p to hit a band is forbidden (their PREREG says so, and I agree) | open |
| A10 | The "flagged in the source" gate lets through 26 of 28 deviations >2 % in Kim 2008 table 2, largest unflagged 20 % | every consumed mean-value row gets its own deviation check against its nominal value, independent of the source's flag; no cell may consume the row without it | open |
| A11 | `sufficient_age_upper_s` in EXTERNAL_SOLVER r25 is a LOWER bound (pass for age ≥ 1088,1966 s), and the bracket width 0,00576 s is bisection resolution 1509,7/2¹⁸, not measured sharpness | field renamed to `sufficient_age_lower_s` with direction as its own field, and the bracket's provenance written into the artifact | EXTERNAL_SOLVER r27 |
| A12 | The measured serum creatinine increase in PMID 32989831 is not in CONSUMABLE.json — it is the independent reference that would test the prediction that creatinine secretion is blocked 61–66 % | number retrieved from the paper and compared with what 61–66 % block implies for serum creatinine | the swarm |
| A13 | CLOSED 2026-10-04: CORNEA_SHAPE r24 recalculated `same89` in the workspace (`same89_recomputed_predictions = True`) and the fields carry predictor in the name | the two chains now meet within 1,24 % in magnitude and 1,80 % in vector (0,2808/0,4985 against 0,2843/0,4897 on the same 69 eyes); our chain reports both metrics because the arm order reverses between them | done |
| A14 | `data/msk_smoketest/` in the readable repo carries static optimization for stoop lifting WITH and WITHOUT a box, readable as text without opensim | pair consumed: `tasks/assembly/load_redistribution_decision.py` shows that uniform load scaling explains 70,2 % of 2 523,9 N redistribution, thus 30 % structural; remaining: more poses, because the file carries ONE pose per arm | open |
| A15 | The population steadfastness orientation in `toric_decision.py` (−0,30 D at 90°) is a constant without a citation, and the gain is 0,3401 D at 90° against 0,0828 D at 0° on the same 69 eyes | orientation bound to a cited nomogram (Koch/Baylor or equivalent) with DOI, or the headline permanently quoted as an interval over both conventions | open, blocks every citation of 0,62 → 0,28 D |
| A16 | The toric margin's certification rests on rho = −0,796: cancellation removes 51,8 % of the error, and full repair of one stage makes the system 1,65× and 1,26× WORSE respectively | margin never cited as a bound that survives improvement; gate renamed or made two-sided, because "upper CI below 0,15" lets strong anticorrelation through under the name decorrelation | open, blocks every citation of the toric system margin |
| A17 | STANDING GATE: every favourable aggregate or margin that may depend on cancellation between components must carry a repair response — repair one component and report whether the aggregate grows | three cases tonight (toric margin ρ = −0,796 becomes 1,65× worse; swarm job 1,90–2,74× underestimation; LASER_SURGERY 93,51 % strict underestimation) make it a rule and not a one-off check | standing |
| A18 | A17's standing gate partly measures the decomposition, not the system: the proof lane showed that the same total error can give different repair effects by only changing the stored intermediate predictor, and that no finite observation universally separates physical compensation from a constructed decomposition | gate reformulated as WARNING and never as proof; the conditional test in `results/PROOF_LANE_CANCELLATION/` runs before any margin is called physical | open |

