# Next round — continue from R3

Everything PENDING_INDEPENDENT_REVIEW. The parent target is OPEN; innovation gate FAIL, strongest equally informed control TIE. No new FE/LP is needed before a real work port exists.

## Read first

`CHECKPOINT_R3.json`, the latest section in `RESULTS.md`, `SOURCE_REVIEW_R3.md`, `SOURCE_TARGETS_R3_v2.json`, `r3/metrology_results.json`, `VERIFICATION_R3.json`, `NEXT_WORK_PORT_R3.json` and frozen `r3/BINDINGS_PORTS_R3_SNAPSHOT.json`. Read current steering and other lanes' latest ports; their empirical h_pz and skin intercept were null in R3. R1/R2's executable layer/history operators and all misses remain.

## Next key operation

Change the observation operation from normal load/damage depth to **the work's two force components and measured contact/crack kinematics**. The published NBS series is actual skin cutting with varied geometry, but Fn has no effect in direct work balance when vn=0. Acquire an original that has simultaneous Ft, Fn and displacement, real new crack area and actual local radius. Prioritise the Doran original through a legitimate author/institution host, later instrumented skin-cutting studies or available data files; never use secondarily cited forces as primary curves. DOI/full text remaining in `SOURCE_MANIFEST_R3.json` and saved errors are starting points, not evidence.

1. Freeze the source's numerical data and their joint uncertainty before a new prediction. Fill `NEXT_WORK_PORT_R3.json` in a **new** file. Leave missing fields null. Matched material, body site/prestretch/direction, first cut versus later passage and speed must be explicitly stated. Distinguish a real cut from a crushed track. Residual damage depth is not automatically loaded penetration or new crack area.
2. Integrate `Wext=∫Ft dx+∫Fn dz`. Keep a clear interval for `Geff=(Wext−ΔPsi−Dcontact)/ΔA`. Use a shared reference/material uncertainty for all geometries, not a separately fitted Γ for each observable. A full equally informed cohesive contact model is the control; count its acquisitions, solutions and fallback.
3. If chemical Γ0 is requested, independent full first-pass pre-cleavage work is needed: `Γ0=Geff−B_before`. R2's Γ0+B0 null direction does not disappear with more radii without a mechanism that identifies B0. Synthetic 150 or the scissors value 2500 J/m² must not become measured Γ0.
4. Freeze one training geometry and another radius held out from all parameter fits. Predict measured work/force, then a new gap through the preserved layer operator if the ports match. Run the actual control only when the data makes this meaningful. No history-free/10× gain without full cost and physical error budget.

## If only NBS data remains available

Preserve R3's negative results; do not make a new fit and call it prediction. A->other geometries misses by up to91,12%; B->K width model gives at least87% depth against40±2%. B+F12lb gives K64% under the RF25,4µm scenario, with shared rounding sensitivity62,5–65,5%. Actual RF is indeterminate. The23,787µm inverted from K is diagnostic, not a new test hit. Necessary local RF error is approximately0,503µm for an8pp depth budget **if** this closure and RK were exact and other errors zero; no achieved measurement precision is claimed.

Furthermore, no constant F offset can match12/16/20lb with the declared resolution scenarios:8,000–8,091 /9,143–9,524 /14,857–15,238lb. This is a failed closure test that is not solved by changing RF. Moving forward requires variation/contact kinematics along the actual cutting sites and synchronous force–motion; a load-dependent fitted offset does not make the model predictive.

## The BINDINGS coupling

Request deltaV/DeltaA for an actually advancing fracture front with full first-pass energy, anchor detachment/fibre rupture/delamination and the wedge's imposed contact work. Microcell contour length is not measured process-zone depth, partial hysteresis is not full fracture energy. The same radius can give different work. Do not transfer a porcine microcell scenario to human autopsy skin as a shared posterior.

## Preservation and running

All R3 results are written exclusively. `geometry_r3.py --out <ny path i lane>` can be reproduced to a new directory; metrology_r3.py and verify_r3.py refuse existing default files. Reproduction requires a new output destination in a copy, no overwriting. Run with OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2. `finalize_r3.py` must not be rerun: it rejects existing artifacts/sections.

Change only the lane's own results or an allowed shared-data directory. No subagents, graph mutations, cloud or product code. Local GRAPH_FEEDBACK_R3.json awaits the coordinator's review.
