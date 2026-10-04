# Next load-bearing operation after round 1

Status: PENDING_CONSTRUCTION, PENDING_INDEPENDENT_REVIEW. Same parents GRAPH_INVERSE/MULTIPHYSICS, BT-CTX-SURG-INCISION, Q033/Q036/Q049. The round's broad capability/cost gate FAIL; all equally informed method controls TIE. Retain all raw negative results. Read CHECKPOINT_R1.json, RESULTS.md, DECOMPOSITION_R1.json, COST_ACCOUNTING_R1.json and r1/centered_residual/summary.json first.

Capability → a new incision tip and supplied fracture history in nonlinear surroundings without global state, full solver or new global training acquisition per query.
Conflict → local topology change, non-local compliance, cohesive memory and max-gap budget0,1µm.
Obstacle → the new tip positions' cusp lies outside the17D global basis; even κ±5% gives0,298–0,649µm error, and the residual bound is too wide.
Changed operation → compiled support family with shared deformations plus local Green/tip responses, and a signed history-residual sector.

## Start immediately

1. Read PREREG_R1_CENTERED_RESIDUAL.json. `compiled_query_r1.py` imports only NumPy/stdlib and reads `r1/centered_residual/operator_tier2.npz`. `--length-mm 9` is conditionally accepted; `--length-mm 8` requires fallback. There is no solver/globalstate dependence in this warm consumer. Full fallback has actually been run in `centered_residual_r1.py`, not merely proposed.
2. Use exactly the same nx64×ny32, β20, strain0,10 and synthetic materials as preferably the first instrument. Retain the missed lengths8/16/24/32mm and κ±5% at lengths1/9/17/25mm. Frozen gates: gap≤1e−7m, cert≤1e−7m, no warm global solve or state recreation. Full evolving cohesive history and biological accuracy remain UNKNOWN until separate sectors close.
3. Acquire, **without validation ground truth**, nonlinear intact prestretch and the tangent's seam-Green columns `q_j=Jref^{-1} B e_j`. They have been acquired temporarily for linear compliance/cert purchases but the global matrix is not stored in the portable package; new explicit acquisition is required and counted. Save it in this lane/shared-data, never in private sources.
4. Build the trial family as `U=Uref+Pshared z+q_tip eta` or two/four supported tip/wake columns. The first Pshared may reuse the five original training supports0/4/12/20/28mm, with documented acquisitions. The new q column must not be a full solution at a missed heldout. This tests whether a physically local deformation of the basis replaces acquisition of another entire state direction.
5. Compile the residual at every possible mesh tip, preserve `r=Cv+Bs` and `2sᵀDv`. No new exterior phase may be disguised as coefficient reading. A small d per support can avoid a single full m-port basis, whose cubic monomials and Gram grow O(d³)/O(d⁶). Measure actual package size/setup; do not call a linear,40mm wake constant work.
6. Check polynomial identity against the full RHS as **validation**, but use only the compiled package in query. Start with lengths8 and32mm. For a larger remainder: distinguish residual at the tip from distributed nonlinear strain/older wake. Add only a physically supported direction when the residual shows which missing field it represents; execute the change in the same round. A blindly larger global ROM is not the next main operation.
7. For supplied κ, the full bulk is convex and the residual theorem applies. Letting κ evolve requires explicit contraction/stable branch or a branch/variation enclosure. A negative cohesive tangent during softening must not be slipped into the same monotonicity theorem. Test supplied κ first, then actual history with recorded maximum opening and dissipation; no new cell-injury claims.
8. Give conventional local enrichment/Green/residual-Gram control the same acquisitions, support and exact algebra. Also measure a reusing sparse full Newton with previous solution/tangent. The difference between query with acceptance and query with fallback must be visible; all renewals and failures count. Report TIE if the method is the same. A new support contract that works is useful capability but not automatically performance novelty.

## Physical port that needs acquisition

The Goda gel test failed for a single tear/cutting energy but the split into Γ0 and Wf gave10,03% slope error. It is a known published mechanism, not skin calibration. Skin cross-prediction requires measured cleavage intercept/friction, cohesive peak or failure-strain, prestretch and an actual gap/strain map on the same sample. The molecular modules' Tm/moduli do not identify these ports. If open matched primary data can be found, freeze them before new fitting. The in-vivo/ex-vivo prestretch reference in WORLD_REFERENCES_R1.json is a separate observation and has not been used for fitting in round1.

Bleeding is UNKNOWN: an oriented vessel map, radii, upstream pressure/conductance and tissue-coupled hemostasis are needed. OMAG area fractions or Ivy normal intervals are insufficient. Q033 exudate remains fluid flow, not blood. Q036 coupling is synthetic and one-way; the model has no edema/perfusion→mechanics feedback yet.

## Reproduction and resources

COMMANDS_R1.md has actual commands and preserved versions. Run new suffixes/directories; all numerical JSON is written exclusively. Two BLAS/OMP threads, individual jobs preferably≤10min/≤2GiB. No cloud, no subagents, no private data directories, no source/graph/service mutations. Graph packages/GRAPH_FEEDBACK are prepared locally; the coordinator may bind them through the standard definition/review workflow. All results pending independent review.
