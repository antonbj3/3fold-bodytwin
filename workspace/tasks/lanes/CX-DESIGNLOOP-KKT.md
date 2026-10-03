# CX-DESIGNLOOP-KKT — (brief = results/CX-DESIGNLOOP/NEXT_NODE.md, dispatched by the coordinator 02:20)
# Proposed next node: CX-DESIGNLOOP-KKT-CAPACITY

Status: proposed follow-up, not dispatched or admitted in the graph. The current working rank/packet for BT-IM2-JOINTFIT calls for force/contact inconsistency diagnosis but has no scientific dependency edge or specific coverage for this femur-design KKT cusp.

Question: Is the 140/141 KKT failure at lengthening −0.3421068643 mm a numerical N43 active-set/scaling failure or an actual equilibrium/capacity boundary?

Frozen comparator: CX-DESIGNLOOP results.json, the full 141-frame z001 profile at −0.350 mm, and PREREG.md. Rebuild the failing shape and identify its exact frame, primal residual, dual residual and active support. Verify with an independent high-precision or differently scaled solve and a per-frame LP feasibility/capacity check at all 141 frames. Probe both sides of the cusp without smoothing or replacing failed points. A true boundary requires a separating infeasibility certificate; an algorithmic failure requires a passing independent solve with the same A,b,N and a documented residual.

Falsifier: the independent solve also fails at the same frame with a certified infeasibility witness, or the apparently accepted neighbor fails a 141-frame capacity check. Keep the earlier 140/141 status until a full independent rerun changes it. The broad CT↔motion registration and CART-KNEE contact gap remain separate prerequisites before any clinical design claim.

Deliver RESULTS.md starting with `# CX-DESIGNLOOP-KKT` in results/CX-DESIGNLOOP-KKT/, results.json and pytest. Write PREREG.md + sha256 first. Run under bigmem.lock with 2 threads. lane runner has full permissions in the workspace; ~/projects/bodytwin is read-only. Internal data stays local.
