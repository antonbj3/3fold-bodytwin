# Ambitious coupled history inverse — active evidence checkpoint

**Review:** PENDING_INDEPENDENT_REVIEW. **Innovation gate:** not passed.
The current attack continues; this file preserves completed evidence.

## Capability and obstruction

The desired consumer edits future geometry, a physical coupling and a dosing
schedule, then makes robust decisions under many local exposure constraints
without replaying every region and history for every design. Ordinary sparse
BDF and source-aware context ROM already solve individual design queries well.
Increasing region count or wrapping those methods in an inverse does not pass
the innovation gate.

The shared benchmark contains native Dalla Man organ feedback, nonlinear
finite-site binding, reversible washout, gate history and neighboring exchange.
Its source rows couple permeability, capacity, binding rate and memory via one
synthetic μ. Added biological closures and empirical calibration are UNKNOWN.
Finite μ rows are scenarios, not a calibrated uncertainty envelope. Every
structural certificate remains conditional on numerical ODE/float remainder,
which is not rigorously bounded here.

## Executed changes and outcomes

| Operation | Actual outcome | Barrier / executed next change |
|---|---|---|
| 36-design robust inverse, N32/128/512 | 3 feasible of36; changed lattice empty | Finer shared-dose probe finds0.25U feasible: lattice resolution, not impossibility |
| Conserved inventory/native-potential error coordinates | Slower and looser than component control | Preserve negative; continuous inverse-domain operation |
| Ordered two-pulse fixed-total source manifold | Early transfer domain5.76× triangle bound; late1.01× | Endogenous history destroys simple cancellation; actual material/history anchors |
| Two-phase/context native histories | Mean-only history wrongly rejects all12 material designs | Preserve phase/context histories, not marginal labels |
| Local Markov + all-local guard extrema | Cheap boundary Newton extraction improves measured setup | Two proof repairs required; previous v3/v4 analytic qualification withdrawn |
| Continuous local geometry domain | 15/23 zero-new-ODE queries,8 fallback | Cost0.982× strongest same-tube control: innovationFAIL |
| Capacity-release inverse | Exact event preserves mass and prior error | Conventional exact-event/controlTIE |
| Reachable support0/8 historical inverse | Correct κ=.30; means-only cache κ=.05,5falsefeasible | Source/history-support ROMTIE; execute conditional deletion |
| History deletion + new chord + continuous domain | Aged far history:68states actually removed; source3×two domains accepted at0.1% goal | Conventional component radii tighter; fresh history's context-floor destroys usefulness; execute signed moment/capsule pivot |

## Latest composed operation: actual deletion and coupling

`FORGETTING_DOMAIN_GATE.json` predates the new runs. N512, sourceμ=.8/1/1.2,
actual reachable historical support0 or64, ages0/120min, a new chord0↔16 or
0↔64, and future local perfusion/capacity rectangles±0.0001/±0.0005 were used.
Future SC3U at2min is a synthetic forcing code gate. One reduced future anchor
per source/history/coupling is acquired; domain certificates solve no ODE per
new design. The query is future exposure at node0, with a new clock. Prior
cumulative exposure is not silently erased from cumulative queries.

The new chord is consequential: for fresh support64 its endpoint changes
future node0 exposure by9.62–9.79mg/kg·min. When it touches historical support64,
that endpoint remains active and the old collar cannot be discarded. For the
far chord0↔16 after120min, deleting the old collar removes68states and passes
both declared rectangles at a0.1% relative exposure goal. Sharp0.001% goals
remain unresolved. This is a local query demonstration, not closure of the
original arbitrary all-region constraint capability.

`FORGETTING_DOMAIN_VALIDATION.json` contains120 independent full-network
corner/center forecasts: zero observed enclosure misses or false accepted
decisions. Corner replay validates the construction; it does not prove a
continuum by sampling. Past snapshots, old/dropped states and exact descriptors
are saved in `forgetting_inputs/`; initial input acquisition includes full-N
projection/witness work and is not called O(active).

The first composed run costs5.93CPU-s and~167MiB peak RSS;120 heldout forecasts
cost20.18CPU-s. Prefix aging costs0.381CPU-s. The previous source acquisition,
archive/history construction, witness extraction, anchor ODE, certificate and
fallback are separate recorded costs. The repeated run adds state exports;
its cost is retained, not hidden. No practical gain against the strongest
control has been established: aged component radii are~0.17–0.18 versus
candidate~0.31mg/kg·min at±0.0001. Their bootstrap is explicitly a separate
qualification requirement; numerical radii are not promoted to certified
control acceptance without it.

## Current load-bearing pivot

Fresh deleted history has tiny distant observed influence (~1e-5mg/kg·min),
but an absolute context floor produces61–105mg/kg·min bounds and global tube
failure. Deletion exactly preserves signed free/bound/gate means per context.
The next construction retains signed defect moments and temporal influence
capsules, so native mean forcing is charged for nonlinear covariance/Jensen
effects instead of an invented net inventory perturbation. Body and graph
workers own the mechanism and proof; this consumer supplies exact original,
retained and dropped state inputs and will use the resulting operation in a
whole-domain inverse. Identical moments, invariants and witnesses remain
available to the strongest conventional comparator.

## Whole-domain inverse and matched spatial control now executed

`SIGNED_WHOLE_DOMAIN_RESULTS.json` and `FINAL_DOMAIN_DECISIONS.json` close both
continuous rectangles for all three sources under one common0.1% exposure
cap/floor interval. Sharp0.001% goals remain UNKNOWN. No new ODE is solved per
design. History/native and geometry sectors use one widened physical tube;
their native/free/gate error bounds are summed before admission. A first run
correctly refused admission because the total gate error exceeded its0.005
enclosure margin. The saved repair recomputes the actual Jacobian/Hessian and
positive supersolution with a0.01 regional gate enclosure, while preserving
native margins and the physical exposure error budget.

`sparse_hessian_control.py` implements an actual different spatial comparator:
the same analytic signed-mean/RMS/native sector, followed by a full3N sparse
positive reaction-diffusion comparison and `expm_multiply`. Initial physical
errors are taken after preserving the query collar; signed context residuals
are bounded at their individual interfaces and compensation regions. The
control uses the same past inputs and future BDF anchor, and may discard old
history after constructing its own sufficient error vector. Its radius is
**2.164–2.172**, tighter than the capsule radius **2.206–2.215mg/kg·min**. Both
close the declared whole boxes. A first scalar implementation took~0.06s;
the saved optimized implementation takes~0.0053s for its spatial comparison.
Identical radii before/after optimization prevent a weak-code cost comparison.

The capsule query costs~0.001s, but the retained-context ROM's already acquired
domain comparison costs~0.00035s. Including actual past acquisition, projection,
anchor, contrast and failed RMS-discovery attempts, candidate/retained-ROM cold
ratios are **0.54–0.67**, failing the0.5 gate. The shared120min prefix acquisition
is an additional common cost; adding it pushes these ratios toward1. No10×
full-cost gain or new method class is established. Warm capsule evaluation has
actually been called without the discarded old model/state.

Positive scope is a query and coupling announced while selecting protected
collars. It is not arbitrary new endpoints announced after deletion. The body
peer's actual post-deletion0↔64 test fails the0.1% budget by0.141–0.144%; initial
mean conservation and tiny organ error do not rescue that local exposure goal.

## Executed heterogeneous endpoint pivot and strong compact baseline

The next source uses20 genuinely heterogeneous, overlapping past pulses under
the same endogenous native feedback. Inputs are frozen in
`heterogeneous_history_inputs/`. The peer's325 Haar coefficients are three
component vectors each; their required storage exceeds a fair compact retained
stock baseline. A translated single-pulse dictionary is also an eligible
conventional baseline and is not claimed as the next innovation.

`quantized_history_control.py` executes a strong conventional control on those
actual sources: bitpacked free/bound/gate stocks, exact context means and
covariance, native/ledger state and explicit reconstruction error bounds. It
resets only future-window exposure. Exact source-bound binary packet sizes are2543/2577/2611bytes. Reconstruction
errors are within the declared endpoint tolerances: free≤0.00260mg/kg,
bound≤1.66e-5mg/kg, gate≤1.33e-7. A new endpoint64 is retrieved after encoding
without reading the old/full state. Independent full future replay for the new
0↔64 chord has relative exposure error≤2.99e-6. This is numerical validation,
not a rigorous future certificate. The physical past-decoder construction must
beat this control, rather than an8.8KiB float64 baseline that stores old E.

## Provenance and limitations

`SOURCE_USAGE.json`, `INPUT_MANIFEST.json`, `PARENT_INPUT_MANIFEST.json` and
per-copy manifests pin actual W01/W07/W09/500-M423 mechanisms and raw parents.
The native NFkB ancestry of the old barrier is distinct from the current
glucose/endocrine forcing. No NFkB↔glucose biological equivalence is asserted.
`COMMANDS.md` preserves executed probes and failures. Previous failed gates,
old code and numerical outputs remain on disk. No canonical scientific status,
source graph, cloud policy or service is changed by this lane.

## Late-endpoint whole-domain strong control and next proof operation

`QUANTIZED_LATE_INVERSE_RESULTS.json` now closes the conventional packet control
end to end: late0↔64 chord, actual coupled future, packet epsilon, ROM residual
and continuous±0.0005 local perfusion/capacity box under one component tube.
All three joint source boxes meet the0.1% exposure goal; radii1.70064–1.82506
mg/kg·min.15 independent full future forecasts in
`QUANTIZED_LATE_VALIDATION.json` show zero misses or false accepted decisions.
The packet is the only oldstock input at future-query time. Rigorous numerical
ODE/float remainder remains UNKNOWN. Read/model~0.0065s, anchor0.178–0.195s,
certificate~0.0006s; total three-source caller0.576CPU-s.

The counterfactual-past branch has a frozen four-dimensional prospective family
in `COUNTERPAST_LATE_4D_PREREG.json`; it has not yet closed that domain.
An executed prerequisite changes the transport operation: protected-centre
shortcuts have zero generator contribution to a nearest-centre cosh potential.
48 exact stopped CTMC cases and3 adverse scope checks PASS. For rate4/T30/r9,
newbound0.0016892 replaces oldresetbound0.20439; exactexit0.00033103. Analytic
call~7µs versus exactmatrixexponential~116µs plus~450µs generator setup.
The conventional control may use the same supersolution. The first unnormalized
floating residual test failed due to scaling and is preserved; normalizing by
the boundary potential repaired conditioning without changing tolerance.
This mathematical prerequisite is not a completed innovation-round gain.

## Strong changed-past/coupling family comparator actually executed

`COUNTERPAST_ATLAS_CONTROL_RESULTS.json` acquires true native-endogenous past
fractions0.155/0.2325/0.31 at64 and late chord rates0.08/1/4. The same9 samples
fit two eligible tensor controls. Four independent mixed heldouts give max
free136 error0.1182 for raw rate (fails0.02), versus0.004504 for physical
transfer `1-exp(-2d*1min)` (passes sampled goal). Physical-basis exposure error
0.39574mg/kg·min meets the sampled0.1% goal. Whole-domain interpolation error
is UNKNOWN, m0/cap0 are fixed and no4D decision admission follows.

The actual strong family control has no new-query ODE, setup2.438CPU-s and
median query33microseconds. Total4.558CPU-s includes validation; an earlier
~5s final-write failure is preserved separately. These are known controls,
not a completed innovation. The decisive missing operation is a chronological
endogenous family witness with a nonlinear remainder, rather than final moments
or unchanged native ports.

The matched1e-5 precision control (`COUNTERPAST_ATLAS_PRECISION_RESULTS.json`)
reduces setup to1.435CPU-s but fails a heldout sharp local query: interpolation
error0.025567mg/kg versus0.02, while direct same-precision error is0.002802.
This acquisition-noise amplification is a concrete whole-family proof/cost
barrier. Strict and fast evidence remain separate; no4D corners were added.
