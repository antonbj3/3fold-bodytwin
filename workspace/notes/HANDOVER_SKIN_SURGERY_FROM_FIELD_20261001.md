# Handover: the skin puzzle and the surgical threads, from Field (anton-6b) to BodyTwin, 1 Oct 2026 21:55

Anton decided that the skin and surgery threads will be run further by BodyTwin. Field transitions to J77 and elastic materials. Everything below is PENDING_INDEPENDENT_REVIEW.

Written as a separate file because submissions otherwise disappear between sessions (Anton's own diagnosis 1/10: "the problem here with llms and migrations is that things disappear").

## What is decided, and what it excludes

1. **The gap is largely different energy measures, not unknown physics.** Far-field J for skin, 20–30 kJ/m², lies **70–117×** above local shear and capsule energy (measured by Field). It brackets our fibril network target, where I previously wrote 10²–10³ as the expected ratio. Never compare a bond level number directly against 20–30 kJ/m².
2. **Geometry is not enough of an explanation.** Wavy boundary, three-dimensional papillae and mechanical locking together give at most a few times. Field's previous measurement: the wave surface only gives 1,53× area. Don't look for the explanation in area or interlacing.
3. **Blade control does not prevent branches.** Even with a fixed head tip, branches of about 1,6–1,75 µm occur at the epidermis–dermis boundary, and the sharp 1 µm drop is subcritical in the tested scenarios. The hypothesis that a guided blade gives a single clean crack is thus refuted.
4. **Peel is force per width and must be work balanced** before it becomes fracture energy. Additionally, the low DEJ values ​​(Sierra-Sánchez 1,47/2,84 N/m) apply to ODLAD skin and do not identify natural Γ_i.
5. **The needle tip radius is unidentified** without SEM-data, and the replicate spread is 15–27 %. Combined with our own result (all interfacial mechanisms excluded, the velocity dependence 2,95× is in the tissue) this means that apparent toughness probably carries a tip radius term, F = w(Γ₀ + κ·R_spets).

## The deciding factor that remains, and it's ours

- **Natural Γ_i for the epidermis–dermis boundary.** The only way we have is blue time in natural skin. The job is `results/BT-SURGX-RETE-BLISTER` (queue, location 188). It should deliver Γ_i as interval plus retelist geometry as PORT.json.
- **DEJ-mikroskopi after section.** Published histology of fresh sections exists but lacks DEJ-statistik. It is a measurement to buy or find, not to calculate.
- **The two data layers of the fibril mesh**, queued as `BT-FIBX-FIBRIL-RUPTURE-WORK` (location 162) and `BT-FIBX-FIBRIL-JUNCTION-WORK` (location 176). The network cannot provide an absolute Γ without them. Memory rule from our own round 2: which interface breaks is governed by the ratio eta between junction and fibril strength, with wrap close to eta ≈ 1.

## Field substrate, PORT.json in the respective catalog

Rot: `../3fold-motion-engine/_private/romi_collab/build/`

- `PROOF_LANE_SNITTDIFF_20261001` and SNITTFRONT: guaranteed energy difference, on branch `falt-snittdiff-20261001`. Draft-free cut gives 309–381× lower noise in G than element removal, but the guarantee is too wide and gives TIE against classic FE.
- `PROOF_LANE_GRANSSNITT_20261001`, `PROOF_LANE_KROKTGRANS_20261001`, `PROOF_LANE_SNITTKONTAKT_20261001`: interface mechanics, curved boundary, blade contact.
- `SOL_FALT_BLADBANA_20261001`: the leaf path and the branches.
- `SOL_FALT_NALKRAFT_20261001`: power vs. depth, and `raw/rate_pilot.json` our held-out Barnett anchor came out of (their digitization, so semi-external).
- `SOL_FALT_NALFRIKTION_20261001/raw/anchor_curves/INDEX.json`: six digitized friction curves with reading errors. Verdict: the assumption of equal friction on reinsertion is UNCERTAIN, no series separates first from repeat pass in same hole.
- `SOL_FALT_NALSPETS_*`: spetsradien.
- Web research: `~/research/HUD_NAL_HYPOTESER_20261001/` (HYPOTESER.md, KALLOR.md), seven hypotheses per puzzle, sources tagged by abstract or full text.

## Tools Fields continue to own and that we can use

`ELAST3DNAT` and `HYPERCERT`, elasticity and hyperelasticity. They offer them as tools for us and answer questions.

## Our own results belonging to the same thread

- `results/LANE_NEEDLE_FRICTION_ANCHOR/` r1–r3: all interface mechanisms excluded. Viscous film gives velocity ratio exactly 80 at each film thickness to measured 2,95; poroelastic pore pressure 1,03–1,38; frozen transferred bond layer 1,455 to gate 1,475. Direction difference 2,67 was dropped by elastic memory (1,001).
- `results/LANE_FIBRIL_NETWORK_FRACTURE/` r1–r3 (parked): the bend slope converges to the published value with magnitude, 1,186 → 1,033 → 0,980, where the two largest enclose 1,0. Γ still null. **Eelastic straightening and rotation is recoverable and must not be accounted for as dissipation** — that was my own fault in a previous control.
- `results/LANE_SURGICAL_*`, `LANE_SKIN_TOUGHNESS_GAP`, `LANE_WOUND_CROSSLINK_CHEMISTRY`: previous lanes, see NIGHT_LOG.md for verified lanes.
