# Steering LANE_SURGICAL_INCISION — round 4 (coordinator, 1/10 02:15)

Review of round 3: you found NBSIR 73-262 (normal load versus damage depth — 91 % error in force, not cutting work) and — more importantly — needle penetration where a held-out geometry was predicted with 11 % central error (±38 % digitization; verified r3/geometry/needle_geometry_results.json), effective toughness ≈ 1,8 kJ/m². This is the first unfitted success against skin tonight.

## Critical operation: needle series as an edge-radius series

Needle insertion into skin is extensively published with **series**: different diameters/gauges, bevels, tip geometries and speeds, often with force–depth curves separated into puncture, cutting and friction. A diameter series is exactly the radius series that separates Γ (fracture) from friction/deformation (Shergold & Fleck-type scaling; Barnett 2016, which you already have). Proceed as follows:

1. Collect at least two independent published needle series in skin (preferably human or pig) with numbers in tables/figures; freeze them.
2. Calibrate on a subset (one gauge or one study), predict the rest (other gauges, other studies) without refitting. Report error per held-out case.
3. Separate Γ_pierce and friction; compare Γ_pierce with the Γ_cut scenario (150–380 J/m²) and with the SKIN_TOUGHNESS_GAP lane's mechanism table when available.
4. Close the lane with final PORTS.json for RESPONSE (damage-zone width per tool, Γ_pierce/Γ_cut with uncertainty) and a measurement specification for scalpel data (Ft/Fn, edge radius, crack area).

Inga nya syntetiska FE-/LP-maskiner.
