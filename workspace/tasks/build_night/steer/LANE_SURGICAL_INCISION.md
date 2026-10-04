# Steering LANE_SURGICAL_INCISION — round 4 (the coordinator, 1/10 02:15)

Review of round 3: ni hittade NBSIR 73-262 (normal last mot skadedjup — fel 91 % in force, not cutting work) and — more importantly — needle penetration where a hold-out geometry was predicted with 11 % centralt fel (±38 % digitalisering; verifierat r3/geometry/needle_geometry_results.json), effektiv seghet ≈ 1,8 kJ/m². It's the first unbeaten success against skin tonight.

## Bearing operation: needle series as egg-radius series

Needle insertion into skin is abundantly published with **series**: various diameter/gauge, fasning/bevel, tip geometry and speed, often with force–depth curves divided into puncture, cutting and friction. A diameter series is exactly the radius series that separates Γ (break) from friktion/deformation (Shergold & Fleck type scaling; Barnett 2016 which you already have). Do this:

1. Collect at least two independently published needle series in skin (preferably human or pig) with numbers in tabell/figur; freeze them.
2. Calibrate on a subset (a gauge or a study), predict the rest (other gauges, other studies) without refitting. Report errors per hold-out case.
3. Separate Γ_pierce and friction; compare Γ_pierce with the Γ_cut scenario (150–380 J/m²) and with the SKIN_TOUGHNESS_GAP lane's mechanism table when available.
4. End the lane with final PORTS.json for RESPONSE (damage zone width per tool, Γ_pierce/Γ_cut with uncertainty) and a measurement specification for scalpel data (Ft/Fn, edge radius, crack area).

Inga nya syntetiska FE-/LP-maskiner.
