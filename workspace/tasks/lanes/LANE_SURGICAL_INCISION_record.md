# LANE_SURGICAL_INCISION

Anton's direction 30/9: BodyTwin should provide high resolution where a scalpel opens the body — materials, cells, tissues and bonds. This is a new lane in tonight's pursuit of breakthroughs; results folder `results/LANE_SURGICAL_INCISION/`.

## Desired capability

A multiscale incision operator: given blade geometry (edge radius, angle, speed, depth/path) and a local tissue stack (epidermis, dermis with collagen network and tension lines, subcutis, fascia, muscle, vessels), predict cutting force over depth, damage-zone width (cells, collagen, vessels), wound gap, bleeding flow until hemostasis and what determines the initial state of healing — with high resolution only where the blade is and a coarse model in the rest of the body, without global recomputation when the incision is extended.

## First principles to break down into

- Bonds: collagen triple helix, crosslinks, fibril sliding; fracture energy per bond → fibril → fiber → network → tissue fracture toughness (J, G_c). Thermal denaturation during electrosurgery.
- Cutting mechanics: work = new fracture surface × toughness + friction + far-field deformation (Atkins/Mai-type cutting mechanics, cohesive zones); the role of edge radius; difference between cutting/tearing.
- Prestress and anisotropy: in vivo prestress and fiber direction (Langer's lines) determine the gap; compare with Q033's elliptical opening model and the cadaver anchor 4,2–6,0 N closure force.
- Cells: damage zone from strain/temperature, cell membrane rupture, ischemic margin (Q036's state vector).
- Vessels: vessel density and caliber per layer, Poiseuille at severed ends, platelets under shear, thrombin generation.

## Existing code to build on (read, extend, do not redo)

- `results/BT-HX-Q033/` (local incision: gap, closure force, Darcy, evaporation), `results/BT-HX-Q036/` (ischemia/reperfusion), `results/BT-HX-Q049/` (thin layers), `tasks/free48/sources/SURG_*/`.
- Private code, read-only: `~/projects/bodytwin/scripts/msk/{skin_pulp_mechanics.py, collagen_triple_helix_thermal_stability.py, tendon_collagen_hierarchical_mechanics.py, myofascial_transmission.py, coagulation_hemostasis.py, platelet_hemostasis.py, wound_healing_cascade.py, acute_phase_inflammation.py, bone_fracture_toughness_lefm.py, skin_barrier_tewl.py}`.
- Graph nodes (context): SKIN-SUBCUTIS-DECOMPOSITION, SKIN-DERMAL-STIFFNESS-COLLAGEN-COMPOSITION, MSK-FASCIA-NETWORK, MODEL-TENDON-COLLAGEN-MECHANICS, MODEL-COLLAGEN-TRIPLE-HELIX-THERMAL-STABILITY, MODEL-COAGULATION-*, MODEL-PLATELET-ACTIVATION-AGGREGATION, WOUND-HEALING-CASCADE. Graph goal for binding: `BT-CTX-SURG-INCISION` (kind define/review).

## References to the world

Look up published measurements of cutting force/cutting energy in skin and soft tissue (scalpel, needle, cutting mechanics), skin fracture toughness, in vivo prestress and gap, vessel density per skin layer and bleeding time. Preregister which numbers the model should match before running it. Do not invent measurement data.

## Data boundary


## Starkaste kontroll

Uniform fine FE with a cohesive zone/phase-field fracture (or XFEM) with the same material data and the same prestress, full cost. The innovation must show the same error at substantially lower cost, or a prediction the control cannot provide (e.g. damage zone + bleeding + gap from the same bond parameters).
