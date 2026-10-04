# LANE_LASER_SURGERY — laser-tissue interaction, and the non-contacting robot

Results directory `results/LANE_LASER_SURGERY/`.

## Why the lane exists, and why it is cheap
The operator's seed: laser surgery in various forms — eye, skin treatments — and "there we already have very
a lot from the lithography work". The same lane also gets to work with medical robots.

The seed bears, and for three reasons that are sated tonight rather than assumed.

**1. We have a LIVING optical defect that the laser literature is densest on.** The bleeding lane tries
determine if an operation field is obscured and does not proceed: `held_roggan_optics = FAIL`, 6 of 10
anchor consistencies, maximum relative error 0,304, and the repair attempts went the wrong way (to remove
packaging raised the error to 0,718). The missing physics is absorption and distribution in blood and
tissue as a function of wavelength. It is exactly the size family that laser surgery accurately measures,
because its entire dosage rests on it. A lane that fetches those constants against published reference
lifts a gate that is standing still.

**2. We already have a receiver built.** `tasks/free48/sources/SURG_COLLAGEN/` emits
**R_thermal ∈ [800/177, 80000/503]**, thus 4,5198 to 159,0457 — the work to cut native tissue
divided by the work of cutting the same tissue above the collagen's melting temperature. Γ_fibril is deleted in the quota,
which is why it can be emitted while the absolute number is null. The cell was built for a warm
Blade. **A laser is the canonical case of a tool that crosses that boundary**, and the anchor line
gave the quota a source-backed primary set with executable consumer tonight. So this is one EDGE to one
existing greatness, not a new branch.

**3. The observability thesis also applies here.** Measured over 5 984 swarm reports: 13,3 % make a
explicit claim of structural unidentifiability against 0,8 % complaining of parameter uncertainty.
The value of a new track is in if it carries a OBSERVABLE we are missing. Laser surgery carries several: they
the optical coefficients, thermal relaxation time, and ablation threshold in fluence — all with unity and
all published against reference.

## What actually exists, where it lives, and what is NOT transferred
The material comes in two levels, and the difference determines what you get to build on.

**Experiment script** in `source_repository/scripts/physics_exp/` — mirrored from
cad-to-simulation lane repona (mainly lane I, some lane A) and also into `mechanism`:
`euv_reflectivity_real_roughness_cell.py`, `etendue_budget_cell.py`, `laser_speckle.py`,
`slm_maskless_writer_map.py`, `dual_triplet_cell.py`, `i_asml_stack_apriori_cert_map_cpc_gated.py`.
Broader: 47 files mention lithography, 124 EUV, 245 laser.

**Promoted engine modules** in `../3fold-physics/src/physics_engine/`
— these are packaged modules, not experiments, and are what you should primarily focus on:
- `scattering/` — **the most important one for us.** Mie scattering in four files
  (`p18_mie_scattering_render_match.py`, `dielectric_sphere_permittivity_recovery_mie.py`,
  `p18_dielectric_sphere_mie_forced.py`, `p18_mie_sphere_rcs_render_match.py`) plus Fresnel.
  **Tissue optics IS Mie scattering** — cells and organelles as scatterers — plus absorption. It is
  exactly the physics bleeding lane's failing optics gate needs for a reduced
  spreading coefficient, and it is available as an engine module with render-match against reference.
- `wave_optics/` — 31 files, among them `airy_diffraction_limit.py` for spot size and
  `wave_optics_cell.py` with GPU variant.
- `litho/` — 14 files: aerial image scalar with GPU variant, `litho_depth_of_focus_cert.py`,
  `litho_euv_stochastic_ler_cert.py`, the etch level set. What is transferred is dose distribution over one
  field and depth of field as a TOLERANCE instead of a value.
- `process/laser_keyhole_absorptance.py` — absorbed share at a keyhole, i.e. the connection from
  incident energy to deposited energy.
- `ray_optics/` — lens design, achromat, thermo-optical tuning.

**A TRAP, and I almost walked into it myself:** `wave_optics/laser_threshold.py` is about
**laser cavity threshold** — rate equations, gain clamping, threshold pump 1/(B·τ_c·τ_2),
Siegman as a reference. It is NOT ablation threshold in tissue. The names are similar and the physics are
another. Do not use it for tissue.

**THE GAP, and my own correction of it:** I first wrote that the thermals are not covered because
`thermal/` only has four files about Kelvin–Helmholtz and a flame. It was an absence claim
delimited to a DIRECTORY, and the module is located elsewhere:
`process/moving_heat_history_v1.py` is a moving heat source with history, that is, precisely
The Rosenthal machinery for a SCANNING beam, and `chains/process_heat_horizon_probe.py` are available
next. Use them. That's the fourth time tonight I've pulled a fake absence from the wrong search location, so
check yourself in all four pools before saying something is missing (see COMMON.md).

What does not transfer at all: resist chemistry and mask design. Say which files you are writing off and why.

## DENTAL, operator follow-up question — and dental already has more than this lane
Measured in dental's own workspace `local_path`: **315 files mention
laser, 36 mentions Er:YAG, 1 202 mentions ablation, 2 218 dentin, 4 330 pulp.** So do not duplicate
their work. The dental task of this lane is the SHARED QUANTITIES, and they are three, for the three
the tissues have three different mechanisms:

1. **Enamel, ~96 % mineral.** The ablation is water-mediated — the absorption is in the water and
   the mechanism is microexplosion rather than melting. Er:YAG at 2,94 µm lies on the water's
   absorption peak, which is the whole reason that wavelength is used clinically. Shared greatness with us:
   absorption coefficient versus wavelength, the same magnitude the bleeding lane's optic gate falls on.
2. **Dentin, about one fifth collagen by weight.** Here we have a FINISHED RECEIVER:
   `SURG_COLLAGEN` emits R_thermal ∈ [800/177, 80000/503], the work of cutting native tissue divided by
   the work above the collagen's melting temperature. The collagen matrix of the dentin crosses the same boundary during a
   laser pulse. **There is an edge between our cell and dental's domain, and it is drawable today.** But
   read the warning above: coagulation at ~60 °C is time dependent while Tm is an equilibrium temperature, and
   the collagen cell has no time variable.
3. **The pulp — and here sits the binding term.** The clinical ceiling is not the ablation but
   the temperature increase inside the pulp. Get the threshold against published reference with device; treat
   every number you remember as a candidate to verify, not as known. It is one
   **heat diffusion problem under a scanning source**, thus `process/moving_heat_history_v1.py`,
   and it is the same physics as the thermal relaxation time in the skin.

**The cross against dental's own results, coming in tonight:** their lane X26 graded all eight
overall dental physical decisions such as CLASS 2, cross-size floors, each with a named
binding magnitude — among them signed recorded scanner error, contact force distribution, and
scanner attenuation plus leg lag. **Question worth asking:** if the temperature rise of the pulp is it
binding quantity for laser dentistry, it is CLASS 1, i.e. precision limited and down payment
with resolution, in a world where everything else dental is CLASS 2? A non-contact tool removes
the contact force distribution, which is one of their eight binding quantities. If the answer is yes, it is laser
only dental decision path that can be certified by buying resolution — and that would be one
heavy find. Run the classifier on it, and **coordinate with dentallan instead of guessing**;
their table is in `results/LANE_X26_DECIDABILITY/DECIDABILITY_TABLE.csv`.

## The physical dividing lines that determine everything, and which must be declared per question
The laser-tissue interaction is sorted by two quotas, and a lane that does not declare them will
mix up four different physicists:
- **pulse duration versus thermal relaxation time** for the target structure. Shorter ⇒ the heat stays in the target
  (selective photothermolysis, skin treatments). Longer ⇒ the heat diffuses out and damages the neighbor.
- **fluence against ablation threshold**. Below ⇒ coagulation, above ⇒ material removal.
Four regimens follow: photochemical, photothermal, photoablative, and plasma-induced/photodisruptive. Each
question must name its regime, because a constant measured in one regime does not apply in another.

## Do like this
1. **Get the optical coefficients against published reference** for blood and the skin layers
   `SURG_HEMOSTASIS_model.py` declares (epidermis, papillary dermis, reticular dermis, subcutis):
   absorption coefficient, scattering coefficient, anisotropy and thus optical penetration depth per
   wavelength, with unit and provenance class. Deliver them in a form the bleeding line's optic gate can
   consume — say which of its ten anchor consistencies each addresses.
2. **Connect R_thermal to a laser.** The collagen cell declares the melting temperature 310,15–313,15 K.
   Calculate the fluence and pulse length that bring a given tissue volume over that limit, and say what
   R_thermal then predicts for the cutting work. It is a falsifiable prediction against published ones
   laser cutting powers, and it is the lane's strongest delivery if it can be done.
   **Warning already posted:** an entry in the surgical source suggests that tissue coagulation
   occurs at ~60 °C and is TIME DEPENDENT (CEM43 dose), while Tm is a molecular equilibrium temp.
   The collagen cell has no time variable. If true, R_thermal is underspecified by one
   laser, and **showing it is as good a result as connecting it**.
3. **The eye and the skin are two different questions — choose ONE per round.** The eye: corneal ablation has
   published thresholds and a clinical reference in refractive outcome. The skin: selective photothermolysis has
   published relaxation times per target structure. Tell me which one you chose and why.
4. **The non-contacting robot, and here is a structural difference worth noting.** Four
   lanes has built on the robot's binding term being registration time, 150 s per registration, and that
   1 192 updates fit in one procedure. A laser is **non-contact**, removing the entire
   the contact force branch — including the DEJ peel requirement ≤0,1 mN which is the binding mechanical requirement for a
   contacting tools. Question: what becomes the binding term when the contact force disappears? Candidates
   is line of sight, thermal accumulation between pulses, and eye movement or tissue movement between
   registration and shots. **NOTE a fresh uncertainty:** the bleeding lane found
   `published_timing_excludes_spatial_registration = True`, so the 150 s may be incorrectly delimited.
   Lane `LANE_REGISTRATION_SCOPE` is working on it — don't build on the speech without reading its
   port first.

## Strongest control and falsifier
- **Check:** for the optics, the bleeding lane's current optics model that supports 6 of 10. The profit is how
  many of the ten anchor consistencies that pass after your round. For R_thermal: published
  laser cutting powers, not our own speech.
- **Falsifier 1:** if the published optical coefficients of tissue are so wavelength and
  condition that no set can handle more than 6 of 10, the occluded area is not
  certifiable from optical data and it must be stated outright as an acquisition record.
- **Forger 2:** if R_thermal's prediction is outside of published laser cutting powers and more
  than an order of magnitude, the ratio is not transferable to thermal denaturation by laser, and
  the collagen chain is underspecified. It is a perfectly good and important outcome.
- **Falsifier 3 for the robot part:** if the contact force branch removal does not change which term
  binds, non-contact is irrelevant to the robot calculation and should be reported.
- **Forbidden:** to use a laser constant measured in one regime for a question in another without saying
  the; to use a supplier figure as a measurement; to build a new cell before an edge to one
  existing greatness is drawn; to touch `~/projects/bodytwin` other than as read view.

## Delivery
`PORT.json`: the optical coefficients per layer and wavelength with unit and provenance class, which
of the bleeding lane's ten anchor consistencies they address, R_thermal's laser prediction against published
reference, and the binding term for a non-contacting tool. Plus the list of lithography files you
wrote off, with reasons.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
