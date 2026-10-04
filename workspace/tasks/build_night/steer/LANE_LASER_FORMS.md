# LANE_LASER_FORMS — one dose-to-response chain, instantiated for every clinical laser regime

Result directory `results/LANE_LASER_FORMS/`.

## Why one lane and not five
The operator wants laser covered in all its forms. The forms differ in which mechanism dominates, not in
the chain: **delivered dose → absorbed energy distribution → tissue state → irreversible change**. So
build the chain once and instantiate it, because that is also the only way a disagreement in one regime
can inform another.

The five regimes, with what distinguishes them physically:

| regime | dominant mechanism | distinguishing parameter | clinical use |
|---|---|---|---|
| ablative | vaporisation above a fluence threshold | fluence J/cm² | corneal reshaping, skin resurfacing, hard-tissue cutting |
| coagulative | protein denaturation over a thermal history | temperature-time | retinal photocoagulation, vessel sealing |
| photothermal selective | absorber-selective heating with pulse shorter than thermal relaxation | pulse duration vs relaxation time | pigment and vessel targeting |
| photochemical | absorbed photons drive a reaction, little heating | wavelength and total dose | photodynamic therapy, crosslinking |
| photodisruptive | plasma formation and mechanical disruption, ultrashort pulse | irradiance W/cm² | intrastromal cutting, capsulotomy |

**The distinguishing parameter is the point.** The same delivered joule means different things in each
regime, which is exactly why a single summary number has failed us twice already in the laser work.

## What the existing lanes already established, and must not be redone
- A published ablation threshold of **1.15 J/cm²** was acquired, defined as the x-intercept of a mass-loss
  regression with the mass corrected for evaporation. Our chain sat **26.5 % above it** with a slope
  error of −25.0 %, and the lane's own diagnosis was a mismatch between the source convention and the
  observation definition. Read `results/LANE_LASER_SURGERY/` before touching the ablative row.
- **Delivered energy is not sufficient to determine removed mass**: two cases with identity error exactly
  0.0 J gave 0.0 against 1.255e-9 kg.
- **Surface temperature history is not sufficient for the state at depth**: identity error 0.0 K, gap
  0.3806 K, with numerics three orders cleaner at 4.33e-05 K.
- **CEM43 with time and peak temperature also fails**: 7 counterexamples at identity error 0.0 giving
  5.04 percentage points difference. So the standard clinical thermal dose measure does not determine the
  response in our chain either. `LANE_DEPTH_TEMPERATURE` owns the depth field — coordinate, do not
  duplicate.

## Do this
1. **Write the chain once with the regime as a parameter**, so the ablative and coagulative rows share
   the absorption and conduction steps and differ only in the irreversibility criterion. Report which
   steps are genuinely shared and which are not; a chain that needs a different absorption model per
   regime is a weaker claim than one that does not, and both are worth knowing.
2. **Give every regime a published threshold with its measurement definition**, from the literature
   collection. The optics file carries 236 records with units and DOIs, among them `mu_a` = 3.2665 cm⁻¹
   for whole blood at 630 nm at haematocrit 0.45, and the thermal-confound file carries 25 records with
   13 PMID-asserted against the XML. **Read `tasks/build_night/LIT_REFS_FILTER.md` first** — two files are
   excluded entirely, the rest are read with a row filter, and you must report how many records the
   filter removed.
3. **Then run the sufficiency test per regime**, since it has now failed in three different forms here.
   For each regime: is its conventional dose measure sufficient to determine the irreversible change?
   Construct two states with the identical measure to machine precision and report the downstream
   difference. Five regimes, five answers, each with its identity error.
4. **And report the wavelength dependence where it is the distinguishing parameter.** The corneal
   transparency data gives a scattering exponent of −4.8765; the optics file gives absorption against
   wavelength. A regime whose selectivity rests on wavelength cannot be modelled with a single
   absorption coefficient, and if ours does, say so.

## Control and falsifier
- **Control:** per regime, the published threshold with its own measurement definition. This is a
  capability claim only where an external facit exists; where it does not, it is an information link and
  must be labelled so.
- **Falsifier:** if the five regimes need five unrelated chains, the unification fails and the honest
  deliverable is five separate models plus the reason they do not share structure. Report that as the
  headline rather than forcing a common form.
- **Forbidden:** tuning any parameter to hit a published threshold — declare and freeze every free
  parameter before the comparison, because a constant chosen to produce a target is exactly the defect
  found in our front-velocity constant; reporting a dose measure as sufficient without its identity
  error; duplicating the depth-temperature field.

## Delivery
`PORT.json` with the chain, a row per regime carrying the published threshold and its definition, our
prediction and the deviation, the sufficiency verdict with identity error, and the count of filtered
records.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
