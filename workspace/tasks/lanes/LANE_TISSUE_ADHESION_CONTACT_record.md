# LANE_TISSUE_ADHESION_CONTACT — tissue contact as adhesion, not as Coulomb friction

Result directory `results/LANE_TISSUE_ADHESION_CONTACT/`. Tissue coupled to geometry as support for
medical robots.

## The measured state
Two independent findings tonight point to the same thing, and both are outside the eye:

1. **The needle.** Edge `T-E18` now carries an identifiability theorem: with R friction, S other resistance and C
   cutting force, `F₁ − F₂ = C₁` holds exactly if and only if `(R₁−R₂) + (S₁−S₂) − C₂ = 0`. From two
   total-force records only `C₁ + (R₁−R₂)` is identified; the friction difference is not separately
   identifiable. Witnessed with 13 104 indistinguishable pairs at observation error exactly zero. Simultaneous
   tip-force measurement solves it — and it exists: Fukushima and Naemura, *ROBOMECH Journal* 1:14
   (2014), doi `10.1186/s40648-014-0014-7`, coaxial needle measuring total force and tip force separately in
   real time.
2. **The mechanism correction.** PMID 41522820, doi `10.1038/s41524-025-01869-y`: in soft materials
   tangential stress is dominated by **adhesion and damping**, and Coulomb friction is negligible at
   low contact pressures. Also PMID 33445344, doi `10.1021/acsbiomaterials.8b00387`: adhesion force
   between tissue and bipolar forceps against tip temperature, **in newtons**.

And the external constant circulating, µ = 0,4, lies 1,69–1,88 σ above the µ it itself cites
(0,295 ± 0,056 static, 0,255 ± 0,086 dynamic, Urrea et al., *JMBBM* 56:98–105, 2016). It is not
ours — no BodyTwin model carries a friction coefficient — but it is what everyone else calculates with.

## The obstacle
Our entire contact description is a friction coefficient that does not exist in our models and that in
the literature is the wrong quantity at the pressures an instrument actually works at. What is needed is a
tangential load law where adhesion is the carrying term and temperature is an input, not a
disturbance.

## The operation
Build the law and make it a **decision**: given tip temperature and contact geometry, what
tangential load should a robot-held instrument tip plan for? Return the load, not the description
of the load. Adhesion against temperature is measured in newtons in PMID 33445344 — use it as the reference, not as
illustration.

## The strongest control
Coulomb with µ = 0,295 ± 0,056, thus the best measured friction coefficient in the literature, run against
the same reference. If the adhesion law does not beat it, it should not go in.

## Falsifier
If the adhesion term does not dominate at the contact pressures an instrument tip actually exerts — calculate
that pressure from tip geometry and planned load and state it in kPa — then PMID 41522820 does not apply to our
case, and then the edge should say Coulomb suffices in our pressure range.

## Regler
`PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material, inget som
names the collaborator. Each source with PMID eller DOI, volume and pages.
