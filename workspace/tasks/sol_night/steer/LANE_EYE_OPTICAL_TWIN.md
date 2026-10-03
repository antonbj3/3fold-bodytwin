# LANE_EYE_OPTICAL_TWIN — the eye as an optical unit, all the way to the wavefront

Resultatmapp `results/LANE_EYE_OPTICAL_TWIN/`.

## The operator’s seed
"Hyperrealistic digital twin of an eye, just like everything else." And: the eye as an **optical
unit**, not as tissue. He remembered there was previous work but not where. There is, and
it is better than expected — but it is disconnected from the workspace, and it itself names what is
missing.

## THE PREVIOUS WORK, verified and with path
### 1. The cornea’s refractive power is MEASURED against open data
`source_repository/data/MECHANISM_ANCHOR_GRAPH.json`, node 13,
`EYE-OPTICS-FORWARD-MODEL`. Hidden state: *"eye refractive power/wavefront (what the retina
receives)"*. Leg: corneal curve raytrace with imposed OPD, on a NORMAL subset of 889 eyes
with keratoconus excluded. Anchor: NHANES. Verdict: **the refractive power leg MEASURED on open NHANES.**
Evidence: `reports/probes/mt_eye_corneal_power.json`.

The build script is in
`data/body_twin/agent_scratch_preserved/extern__mt_cornea__build_json.py` and carries the numbers:
- the cornea’s net power ~43 D **coincides with NHANES n = 12 685**, mean 43,57 D against
  the cell’s quoted 43,05 D, z = −0,33, i.e. within 1 SD = 1,58 D
- **the tear film’s telescoping:** air-to-tear gives 43,38 D = 89,4 % of the anterior power, while
  tear-to-stroma gives 5,16 D = 10,6 %. Stroma-to-aqueous humor (−5,88 D classically) almost cancels
  tear-stroma, so **the air-tear transition drives the net power almost alone.**
- three openly reported convention gaps, and they are the lane’s work: the posterior surface gives −4,97 D with
  clinical keratometric index against −5,88 D in Gullstrand, i.e. 15,5 % gap; total power naively
  43,05 + 19,11 = 62,16 D is WRONG because the separation term is missing, two-surface with realistic
  anterior chamber depth gives 60,0 D matching "reduced eye", but Gullstrand-exact 58,64 D is
  another 1,36 D or 2,3 % lower and that is **unexplained**; principal plane H′ was reproduced
  to 1,217 and 1,233 mm but the web-verified Gullstrand value 1,602 mm shows that both
  **underestimate by 24 %**, and the self-consistent axial length is ~24,39 mm, not rounded
  24,0 mm.

### 2. The cornea’s TRANSPARENCY is modeled down to the fibril lattice
`source_repository/data/corneal_transparency/corneal_transparency_results.json`.
I have independently recalculated it and it is internally exact: fibril radius 14,0 nm and area fraction 0,28
give 2D density 4,5472840883398674e-4 per nm² and hexagonal spacing 50,3916571460783 nm — both
reproduce their numbers to the last digit. Refractive index 1,411 in the fibril against 1,365 in the matrix,
Δn = 0,046. Thickness 500 µm. **Wavelength exponent −4,876**, i.e. steeper than Rayleigh’s −4.
And the falsifier is the mechanism: a POISSON-distributed fibril placement would scatter
**20,2 to 59,5 times more, mean 38,1**, than the physiological lattice. Transparency thus comes
from ORDER, not smallness. Validation: 7,93 % median relative error at high disorder, and
the control with a perfect lattice is nonzero because of a declared sidelobe artifact from a finite
patch — not physics.

### 3. What is MISSING, and the graph says it itself
Node 13 carries an explicit limit: *"HONEST_LIMIT: this is SPHEROCYLINDRICAL power, NOT the full
OPD/wavefront raytrace (inherited raytracer is projector-only, 0 eye code)."* Thus
refractive power exists as a spherocylinder, but **no wavefront**. And the plan already exists as node 127,
`wave-eye-forward`: feed Navarro 1985, Liou–Brennan 1997 and Gullstrand’s surface prescriptions — R, T and n per
surface from Atchison & Smith’s tables — into the existing wave optics engine.

### 4. The retina exists as its own track
Node 173, `NEU-RETINA-PHOTORECEPTOR`: foveal angular Nyquist bound rederived from Curcio’s
cone spacing map, 199 000 cones/mm² at the peak. Plus four open contradiction nodes to extract questions from:
tear film optics, lens accommodation, pupil reflex as aperture, and cone Nyquist excess.
Node 14 carries an honest negative that must be respected: **no open dataset with retinal image and
brain MRI in the same individual exists** — converging zero across five archives. Do not search for it.

### 5. The engine modules the wavefront should be built on
`the public staging tree/3fold-physics/src/physics_engine/`: `wave_optics/` med 31 filer,
among them `p20_jwst_zernike.py` (Zernike aberrations, built for telescopes but the basis is the same),
`airy_diffraction_limit.py`, `optics_sphere_psf_diffraction_lambda_scaling_real.py`,
`metalens_psf_imaging_cert.py`; and `ray_optics/` with lens design, achromat and `optics_departure.py`.
The tools thus exist — they have been aimed at telescopes and metalenses, never at an eye.

## Do this
1. **Close one of the three convention gaps, and choose the unexplained one.** The 1,36 D between
   the two-surface model with realistic chamber depth (60,0 D) and Gullstrand-exact (58,64 D) are recorded as
   unexplained. 2,3 % is small and therefore dangerous: it looks like rounding and may be a
   missing term. Calculate through with thick lens formalism, all separation terms written out, and say
   which term is missing or show that the gap is a convention difference in where the principal planes are placed.
2. **Build the wavefront, as node 127 already specified.** Feed the Navarro or Liou–Brennan prescription
   through `ray_optics` and `wave_optics` and deliver **OPD across the pupil plus Zernike coefficients**,
   not just spherocylindrical power. Control: the published Zernike spectra for normal eyes, and
   for corneal power the NHANES anchor that is already MEASURED. Report per Zernike term, never an
   aggregate RMS number — an aggregate hides which aberration is wrong.
3. **Connect the transparency cell to the wavefront.** It calculates scattering from lattice order with exponent
   −4,876; the wavefront calculates phase. Scattering and phase are two channels from the same structure. Minimum
   connection: how much of the retinal PSF’s energy lands in the scattering tail instead of the
   core, as a function of lattice order? That number is **Strehl-like and is the quantity
   that actually measures visual quality** — and it is the quantity refractive surgery changes.
4. **LASIK as the first application, and this is where the chain becomes decidable.** An ablation profile
   changes the cornea’s anterior radius and thus power. But the air-tear transition carries 89,4 % of anterior
   power — so what happens to the tear film over an ablated surface? And ablation removes stroma,
   i.e. fibril lattice, so the transparency cell predicts a scattering change: **haze after
   refractive surgery IS loss of lattice order**, and the Poisson-to-physiological ratio 20–59× is
   the upper ceiling on how much scattering total disorder would give. Calculate what a realistic
   degree of disorder gives, and compare against published haze after PRK or LASIK.
5. **Previous work first, never rebuild.** Everything above is in the old project and nothing in
   the workspace calls `corneal_transparency` — it appears only in a search file and a roadmap. Search
   all four pools (see `tasks/build_night/COMMON.md`) before every absence claim. I have drawn
   four false absences in one evening, all from the wrong search location or search term.

## Strongest control and falsifiers
- **Control for refractive power:** NHANES, n = 12 685, mean 43,57 D with SD 1,58 D. For the wavefront:
  published Zernike spectra. For scattering: the transparency cell’s own 7,93 % validation at high
  disorder. Never our own number as control.
- **Falsifier 1:** if the 1,36 D cannot be derived from any term or from a
  convention difference, the two-surface model is wrong and not just coarse — say so.
- **Falsifier 2:** if the wavefront from published surface prescriptions gives Zernike coefficients outside the
  published spectra for normal eyes by more than their spread, the connection between prescription and
  our wave optics engine is wrong, and that is a finding about the engine and not the eye.
- **Falsifier 3:** if the transparency cell’s scattering change at realistic postoperative
  disorder lies far below published haze, lattice order does not carry haze and the mechanism is
  another — a fully valid outcome.
- **Prohibited:** reporting aggregate wavefront RMS instead of per Zernike term; searching
  for the fenced retina-plus-brain-MRI dataset; rebuilding the transparency cell; writing
  i `~/projects/bodytwin`.

## Delivery
`PORT.json`: OPD and Zernike coefficients per term with unit, the 1,36 D explained or
classified, scattering fraction in the PSF tail as a function of lattice order, and the LASIK profile’s
predicted power and scattering change against published ground truth.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.

# New round (coordinator, 3/10 00:55)
You closed the gap that's been unexplained since July: 58,636095477774106 D reproducerar
Gullstrand exactly, so 1,36 D was a formalism error and not a convention difference. Write it as a
correction of the old work, with both formalisms side by side.

Obstacles: three named failures, and one of them is measurablely insane.
zero absorption transmission lower limit 0,9598 Towards the direction 0,80, that is 20 % Too loud.
zero absorption limit SKA be above a measurement containing absorption, so the question is not whether it
are high without if they 20 % is reasonable absorption or if the spread is absent.
which is required to go from 0,9598 till 0,80, and compare with published absorption in the cornea at
450 nm. If the requirement is excessively high, a spreading term is missing.

Changed operation: take the distribution tail share 1,665e-4 to a Strehl-like quantity and connect
den till gitterordningen i corneal_transparency, som har oordningsfaktorn 20,2 till 59,5The question is:
how much the tail grows per percent lost grating order. It's the chain of refractive surgery
changes, and it's the lan's strongest possible delivery.

And the signs: clinical spherical aberration and haze in clean material removal have both FEL TECKEN
with you. An incorrect character is not a matter of accuracy. Find which term changes characters, for if
Haze becomes negative when removing material, the model says that ablation makes the cornea clearer.

Falsifierare: om C40 = 0,1895 µm is outside published Zernike spectra for normal eyes with
more than their spread, is the connection between the surface recipe and our wave optic engine failure, and it is a finding
if the engine and not about the eye. Report by term, never an aggregated RMS.

# New round (the coordinator, 3/10 02:20)
You went from assumptions to assumptions. 1 455 measured fibrillation centres, and the transmission held is:
again between your two scenarios. That pattern has now been repeated two times in a row with different
methods, which makes the circle a result and not a coincidence.

Prioritet 1, and it is cheap and crucial: you counted that the required absorption is 3,64 till 4,03 per
cm to unite the transmission. Look up published corneal absorption at the same wavelength and
compare. Located the published value in the interval is the missing term absorption and the question is
Closed. If it is far below it is the missing spread, and then it is the structure and not the material.
It's a yes-or-no test against external reference observations and it should be done before anything else.

Prioritet 2: a one-percent tail changes the result by up to 30,6 procent. Det betyder att din
Conclusion depends on the end of the distribution, where you have the least data. Report how many of them 1 455
centres lying in that tail — if there are about ten points, sensitivity is a warning of:
sample size and not physics.

And hold on to what you've won: full-eyed strength 58,636 D who closed the old 1,36 the D-gap; and
spread of the tail in the retinal point spread. The two signs came out wrong, clinical
spherical aberration and haze in material removal, is still unanswered.

Falsifiers: if published absorption does not explain the gap or is excluded by it, then
The transmission cannot be determined from the available data and must be provided as an acquisition item.

# New Round (koordinatorn, 3/10 04:05)
Pressure transfer fell unequivocally: foretold 6,87 till 21,13 mmHg towards the side 84,0, i.e. four more
twelve times too low, with the numerical calculation verified (3,2e-16 full field error, zero non-overlapping error over 24 fall).
It is physics and not calculation, and it is a useful negative: the fixed charge is sufficient
not as an explanation for the swelling pressure.

Changed operation: say which term is missing, with quantity and unit. Candidates to try and not
‘Assume' means the elastic counter-pressure of the collagen net, which carries pressure without being unsmotic; ion activities which:
are not ideal at that charge density; and bound to free water. Count for each which
the size it must have to close a factor four, just as the metabolite chain did with its
clearance. Anyone requiring an unreasonable value is precipitated.

And a warning from your own speech: the thermal spread at 25 nm grows from 3,18 till 13,29 when
the charge is reduced to one percent, i.e. 4,2 Conclusions at low charge are therefore very
It's safer than drunk, and it's gonna stand up.

Keep and don't leave: full eye strength 58,636 D who closed the old 1,36 D-glappet, inringningen av
The transmission between the lattice and fluid model in two independent rounds, and that the independent
human absorption coefficient at 450 nm is not available — It is now an acquisition position.

Falsifier: if none of the three candidate terms can carry a factor four within published
scales, the swelling pressure is not derived from this model class and it should be said out.

# NY RIKTNING (operator 3/10 09:15): THE OBJECTIVE IS ATT SIMULERA OGONOPERATIONS
Anton: "We should be able to simulate surgery of the eye."That's the lan's goal from now on, and
It changes the priority. Everything you built is kept — Most of it is just what is needed.

VAD VI HAR, and it is the pieces that determine UTFALLET:
1. Transparens ur gitterordning, med Poisson-mot-fysiologiskt-kvoten 20,2 till 59,5 as a ceiling for what
   Total disorder costs in spread. It IS haze after refractive surgery.
2. The swelling pressure, which went from 4-12x too low to contain held 84,0 mmHg med 219 av 324 fall
   inom 25 It determines post-shrunk and ectasirisk.
3. The wave front machine per Zernike term, plus the dispersal tail's share in the retinal point spread.
   It's optical quality, that is, what the patient sees.
4Geometry on a nanometer scale with position errors down to 1,5e-06 nm, and four representations where:
   valet gav 676 twice the better limit.

VAD SOM FATTAS, and it determines the incision GEOMETRI:
A. En 3D-ASYMMETRISK corneal surface. Breaking power today comes from the surface recipe — radier, tjocklekar,
   index per yta — that is, rotationally symmetric. A real cornea is not, and astigmatism and
   higher order aberrations IS the deviation from symmetry. That's exactly the refractive surgery
   treats, so without that surface, no surgery can be simulated.
B. FOTODISRUPTIONA femtosecond laser cuts a tab INNE in transparent tissue, which works exactly
   Because the breakthrough is driven by intensity and not by absorption.
   The Laserlan has the thermal regimes; this does not have the optical breakthrough regime.

AMENDED OPERATION, i ordning:
1. Build the asymmetrical surface as a Zernike expansion of HÖJD, not of wave front.
   is a change in height; a change in height results in a change in the wave front through the refractive index difference.
   It is the shortest route from intervention to optical effect and the re-use of the machine you have.
2Chain the three pieces together for ETT ingrepp: en ablationsprofil ger (a) new height area and therefore new
   wave front per term; (b) removed stroma volume and thereby changed swelling and post-shrunk,
   (c) disturbed grating order and thus spread in the tail of the point distribution. Report the three
   The outcomes for the same profile. It's the simulation, and the three parts already exist.
3Tell me what's going on. NOT can be simulated and why, in a row each. The font requires photodisruption;
   The healing over weeks requires remodeling that no lane has.

STARKASTE KONTROLL: publicerade utfall efter refraktiv kirurgi — refractive accuracy, haze and
higher order aberrations before and after. Compare by term, never as an aggregate RMS.
FALSIFIERARE: if the height-to-wave front chain does not reproduce the published force change for a
known ablation depth profile, the clutch is wrong and no operation can be simulated until it is fixed.
PROHIBITED: simulating an operation on a rotational cornea and calling it an operation.
