# Styrning r13 — Close the circuit back to the wave front, it is open since r7

## Obstacle, measured in lanen's own files
the lane owns real optics and it should be said first: `PORT_R7_V2.json` berries **total brytkraft 58,636 D med
28 enskilda Zernike-termer** (`ZERNIKE_TERMS_R1.md`), OPD sparad som `.npz`, and a
**Dispersal tail at an arc minute at: 1,6649e-4**.

But that line of tail bears its own warning: *"R1 definition only, **not converted from this axial phase
model**"*. And I searched the whole `r12/`: **noll filer** mentions wave front, Zernike, cornea, stroma,
breaking or spreading. 8 till 12 has worked with charging inventory and
Poisson–Boltzmann-intervall — physically relevant as a basis for hydration and swelling pressure; and
calculationally clean (`max_charge_ledger_error = 6,62e-16`) — men **none of the five rounds have:
producerat en optisk konsekvens.** The circuit from charging to wave front has been open since r7.

## The operation this round: an entire chain, no new branch
The operator's goal is to be able to: **simulate eye surgery**An operation changes the thickness, curvature
and hydration; the twin must tell what happens to the vision. EN komplett kedja denna
rounding, even if each joint becomes coarse:

1. **laddningsinventarium Q_t → hydrering/svoperating pressure** (det du redan byggt i r8–r12),
2. **→ change in stroma thickness and refractive index**,
3. **→ changed Zernike coefficients** i de 28 termer du redan har,
4. **→ a surgically meaningful reading measure**: break force in diopters and strehl numbers; or
   wave front RMS.

Report one number per line and indicate which line is weakest. A chain with a coarse line is worth more than
a fifth charging interval, as it is the first to be tested against clinical data.

5. **Convert the distribution tail or tell me why it cannot be converted.** 1,6649e-4 applies to:
   R1-the definition and not the axial phase model. Either the conversion is done now, or it is written
   det som en namngiven brist — but it must not remain as a number without its validity.

## Strongest control and falsifiers
- **Kontroll:** a rotational symmetry of the cornea with the same mean thickness. It is prohibited to simulate a
  surgery on a rotational cornea as if it were patient geometry, but which KONTROLL is it right:
  profit is the part of the wave front change that the symmetry cannot reflect, i.e. the non-symmetric
  Zernike-termerna.
- **Falsifierare:** if the whole charging chain changes the breaking force less than the measurement uncertainty of a clinical
  autorefraktor, cirka 0,25 D, the charge level is not load-bearing for optics — and then:8–r12
  is reported as evidence for swelling, not vision. That would be a clear and useful negative
  resultat.
- **Forbidden:** a sixth batch of charge intervals without optical reading; to report a
  Zernike's without specifying the pupil diameter; to treat the human
  450 nm absorption as known — din egen r4-fil har den som UNKNOWN.

## APPENDIX — the conversion you lack is already counted in the source project
`source_repository/data/corneal_transparency/corneal_transparency_results.json` (readable;
skrivskyddad) bearing the exact relationship disturbed lattice order → spridning: `C_PREF` = 16,3511, a = 14,0 nm,
φ_areal = 0,28, d_hex = 50,3917 nm, n_fibril = 1,411, n_matrix = 1,365, Δn = 0,046, L = 500 000 nm,
**oordningens pris `falsifier_ratio_poisson_over_physio` = 20,2466 till 59,4856 med medel 38,1381**; and
**wavelength exponent −4,8765**. Valideringen anger median relativt fel 7,9271 % in case of major disorder and
declares that the perfect guitar gives zero response as side lobe artifact from finite patch.

Dina egna r7-standard values (d = 53 nm, r = 14,73 nm) located close d_hex and a, so the scales are
comparable. **Try your spreading tail 1,6649e-4 at a arc minute is compatible with this
Disorder ratio, and indicate the margin.** The domain is our own previous account and therefore nothing external
facit — it is input and comparison point. the lane LANE_DOMAIN_DATA_TO_CELLS bygger cellen; din uppgift
is the compliance statement and the optical reading.
