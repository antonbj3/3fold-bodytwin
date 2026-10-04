# Filter for data/tissue_lit_refs — applies to: ALLA Lanes, without exception

Operator's instructions 2026-10-03: the private target tissue should not be used at all.
sex and swarm without work itself: don't read, don't consume, don't quote, don't count.

## Two files used NOT
Two of the ten files directly concern the excluded target tissue — 28 av 47 respektive 31 av 54 poster
affected, and the rest of the other is apparatus dosimetry; both are off-repo and are not used.
File names are not here: they carry the terms in themselves.

## Other eight files are used with line filters
Skip each entry whose text anywhere matches the block list, case insensitive, plus
`device dosimetry`. The list is outside of the repot and read via `tasks/assembly/excluded_terms.py`;
it can't read it all rejected.

Affected entries per file, measured 2026-10-03: crosslink_release_routes 17 av 93, nir_photobiomodulation_optics
62 av 236, pbm_dose_shape_wavelength 4 av 44, pbm_thermal_confound 1 av 25, stretch_duration_dose 6 av 46,
synthesis_capacity 1 av 38Two files are completely clean: `bending_vs_stretch_mode` 62 items and
`pbm_load_interaction` 26 poster.

## What remains after the filter
**549 av 671 poster.** The emphasis remains where it is needed: skin, dermis, layer stack with fascia,
tendons, cortical bone and whole blood without the anatomical context.
samlingen — laserlane and fibrillane — almost nothing is lost: the thermic file has an affected entry of 25
and bend/stretch file zero of 62.

Present in each report how many records the filter removed, so that a digit never rests on a disfellowshipped
post.
