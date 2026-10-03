# Filter for data/tissue_lit_refs — applies to: ALLA Lanes, without exception

Operator's instructions 2026-10-03: the private target tissue should not be used at all.
sex and swarm without work itself: don't read, don't consume, don't quote, don't count.

## Two files used NOT
- `measurements.tunica_turnover_gap.jsonl` — 28 av 47 items concerned
- `measurements.pemf_dose_matrix.jsonl` — 31 av 54 items concerned, the rest are apparatus dosimetry

## Other eight files are used with line filters
Skip each entry whose text anywhere matches, case insensitive:
`excluded_category`, `excluded_category`, `excluded_category`, `sinusoid`, `excluded_category`/`excluded_category`, `excluded_category`, `excluded_category`,
`device dosimetry`.

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
