# LANE_MENISCUS_LOAD_LABEL — which load the meniscus pressure values carry

Results directory `results/LANE_MENISCUS_LOAD_LABEL/`. Decision **outside the eye**.

## Measured state
Eight edges in the network are `UNKNOWN` with note `load_label_unsupported`, set by the coordinator
2026-10-03, because they carry a pressure and area but no load to which the two numbers can be
calculated back. The sources are `MECHANISM_MENISCUS_LOAD*` in the readable repo.

I calculated the products myself today (2026-10-04):

| source | area | peak pressure | product | against stated load 1000 N |
|---|---|---|---|---|
| intact meniscus | 1150 mm² | 3 MPa | **3450 N** | 3,45× |
| total meniscectomy | 520 mm² | 6 MPa | **3120 N** | 3,12× |
| Rivarola, FE | 110 ± 8 mm² | 1,2 ± 0,2 MPa | **132 N** | **0,132×** |

## The obstacle, exactly
The two experimental rows give products 3,12× and 3,45× above the stated load.
That is not an error: peak pressure times area **should** exceed load, because load is
mean pressure times area. The peak/mean ratio 3,1–3,5 is physically reasonable for a contact distribution,
and the two rows are mutually consistent to 1,106.

The FE row gives 0,132×. That ratio requires a **peak pressure below mean pressure**, which is impossible.
Thus at least one of three numbers in the FE row is not what the label says: either 110 mm² is not the entire
contact area, or 1,2 MPa not the peak, or the load not 1000 N. The span between the FE row and the
experimental ones is **a factor of 26**, so they cannot all apply at the same load.

## The operation
Determine which of the three numbers in the FE row does not match its label, from the source itself, and write
an executable file in `tasks/assembly/` taking a stated axial load and **choosing** a contact pressure —
not describing the choice. The peak/mean ratio should be a quantity the decision reads, not a constant
hidden in the code.

## The strongest control
Mean pressure = load / area, without any peak/mean ratio at all. If the decision does not beat that control against
the source’s own pairs, the ratio is not load-bearing and should be removed.

## Falsifier
If the peak/mean ratio from the two experimental rows (3,12 and 3,45) does not also reproduce
Baratz’s partial row (−10 % area, +65 % peak pressure) within the source’s own spread, the ratio is not a
property of contact but an artifact of two points, and the edges should remain UNKNOWN with
that reason instead.

## Rules
Status `PENDING_INDEPENDENT_REVIEW`, no statement about biological validation. No material outside scope,
nothing naming an individual person. End with RESULTS.md, WORK_STATUS.json, NEXT_ROUND.md and
`results/LANE_MENISCUS_LOAD_LABEL/night_rounds/r<N>.json`.
