# Consumable findings from the 130 net-naming swarm reports

Demand: `bodytwin.tissue_constraint_net`, 91 OPEN + 43 UNKNOWN of 158 edges. All 130 reports are
`BT-NET-<edge>` jobs queued from that net, so demand maps 1:1 by job name. `PROPOSALS.json`: 11 contradicts,
1 fills_gap, 3 confirms. Nothing written to the net.

## The strongest, number first (ten edges, eight entries)

1. **9.0909 MPa** implied mean (1000 N / 110 mm²) against a stored peak of 1.2 MPa — peak 7.58× below mean, impossible on one patch. HARVEST-E0082, Rivarola FE row, PMID 41539441.
2. **28.5714 MPa** (1000 N / 35 mm²) vs peak 3.3 MPa, ratio 0.1155 — the most extreme. HARVEST-E0086. Same defect:
   **10.5263 MPa** vs 1.4 MPa, ratio 0.1330 (E0088, 95 mm²); **12.5000 MPa** vs 2.1 MPa, ratio 0.1680 (E0084, 80 mm²).
3. **2.1 MPa** tear peak sits *below* the 3.0 MPa intact peak two rows up, reversing the ordering the edge family exists to assert; an independent 1000 N cadaveric series gives 4.5 ± 1.4 MPa. HARVEST-E0085.
4. **0.4611 mm²** — the ALL cross-section implied by `linear_stiffness = 350 N` read as E·A at E = 759 MPa (PMID 1400518). As N/mm the same 350 is coherent: E·A/L = 345 N/mm at A = 50 mm², L = 110 mm. HARVEST-E0057.
5. **4.00 N/mm** — ISL `linear_stiffness = 60 N` over a 15 mm reference length, against measured 42 and 51 N/mm (PMID 26726784): 10.5–12.8× low. HARVEST-E0060.
6. **60 N** measured in the PLL at 7.1 N·m flexion (PMID 8478347) where the net holds 0.000 N, and 0.525 N as its maximum anywhere — 114× . HARVEST-E0013.
7. **4 N** instrument detection floor (PMID 8478347) against the net's 0.286 N; ISL failure force 162 N at L1/2 (PMID 31903363) is 567× the stored value. HARVEST-E0015.
8. **130 N** measured in the ALL in extension against the sweep maximum 13.014 N at +50° — 9.99×. HARVEST-E0042.

Also: **60 N** SSL at 6.9 N·m vs 6.422 N (E0016); the 4–162 N force scale as `fills_gap` for all 50 ligament edges;
three `confirms` — the two Fukubayashi rows (peak/mean 3.45, 3.12, PMID 6894212) are coherent, and LF `200 N` is
admissible as E·A (96–462 N), the *opposite* verdict from ALL and ISL, not to be averaged with them.

## What I could not verify, and why

- **No network egress** (PubMed eutils returned empty), so "the source bears the number" was closed throughout; every
  verification above is recomputation from the raw detail-layer file the edge points at, in
  `source_documents`. Published numbers are marked as quoted from the report, not confirmed.
- **HARVEST-E0096** (spectrin 10 nm vs 7.5 nm) — the finding is that the two numbers are different objects, an
  in-situ inferred anchor length against a WLC tetramer model parameter. Nothing recomputable locally.
- **HARVEST-E0092** (coronary gain) — 0.46 vs −0.20 is over 60–100 mmHg, two later series over 120–60 at 20 mmHg steps.
- **HARVEST-E0051 / E0055 / E0056 and siblings** — they refute a zero by calling +90° flexion; the raw file says
  the opposite (`MECHANISM_SPINE_LIGAMENTS.md` L144: positive = extension). Built on an inverted axis.
- **HARVEST-E0071** (SERCA/NCX) — three sources, three observables; I cannot tell whether 76.6/23.4 and 74/23/3 are
  the same partition, so I will not call it either way.
- **32 of the 130** (`BT-NET-H-E1..H-E34`) report a lane outcome key, not a tissue quantity; several carry no units.
- **Quantity identity**: one meniscus column header covers FE `peak_stress_MPa`, cadaveric peak
  intra-compartmental pressure and mean F/A; and `BT-NET-HARVEST-E0085` attributes the 1000 N triple to Zhang 2015
  while the raw evidence file attributes it to Fukubayashi & Kurosawa 1980 (PMID 6894212) — the raw file wins.
- None of the 14 target variables carries a `plausible_range`, so the range gate ran against a range stated per
  proposal. `review_state: PENDING_INDEPENDENT_REVIEW` on every entry; `exclusion_filter.py` over this directory: 5713 ADMIT, 0 non-admit.
