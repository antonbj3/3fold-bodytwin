# PROOF_LANE_TIP_FROM_SHAFT — what can be known about the tip when only the shaft is instrumented?

## Why the question looks just like the one you have already solved
You proved today for the needle: with R friction, S other resistance and C cutting force,
`F₁ − F₂ = C₁` holds **exactly if and only if** `(R₁−R₂) + (S₁−S₂) − C₂ = 0`. From two total-force records
only `C₁ + (R₁−R₂)` is identified; neither the friction difference nor its ratio is separately
identifiable. You witnessed it with 13 104 indistinguishable pairs at observation error exactly zero. And you
found the source correction: Fukushima and Naemura, *ROBOMECH Journal* 1:14, pp. 1–8 (2014),
doi `10.1186/s40648-014-0014-7`, measures total insertion force and tip force **separately in real time**
with a coaxial needle.

The same structure applies to a robot-held surgical tip, but with a tighter constraint: **a robot has
sensors in the shaft or joint, not in the tip.** The coaxial needle is a laboratory instrument.

## The measured state of the tissue side
- PMID 41522820, doi `10.1038/s41524-025-01869-y`: in soft materials tangential stress is dominated by
  **adhesion and damping**; Coulomb friction is negligible at low contact pressures.
- PMID 33445344, doi `10.1021/acsbiomaterials.8b00387`: tissue adhesion force against bipolar forceps as
  a function of **tip temperature, in newtons**.
- The constant circulating in the literature, µ = 0,4 in an ocular trocar FE study, lies 1,69–1,88 σ
  above the µ it cites (0,295 ± 0,056 static, 0,255 ± 0,086 dynamic, Urrea et al., *JMBBM*
  56:98–105, 2016). It is **not ours** — no BodyTwin model carries a friction coefficient.

## The question
Determine exactly what is identifiable about the tip's tangential load from **shaft-side measurements
alone**, when the tangential load is dominated by adhesion and depends on temperature.

Theorems, not reasoning, and one measurable quantity per condition. Specifically:

- which combination of adhesion, damping and other resistance the shaft signal determines, and
  which remains a free direction — the same form as `C₁ + (R₁−R₂)` in the needle theorem;
- whether **temperature** as an independent input breaks the degeneracy, given that adhesion is measured against
  temperature in newtons, or whether it merely changes the parametrisation;
- which **one** additional shaft-side observable suffices — speed variation, damping phase,
  oscillation at a known frequency, insertion depth — and which provably does not suffice;
- and whether there is a quantity identifiable **independently of** the tissue model, hence robust to
  the adhesion law being wrong.

## The strongest control
Coulomb with µ = 0,295 ± 0,056, hence the best measured coefficient in the literature, as the
control model. State what the shaft signal determines under it, so that the difference from the adhesion case is
a theorem and not a preference.

## Falsifier
Build the counterexample search yourself: construct pairs of tip states that give an identical shaft signal and
report the count. If you find none for a quantity you call unidentifiable, the claim is wrong.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no statement of biological validation, nothing that
names an individual person. Every published source with PMID or DOI, volume and pages.
