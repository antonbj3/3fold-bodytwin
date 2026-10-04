# Measurements that decide something — what tonight's lanes actually delivered

Written 2026-10-02 in response to Anton's question *"what did you do with the work the lanes did?"*. Fourteen lanes have written a PORT and eight have written measurement specifications. This is what that work comes down to: a ranked list of measurements, each of which **reverses a conclusion**, rather than improving a number.

Each row: what to measure, what it decides, and what happens if it is not done. Everything PENDING_INDEPENDENT_REVIEW. The numbers are verified by the coordinator against the lane's raw-data file where it says "verified".

---

## 1. Cread at needle–tissue contact — CAN BE ORDERED, n=7
**Measure:** Cread, the contact readout ratio, on porcine skin at the speeds and geometry specified by `MEASUREMENT_SPEC_R7`.
**Decides:** which of two remaining mechanisms carries needle force's speed dependence — poroelastic pore pressure with drained support, or specific frozen transferred fast rebinding. Predicted separation **1,0 versus 1,4552**, that is, 45,5 %.
**Replicates:** `nominal_direct_Cread_CV27_min_n = 7` at 27 % CV; 18 if four independent ports are required. Verified: one-sample n ≈ 4 at 95 % and 80 % power, so 7 has margin.
**Without it:** nine rounds remain with two candidates and neither can be chosen.

## 2. Mode-resolved separation work for the dermo-epidermal boundary — NO MEASUREMENT EXISTS
**Measure:** separation work in opening, shear and mixed mode on **natural** skin, with the work balance written out (force per width is not energy per area).
**Decides:** the **sign** of our age prediction for blister times. The ratio k = Γ_II/Γ_I is not uncertain but **entirely unbounded** — the same bond inventory permits k from 1 to 10⁶.
**Three derivation routes are closed:** measurement (none published exists, `negative_result: true`), geometry (corrugation, interweaving and friction all add shear dissipation and can only increase Γ_II), bond inventory (unbounded across six orders of magnitude).
**Without it:** the blister mechanism cannot be decided at all. It is the list's only row where we exhausted every alternative to measuring.

## 3. Distribution of rete-ridge geometry — LOCATORS CONFIRMED
**Measure:** wavelength and amplitude per body region and age, with section-plane bias handled (a wavelength from a 2D section is a chord through a 3D surface).
**Decides:** which side of the crossover real skin lies on. The crossovers are 0,148 / 0,227 / 0,288 / 0,334 / 0,368 for 100–300 µm.
**Status:** both candidate DOIs I sent labeled unverified **resolved and matched**.
**Without it:** the mechanism is regime-dependent in a parameter we do not know.

## 4. Actual clinical update rate and number of updates per procedure — THE ROBOT'S BINDING TERM
**Measure:** actual intraoperative cadence, and inter-image times that would make a target speed computable.
**Decides:** whether computation time is a limit at all. At 7,25 s per update, a procedure of the published mean duration fits **1192 updates** — computation therefore does not bind at clinical cadence, and what binds is unreported.
**Verified:** 144,1 min × 60 / 7,25 = 1192,5. And 1 s is not reachable: both terms 20× faster give 7,25 s; 1 s requires 145× on both.
**Without it:** the robot branch optimizes a term that is not the limit. Four rounds showed that it is not.

## 5. Joint distribution of measured GFR and tubular secretion capacity — SAME INDIVIDUALS
**Measure:** both quantities in the same subjects, with a reference method for GFR (not an estimation equation), preserving the correlation.
**Decides:** whether our only anchor that changes sign across the renal-function range does so in a **patient** or only in a fixture. The boundary 56,67 mL/min currently rests on a synthetic low-secretory fixture, and the nominal fixture does not reverse at all.
**Without it:** 27 proved signs are conditional on an unchanged patient, and 25 of 27 do not even have the axis as a port in the model (`UNTESTED_AXIS_ABSENT`).

## 6. Raw degrees of freedom and covariance between amount and activity in the same preparation — A REPORTING REQUEST
**Measure:** no new experiment. Request the numbers already computed but not printed, from supplements or the authors.
**Decides:** whether the enzyme residual exists. At 1,96 SE it is 1,0437–1,4742, but at 3 SE and Student df=3 the interval encloses **1,0**. Only 2,1 % extra symmetric error is needed for it to disappear.
**Without it:** we cannot say whether a gap exists. The main gap 2,9–9,1× is already resolved into a 2,17–7,76× cohort artifact plus 1,18× assay normalization.

## 7. Absolute membrane flux AND free species in the SAME preparation — the ratio alone is insufficient
**CORRECTED 2/10 22:12.** The entry previously read "intracellular accumulation ratio — the only non-flow-limited route". That is wrong, and the lane proved it: scaling a,b → λa,λb preserves steady state while flux scales with λ, so **a ratio is invariant under precisely the transformation that changes the answer** and cannot bound turnover alone.
**Measure:** absolute stationary membrane flux and the free species inside, in the membrane and outside, in the same preparation, with independent support for species, stoichiometry, membrane potential and pH. More expensive than the ratio, but it is what actually decides.
**Decides:** transporter capacity, which is **structurally unidentifiable from clearance**. The certificate: metformin has extraction ratio 0,698 and the capacity interval is unbounded above; TMAO 0,368 and creatinine 0,200 succeed. My hand calculation gave 0,686 and 0,362 — independent agreement within 2 %.
**And the defensible finite bound failed:** `one_SD_empirical_support_defensible = False`, D_2SD = −16,75.
**Without it:** the 243× gap can never be decided, because the quantity we measured is saturated by flow.

## 8. Fibril fracture work and interfibrillar junction work — TWO DATA LAWS
**Measure:** the full force–opening curve to failure on a single fibril, and separation work for an interfibrillar junction.
**Decides:** Γ per layer, which is null today. The net's bending slope already converges toward the published value (1,186 → 1,033 → 0,980, the two largest enclose 1,0), but absolute Γ requires these two.
**Without it:** the recruitment mechanism can only be tested as a ratio. Against the corrected target 1,2978–2,1836, 7 of 8 realizations lie below the floor, mean 1,2090, max 1,9629 inside — too weak but not ruled out.

## 9. DEJ microscopy after incision
**Measure:** histology of fresh incisions with DEJ statistics. Published histology exists but lacks the statistics.
**Decides:** which interface actually fails in an incision.

---

## What is NOT on the list, and why

- **More modeling of what we already measured.** Five lanes closed tonight on questions that were wrongly posed, not difficult. None needed a measurement — they needed the question to be reformulated.
- **More cells before the edges exist.** 0 of 202 robot ports carry a computed number; 3 of 903 catalog pairs are decidable. Dental has the evidence: their coverage push held only 44 % under independent audit.
- **The surgeon's tremor.** It is already measured: 156 µm RMS against a 3,5 µm membrane, ratio 45 and ~1 500 at the fovea's thinnest point. At that scale the limit is the hand, not skill — and that is an argument for robotics requiring no new measurement.
