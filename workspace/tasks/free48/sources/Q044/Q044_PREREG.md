# BT-HX-Q044 — preregistration

Status: frozen before the first numerical model run. This is a synthetic, first-principles feasibility model; it contains no internal or measured EEG data.

## Builds on

- `inputs/QUESTION.md`: Q044, K04/K06/K07, with the requested outputs electrode potential and source distortion for a frozen neural source.
- `inputs/NIGHT_PREAMBLE.md`: the separation of source, prediction, and verification, the synthetic/no-measurement rule, and the requirement to state the next resolution step.
- The interrupted session was checked first. This directory contained only `BRIEF.md`, the two input files, `ALLOW_WEB`, and `agent.log`; there was no prior model or result to extend.

## Not redone

No individual MRI/CT anatomy, no subject-specific skull thickness, no measured tissue impedances, no measured EEG, no clinical localization result, and no claim of human performance is used. The reference values below are literature anchors, not data used to fit this model.

## Question and mechanistic hypothesis

At quasi-static EEG frequencies, Maxwell's volume-conduction equation is solved in a nested spherical volume conductor. A frozen tangential current dipole is placed in the brain and its scalp potential is measured through a finite set of 10-10-like electrodes. The primary question is how much the referenced electrode vector changes under ±50% changes in skull conductivity, skull thickness, and CSF thickness, and how much a mismatched inverse model displaces the recovered source depth.

The preregistered mechanistic prediction is that the low-conductivity skull is the dominant control for the scalp topography, while CSF thickness and scalp position provide secondary geometric controls. This is a prediction about model output, not a measurement claim.

## Frozen outcomes and criteria

The primary scalar is the normalized topographical difference

`RDM = ||u_test/||u_test||_2 - u_ref/||u_ref||_2||_2`,

with `u` the electrode vector after the selected reference subtraction. The secondary scalars are the magnitude ratio `MAG = ||u_test||_2/||u_ref||_2` and the absolute source-depth error from a reference lead-field fit.

The following criterion is frozen before running the model:

1. For the ±50% skull-conductivity perturbations, the larger absolute RDM is at least `0.05`.
2. The largest skull-conductivity RDM is greater than the largest RDM from the corresponding ±50% CSF-thickness perturbation.
3. At least one skull-conductivity perturbation produces an absolute inverse source-depth error of at least `1 mm` in the fixed-source mismatch test.
4. The analytical homogeneous-sphere limiting test, finite-valued output, and unit check pass.

A failure of any criterion is reported as a failed preregistered prediction; it is not repaired by changing the threshold after seeing the output. The inverse result is a controlled synthetic sensitivity result, not a validation of an inverse solver.

## Reference value (published primary source)

Vorwerk, J., Wolters, C. H., & Baumgarten, D. (2024), *Global sensitivity of EEG source analysis to tissue conductivity uncertainties*, Frontiers in Human Neuroscience 18:1335212. DOI: `10.3389/fnhum.2024.1335212`. PMCID: `PMC10963400`.

- Table 1, **Tissue conductivity intervals (mS/m)**: standard values are skin `430 mS/m`, skull `10.0 mS/m`, CSF `1,790.0 mS/m`, gray matter `330.0 mS/m`, and white matter `140.0 mS/m`; the skull interval is `1.6–33.0 mS/m`. The values used here are the Table 1 standards, converted to S/m, with the brain represented by separate white/gray shells.
- Figure 2, **median first-order and full Sobol indices versus source depth**: the text reports a median skull-conductivity first-order index of about `60%` up to `35 mm` and a full skull-conductivity index above `80%` for depths below `30 mm`.
- Figure 7, **average localization error**: the text reports errors up to `10 mm` in sulci and below `5 mm` for superficial sources.
- Figure 9, **conductivity/localization correlations**: the text reports a skull-conductivity correlation of about `0.6`, with other tissue correlations below `0.1` in the relevant regime.

Those reported Sobol and localization values are the literature reference. They are not numerically equated to this reduced spherical model's RDM; the frozen RDM threshold above is an independent mechanistic test.

## Model-frozen setup

- Geometry: five concentric conductive compartments, with radii (m) `0.050` white matter, `0.080` gray matter, `0.083` CSF inner boundary, `0.090` skull outer boundary, and `0.095` scalp outer boundary. The inner three radii and the two layer thicknesses are explicit assumptions for a first runnable model, not anatomical measurements.
- Conductivity (S/m): white matter `0.140`, gray matter `0.330`, CSF `1.790`, skull `0.010`, scalp `0.430`, from Vorwerk et al. Table 1.
- Source: one frozen tangential finite-difference current dipole, center `20 mm` below the CSF/skull interface, pair current amplitude `1 nA`, pair separation `6 mm`; the pair moment is `6 nA m` and has fixed x orientation. The finite pair is used instead of an unbounded singular dipole and is part of the frozen synthetic setup.
- Electrodes: a deterministic 10-10-like 21-electrode set, represented on the scalp surface. Raw potential is calculated first, then the selected reference is subtracted.
- Reference solver: average-reference by default; linked bilateral mastoid reference is implemented as a placebo/robustness alternative. The inverse solver is a fixed-orientation least-squares fit over a radial candidate-depth grid, using a lead field generated with the reference parameters.
- Numerical boundary condition: quasi-static potential is zero at infinity. The spherical-harmonic transfer solution is evaluated at the scalp boundary, so no artificial grounded scalp boundary is used.

## Perturbation protocol

The controlling parameters are varied one at a time by exactly ±50%, with all other parameters held fixed. The source center and moment, electrode labels, and reference definition are held fixed. Skull thickness changes the outer skull radius while preserving the inner skull and scalp radii. CSF-thickness changes the CSF outer boundary while preserving the inner brain and the outer skull/scalp radii; the frozen source center is retained to isolate a tissue-boundary effect.

## Error definition

An error is a non-finite potential, a failed unit or analytical-limit check, a reference shift that changes referenced values by more than numerical tolerance, a malformed perturbation, or a reported result without its input parameters and source/derivation label. A model failure invalidates the sensitivity result. No measured-data claim is permitted.

## Integrity

`PREREG.sha256` is the SHA-256 digest of the exact UTF-8 bytes of this file and is generated before `model.py` is executed.
