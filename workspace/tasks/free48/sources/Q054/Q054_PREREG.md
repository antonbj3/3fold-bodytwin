# BT-HX-Q054 — preregistration

## Hypothesis and predicted quantity

At the same applied normal load \(P\) a reduction of the reduced radius of curvature \(R^*\), a positive radial form deviation and a positive initial clearance should move the load toward smaller area and/or give higher local pressure. For a load-controlled analysis, only a homogeneous increase in clearance \(c\) should not change the pressure field; it should only increase indentation. At the same preset indentation, however, increased clearance should reduce contact area and load.

The predicted quantity is the difference between nominal and as-built:

\[
\Delta A/A_0=(A_\mathrm{asbuilt}-A_0)/A_0,
\quad
\Delta \bar p/\bar p_0=(\bar p_\mathrm{asbuilt}-\bar p_0)/\bar p_0.
\]

## Reference value

Verified primary source: Anderson AE, Ellis BJ, Maas SA, Peters CL, Weiss JA. *Validation of Finite Element Predictions of Cartilage Contact Pressure in the Human Hip Joint*. Journal of Biomechanical Engineering. 2008;130(5):051008. DOI: `10.1115/1.2953472`. Results and Figure 4 report experimental contact area `321.9–425.1 mm²`, experimental mean pressure `4.4–5.0 MPa` and peak pressure `10.0 MPa` (the film's detection limit). Source: the publication's abstract, Results and Figure 4; DOI and PMCID `PMC2840996` looked up before this run.

The source is used as a size and plausibility reference. It is not a validation of the synthetic geometry in this assignment.

## Frozen assumptions for the first run

- Model: axial elastic contact of a local conforming ball/socket, with a compressible layer and an elastic foundation.
- Geometry: \(g(r)=c+r^2/(2R^*)+e(r)\), where \(e(r)=e_a\sin(2\pi r/\lambda)\).
- \(R^*=(1/R_h-1/R_s)^{-1}\); \(R_h\) and \(R_s\) are local synthetic CAD-like radii.
- Nominal scenario: \(R_h=25.0\) mm, \(R_s=50.0\) mm, \(c=0\), \(e_a=0\), \(P=1800\) N.
- As-built scenario: \(R_h=24.75\) mm, \(R_s=50.25\) mm, \(c=0.05\) mm, \(e_a=0.04\) mm, \(\lambda=8\) mm. This is a synthetic perturbation, not measured data.
- Layer parameters: \(G=6.8\) MPa and \(\nu=0.495\) from the reference study; \(h=1.55\) mm as a midpoint of the femoral and acetabular thickness intervals reported there. \(E_c=2G(1+\nu)\), and the simplified two-layer contact modulus \(E_c^*=E_c/[2(1-\nu^2)]\) is used only to define \(k=E_c^*/h\).
- \(P=1800\) N is a synthetic load scenario chosen to place the nominal model in the source's order of magnitude; it is not a patient load.
- No weighting, friction, fluid pumping, ligaments or foreign bodies are included.

## Frysta kriterier

1. Load equilibrium: \(|\sum_i p_iA_i-P|/P<10^{-10}\).
2. The unit check should pass all declared dimensions.
3. Parabolic null case (\(c=0,e_a=0\)) should match the analytical solution for contact area, peak pressure and mean pressure within 1 percent at a finer radial grid.
4. Load-controlled placebo: only increased homogeneous clearance should give unchanged \(p(r)\), area and mean pressure within \(10^{-9}\), and increased indentation by exactly the clearance increase.
5. The nominal output should be within the source's published intervals for area and mean pressure; this is a size requirement, not external validation.
6. As-built should have smaller area and greater mean pressure than nominal in the chosen load-controlled scenario. The form deviation should give a nonidentical radial load distribution.
7. Sensitivity is run for `reduced_radius`, `layer_modulus` and `form_error_amplitude` with factors `0.5, 1.0, 1.5`.

## Countertests and error definition

- The null model is the parabolic analytical solution.
- Placebo is a homogeneous shift of the gap without change of shape or load.
- Negative pressure, nonmonotonic load–indentation, broken load equilibrium, failed unit check or analytical error count as errors.
- Synthetic CAD/shape values must never be reported as measured data.
- Underlying fluid processes, three-dimensional contact and shear do not count as validated.

## Frysning

This preregistration is written before the first run of `model.py` and `test_model.py`.
