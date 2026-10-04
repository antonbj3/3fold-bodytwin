BT-HX-Q054

## Builds on
The working directory contained no partial results from a previous interrupted session. `PREREG.md` froze the hypothesis, criteria and synthetic scenarios before the run. [Source] Anderson et al. (2008), DOI `10.1115/1.2953472`, Figure 4/Results: 321,9–425,1 mm² and 4,4–5,0 MPa. These values are a magnitude reference, not validation.

## Results
[Derivation] `model.py` uses \(R^*=(1/R_h-1/R_s)^{-1}\), the gap \(g(r)=c+r^2/(2R^*)+e(r)\), \(p(r)=k\max(u-g(r),0)\), \(k=E_\mathrm{contact}/h\), and load equilibrium \(P=\int p\,dA\). Nominally at synthetic \(P=1800\) N: 360,50 mm², 4,993 MPa mean pressure, 9,961 MPa peak pressure. As-built: 346,36 mm² (-3,92 %), 5,197 MPa (+4,08 %), 10,015 MPa (+0,54 %), and 6,06 % greater indentation. The pressure field RMSE is 0,104 MPa (`results.json`).

At the same indentation as nominal, as-built gives 325,15 mm² (-9,80 %) and a load 0,887 times the nominal load. [Derivation] Sensitivity ±50 %: \(R^*\) gives area -28,6/+26,5 % and mean pressure +40,0/-21,0 %; layer modulus gives area +48,3/-19,7 % and mean pressure -32,6/+24,6 %; shape amplitude gives area ±1,36 % and mean pressure ±1,4 %. [Hypothesis] Geometry and layer stiffness govern more of the radial load distribution than this sinusoidal error profile.

## Not redone
No internal data, measured data or external musculoskeletal solver runs have been used; the values for CAD, load and error profile are explicitly synthetic assumptions.

## Validation and next steps
All frozen checks pass (`test_model.py`: 5/5; load balance ≤4e-16; analytical error <0,07 %). The model is axial and uses an elastic Winkler-layer approximation; poroelasticity, fluid, friction, shear and 2D-COP are missing. Synthetic CAD values are not measured data. Next steps are registered CAD/as-built deviation, independent load–indentation–pressure testing and then a layer/FE model with angular and shear verification.