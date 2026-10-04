# MITOSTRESS — first-principles research cluster

Status: proposed work, PENDING_INDEPENDENT_REVIEW. Equations below are model specifications requiring closure and units.

## M01: ATP supply and stress demand

- Constituents: ATP; ADP; Pi; oxygen; glycolysis; respiratory complexes.
- Equation skeleton: `dN_ATP/dt = J_oxphos + J_glycolysis - J_maintenance - J_response`.
- Observable: ATP flux mol/s; oxygen uptake mol/s.
- Discriminating test: Compare fixed ATP supply with demand feedback at matched total substrate; reject a benefit if compensation by glycolysis explains it.

## M02: Steroidogenesis across mitochondria and ER

- Constituents: cholesterol; StAR; CYP11A1; ER intermediates; CYP11B1; NADPH.
- Equation skeleton: `dn/dt = S v(n,enzyme,NADPH) + transport; cortisol output is a boundary flux`.
- Observable: intermediate amount mol; cortisol secretion mol/s.
- Discriminating test: Test transport-limited against enzyme-limited secretion; identify an intermediate measurement that separates them.

## M03: Catecholamine turnover and local redox

- Constituents: catecholamine pool; outer-membrane MAO; peroxide; clearance.
- Equation skeleton: `dN_CA/dt = synthesis_external - v_MAO - clearance; dN_H2O2/dt = yield*v_MAO - detox`.
- Observable: catecholamine and peroxide mol/s.
- Discriminating test: Separate gland/vesicle synthesis from mitochondrial degradation and compare with lumped hormone clearance.

## M04: HPA feedback with mitochondrial capacity

- Constituents: CRH; ACTH; free cortisol; bound cortisol; receptor; delay.
- Equation skeleton: `dx/dt = f(x(t),x(t-tau),u,capacity); secretion port from M02`.
- Observable: cortisol concentration mol/m3; lag s.
- Discriminating test: Preserve Q023 failed onset gate; determine whether added capacity is identifiable separately from feedback delay.

## M05: Tissue receptor occupancy

- Constituents: free hormone; receptor; bound receptor; trafficking; tissue exposure.
- Equation skeleton: `dB/dt = kon*C*(Rtot-B) - koff*B - kint*B`.
- Observable: bound receptor mol/m3; downstream response.
- Discriminating test: Compare kinetic occupancy with instantaneous Hill response during short pulses; expression is not affinity or gain.

## M06: Adaptation, damage and recovery

- Constituents: mitochondrial capacity; damage; repair; mitophagy; repeated stress.
- Equation skeleton: `dQ/dt = biogenesis - loss(Q,D); dD/dt = injury(stress,Q)-repair(D)`.
- Observable: capacity mol ATP/s; recovery time s.
- Discriminating test: Find regimes where repeated stress improves, exhausts or leaves capacity unchanged; do not hard-code acute good/chronic bad.

## M07: Crista geometry and substrate limitation

- Constituents: crista junction; ADP; ANT; ATP synthase; diffusion surface.
- Equation skeleton: `partial_t c = div(D grad c)-r; -D grad(c).n = membrane flux`.
- Observable: ATP mol/s; transport time s.
- Discriminating test: Reuse Q077 failure of 10% flux gain and segmentation uncertainty; compare equal area/volume with equal enzyme count.

## M08: Calcium and ER contact energetics

- Constituents: ER calcium; mitochondrial calcium; MCU; NCLX; membrane potential.
- Equation skeleton: `dN_Ca,m/dt = J_MCU-J_NCLX; proton/charge balance coupled to ATP`.
- Observable: calcium mol; voltage V; ATP mol/s.
- Discriminating test: Compare transient calcium-driven gain with depolarization cost under conserved calcium and matched oxygen.

## M09: Redox buffering and ROS locality

- Constituents: NADH; NADPH; glutathione; ROS; detoxification enzymes.
- Equation skeleton: `dn/dt = S v; constrain reducing-equivalent and cofactor pools`.
- Observable: redox ratios; peroxide mol/m3.
- Discriminating test: Ask whether a redox observation separates ETC limitation from antioxidant limitation at equal ATP.

## M10: FGF21, GDF15 and circulating signals

- Constituents: cell stress; secretion; blood; clearance; target receptor.
- Equation skeleton: `dN_blood/dt = sum(J_secretion)-CL*C; target-specific occupancy`.
- Observable: plasma mol/m3; secretion mol/s.
- Discriminating test: Separate release from clearance before interpreting elevated signal as mitochondrial damage.

## M11: mtDNA release and innate sensing

- Constituents: damaged mitochondria; cytosolic DNA; extracellular DNA; cGAS-STING; nucleases.
- Equation skeleton: `dN_DNA,k/dt = release_k + transport_k - degradation_k; sensing as explicit closure`.
- Observable: DNA copies/volume; signaling response.
- Discriminating test: Distinguish vesicular/free and cytosolic/circulating DNA; test alternative clearance explanation.

## M12: Immune metabolic reprogramming

- Constituents: T cells; macrophages; glucose; lactate; oxygen; mitochondria.
- Equation skeleton: `population dynamics coupled to stoichiometric ATP/carbon budgets`.
- Observable: cell count; uptake mol/s; killing cells/s.
- Discriminating test: At matched viable-cell number test metabolic limits on tumor killing versus inflammatory output.

## M13: Circadian and lifestyle inputs

- Constituents: light/sleep timing; meal substrate; activity; HPA phase; tissue clocks.
- Equation skeleton: `phase oscillator or measured drive -> exposure -> compartment balances`.
- Observable: phase rad; energy J; hormone time series.
- Discriminating test: Separate time-of-day and feeding confounding from stress-induced mitochondrial response.

## M14: Sex, age and cell-type heterogeneity

- Constituents: receptor distribution; endocrine state; mitochondrial abundance; tissue composition.
- Equation skeleton: `hierarchical parameters by cell/tissue; observation mixture y=H(x,composition)`.
- Observable: flux per cell and per tissue mass.
- Discriminating test: Compare within-cell-type change with mixture shift; prohibit a universal female-better multiplier.

## M15: Whole-body energy claim audit

- Constituents: basal expenditure; net stress increment; duration; activity; body mass.
- Equation skeleton: `E_increment = integral(P_stress-P_matched_rest) dt`.
- Observable: J and J/kg over specified duration.
- Discriminating test: Trace the mile-running comparison to direct measurements; return unsupported if only fibroblast or inferred expenditure exists.

## M16: Best experiment and bounded resolution

- Constituents: state; observation; noise; correlated errors; candidate measurement.
- Equation skeleton: `y=h(x,theta)+noise; propagate discretization and parameter error separately`.
- Observable: prediction interval; information gain at fixed cost.
- Discriminating test: Find the cheapest independent measurement distinguishing a mitochondrial mechanism from endocrine-only baseline; preserve source covariance.

