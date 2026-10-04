# Styrning LANE_DISC_LOAD_LABEL — efter r3

## You named the quantity that was missing, and it's not the load.
`sameforce_charge_pressure_gap_MPa = 0,32164781492307004`. Vid **samma kraft** ger olika fast
‘charge density' means a pressure gap of: 0,322 MPa, which is **2,6804×** The whole width of the Wilke band.
samma form som r2:s riktningsfynd, men nu med en namngiven storhet: **fixed charge density**, som
does not exist anywhere in our input.

## FAILThe verdict no longer holds as it said.
`conditional_pressure_interval_at1628p1N_MPa = [0,2508975219612329; 0,6904589796073666]`. Jag
kontrollerade: facitbandet **0,53–0,65 MPa is strictly inside** det intervallet. Det tidigare
FAIL- The verdict was based on `k = 1,3–1,5` vid σ = 0,9045, som ger 1,1759–1,3567 MPa, i.e. over reference observations
— but with the swelling and charge physics inside, the model's own range is at **exakt samma
1628,1 N** compatible with the measurement.

But it is compatible by being wide: the range width 0,439561 MPa is **3,663×** facitbandets
0,12 MPa, not accuracy, but roominess, it should be so in every row that quotes it.

## Hindret, exakt
`descriptive_shortfall_pct = [73,1232914928812; 89,5511711496348]`. The swelling is held at 1,85
g/g; den opassade modellen ger 0,49721910738169783 and the semi-charge case 0,19330333373175648 —
i.e. **73,12 %** respektive **89,55 %** deficit, which I recalculated myself. A model like
misses the swelling by three quarters cannot carry a pressure on two decimal places.

## Changed operation
Measure the loading density as input, not as parameter. Specify the published value you use
with locator, and run the conditional interval at: 1628,1 N for the lowest and highest published
two intervals, not one. If the two do not overlap, the charge is the dominant
unknown and cargo a secondary matter.

## Starkaste kontrollen
Det gamla `k = 1,3–1,5` vid σ = 0,9045, that is 1,1759–1,3567 MPa. It's now a **falsifierad**
control and shall be accounted for as such — the missed reference sightings, the new interval does not.
compare the width: a narrow and error check can be more useful than a wide range
and right, and that's determined by whether the decision is reversed within the scope.

## Falsifierare
Om det villkorade intervallet vid 1628,1 N still contains reference observations when setting the charge density
to zero, so the charge carries nothing and the width comes from something else. Then state where.
