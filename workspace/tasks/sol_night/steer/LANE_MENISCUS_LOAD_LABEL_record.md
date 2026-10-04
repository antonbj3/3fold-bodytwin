# Styrning LANE_MENISCUS_LOAD_LABEL — efter r5

## You closed the κ question, and the answer was the label "peak pressure"
With real nodal areas and **mean pressure** from an executed FE — 107 446 corner tetrahedra, 57 921 free
degrees of freedom, nonlinear residual 8,561430045960041e−07 — the force closes:

| kondyl | area | medeltryck | produkt |
|---|---|---|---|
| 1 | 199,01469891419063 mm² | 1,2885913888072562 MPa | **256,4486 N** |
| 2 | 252,57164890113802 mm² | 1,0695118722455845 MPa | **270,1284 N** |
| | | **summa** | **526,577 N** |

Against the imposed load **500 N** (Leeds DOI `10.5518/981`, `Full_Model_Output_Data.csv`, case
K2-SG-MN-FX) this is **+5,32 %**, which is within your own area-relative errors 4,80 % and 4,44 %.

Thus: `κ ≥ 1` never applied to the triple we had. The inequality requires **peak pressure**; the four old
rows carried a peak pressure against a **subregion's** area against a **total** load. With mean pressure and the full
nodal area, `p̄·A` is identical to the load on that surface, and it closes. The incorrect label was
"peak pressure".

## And the share turned back in your favor
`finest_share_error_pp = 1,2674576562528372` against equal sharing's 2,566389219183513 — I calculated
the ratio: the FE model is now **2,0248×** better than equal sharing, and 5,6007× better than r4's variable
model. In r4 equal sharing beat you; now it does not.

## But the control still beats you on one measure, and that must not be lost
The output-fed mean reading `500/(209,04+264,32) = 1,05627852 MPa` has per-condyle errors 8,150 %
and 5,628 % — and **lower maximum mean error than your new model**. And the contrast is still wrong
by almost half: `0,21907951656167168` against held `0,15`, relative error **0,46053011041114544`,
which I checked against your own numbers.

Thus: you win on the distribution, the control wins on the maximum mean error. Both must appear on
the same line every time this is cited.

## Changed operation
Attack the contrast. It is the only quantity where you are 46 % wrong while everything else is below 5 %.
Determine whether the error is in the pressure distribution within a condyle or in the distribution between the condyles —
they can be distinguished, because the share is already correct to 1,27 pp.

## Strongest control
The output-fed mean reading, 1,05627852 MPa. Every improvement in the contrast must be scored against
its maximum mean error, not against equal sharing.

## Falsifier
If the contrast error does not fall below 0,15 relative when the pressure distribution within the condyle is refined, then
the contrast is supported by the material model rather than the geometry, and the edge must then say that the contrast is
material-dependent and not geometry-dependent.
