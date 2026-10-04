# Control LANE_MENISCUS_LOAD_LABEL — after r5

## You closed the κ question, and the answer was the label "top pressure"
With real nodal areas and **mean pressure** from a run FE — 107 446 heagle tetrahedra, 57 921 free
degrees of freedom, non-linear residual 8,561430045960041e−07 — closes the force:

| kondyl | area | medeltryck | produkt |
|---|---|---|---|
| 1 | 199,01469891419063 mm² | 1,2885913888072562 MPa | **256,4486 N** |
| 2 | 252,57164890113802 mm² | 1,0695118722455845 MPa | **270,1284 N** |
| | | **summa** | **526,577 N** |

Against the imposed load **500 N** (Leeds DOI `10.5518/981`, `Full_Model_Output_Data.csv`, case
K2-SG-MN-FX) it is **+5,32 %**, which is within your own area relative errors 4,80 % and 4,44 %.

Thus: `κ ≥ 1` never applied to the triple we had. The difference requires **top pressure**; the four old ones
the rows carried a peak pressure against a **sub-region** area against a **total** load. With medium pressure and whole
the node area is `p̄·A` identical to the load on that surface, and it closes. The label that was wrong was
"topptryck".

## And the odds turned back in your favor
`finest_share_error_pp = 1,2674576562528372` against the equal division's 2,566389219183513 — I counted
the quotient: FE model is now **2,0248×** better than equal division, and 5,6007× better than r4's variable
model. In r4 the tie struck you; now it doesn't.

## But the check still hits you on a measure, and it must not be dropped
The output fed mean reading `500/(209,04+264,32) = 1,05627852 MPa` has per-condyle error 8,150 %
and 5,628 % — and **lower maximum mean error than your new model**. And the contrast is still wrong
by almost half: `0,21907951656167168` towards direction `0,15`, relative error **0,46053011041114544`,
which I checked against your own numbers.

So: you win on the distribution, the control wins on the maximum mean error. Both must be in
same line each time this is quoted.

## Changed operation
Attack the contrast. It is the only magnitude where you are 46 % wrong while everything else is below 5 %.
Determine if the error is in the pressure distribution within a condyle or in the distribution between the condyles —
it is distinguishable, because the proportion is already right to 1,27 pp.

## Strongest control
The output-fed mean reading, 1,05627852 MPa. Any improvement in contrast should be scored against
its maximum mean error, not against the equal division.

## Forgers
If the contrast error does not drop below 0,15 relatively when the pressure distribution within the condyle is refined, then
carry the contrast not by the geometry but by the material model, and then the edge should say that the contrast is
material dependent and not geometry dependent.
