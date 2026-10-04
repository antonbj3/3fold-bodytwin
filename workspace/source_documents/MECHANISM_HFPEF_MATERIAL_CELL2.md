# HFpEF cell 2 — the material residual: can titin and collagen be separated?

Cell: `data/hfpef_diastolic_v1/`. Pre-registration: `PREREG_CELL2.md` (written and committed before
any cell-2 anchor value was in hand). Code: `scripts/tissuetwin/hfpef_cell2.py`.
Evidence: `data/hfpef_diastolic_v1/cell2_evidence.json`. Anchors: `anchors/cell2_*.json` +
`anchors/cell2_anchors.md`. Cell 1: `docs/MECHANISM_HFPEF_DIASTOLIC.md`.

## The tension cell 1 created, and what the literature says about it
Cell 1 removed chamber geometry and left a passive material stiffening of **3.12x** (preserved-EF)
and **3.51x** (reduced-EF) against the same healthy reference — ratio **0.89**, i.e. within what
that anchor resolves, the two diseases stiffen the tissue *alike*. The pre-registration set that
against the field's best-known molecular story — titin shifting toward the **stiff** short isoform
in preserved-EF disease and toward the **compliant** long isoform in dilated disease — and named
three ways they could be reconciled, before looking.

## G4 — the titin dichotomy is half-refuted, counted in labs that MEASURED it

| group | table-grade verdict | labs |
|---|---|---|
| preserved-EF | **NO SHIFT** | Charleston/MUSC (table-grade); Amsterdam/Paulus agrees, citation-only |
| reduced-EF | **no table-grade record exists** | Muenster/Linke and Houston/Granzier, both citation-only |

The preserved-EF half of the textbook claim is **not supported by either lab that measured it
directly**: N2BA/N2B 0.53±0.14 (control) vs 0.59±0.14 (HFpEF), not statistically different, and
"no conspicuous changes in sarcomeric protein composition" from the second lab. The reduced-EF half
is repeated by two labs but **every one of its records is abstract-only**, so this cell does not
grade it. The mechanism those same labs point to for preserved-EF is titin **phosphorylation**, not
isoform composition — a different quantity, not captured by an isoform ratio at all.

Counting rule, applied deliberately: repetition count and measurement count are kept apart. This
claim is one of the most-cited statements in the field, which is exactly the condition under which
they diverge.

## ★ How large a shift could G4 actually have seen?

The refutation above rests on a table-grade null: Zile's N2BA/N2B is 0.53 ± 0.14 in controls against
0.50 ± 0.13 in HTN(−)HFpEF, and the paper states there were no differences between the three groups.
**A null is not a measurement until it is bounded** — "not significantly different" says a
difference was not detected, not how large one would have had to be to show up. Nobody here had
asked, and the numbers to ask with were in the record's own prose.

| group-size reading | comparison | difference | \|t\| | smallest detectable difference |
|---|---|---|---|---|
| whole cohort, 17 / 31 / 22 | control vs HTN(−)HFpEF | +0.030 | 0.73 | 0.115 = **22 %** of control |
| | control vs HTN(+)HFpEF | −0.060 | 1.33 | 0.127 = 24 % |
| assay subset, 9 / 10 / 6 | control vs HTN(−)HFpEF | +0.030 | 0.48 | 0.174 = **33 %** |
| | control vs HTN(+)HFpEF | −0.060 | 0.81 | 0.207 = 39 % |

**The null is real** — no arm reaches |t| = 2 under either reading, so G4's refutation is not an
artefact of how the paper worded it. **But it excludes only a large shift.** A change of a fifth to a
third of the control ratio is not excluded, and the textbook claim states no magnitude that could be
compared against it.

The anchor says plainly that the isoform-assay group sizes are not itemised — the whole cohort is
17/31/22 while the assay subset "overlaps the stiffness subset n = 9/10/6". Both are computed and
**the larger bound is the one quoted**; taking the tabulated cohort sizes would be choosing the
flattering denominator.

## G3 — the split cannot be formed: PASS by abstention (pre-registered)
The titin/collagen split is strain-dependent, so it can only be formed where an isolated-titin AND
an isolated-collagen stiffness number both exist **at a stated sarcomere length**, for **both**
disease groups. Coverage on disk:

| family | preserved-EF | reduced-EF |
|---|---|---|
| A titin isoform | 2 records / 2 labs / 1 table | 3 / 2 / **0 table** |
| B titin stiffness isolated | 3 / 2 / 2 table, 2 with SL | 3 / 3 / **0 table**, 1 with SL |
| C collagen composition | 3 / 3 / 1 table | 3 / 3 / 2 table |
| D collagen stiffness isolated | **1 / 1 / 1 table, 1 with SL** | **1 / 1 / 0 table, 0 with SL** |

**Family D for reduced-EF does not exist.** One number in the entire searched literature isolates
the collagen contribution at a stated sarcomere length, and it is preserved-EF only. So the
pre-registered verdict stands as written before the harvest: *the literature does not support the separation* —
recorded as a PASS by abstention, which is precisely what the brief asked for if it came to that.

## The independence contract, enforced rather than asserted
Six papers span more than one anchor family and are therefore **n_eff=1** for those pairings:

    Zile 2015        A + B + C + D   <-- all four families, same hearts, preserved-EF only
    Borbely 2005     A + B + C
    Makarenko 2004   A + B + D
    van Heerebeek 2006  B + C
    Nagueh 2004      A + B
    Borbely 2009     A + B

Zile 2015 is the extreme case the pre-registration named in advance: it is the only source with a
true with/without-titin extraction *and* an ECM-isolated stiffness number *and* an isoform ratio
*and* a collagen content, all in the same hearts — so any cross-family statement built on it is one
path, not four. That is recorded, not laundered.

## What this leaves, and what would decide it
- The one true head-to-head on collagen (van Heerebeek 2006, systolic vs diastolic heart failure)
  reports collagen **equally elevated in both** — which is the direction cell 1's residual ratio of
  0.89 would need. It is citation-grade only (the full text is behind a 403), so it corroborates
  without deciding.
- Reading the two cells together, the hypothesis the data is *consistent* with — and cannot
  establish — is the pre-registered H3: in reduced-EF the titin shift is toward compliance while
  collagen rises, in preserved-EF titin composition is unchanged (phosphorylation instead) while
  collagen rises, and the tissue-level sums land close together. Stated as a hypothesis, because
  G3 abstained and the offsetting weights are exactly what cannot be formed.
- **The single measurement that would decide it:** a with/without-titin (or decellularisation)
  passive-stiffness measurement on **reduced-EF human myocardium at a stated sarcomere length**.
  That one record fills the empty D:reduced-EF cell and makes R computable, and with it G2.

## ★ The collagen elevation IS bigger than its own noise — and the numbers were in the prose

The register's detectability inventory asks whether a group difference exceeds the dispersion behind
it. When it was built, **every anchor family this cell reads appeared to state no dispersion at
all** — five families, zero. That looked like a hard limit of the source literature.

It was not. Three records carried mean ± SD or median [IQR] **with group sizes**, written into the
`value` prose where no instrument could reach them. Structured out of the records' own text — no new
source, no invented number — they give:

| study | disease | control | ratio | n | two-sample \|t\| |
|---|---|---|---|---|---|
| Borbely | 7.5 ± 4.0 % | 3.8 ± 2.0 % | 1.974x | 12 / 8 | **2.73** |
| Mohammed | 9.6 % [6.8–13.5] | 7.1 % [5.1–9.0] | 1.352x | 124 / 104 | **4.73** |
| Hahn | 6.27 % [3.12–8.64] | 2.77 % [2.45–3.59] | 2.264x | 78 / 13 | **6.74** |

**All three resolvable, in three independent labs, by two different summary statistics.** The
collagen elevation behind G5's ceiling is real, not noise — a thing this cell asserted for many
ticks and could not show until the dispersions were lifted out of prose.

What is still absent: **Zile's record**, which supplies the 2.300× ceiling itself, states a percent
increase with no dispersion, so the ceiling's own record remains untestable.

**G4's titin verdict is a different case, and a better one.** The prose-dispersion scan finds
Zile's and Borbely's titin-isoform records carrying ± values, and four more in the titin-stiffness
family — six records whose dispersions sit in prose. That verdict is not untestable; it is
*unstructured*, and structuring it is a bounded piece of work rather than a missing measurement. The
inventory's reach across the whole project rose when IQRs and these three blocks became readable —
and still reaches under a third of anchor records. The live count is the generated
`dispersion_reach` table in `MECHANISM_HFPEF_DIASTOLIC.md`; earlier figures quoted here were three
harvests old and are not restated.

## G5 — the weaker question G3's abstention does NOT block

G3 abstained because the titin/collagen **split** cannot be formed: family D (isolated collagen
stiffness at a stated sarcomere length) is empty for reduced-EF. But a bound survives that
abstention. Under parallel elements, if collagen carries a fraction `w` of the passive stress and
its content rises by a factor `c`, the tissue stiffens by `(1-w) + w*c` — which is **bounded above
by `c`, at `w = 1`**. So a measured content ratio can be tested against the required stiffening
*without knowing `w` at all*, and the test is maximally generous to the collagen hypothesis.

| preserved-EF collagen content ratio | value | grade | lab |
|---|---|---|---|
| Zile, HTN(+)HFpEF CVF +130% | **2.300x** | table | Charleston/MUSC |
| **Hahn**, HFpEF 6.27% vs donor 2.77%, n=78 vs 13, P=0.003 | **2.264x** | table | Johns Hopkins |
| Borbely, 7.5±4.0% vs 3.8±2.0% | 1.974x | citation_only | Amsterdam |
| **Mohammed**, 9.6% vs 7.1%, n=124 vs 104, P<0.0001 | 1.352x | table | Mayo (autopsy) |
| Kasner, "significant increase", no number | — | citation_only | Berlin |

| | |
|---|---|
| required tissue stiffening (cell 1, held out) | **3.020x** |
| ceiling on tissue stiffening at w=1 | **2.300x** |
| **shortfall** | **1.313x** |

Verdict `CONTENT_ALONE_CANNOT_DELIVER_THE_RESIDUAL`, and the ceiling is the most generous reading
available: the largest ratio stated by any preserved-EF record, in any form, at table grade.

**A second harvest was run specifically to break this gate.** The brief said so outright: if one
well-measured preserved-EF record shows a content ratio above the required stiffening, G5 flips. It returned two new
records, both table-grade, from two labs not previously in the family — and **neither exceeds the
existing ceiling**. The gate is not flipped, and the attempt is on the record rather than the
conclusion being left to rest on a single unchallenged max().

What did change is how much weight the ceiling carries. It was one number from one lab; it is now
**two independent labs landing 1.6% apart** — Zile's 2.300x (Charleston/MUSC, surgical biopsy) and
Hahn's 2.264x (Johns Hopkins, endomyocardial biopsy, n=78 vs 13 donor controls, P=0.003, verified
directly against PMC7604801 Table 3 rather than taken on the harvest's word). Four labs now supply a
number; three are table-grade.

The harvest also **rejected a fabrication that would have flipped the gate**: a search synthesis
asserted a "~5-fold" HFpEF collagen increase attributed to Zile. Fetching the putative source found
no such statement, no 5-fold figure and no such attribution. Logged under "checked and rejected" —
a number that would have changed a verdict is exactly the kind that must be run to ground.

### The ceiling is a ratio, and a ratio is only as good as its denominator

Pass 2 also put absolute percentages on disk for the first time, and they show that **the control
groups disagree with each other more than some studies' disease effects**:

| study | diseased | control | control type | same-study ratio |
|---|---|---|---|---|
| Hahn | 6.27 % | **2.77 %** | unused donor hearts | 2.264x |
| Mohammed | 9.60 % | **7.10 %** | autopsy, noncardiac death | 1.352x |
| Borbely | 7.50 % | 3.80 % | — | 1.974x |

The controls span **2.563x** among themselves, and Mohammed's *control* (7.10 %) is higher than
Hahn's *diseased* myocardium (6.27 %). Absolute collagen area is not commensurable across these
designs — fixation, sampling site, post-mortem interval and morphometry threshold all move it.

That matters because G5's ceiling is a ratio. Permitting cross-study pairing raises the maximum from
**2.264x to 3.466x** (Mohammed's diseased over Hahn's control), which **clears the 3.020x that would
flip the gate** — 1 of 6 cross pairings does it. G5 uses same-study ratios only, which is the correct
discipline and the same one the splice falsifier enforced for the stiffening exponent, and this does
not propose relaxing it. But the gate's survival is now known to be a property of that discipline
rather than of the collagen literature, and that is a measured number instead of an assumption.

**The first version of this gate got the number wrong, and its own falsifier caught it.** G5's
reader accepted one phrasing -- `A ± s vs B` paired means -- so it saw Borbely's 1.974x and was
blind to Zile's `increased ~130%`, a **table-grade** preserved-EF positive worth 2.300x. An unparsed
record can only *lower* a `max()`, so the blind spot biased the ceiling toward the verdict the gate
reached. Two consequences, both recorded rather than quietly patched:

- the ceiling was understated by **1.165x** and the shortfall moves 1.569x → **1.313x**;
- the prose here previously said *"the only table-grade preserved-EF collagen record is a NULL"*.
  That was wrong. Zile carries **two arms in one record** — HTN(−)HFpEF is a null, HTN(+)HFpEF is a
  positive — and collapsing a two-arm record to one flag erased half of it. The gate now records
  `has_null_arm` and `has_positive_arm` separately.

The verdict survives all three ceilings the falsifier could construct (as-run 1.974x, all-forms
2.300x, table-grade-only 2.300x), because all three sit below the required stiffening. The number moved; the
conclusion did not. `scripts/tissuetwin/hfpef_cell2_g5_falsifier.py` keeps its own independent
reader so the two implementations can disagree.

What this does *not* touch: stiffness **per unit** collagen. Crosslinking, fibre orientation and
network percolation all raise that, and none of them is a content measurement. So the residual
needs either super-linear stiffening per unit collagen, or a titin contribution — and cell 2
refuted the titin **isoform** route for preserved-EF while cell 5 found the **phosphorylation**
route non-discriminating and never independently replicated. The budget is not closed by anything
currently measured; G5 narrows what could close it.

## What the residual's own interval does to G5

Cell 1 reports the material residual as an **interval**, 1.499 … 2.941 … 4.554, and G5 — like every
other downstream gate in this cell family — read the median and discarded the rest. The
interval-propagation falsifier (`hfpef_interval_propagation_falsifier.py`) recomputes G5 at both
ends:

| residual taken at | shortfall vs the 2.300x ceiling |
|---|---|
| low, 1.499x | **0.652x** — BELOW ONE: content would suffice |
| median, 3.020x | 1.313x |
| high, 4.554x | 1.980x |

**The verdict no longer survives the whole interval, and that changed with better cohort coverage,
not with new collagen data.** When cell 1's absolute run grew from 4x4 to 5x5 cohorts — two records
had been excluded by the spelling of a JSON key — the residual's low end fell from 2.516x to
**1.499x**, below the 2.300x ceiling. At that end the verdict reads `CONTENT_COULD_SUFFICE`.

So G5 holds **at the median and above** and fails at the bottom of its own input interval. The honest
statement is that collagen content cannot deliver the residual as measured, and that this depends on
the residual not sitting near the low end of cell 1's spread. Quoting 1.313x alone is quoting the
median of an interval as though it were the quantity — which is exactly what the
interval-propagation falsifier was built to catch, and it caught it on its own project.

## A fail-open gate of my own, caught and fixed
The first version of G4's classifier tested for a list of null phrases and **fell through to
"shift"** for anything else. Zile 2015 — whose text reads "NOT statistically different" and "There
were no differences in the N2BA/N2B titin ratio between the three groups" — was therefore counted
as *evidence for* the very dichotomy it refutes. A default that lands on the expected answer is a
fail-open gate. Both directions must now be positively evidenced and anything else is recorded as
UNCLASSIFIABLE, counted in neither arm.
