# LANE_SUFFICIENCY_SWEEP — ask the sufficiency question systematically, because it almost never gets asked

Resultatmapp `results/LANE_SUFFICIENCY_SWEEP/`.

## The measurement that motivates the lane
On the night of 2–3 October, sufficiency fell into **six independent subsystems** using six different mechanisms:

| subsystem | the magnitude that was not enough | proof |
|---|---|---|
| hemostasis | scalar plug coverage θ(t) | 21 plugs with identical θ to 4,37e-16, 61,53 % downstream difference |
| bleed optics | flow counting | flow-identical coatings provide oppositely directed sensor bias |
| tissue history | physical number | order number required; 192 pair with identical current position, 192 counterexample |
| laser dosage | delivered energy | 1 212 of 1 212 chronology tests fall, the energy equal to 6,01e-16 |
| eye assay | measured signal | 6 pair with signal error **exact 0,0** and 69,6 % charge difference |
| zero flow | average power | 63 180 of 100 000 violations against **analytical** 63,2121 % |

The last line is the most important: the proportion matches a closed expression to 0,00032, so it is
structure and not sampling noise.

**And then I measured the distribution before making it a system statement, which changed
the conclusion.** Strict matching over **6 312 swarm reports** yields only **33 unambiguous
sufficiency precipitations, thus 0,52 %** — but they span **18 distinct families**, with nine in one
family and sex in every five others. Loose matching gave 279, thus 4,4 %, and was noise; factor 8,5 i
overestimate.

**The conclusion is therefore not that scalars are usually inadequate. It is that THE QUESTION ALMOST ALDRIG
SET — 33 reports of 6 312 — and that it falls when set.** It makes a systematic
review to the highest expected return per run identified by the night.

## Do like this
1. **Build the test as a routine, not as an argument.** For a quantity S summing a
   state: construct two states with **identical S** and measure if a downstream quantity differs
   one. Identical means to machine precision and not approximately — the six cases above have
   identity mistake between 0,0 and 6e-16, and that is what makes them impossible to explain away as
   tolerance gap.
2. **Run it over the 43 cells in `tasks/free48/sources/`.** Each cell has state variables and
   reported output. For each reported output: which summary is used, and there are two
   state with the same summary but different output? Report per cell: tried, failed, held,
   or could not be tested — and why in the last case.
3. **Measure the distribution, not the number.** Three numbers: how many summaries that FALL, how many that
   HELL, and how many could not be tested. An individual case is nothing. And specify the coverage: if
   you only reach twenty summaries out of a hundred, say twenty and which ones.
4. **For each precipitation: name the MINSTA extension that will suffice.** That's what makes the find
   useful instead of just disturbing. Historielanen found that an ordinal number plus two local ones
   mother suffices; the homework that (θ, ε) is enough, i.e. two scalars and not the whole path; the laser lane
   that **a single functional** reproduces all 1 212 rejections. The pattern so far is that
   the extension is small — two states or a functional — and if it holds generally it is
   the most useful single result of the night.
5. **And test if a fall is transferable.** If the same summary form falls in two cells with
   
   

## Strongest control and falsifier
- **Control:** not to ask the question, i.e. the current state, where it is asked in 0,52 % of the jobs. The profit is
  
  
- **Falsifier:** if sufficiency HOLDS for most summaries when systematic
  
  
  

  
  

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.
