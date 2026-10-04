# LANE_BLEEDING_VISIBLE — make the bleed SEDD, not just calculated

Resultatmapp `results/LANE_BLEEDING_VISIBLE/`.

## The question Anton asked
"can you actually see bleeding in the simulation, in digital twin context". The answer today is no, and
it is important to know exactly where the line is. `tasks/free48/sources/SURG_HEMOSTASIS/` now has one
cell giving a **source term**: 0,0245–0,0979 mL/min via the perfusion path of a 10 cm section 2 mm
deep, and 0,31019 mL/min via the vessel counting pathway. It's a number per minute, not blood moving.
No free surface, no puddle, no obscured surgical field, no clot that stops the flow locally.

## What is missing, in order, and what we already have
1. **Source term on a surface** — EXISTS since today, but its two paths differ **3–13×** and both
   rests on text book value (STANDARD). It is the weakest link and should not be hidden.
2. **Free-surface flow over the geometry** under gravity and wetting, so the blood flows and collects
   where the geometry says. Here is actual machine: motor engine GPU particle/MPM layer. That's it
   only part of the chain that does not need to be built from scratch.
3. **A law from fibrin to rheology** — MISSING. `coagulation_hemostasis.py` in the same folder gives
   fibrin over time from a reduced 15 species model, but nothing translates fibrin concentration
   to viscosity or yield point, and without that law the flow never stops by itself.
4. **What the observer sees** — hidden field fraction as a function of time, thus a quantity and not one
   image.

## Why it's worth building, and the reason is NOT realism
Obscured field is the likely driver of **re-registration**, and re-registration is the bot branch's
binding term: 150 s per registration, of which 124 s solution and 21 s resampling, and at 7,25 s per
update accommodates 1192 updates in an intervention of published average length. The calculation is binding
so not. If blood covers landmarks, registration must be redone, and **that's the only way
we found through which bleeding at all enters the robot's binding limit**. A visible pool of blood which
just looks real doesn't change any decision. A hidden field share that forces a 150-second
re-registration changes the entire branch's calculation.

## Do like this
1. Read the cell and its `results.json` first. Use its source term; don't build a new one.
2. **Determine if free-surface flow is needed at all for the question.** At 0,02–0,3 mL/min over a 10 cm section
   is the film thickness perhaps below that which obscures something. Calculate it FIRST: hidden area per minute out
   the source term, the wetting angle and the slope of the field. If the answer is that the field is not obscured by dermal
   depth, the entire rendering chain is unmotivated at that depth, and DET is the result.
3. First, if the area count says the field is obscured: specify what the free-surface solver must handle, with
   cell count and time step, and tell which of the engine's existing solvers is sufficient.
4. **Name the missing fibrin-to-rheology law** as an acquisition record with unit and how it
   would be measured, and put it in the `notes/ACQUISITION_LIST.md` format. Don't make up a law.

## Control and falsifier
- **Control:** not to simulate the bleeding at all and instead assume the field freely. The win must
  be a changed number of re-registrations per procedure, not an image.
- **Falsifier:** if hidden area per minute at dermal depth is below the threshold to hide a
  landmark, the bleed does not bind the robot, and the chain should not be built. Report it plainly.
- **Forbidden:** to render anything; to smooth the two paths 3–13× to an average value; to
  use a textbook value as a measurement.

## Delivery
`PORT.json` with hidden area per minute as a function of cut depth, the threshold for hiding a
landmark, and the number of re-registrations per procedure with and without bleeding. Plus the acquisition item for
fibrin-to-rheology.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
