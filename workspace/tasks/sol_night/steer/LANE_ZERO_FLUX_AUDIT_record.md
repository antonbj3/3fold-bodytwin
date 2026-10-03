# Control r15 — your folding of the imported border is the result of the round; make it a rule

## What you did, and why it was right
I fed in a swarm finding that two pinched surfaces with no sink gives q(oo) = 0 exactly. You fell its
range with **32 of 32 finite-volume-reference at maximum relative error 1,188e-12**: at OLIKA
maintained boundary pressure is the stationary flow **q_ss = −K·dP/L and not zero**, with
|q_ss| = **3,199728e-5 m/s**, so 0,19198368 mL/min over 1 cm². And you specified the correct range for
the null result: **closed finite reservoirs** gives q(oo) = 0, with 32 local ODE reference matched to
7,211e-12. Plus the unit points: K/L is m/(Pa·s), q is m/s, and tracer permeability in m/s is another
greatness. I have corrected my own accounting accordingly.

## The operation this round
1. **Make the two-regime rule a gate that other cells can run.** Every flow cell of ours has one
   boundary conditions that are either a sustained pressure ridge or a closed finite reservoir. Calculate how
   many of our flow cells that have which, and mark those where a steady-state final value is used for
   to identify a conductance — for there the identification is invalid according to your own line.
2. **And keep the sharp corollary:** a stationary final value does not identify G, but a
   time-resolved transient with known compliance can. It is a measurement protocol requirement and must be written as one
   such, not as an observation.
3. **The external reference you acquired must be counted.** The full text with DOI 10.1039/D2SM01356H gave
   128 + 32 independent local material ODE reference matched to 8,97e-13. State explicitly that the source's
   5–10 wt% agarose and five hour hold **doesn't** transfer to our tissue — you wrote it, keep it.
4. Do not proceed with the agar batch energy budget until point 1 is complete. It is an acquisition item.

## Control and falsifier
- **Control:** the imported border that you already fell. For point 1 the control is the current position, there
  the boundary condition is not declared per cell.
- **Falsifier:** if all our flowcells have already declared boundary conditions, the rule is unnecessary as
  gate — say it and count it as a reassuring result.
- **Forbidden:** to reuse source agar parameters as ours; to report a stationary value which
  a conductance identification.
