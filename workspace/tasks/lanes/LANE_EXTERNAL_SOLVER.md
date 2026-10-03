# LANE_EXTERNAL_SOLVER — the interoperability seed, now that a real external solver exists here

Result directory `results/LANE_EXTERNAL_SOLVER/`.

## Why this lane exists now and could not before
The operator's first seed is that external biology models should be connectable so that more questions
become askable. Its blocker has been that interoperability was never actually tested: one solver
connection was rejected at 0 % valid frames, and FEBio was closed without ever being run. I checked why:
**FEBio is not installed on this machine and no finite-element solver was available at all** — the
inventory found gmsh, trimesh and scipy, which mesh and manipulate but do not solve.

That is now fixed. `scikit-fem` 12.0.2 is installed and verified by assembling a linear-elasticity
stiffness matrix on a refined tetrahedral mesh: 345 degrees of freedom, 10 383 non-zeros. It is a pure
Python finite-element solver with no system dependencies, developed independently of this project, which
is exactly what an interoperability test requires — **an independent implementation, not a second copy
of our own assumptions.**

## Do this
1. **Pick one quantity that one of our cells already computes and that a finite-element solver can also
   compute.** Tissue stiffness response is the natural candidate: the literature collection's
   bending-versus-stretch file has 62 records with the unit *ratio (bending/stretch at matched peak
   tensile strain)* and a computed floor of 0.25 at a strain of 0.1 for pure bending with the neutral
   axis at mid-wall. Our own fibril and collagen cells carry the same kind of quantity.
2. **Compute it both ways and report the deviation.** Ours, the external solver's, and the published
   record. Three numbers for one quantity. State the mesh and the element order the external result used,
   and its convergence under one refinement — an external number without its discretisation is not
   better than ours just because it came from elsewhere.
3. **The deliverable is the INTERFACE, not the agreement.** Write down exactly what had to be supplied to
   make the external solver answer our question: geometry representation, material law, boundary
   conditions, units. That list is the interoperability result, because it is what any future external
   model will have to be given. Count the items and name the ones we could not supply from our own cell.
4. **Then say which NEW questions become askable.** The seed is not about reproducing a number we have;
   it is about questions we cannot pose today. Name three concretely, with the quantity and unit each
   would deliver, and mark which of the three the installed solver can already answer.

## Control and falsifier
- **Control:** our own cell's value for the same quantity, at its current resolution. The comparison is
  an information link, not an algorithm contest — neither implementation is a facit for the other, and
  the published record is the only external reference.
- **Falsifier:** if the external solver cannot be given our question without inventing a material law or
  a boundary condition that our cell does not contain, then the interface is the missing piece rather
  than the solver, and the honest result is the list of what is missing. That would be a more useful
  outcome than agreement, and it must be reported as the headline.
- **Forbidden:** calling agreement a validation — an independent implementation of the same assumptions
  agreeing with us is a consistency check, not evidence about tissue; tuning either side to match the
  other; reporting the external number without its mesh and convergence.

## Delivery
`PORT.json` with the three numbers for one quantity, the external solver's mesh, element order and
refinement check, the itemised interface list with the unsuppliable items named, and the three new
questions with units.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
