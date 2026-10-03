# Spec: get the detail layer into the graph

Written 2026-10-03 evening, for tomorrow’s push. You do not need to remember any of this — everything is
below, and the only thing you said you do not want is a narrow push.

## The state, measured and not estimated

Original repo `~/projects/bodytwin` (read-only) against the workspace `~/projects/3fold-workspaces/bodytwin`:

| lager | i originalet | i arbetsytan |
|---|---|---|
| MECHANISM detail documents | **313 unique topics** | **1** |
| python scripts under `scripts/` | 15 056, of which **9 688 carry numbers with units** | — |
| cells | — | 45 |
| data directories under `data/` | 164 | — |
| `scripts/` subdomains | 31 (msk, tissuetwin, kernel, eye, geom_recon, hand_leg, hum_retarget, physics_exp, opt, elec …) | — |

The constraint network has 65 variables and 56 edges. **Zero concern the spine**, and none of the 300
creative briefs can ask about anything in the original repo, because every generator I built tonight reads
only the workspace. That is the whole explanation for why it is not in the graph: nothing has ever pointed
there. No decision, just a blind spot.

## Why it is worth a broad push and not an inventory

Spot checks of three spine documents show completed work of the same kind that took all night to build for
the eye chain:

- **The intervertebral disc** carries external ground truth with a named falsifier: predicted pressure in the disc’s nucleus
  at peak walking ∈ [0,53, 0,65] MPa against Wilke 1999, measured in vivo in a living human.
- **The motion segment** carries an **already executed** falsification: does a model without facet joints reproduce the
  measured facet load? **FAILS.** Thus a proven case where a summary discards a load-bearing
  coordinate — exactly the question I have asked eleven times tonight in other subsystems.
- **The ligaments** carry a constitutive table: five ligaments, force in newtons per degree of angle, −90 to +90.

And the breadth: 15 documents on muscle and motor control, 13 on bone, cartilage and joints, 12 on blood, heart and vessels,
6 on fascia and tendons, 5 each on immunity, nervous system and skin. **Two on the eye** — and those are the two
I have built the night on.

## The order, and why this one

**1. Edges before briefs.** An edge whose text lacks numbers gives a worthless brief. This is measured: the 69
template-bound jobs gave zero consumable quantities while the creative ones with numbers in the brief gave 100 % with DOI
or PMID. So harvest numbers first, generate questions afterward.

**2. External ground truth before everything.** A decision without external ground truth is an opinion. Documents with a published
measurement and an executed falsification are those that can directly become a decision, and should rank ahead of
number density. The intervertebral disc and facet segment are the two known ones; the harvest seeks the rest.

**3. Scripts after documents.** 9 688 scripts carry numbers with units, but a script is an implementation and
a document is a claim. Claims become edges; implementations become cells, and a cell costs
more to verify than an edge.

**4. The generators are redirected last**, when the network carries the nodes. Then they start asking questions about the spine,
fascia and the rest on their own, which is the point of the entire chain.

## What is already underway

A Sol (`bt-sol-harvest`, lane_runner credits, not coordinator quota) reads all 313 documents and should deliver
four files in `results/SOL_MECHANISM_HARVEST/`: `INVENTORY.json` per document, `EDGES.json` with
edge proposals where the text carries the number, `EXTERNAL_FACITS.json` with each published reference, and `TOP_TEN.md`
ranked by external ground truth and executed falsification rather than number density.

## Four traps that cost most tonight, all have struck

1. **The same field name means different things in different files.** Three corpus searches failed on this. A number belongs
   to a quantity only if one can say in one line why it is the same quantity.
2. **Unit agreement is insufficient.** "801 mm" was attached to the eye’s axial length from an incision job. Right
   dimension, wrong quantity. The network now carries `plausible_range` on 14 variables; they need to exist on all that are
   harvested.
3. **A compound unit is not its last component.** `mu_pa_s` is Pa·s, not seconds — 7 511 rows
   became wrong because of this.
4. **A number chosen to hit a target is no measurement.** A cell held a drag resistance whose own
   table said "chosen to give 0.03 mm/h", and that agreement had been read as a result for weeks.

## Excluding harvest, unconditional

Three things do not get in, and it is not a filter issue but a content issue.

**1. Det excluded_category materialet — helt ute.** Not filtered from what is passed on, not rewritten, not
Counted. It is not read, consumed, not quoted and not included in any sum.
det: excluded_category, excluded_category, excluded_category, sinusoid, excluded_category, excluded_category.

**2Everything related to the collaborator — helt ute.** Named people don't appear in anything I write, and that
the material is not harvested.

**3. excluded_category — conditional, and the condition is substantive.** It gets in. ENDAST whether the quantity can be formulated
‘general' means an electromagnetic field with frequency, flow density and exposure time acting on tissue;
with the tissue identified only by tissue type. If this cannot be written without the
original application leaks through, so it stays out. Doubtful cases stay out — it is
It's cheaper to miss an edge than to have to tear you down.

the validation gate is mechanical and runs prior to administration: `tasks/assembly/exclusion_filter.py` overtakes the harvest
output and rejects each entry that matches, with a line about why. It is fail-closed: an entry that does not go
reading is rejected.

## Vad som NOT shall be made:

- Inget skrivs i `~/projects/bodytwin`It's reading mode.
- Nothing is automatically referred in the net. Suggestions, then a read, then the referral. It was automatic
  release that let in 500- the list and they 69 mallbundna jobben.
- Ingen push, ingen mejl, ingen coordinator-attribution i commits. Allt PENDING_INDEPENDENT_REVIEW.
- Internal data never leaves the machine.

Status: PENDING_INDEPENDENT_REVIEW.
