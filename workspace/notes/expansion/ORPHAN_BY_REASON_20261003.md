# The 259-node component has no goal on purpose, not by oversight

Status: PENDING_INDEPENDENT_REVIEW. Written 2026-10-03. Nothing in the archive was modified; the
archive at `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json` is read-only to this work and was last
written 2026-08-14.

## What the component is

Inside the frozen archive there is one genuine dependency component of **259 nodes** — 255 in the
soft-tissue cluster plus 4 musculoskeletal — among 3661 components of which 3651 are singletons. Its
type mix is the shape of an unlock graph rather than a list: 108 THEORY, 63 ENGINEERING, 60 EMPIRICAL,
28 DATA-ANCHOR, with theory and engineering standing on data anchors. Ranking it returns a reach of 238
and 236 for its two hubs while `goals_reached` stays empty, because the file's single GOAL node is
itself a singleton with no dependencies in either direction.

So the structure is not missing. It is **orphaned from the goal**, and that is why a reach of 238 buys
nothing.

## Why it is not being attached

I measured the content before choosing a goal, and the content decided it. The component is not generic
soft tissue: a keyword census over all 255 nodes gives **549 mentions of one restricted target tissue**
against 282 for tendon and ligament, 213 for skin and dermis, 182 for vessels and 66 for cartilage and
bone. Two of the component's five hubs name that tissue directly.

The operator has said that material stays outside this work. Earlier the same day, 122 of 671 literature
records were filtered out for the same reason. Therefore:

1. **It is not attached under the surgical-decision goal.** That would be wrong on content — a lens or
   cylinder decision is not what those nodes stand for — and it would be a route for restricted material
   into an active goal structure.
2. **No new goal is proposed for it either.** Formulating a goal for the component is the same as taking
   the track on, and that is the operator's decision, not a coordinator's.
3. **It stays orphaned, and this file is the record of why.** An absent declaration must not read the
   same as an unexamined one. Without this note, the next session finds a well-formed 259-node component
   with no goal and attaches one mechanically — which is precisely what two sessions nearly did today.

## What was taken, and what it is worth

The **7 musculoskeletal nodes** were copied read-only into
`notes/expansion/bodytwin_msk_frame_20261003.jsonl` after three checks passed: no node depends on
anything outside its own cluster, no node's text mentions the restricted tissue, and the single contact
between the clusters is **incoming** — a restricted node stands on one musculoskeletal node, not the
reverse — so the edge lives on the restricted node's side and does not travel.

They are a **frame, not evidence**: 4 ASSUMED and 3 OPEN, none proven, all EMPIRICAL in type, and 3 of 7
carry any evidence field at all. The 28 data anchors in the component all sit on the restricted side and
none came along. So the generic part arrives without its anchors and must not be cited as support for
anything in the robot or contact track — it is a vocabulary and a dependency shape to build against.

## The rule this leaves

For the 3651 singletons the standing rule holds: a node is attached when something actually comes to
depend on it, never in advance, because a machine-added edge invents a relation nobody measured. For
this one component the rule is different — it is rankable the day it has a goal, so "archive, query but
do not rank" must not be read as a reason to leave it alone forever.
