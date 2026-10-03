# LANE_TISSUE_BEHAVIOUR_EXPERT

The dense expert as CONTENT: how different tissues behave, tied to geometry, in a form a robot can use. Result directory `results/LANE_TISSUE_BEHAVIOUR_EXPERT/`.

## What the seed actually is, after Anton's clarification 2/10

I first read "a hundred surgical experts into a dense expert" as **machinery** — a search algorithm — and the crossing of fourteen agents then correctly answered that the machinery already exists in the graph engine, where all twelve proposed operators point to existing functions.

Anton meant **the content**: knowledge and tricks of the trade about how different tissues behave, perhaps linked to geometry, as support for medical robots later. Surgery, all forms of tissue.

Both readings stand. The machinery exists; **the content does not**, and that is this lane.

## Why the content is harder than it looks, and where our leverage is

Practitioner knowledge is by definition **expert-reported, not measured**. That is not a reason to exclude it — it is the only source for most of how tissue behaves during a procedure — but it imposes a requirement: every record must carry whether it is measured or reported, and if it is reported it must be traced to its roots.

We have the machinery for exactly that. `claim_federation` traces every source to its root sources and calculates an effective independent count, so fifty textbooks citing the same original support the claim as ONE. It has never run on surgical material, and a lane tonight found that it would be the sharpest available task in the entire material. Pinned commit `352c6d3`, branch `research/typed-throws-20261001`, and the engine is now also at `/opt/agents/graph_engine` on both cloud hosts. **Do not change it.**

And we already have the measured side to connect to, per tissue and partly per geometry: skin tear toughness 20–30 kJ/m² and anisotropy along versus across Langer's lines; an interface toughness with upper bound 17,8–80 J/m²; rete ridge wavelengths 100–300 µm with crossover at 0,148–0,368; fibril diameter 82 ± 14 nm, fiber diameter 2,23 ± 0,96 µm, curvature radius 6,56 µm, opening angle 32°; a needle force that rises 2,95× over 80× speed, of which 28–76 % may be contact age rather than speed.

## Do this

1. **Decide the record's form before collecting any.** Minimum per record: tissue, action (cut, pull, coagulate, retract, suture), **geometric condition** (orientation to fiber direction, layer thickness, curvature, distance to a boundary), expected behaviour with quantity and unit if available, source with locator, and `evidence_class` ∈ {MEASURED, PRACTITIONER_REPORTED, GUIDELINE, TEXTBOOK}. A record without a geometric condition is not usable by a robot, because the robot only knows its own geometry.
2. **Geometry is the coupling, not a side issue.** Anton's wording is "perhaps linked to geometry", and that is where the records become useful: a robot knows where it is and how it is positioned, not which tissue lies underneath. Every record should therefore be indexed by something the robot can observe or plan, not by something it must guess.
3. **Trace the reported ones to their roots.** Run `claim_federation` on at least one cluster of records that appear to have broad support, and report effective independent count against nominal. A trick of the trade that fifty sources repeat but that traces to one observation is still worth having — but it must carry its N_eff.
4. **Connect to the measured where available.** For every record: is there a measured quantity that confirms, contradicts or is irrelevant to it? A contradiction between practitioner knowledge and measurement is the most valuable thing you can find, because one of them is wrong in a context where both are used.
5. **Build small and deep, not broad and empty.** Ten records with geometric conditions, evidence class and root tracing are worth more than two hundred rows without them. Another lane today produced 24 identical follow-up proposals; breadth without content is the most common way to look productive.

## Strongest control and falsifiers

Classify the claim according to COMMON.md. This is mainly `information_link`: the control is what a robot planner can do **without** the records, not an equally informed planner. Giving the control the same knowledge guarantees a tie and is the wrong comparison.

- **Falsifier:** if every record that can be formulated with a geometric condition is already a measured quantity we have, the practitioner layer adds nothing and that must be said. If conversely no record can be tied to a geometric condition the robot can observe, the knowledge is not usable by a robot in this form and the deliverable is why.
- **Forbidden:** mixing MEASURED and PRACTITIONER_REPORTED in the same claim without labeling them; counting citations as independent support; writing a record without a geometric condition; changing the graph engine.

## Leverans

`PORT.json` with the records in the decided form, N_eff for the root-traced clusters, and a list of contradictions between practitioner knowledge and measurement. Plus a line about which of the two readings of the seed bore fruit — machinery or content.

Inga interna data, inga patientdata, inga namngivna personer. Allt PENDING_INDEPENDENT_REVIEW.
