# LANE_ROBOT_PORT_CONNECT

Turn the robot branch into a graph instead of a diagram. Results directory `results/LANE_ROBOT_PORT_CONNECT/`.

## Why this lane

Anton 2/10: the robot part should start being built in our graph, the goal is to outperform everything that exists, and the threshold he set is a robot good enough to force EU regulations to be revised.

The branch now exists: **41 nodes and 202 ports**, merged from two independent research waves and validated cleanly against the graph engine's own validator (`validate` ok, 0 errors, 0 warnings). But it is not yet a graph:

- **Only 7 of 202 ports connect to a quantity the project actually calculates, and all 7 are on the tissue side.** Every robot-side port is dangling.
- **80 of 205 port rows lack a value**, and of the 122 that have one, only **13 rely on more than one source**. The typical case is an article read at abstract level.
- Exactly ONE port name is shared between the two waves (`tremor_band`). All four declared dependencies from the sensor side to the mechanics side resolve to nothing or an empty port, and one of them is requested in meters against ports specified in degrees.
- No node can be PROVEN, because no fold-ledger row exists for any claim here; all 41 are OPEN or ASSUMED.

A branch whose ports connect to nothing is a drawing. That is what should be fixed, not the node count.

## What is already decided and should not be redone

- **The binding limit is registration time, not sensing.** At 150 s per registration, acquisition is 0,002 % of the model update chain: 312 → 6417 Hz buys 3 ms out of 150 s, while 150 s → 1 s is 149,5×. The brainstorm's headline gap of 1,0e8 sits on the wrong term, and its denominator is a work cycle sold as a frequency.
- **The error budget against 1,08 mm target accuracy:** brain shift +195 %, pivot deviation 110 %, tremor 1,0 %, optical error 0,1 %. Three of four headline gaps are far from binding; brain shift versus tremor is 193× in leverage.
- **The binding MECHANICAL requirement is DEJ peel force ≤0,1 mN, not tip force.** A tip force of 0,3 mN meets the cutting requirement ≤1 mN with a 3,3× margin but misses the peel requirement by ~3×.
- **All 15 mechanics nodes are regulatorily inert.** What makes MDR rule 22 incoherent is a damage measurement, not a specification: note 1 covers external subsystems, so wire-versus-voice is the route. We already hold the human path's numbers.
- **A large part of our own measurement spec is bench-only by construction** (matched adjacent samples, before/after calibration of the same tool, destructive sectioning). A machine can neither meet nor miss such requirements, and they should not be counted as robot gaps.

## Do this

1. **Read the branch before touching it:** `external_research_path`, `PROVENANCE.md`, `BINDING.md`. Also read how the graph engine actually represents nodes, ports and validity boxes at pinned commit `352c6d3`, branch `research/typed-throws-20261001`, under `src/graph_engine/`. **Do not change the engine** — the graph lane (anton-4d) owns it.
2. **Connect the robot-side ports, one at a time, to a quantity we calculate.** For each port: name the quantity, its unit, where it is calculated, and which number it gives today. A connection requiring a quantity we do not calculate is not a connection — declare it as `REQUIRES_NEW_COMPUTATION` and say exactly what must be calculated first.
3. **Be honest about the direction.** The project's position: tissue physics and physiological response belong to the intelligence engine, and what a robot controls is its own **precision and planning**. A port making the robot responsible for tissue behavior is misplaced; move it and say what it was.
4. **Prioritize by leverage, not by how many ports are connected.** Registration time is the binding term, so a connection concerning it is worth more than twenty concerning tremor. Say which order you chose and why.
5. **Every number gets a provenance class:** INDEPENDENT_MEASUREMENT / PEER_REVIEWED / STANDARD / VENDOR / UNSOURCED. No vendor number is a measurement. Count how many ports after your round rely on more than one source — that number is the lane's actual progress.

## Strongest control and falsifiers

- **Control:** the branch as it stands, thus 41 validating nodes with 7 connected ports. The gain should be the number of robot-side ports that carry a calculated quantity after the round, not a prettier drawing.
- **Falsifier:** if no robot-side port can connect to anything we calculate, the branch is premature, and the result is a list of what we must calculate FIRST. Report it plainly — it is a fully adequate and important outcome, and better to know now than after twenty more nodes.
- **Forbidden:** adding nodes or ports (the branch should be connected, not grow); using a vendor number as a measurement; counting a bench-only requirement as a robot gap; changing the graph engine.

## Delivery

`PORT.json` with each connected port, its calculated quantity and provenance class; a `REQUIRES_NEW_COMPUTATION` list of what must be calculated first; and one number: how many ports now rely on more than one source, against 13 before.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.


## PIN UPDATED (coordinator 2/10 21:05)
Pinned commit moved from `cf9c1e9` to `352c6d3`, three commits later, seven files changed. The reason is that those three contain exactly the corrections we pursued tonight: *Abstain when the downstream cert says nothing, instead of reading silence as a pass*, *Stop a validity filter from passing a point it never checked*, and *Say whether the declaration could separate two discordant reports, and count unconsumed verdicts*. The older pin had let silence pass as approval.
The engine is now also on both cloud hosts under `/opt/agents/graph_engine`, with the quoted path symlinked there, so a job no longer needs to reimplement it. Verified: import through the quoted path works on both.
