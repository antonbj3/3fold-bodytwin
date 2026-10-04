# LANE_CONSTANT_TO_MECHANISM

Anton 1/10: "our cell simulations should be able to calculate this without textbook constants." The swarm's 43 source families (`tasks/free48/sources/`) currently carry many input constants that are actually computable from lower levels. Results folder `results/LANE_CONSTANT_TO_MECHANISM/`.

## Task

1. **Inventory**: go through the source models' parameter tables and code. Classify each numerical quantity as (a) fundamental/measured at a lower level (molecular constant, geometry, physical constant), (b) **emergent quantity entered as a constant** (e.g. an organ or cell level rate, flow, stiffness, half-life that should follow from mechanisms), (c) a purely numerical/regime parameter.
2. **Rank (b)** by leverage: how many swarm jobs/families read it (STATE.json + results directories), and how much of the output it controls (sensitivity).
3. For the 10 highest ranked: write a **mechanism specification** — which lower quantities it should be calculated from, the equations, which of them already exist in a source family (coupling instead of new construction), and an **external data anchor** (locator + number, same system) that the calculated quantity should hit without the anchor being input.
4. Build and fully test **one** of them (calculate, compare against the anchor, report ratio and uncertainty).

Delivery: CONSTANTS_INVENTORY.json, TOP10_MECHANISM_SPECS.md, and the built test. The sources are not changed by the lane. No internal data.
