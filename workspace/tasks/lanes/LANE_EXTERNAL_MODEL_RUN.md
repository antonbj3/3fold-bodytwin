# LANE_EXTERNAL_MODEL_RUN

Have the ONE external model run all the way through our ports and produce a number that compares to a held-out measurement. Results folder `results/LANE_EXTERNAL_MODEL_RUN/`.

## Why this lane

Anton seeded that the models that exist in biology, released from various labs, should be able to be plugged into BodyTwin frictionlessly so that more questions can be asked. Fourteen agents expanded the seed to 324 questions, and the cross against our internal gave an unambiguous answer as to where it actually sits:

**Interoperability is untested, not absent.** The door machinery is available in versions and unit types. Three external models have been run, two of which as reimplementations and one via SBML. The only true solver plug-in was built, unit tested — and rejected by our own gate at **0 % valid frames**. FEBio is executable but has never been executed in the project: two mentions, zero executions.

And it is the same form as four other finds today: the parts are there, the connection does not work. 0 of 202 robot gates carry a calculated quantity. 27 proven signs have 0 consumers. A typed `memory` field is not read by the assertion layer. A finding requirement was not set in the 709 briefs. **The only axis that makes BodyTwin vastly better is that things actually connect, and no lane has worked on it.**

This lane therefore does not run a new model as a demonstration. It takes the connection that already FAILED and turns it into a run that passes the gate, or shows why it doesn't work.

## Do this

1. **Find the connection that was rejected.** Search broadly and in few words: `python3 source_repository/scripts/tool_find.py solver dropin`, and search in `results/CX-SOLVER`, `results/CX-SOLVER2`, `results/CX-SOLVER-ROBUST`. These contain INTERNA data — read them locally, never pass them on, and don't quote internal identifiers in anything that might leave the machine. If you don't find the rejected connection: say so with the searches you ran, and choose FEBio instead, which is executable and has never been executed here.
2. **Read the gate that rejected before touching the model.** 0 % valid frames is a gate failure, not a model error. Determine which criterion was dropped and if it was dropped for the right reason. A gate that rejects a correct model is a bug in the gate; one who rejects a faulty model did his job. **Say what it was before you cook something.**
3. **Run the model end-to-end through the ports.** Input via our typed ports, output as a number with unit. No manual intermediate steps, because a link that requires a human to move a file is not a link.
4. **Compare against a held-out measurement.** Choose an anchor from `results/LANE_EXTERNAL_FACIT_HUNT/FACIT_INDEX_v2.json` where the model can actually predict the quantity. The anchor is held-out data and must never be input data.
5. **Report the cost of the coupling honestly:** how many lines of glue, how many unit conversions, what assumptions needed to be added. "Friktionsfritt" is the word of the seed and the cost is what determines if it is correct.

## Strongest control and falsifier

- **Control:** our own reimplementation of the same mechanism. The win isn't a better number — it's that a model we didn't write ran through our gates and gave a comparable number. If the reimplementation is better, that is a valid result and should be said.
- **Forger:** if no external model can run end-to-end without manual steps, "drop-in" is not a feature we have, and the delivery is the exact list of what is missing in the port contract. It's a good and important outcome — better than another reimplementation.
- **Forbidden:** to reimplement the model and call it a plug-in (it's been done twice); bypassing the gate instead of determining whether it was right; to let an anchor become input; to send internal data somewhere.

## Leverans

A run, a number, a comparison against a hold-out measurement, and the coupling's cost in lines and assumptions. `PORT.json` with what the port contract lacked. If the outcome is negative: the list of what is missing, so that the next lane can build it.

Everything PENDING_INDEPENDENT_REVIEW. No internal data leaves the machine.
