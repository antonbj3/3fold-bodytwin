# Researcher roles, multimodality and combined interventions

Mandate from Anton, 2026-09-22. This complements `../STARTUP_MESSAGE.md` and broadens its general BodyTwin track.

## Product and research goals

BodyTwin should support different types of medical researchers in posing and investigating questions that are hard to answer with their current tools and partitioned data. Connect observations, mechanisms, models and experiments so the researcher can follow how an answer arises, which assumptions it depends on and which information could change the conclusion.

the collaborator is a concrete use case with documented wishes. Explore more researcher roles and uses at the same time. The role, question and available data should guide the investigation. A shared engine may need different inputs, workflows, visualizations and exports for different researchers.

## Persona research and simulated workflows

Start with a representative sample and expand it according to findings. The following are **hypothetical researcher cases to test**, except the collaborator's documented wishes:

| Researcher role | Example question to investigate |
|---|---|
| the collaborator / biomechanical modeling | How is anatomical shape parameterized, and how can the geometry be connected to one’s own model? |
| Rehabilitation and movement research | Which alternative explanations for a movement pattern can be distinguished with video and other measurements? |
| Tissue and material research | Which composition or structural factors can explain a measured response, and which specimen distinguishes the models? |
| Physiology and systems biology | Can observations at different scales be explained by a coherent mechanism, and where does it break down? |
| Imaging and measurement technology | Which part of an observed change may be due to the instrument or reconstruction? |
| Intervention and treatment research | How can multiple simultaneous or sequential interventions be tested in a model with explicit assumptions? |
| Surgical and prosthetic research | How do anatomy, material and a procedure’s course affect each other, and what can be validated experimentally? |

For each prioritized role: formulate two or three concrete questions, describe current practices and tools, what data exist, which connections are missing and what a useful answer must contain. Support the description of actual practices with primary sources, documentation, examples or actual feedback. Mark assumed needs separately.

Simulate the workflow from question to result: the user’s inputs, refinement of the question, model selection, visible assumptions, computation, comparison and further investigation. Also test incomplete data, conflicting measurements and questions the system cannot yet decide. Review whether the answer is comprehensible for that role and whether the researcher can trace it to data and equations. Simulated people are a way to design and test scenarios; they do not replace feedback from real users.

Let this refine interface needs, for example a question in text, selection in geometry/graph, time series comparison, parameter study or export to existing software. Choose implementation according to tested tasks. Document scenarios in `notes/RESEARCHER_SCENARIOS.md`.

## Modalities and connections between data points

Map relevant routes from imaging, ordinary video/multiple cameras, 3D reconstruction and splats, motion sensors, force measurement, electrophysiology, laboratory measurements and tissue/material specimens to the model. The list is an expansion direction; each first experiment chooses a justified subset.

For each source: state what is actually measured, individual/specimen/region, time, unit, coordinate frame, instrument and protocol, calibration, error model and provenance. Distinguish a measurement from a derived value, a reconstruction, a prior and a synthetic assumption. Data from different people or specimens must not become a claimed common individual through anatomical similarity alone.

Test which combinations add information and which share the same errors or underlying source. When modalities are missing, report alternative explanations and which new observation would distinguish them. A visually convincing reconstruction needs its own tests before its distances, material or movement are used as physical measurement data.

The camera twin is a candidate for representing image formation and its uncertainties: optics, projection, camera movement, timestamps and relevant image processing stages. Map what the existing code actually does. Connect each reconstruction back to raw observation and camera assumptions. See `MECHANISM_MODALITY_SOURCES.md` for local entry points.

## Flera behandlingsfronter

Develop a research contract for **combined and sequential interventions**, including their connection to procedures, rehabilitation, mechanical influence or biological models where a basis exists. No particular treatment or patient is assigned through this task.

An intervention needs a target in the model, an explicit change of parameter, mechanism, control signal or boundary condition, a time course and a basis for that connection. Represent combination, order, duration and relevant feedback. Investigate separate and joint effects with comparisons against the baseline course and each single intervention; do not assume effects merely add.

Look for conflicting goals, alternative mechanisms and uncertainties. Distinguish association in observed data from an established interventional relationship. Simulated combination effects should be tied to the model’s assumptions and validity domain. Show which independent measurement or experiment can test the prediction. Data connection and a graph distance do not suffice as evidence for treatment effect.

## Shared graph and first day

Connect researcher questions to observable quantities, latent states, mechanisms, computation cells, interventions and validation tasks. Preserve the relationships’ meaning and evidence. A workflow node, a physical relationship and a causal hypothesis need separate types/statuses.

The first day should give a broad map of these tracks and a limited number of decisive experiments, together with the collaborator track. Feel free to choose an experiment where several roles benefit from the same building block. Register other questions and dependencies for continued work. Report concretely what could be answered, what improved and what still requires data, model development or user feedback.
