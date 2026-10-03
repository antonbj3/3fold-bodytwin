# CHAIN_PLANNER_FEEDBACK — standing brief for every future dental planner

Written by **DENTAL-CHAIN-PLANNER-20260928** (2026-09-28). Copied automatically into
`inputs/chain_map/` of every `BT-DW48-PLAN-*` job by `tasks/bunny48/controller.py`.
Full detail: `notes/DENTAL_CHAIN_COVERAGE_20260928.json` and
`notes/DENTAL_CHAIN_GAPS_RANKED_20260928.md`.

This is a coverage map of what is MISSING, not admitted science. Every report in this
system is `PRODUCER_OUTPUTS_UNAUDITED`.

## 1. Read this before proposing anything

The purpose of the next proposal is to close a **chain interface** or to **measure**
something, not to derive a deeper version of an inherited model.

## 2. The census that motivates this (re-measure, do not trust)

```
draft mechanism nodes  428   OPEN 386 · ASSUMED 105 · REFUTED 40 · PROVEN 13
producer jobs         1414   target nodes 6 · source directories 8 · outside catalog 0
prosthesis classes      3   implant_supported 1 case · fixed partial · removable 0
binding layer           —   538 draft nodes · 6 result bindings · 1 reviewed node of 102
```

Two consequences you must respect:

* **13 PROVEN of 428.** The system has breadth of hypothesis and almost no depth of
  execution. The bottleneck is not ideas; it is running and measuring.
* **6 target nodes.** Every producer job in the window hit one of
  `DENT-MAT-RESTORATIVE-PSP`, `DENT-GEOM-UNCERTAINTY`, `DENT-MAT-TI64-LPBF`,
  `DENT-MFG-PROCESS-MODEL`, `DENT-PROC-DRILL-THERMAL`, `DENT-MAT-PDL-CONSTITUTIVE`.
  A new mechanism with a new target node is worth more than a better experiment on
  one of these six.

## 3. Six source families that now exist, all OUTSIDE the original eight

Registered in `CATALOG.json` with real source files under `tasks/bunny48/sources/`.
**All six carry `numerical_dispatch_allowed: false`** — jobs against them are routed
to `define`/`review` and are never marked experiment-ready. That is correct, not a defect.

| key | graph target | what it opens |
|---|---|---|
| `LOADCASE` | `DENT-LOAD-MAGNITUDE-SPECTRUM` | the occlusal load case (replaces derived spectra with one acquisition spec) |
| `PROCPROP` | `DENT-MFG-PROCESS-MODEL` | geometry↔**measured** mechanical response: NIST AM-Bench, clinical fractography |
| `REMOVABLE` | `DENT-ACCESS-SOFT-TISSUE-OBSTACLE` | tissue-borne support, access obstacles, mesh-quality error budget |
| `BIORESP` | `DENT-BIO-CRESTAL-BONE-REMODELING` | the biological leg, after three REFUTED drivers |
| `VALIDATE` | `DENT-VAL2-BENCH-PROGRAM` | the 18 specced but unrun pre-registered bench experiments |
| `PRIVACY` | `DENT-DATA-LICENSE-GATE` | dataset access class, consent, external validation |

A proposal whose mechanism is genuinely new should name one of these (or propose a
seventh the same way: real file, real graph node, real dispatch receipt).

## 4. Six branches that are saturated — do not submit another derivation

| branch | why it is saturated | replacement |
|---|---|---|
| inherited fatigue power law / shell eigenstrain | 73 + 141 jobs, one identity, no new input | `PROCPROP` measurement |
| generic Frost mechanostat bands | three already REFUTED | one stimulus-space test, or an empirical risk term |
| post-hoc null-space / identifiability of the six-mode shell | amplitude unmeasured; analysis complete | `PROCPROP` measurement |
| further inherited-arithmetic audits | already 54 % of producer work, 0 bound into the graph | bind the completed `REVIEW-*` results |
| derived occlusal load spectra | 80× epistemic band already admitted; answer to the load question is "no" | `LOADCASE` acquisition spec |
| contact-law re-derivation (Hertz / quasi-plastic) | both contact nodes REFUTED | one shared inverse problem, not two re-derivations |

## 5. The ten ranked gaps, one line each

1. **G01 occlusal load** — no patient-linked bite force exists; 100+ nodes inherit an 80× band.
2. **G02 process-property** — zero measured pairs; AM-Bench is public-domain and unused.
3. **G03 biological** — three REFUTED drivers, no admissible mechanism left.
4. **G04 intake gate** — 42 written proposals silently unadmitted (now fixed; confirm it stayed fixed).
5. **G05 prosthesis class** — removable unrepresented; implant is 1 decision in 1414.
6. **G06 graph hygiene** — the top two ranked gaps are uncreatable `MISSING_DRAFT`s.
7. **G07 metrology loop** — nothing measures an as-built part and updates the design.
8. **G08 review binding** — 1 independently reviewed node of 102; 371 audits bind nothing.
9. **G09 wear/corrosion/aging** — entirely unmodelled; only fatigue exists.
10. **G10 manufacturability** — a post-check, not a search constraint.

## 6. What the next automatic planner must update

1. Re-run the census in §2 from `STATE.json` and `notes/expansion/*.jsonl`. If the
   distinct-target-node count has not risen above 6, the confinement is unchanged and
   §4 is still binding.
2. Check `STATUS.json` → `intake_gate_open` is no longer a field; verify directly that
   unadmitted proposal files are 0 and that the queue is still refilling.
3. For each §4 branch, check whether its **replacement acquisition** has produced a
   result. If it has, the branch reopens; if not, the branch stays paused and you write
   the acquisition specification rather than a derivation.
4. Extend the coverage census to any new node family you create. A new node with no
   `status` is worth less than a binding on an existing one.

## 7. Hard rules for a proposal in this project

* `source_keys` 1–3 catalog keys; `target_id` must match one of them.
* Prefer a key outside the original eight whenever the mechanism is new.
* A job whose prerequisites are open is a **definition or acquisition-planning** task.
  Say so in the body. Do not write it as an experiment.
* Name the changed operator, the primary metric, the strongest matched baseline at
  equal accuracy, the downstream consumer, and the falsifier.
* Never edit a generated source graph, `current/`, `generations/` or an engine pin.
  New work goes to `notes/`, `tasks/`, `results/`.
