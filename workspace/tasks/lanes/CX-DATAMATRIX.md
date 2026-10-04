# CX-DATAMATRIX — a deterministic packet generator: dataset × measurement protocol → Space The swarm packets that PRODUCE DATA the BodyTwin models consume

## Why
- Anton: "it has to generate value, preferably several flies in one hit… it has to pass a higher filter".
- Evidence from today: The swarm adds value when measuring real data with a clear protocol (A914: JW fluoro rollback, isokinetics, CT muscles on 30 persons) and when reproducing others' results (A908: 2 real errors).
- The swarm wastes resources on self-invented builds that hit their own anchor (A914, circular) and on second-order questions.
- So: a systematic measurement matrix over our data, generated mechanically, with no LLM cost per packet.

## The filter — EVERY packet must meet all five (put a `FILTER.md` in each packet with one line per point)
1. **Named consumer:** which BodyTwin model/decision uses the output. Choose from:
   - (a) the knee chain (L1/V0, N1g k(activity), CX-KNEEMERGE);
   - (b) population priors/band (N2b, N50, bodytwin_core.population);
   - (c) null models B24/N1/N1g per activity;
   - (d) the geometry manager (shape parameter distributions, N7c);
   - (e) sparse sensors (virtual IMU → load, A321);
   - (f) what-if (distributions of anatomical parameters);
   - (g) the cartilage/contact layer (Field U384 regime indicators).
2. **Raw data in the packet:** the measurement is made on raw data (markers/GRF/EMG/implant force/CT/mesh/Biodex) copied into `inputs/`. Never on our derived results only.
3. **Several flies in one hit:** each packet produces a standardised TABLE with ≥ 5 measured quantities per unit (trial/person/bone), serving ≥ 2 of the consumers in point 1. Examples:
   - per gait trial: stance/swing times, peak |GRF| in BW, knee flexion at heel strike/peak, peak external knee adduction/flexion moment (ID), N1 and N1g predictions, and implant force where it exists;
   - per femur: length, CCD, AV, head radius, condyle width, mid-shaft cortical thickness.
4. **Non-circular check:** a built-in consistency check that does not reuse the measured value itself. Examples:
   - the ID moment against |GRF| × lever arm;
   - left/right symmetry;
   - repeatability between repeated trials;
   - total force vs the sum of components;
   - scale invariance.
   Report the check's outcome as a number.
5. **Not already measured:** check against `results/*/results.json` and `notes/RESULTS_INDEX.md` (grep for trial/person/quantity). Skip, or measure only what is missing.

## Datasets
Read `tasks/index/DATASETS.json` and `~/projects/bodytwin/docs/DATASET_MAP.md`; find the paths; use internal ones only locally or on OVH.
- Grand Challenge 1–6 (`external_media*Competition-latest.zip`): gait of all types, lunge, stairs, step-up, isokinetics, EMG, eTibia/eKnee.
  - EXCEPTION: `jw_lungef1` eTibia/knee_forces is Field's sealed F-8 facit. It must not be copied into any packet.
- OpenCap LabValidation: markers, GRF, EMG, activities per subject.
- OrthoLoad: hip/knee/spine loads per activity and patient; `scripts/msk/index_orthoload_forces.py`.
- Imperial 35 femur+tibia, VSD CT 30 (`results/CX-MUSCLE-CT2`'s extracts may be reused as raw data), Keast tibia SSM.
- Fukuchi (walking/running), weightlift_grf, OpenBiomechanics.
- NOT LHDL, not dental, not private health/soft-tissue domains.

## Build
1. `tasks/datamatrix.py`: a list of (dataset, unit, protocol). Each protocol is a small, fixed, preregistered measurement recipe: fields, units, the check, the consumer. Generate packet directories `results/BT-DM-<dataset>-<unit>-<protocol>/` with:
   - BRIEF.md (≤ 20 lines, Swedish, following the preamble; NO references to files outside the packet);
   - `inputs/` (raw data + a small reference reader + the protocol's field list);
   - FILTER.md;
   - DATA_SUFFICIENCY.md + a load test that passes before queueing.
2. Protocols (at least 8), e.g.:
   - P-GAIT-KNEE (per gait trial: the fields in point 3);
   - P-ACT-LOAD (OrthoLoad per activity: peak, impulse, duration in BW, for N1g k(activity) and the B24 band);
   - P-FEMUR-SHAPE;
   - P-TIBIA-SHAPE;
   - P-EMG-TIMING (onset/offset of quadriceps/hamstrings/gastrocnemius per trial, vs knee angle);
   - P-ISOKIN (Biodex per person: peak moment per angle and velocity, passive-corrected);
   - P-STAIR-LUNGE (non-gait activities in GC: the knee load pattern + N1g);
   - P-IMU-VIRTUAL (virtual IMU signals from markers per trial, for the sensor models).
3. At least 250 packets. A dry run on 5 per protocol with check_packets 100 % pass. Write `results/CX-DATAMATRIX/SAMPLE.md` containing 2 full BRIEF+FILTER per protocol for the coordinator's review.
4. **Do NOT queue.** Write `results/CX-DATAMATRIX/QUEUE_PROPOSAL.txt` in the format `<A|B|C> swarm BT-DM-…`, spread over A/B/C. The coordinator reviews SAMPLE.md and queues.
5. A reducer `tasks/datamatrix_reduce.py`: gathers the finished packets' tables into `results/DATAMATRIX_TABLES/<protocol>.csv` (one row per unit, with provenance), so that consumers can read them directly.

`results/CX-DATAMATRIX/RESULTS.md` starting with `# CX-DATAMATRIX`: protocols, number of packets per dataset, filter statistics. lane runner has full permissions in the workspace; the raw-data zips are read-only.
