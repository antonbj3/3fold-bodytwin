# CX-SWARMWAVE2 — 60 curated Space The swarm packets (A/B/C) aimed straight at the plan's priorities BT-1–BT-4 (highest value, no busywork)

Anton: "Highest value, everything thought through, no busywork… maximum progress in 8 h." The swarm is free, and Anton wants it maxed with real value. Mechanical packet sources (the factory, "next step" extraction) mostly produced non-actionable packets. You will write 40 thought-through packets instead.


## NEW LEADS since the previous wave (use them first)
- `results/CX-BOOKKEEP3/RESULTS.md` (important findings A1251–A1268, deviating reproductions, three recommended steps).
- `results/CLOUD-*` (coordinator cloud results, wave 1): every finding that can be checked cheaply with our data gets a packet.
- `results/CX-ADJUDICATE` (if done): corrected rows → consequence analyses.
- Field's F-8 (U388/U391: in vivo contact vs eTibia) — we deliver N1g and sync; packets that test BodyTwin's side (the knee angle, GRF sync, N1g for the lunge vs other activities).
- Do not duplicate BT-BW-*/BT-DS-* that already exist (ls results/).

## Read
- `~/research/PLAN_BODYTWIN_FIELD_8H_20260924.md` §10 (BodyTwin's priorities, brainstorm B1–B5) and §14 (BT-1–BT-4).
- `notes/RESULTS_INDEX.md` A359–A369 and the rows they refer to.
- `results/CX-D1PARITY` (running), `results/CX-WHATIF2`, `results/CX-SLACK`, `results/CX-PATHS`, `results/CX-JWGEOM`, `results/CX-GC-STRENGTH*`, `results/CX-STRENGTH-POP`, `results/CX-MUSCLE-CT2`, `results/CX-SPINEARM`, `results/CX-SURGERYFE`, `results/CX-IMUFORCE`, `results/CX-GEOMCERT`, `results/CX-FLUOROLINK`, `results/CX-POPBAND`.

## Packets `results/BT-BW-061` … `BT-BW-130` (queue line `<A|B|C> swarm BT-BW-1xx`, spread evenly over A/B/C; The swarm = SMALL bounded questions, ≤ 60 s of computation)
Distribution (≈ 15 per priority), and also 10 packets that independently reproduce the key numbers of the most important CX results (CX-WHATIF2, CX-STRENGTH-POP, CX-SPINEARM, CX-SURGERYFE, CX-IMUFORCE, CX-POPBAND) with their own code:
- **BT-2, the knee chain:**
  - patella tracking along the femoral groove on the GC bone surfaces (the arm maximum should land around 45°; Krevolin 2004);
  - operating range (L0/slack) from the strength curves, with the patella fixed;
  - GC fluoroscopy as the CT↔marker link (CX-FLUOROLINK's next open part).
- **BT-3, strength from composition:**
  - a DXA → Fmax-per-group law (A359) with sex/age terms;
  - its transfer to the Grand Challenge persons, if their weight/height is available;
- **BT-1/D1 follow-ups:**
  - osteotomy optimum vs uncertainty;
  - the what-if curves' gaps (the 21 points with 140/141 KKT);
  - AV intervals at active-set changes.
- **B1/B5 (bold):**
  - a minimal measurement protocol (OED over DXA/X-ray/landmarks/strength/knee angle, with numbers from the register);
  - L5/S1 disc pressure vs Wilke with an individual lever arm (CX-SPINEARM).

Every packet:
- a question with a number criterion, a null model/counter-test, "Builds on";
- a complete `inputs/` (copy exactly the needed files from results/…, ≤ 50 MB), `DATA_SUFFICIENCY.md` + a load test;
- run `python3 tasks/packetfactory.py check` or `results/CX-SWARMGEN5/check_packets.py`.

A packet that would only find out that data is missing must NOT be built. Instead, find the data (grep the datasets, results/, /media/anton/sdc1-tmp/bodytwin, /mnt/shared_data/bodytwin_work) or choose another question.

Queue only passing packets (`>> tasks/lanes/bt_queue.txt`). `results/CX-SWARMWAVE2/RESULTS.md` starting with `# CX-SWARMWAVE2`: table id → priority (BT-DS in CX-DSWAVE exists: do not duplicate) → question → why it has high value.

Rules: lane runner has full permissions in the workspace. `~/projects/bodytwin` is read-only. restricted model data/the collaborator data may go into local/OVH packets, never to the coordinator cloud or Modal. Dental is paused.
