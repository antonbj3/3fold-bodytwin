# MAP3 — video/athlete mapping across the entire ecosystem (Anton 00:50: "run the mapping right now so it gets done")

Purpose: a complete, substantiated map of all video work concerning humans in motion (athletes, sport, walking, jumping, lifting, hands, humanoids, disassembly), cameras (camera twin, calibration, reflections, rolling shutter) and video corpora — WHERE it exists, WHAT it does, and WHAT CONDITION it is in. Anton himself does not know how well the pipelines work today; that is a main question. The map will be used to plan tonight's BodyTwin work.

## Mandatory entry point (BodyTwin's own tools — use them, free grep is not enough)
From `~/projects/bodytwin` (also works against cad-to-simulation-I code):
1. `python3 scripts/tool_find.py <word>` — run with AT LEAST these words: video, pose, athlete, sprint, jump, forceplate, grf, gait, markerless, keypoint, camera, calib, reflection, rolling, flicker, radar, football, player, highlight, youtube, corpus, clip, hand, humanoid, disassembly, splat, smpl, labvalidation, pose2sim, mediapipe, rtmpose. Note hit count per word.
2. `python3 scripts/capability_index.py --task "<word>"` — gives a STATUS stamp (PROTOTYPE/…); use for condition.
3. `python3 scripts/anchor_graph_tools.py` (see `--help`) and grep in `data/MECHANISM_ANCHOR_GRAPH.json`, `data/ANCHOR_GRAPH.json`, `data/DATA_GRAPH.json`, `data/MECHANISM_DATA_GRAPH.json` for nodes (claim, status, evidence).
4. `bt_memory/` (grep) and `~/.coordinator/projects/-home-anton/memory/` + `~/.coordinator2/` memories (grep) — lessons and where data is stored.
5. Maskinell listning (kan vara under arbete): `~/projects/3fold-workspaces/bodytwin/results/MAP/video/video_scripts_all_projects.tsv`.
6. Data modes: `lsblk`; `/mnt/shared_data`, `/media/anton/*`, `_DISKROLL.md`, datasets flyttade 2026-09-05 till /dev/sda2 (kan vara omonterad — montera NOT, notera bara).

## Condition — how to assess it (without heavy runs)
For every pipeline/cell: latest mtime on code and output, whether result JSON/logs exist and what they say (numbers), whether tests exist, whether it imports modules that exist (`python3 -c "import ast"` check or `python3 -m py_compile`), whether it refers to data that exists on disk today, status stamp in capability_index, graph node status. May run `--help` or a run < 30 s CPU if obviously harmless; no GPU jobs, no downloads, `nice -n 19`, thread ceiling 2. Condition class: {RUN+RESULT (date, number), CODE EXISTS NOT RUN, BROKEN (why), DATA MISSING, UNCLEAR}.

## Output (only under `~/projects/3fold-workspaces/bodytwin/results/MAP3/<ditt-id>/`)
- `inventory.tsv`: pipeline/cell | project | path | what it does (1 line) | input (where on disk, exists?) | latest result (file, date, key number) | condition class | graph node(s) | relevance for BodyTwin (which node/which result it can feed) | found via (tool+word).
- `overlap_bodytwin.tsv`: which of the day's BodyTwin results (V0–V6, G1/G1b/G2, E1/E2, FV1–FV3, M1–M4, HD1, IM1–IM3, BT-V7; see `notes/RESULTS_INDEX.md`) duplicate or should have used something in your inventory.
- Final answer (≤ 300 words): the main pipelines with condition, what is best to build on for BodyTwin tonight, and what is broken/unknown. Do not write SUMMARY.md (the harness stops it); put the summary in the answer.

## Rules
Read-only in all source repos. Never `pgrep -f`, never `nvidia-smi -q`. No own subagents. No mail/publishing/lane_runner. No "does not exist" claims without listed search words with zero hits. Write inventory.tsv continuously.
