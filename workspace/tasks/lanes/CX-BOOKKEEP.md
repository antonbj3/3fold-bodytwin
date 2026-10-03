# CX-BOOKKEEP — book finished The swarm/reserve_worker results in the register (one row each)

Jobs: BT-AN-G10 BT-AN-G11 BT-AN-G12 BT-AN-G13 BT-AN-G16 BT-B143 BT-B152 BT-B153 BT-B154 BT-B161 BT-B162 BT-B163 BT-B164 BT-B165 BT-B166 BT-B167 BT-B168 BT-B169 BT-B170 BT-B171 BT-B172 BT-B173 BT-B174 BT-B175 BT-B176 BT-B177 BT-B178 BT-B179 BT-N59 

For each job:
- Read `results/<id>/BRIEF.md` and `results/<id>/RESULTS.md`, plus `results.json` if needed.
- Write ONE row in exactly the same format as the last rows in `notes/RESULTS_INDEX.md`: | A<n> | <id> (<The_swarm/swarm_worker> <profile if visible in tasks/lanes/bt_queue.txt>) <question in short> | key numbers with verdict (YES/NO/UNKNOWN, numbers from the results) | reservations (not reviewed, packaging deficiency, etc.) | results/<id>/ | A (PREREG) or A or G for audits |
- For audits (BT-AN-G*), give the verdict per audited id: HOLDS / WITH RESERVATIONS / FAILS.
- Numbering starts at A328. A324–A327 are reserved for the coordinator.
- Language: Swedish, as in the register. No embellishment. Only numbers that are in the job's own files. Mark UNKNOWN when the job says so.
- Mark a packaging error when the job says data was missing from the packet.

Deliverables:
- Append the rows to the END of `notes/RESULTS_INDEX.md` (append only; change nothing else).
- Write `results/CX-BOOKKEEP/RESULTS.md` starting with `# CX-BOOKKEEP`, containing the list of booked ids and a short summary (≤ 15 lines) of what stands out: wins, falls, and packaging errors.
