# CX-BOOKKEEP2 — book finished The swarm/reserve_worker results in the register (one row each)

Jobs: BT-AN-G14 BT-AN-G15 BT-AN-G17 BT-AN-G18 BT-AN-G19 BT-AN-G20 BT-AN-G21 BT-AN-G22 BT-AN-G23 BT-B180 BT-B181 BT-B182 BT-B183 BT-B184 BT-B185 BT-B187 BT-B188 BT-B189 BT-B190 BT-B195 BT-MAT-OBS0 BT-MAT-OBS1 BT-MAT-SUM-AUDIT1 BT-N23 BT-N58 BT-SUM-GRID1 

For each job:
- Read `results/<id>/BRIEF.md` and `results/<id>/RESULTS.md`, plus `results.json` if needed.
- Write ONE row in exactly the same format as the last rows in `notes/RESULTS_INDEX.md`: | A<n> | <id> (<The_swarm/swarm_worker> <profile if visible in tasks/lanes/bt_queue.txt>) <question in short> | key numbers with verdict (YES/NO/UNKNOWN, numbers from the results) | reservations (not reviewed, packaging deficiency, etc.) | results/<id>/ | A (PREREG) or A or G for audits |
- For audits (BT-AN-G*), give the verdict per audited id: HOLDS / WITH RESERVATIONS / FAILS.
- Numbering starts at A372. A362–A371 are reserved for the coordinator (running lanes).
- Language: Swedish, as in the register. No embellishment. Only numbers that are in the job's own files. Mark UNKNOWN when the job says so.
- Mark a packaging error when the job says data was missing from the packet.

Deliverables:
- Append the rows to the END of `notes/RESULTS_INDEX.md` (append only; change nothing else).
- Write `results/CX-BOOKKEEP2/RESULTS.md` starting with `# CX-BOOKKEEP2`, containing the list of booked ids and a short summary (≤ 15 lines) of what stands out: wins, falls, and packaging errors.
