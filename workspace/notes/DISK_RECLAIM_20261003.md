# Disk reclaim 2026-10-03 — superseded intermediates in LANE_FIBRIL_NETWORK_FRACTURE/r16

Context: external_mount was at 95 % with 12 GB free and the swarm's job store lives on the same
disk, so filling it stops throughput for every session. This lane held 13.3 GB over r15-r18.

Rule applied: a large r16 file is removed ONLY if its name appears in no file anywhere in the
lane directory. 11 of 40 are cited by r16's own reports and were KEPT, because a file cited by
the round that produced it is still that round's evidence. Reports, ports and logs untouched.

| MB | file |
|---|---|
| 11 | r16/source_1600_1732_v1.npz |
| 21 | r16/dual_r11_0000_0132_v1.npz |
| 277 | r16/DUAL_R11_V1.npz |
| 26 | r16/dual_r11_1412_1572_v1.npz |
| 25 | r16/dual_r13_0128_0288_v1.npz |
| 26 | r16/dual_r11_0452_0612_v1.npz |
| 16 | r16/source_1400_1600_v1.npz |
| 26 | r16/dual_r11_0292_0452_v1.npz |
| 10 | r16/transport_r11_1160_1320_v1.npz |
| 26 | r16/dual_r11_0932_1092_v1.npz |
| 15 | r16/source_0400_0600_v1.npz |
| 15 | r16/source_0600_0800_v1.npz |
| 26 | r16/dual_r11_0772_0932_v1.npz |
| 26 | r16/dual_r12_0108_0268_v1.npz |
| 20 | r16/dual_r13_0000_0128_v1.npz |
| 26 | r16/dual_r11_1252_1412_v1.npz |
| 10 | r16/transport_r11_1000_1160_v1.npz |
| 26 | r16/dual_r11_1092_1252_v1.npz |
| 10 | r16/transport_r11_1320_1480_v1.npz |
| 11 | r16/transport_r11_1480_1640_v1.npz |
| 15 | r16/source_1000_1200_v1.npz |
| 26 | r16/dual_r11_0612_0772_v1.npz |
| 15 | r16/source_0200_0400_v1.npz |
| 15 | r16/source_0000_0200_v1.npz |
| 17 | r16/dual_r12_0000_0108_v1.npz |
| 15 | r16/source_0800_1000_v1.npz |
| 15 | r16/source_1200_1400_v1.npz |
| 26 | r16/dual_r11_0132_0292_v1.npz |
| 26 | r16/dual_r12_0268_0428_v1.npz |

**Total: 29 files, 0.80 GB.** All are .npz working arrays regenerable from the lane's own scripts.
Status: PENDING_INDEPENDENT_REVIEW; the lane's RESULTS.md, ports and round JSON are unchanged.
