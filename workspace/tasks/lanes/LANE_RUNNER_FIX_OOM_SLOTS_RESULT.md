# OOM/slots — completed 2026-09-29

Verified **15.5 minutes without new OOM events** on OVH and UpCloud while queues ran: 2026-09-29 21:10:44 UTC–2026-09-29 21:26:16 UTC, 31 samples. `journalctl -k` contained no new OOM rows; `research.slice/memory.events` was unchanged (`oom_kill`: OVH 36→36, UpCloud 97→97). `dmesg` was also checked. Historical OOM records remain and are not included in the new period.

| Machine | Slice MemoryMax | Reserve slice / host | Ceiling at start | Measured concurrency |
|---|---:|---:|---:|---:|
| OVH | 28 GiB | 2048 / 4096 MiB | 17 agents | 16 |
| UpCloud | 7 GiB | 1536 / 1536 MiB | 3 agents | 2–3 |
| Local | No research.slice | MemAvailable >12 GiB; full PSI avg10 <10 | bt_queue.max=12 | 10–13 BT jobs including another local queue |

The previous session's memory reservation was completed with a shared host check in `oom_slots_runtime/launch_reserved.py` (installed as ubuntu in `/opt/agents`). It reads `systemctl show research.slice` MemoryMax/MemoryCurrent and each direct child cgroup's actual `memory.max`. Reservation is the sum of their full limits, at least current usage, plus remaining memory directly in the slice. New slots are the nonnegative minimum of:

- `(slice_max − reserverat − slice_reserv) / agent_bytes`
- `(slice_max − slice_current − slice_reserv) / agent_bytes`
- `(MemAvailable − host_reserv) / agent_bytes`
- `host_cap − total number of agents`

`BT_AGENT_MEMORY_MIB=1500` is also used at `systemd-run`. An ubuntu-owned `capacity_admission.lock` is held across the final check and registration of each systemd unit. This prevents BodyTwin and Field from reserving the same memory simultaneously. Broken measurement, unlimited cgroup or <2 GiB free disk stops new starts. Existing cached memory makes OVH's practical concurrency about 16, despite host ceiling 17. The OVH queue's cap file was set to 24; the separate memory check governs actual global concurrency.

`bt-queue-ovh.service` and `bt-queue-upcloud.service` were restarted cleanly with different host parameters. Field services use the same check via a workspace-owned adapter and systemd drop-ins; the source project's dispatcher was not changed. Four older Field agents above the UpCloud budget were stopped with preserved packages, verified collection and interruption markers. No SEED agents were stopped.

Local `bt_queue.sh` runs as a third queue via `bt-queue-local.service`, reads the ceiling file continuously and starts `nice 19` jobs outside research.slice. `.local_claim` also protects workers waiting for admission; local and cloud runs cannot take the same job. A separate dispatcher lock does not block the queue producer's lock, and local workers survive restart. The processes' nice=19 and cgroup were verified. This queue stops new starts at 12 combined BT jobs; a separate existing local queue could temporarily increase the total to 13. Memory and PSI gates were retained.

An existing budget timer recreated the policy file as root on UpCloud. Its write command was changed to `sudo -u ubuntu` and the lock stays open during writing; budget parameters were not changed. The swarm ramp also writes as ubuntu and has host ceilings 17/3. Policy/ramp files were restored through ubuntu-written replacement, without changing contents. `/opt/agents` files and locks were checked to be ubuntu-owned. swarm_worker remains paused; quota and configuration were not changed.

**SEED preserved:** all 500 packages remain; every unfinished SEED job is in the queue. One already missing entry, SEED-253, was restored. Earlier results remain. The number of SEED `RESULTS.md` increased from **16 to 21** since the first inventory: BT-FW48-SEED-048, BT-FW48-SEED-049, BT-FW48-SEED-050, BT-FW48-SEED-061, BT-FW48-SEED-071. Result files are run outcomes, not automatically scientifically approved evidence. Provider rate limits/admission interruptions still occur.

Checks: seven regression tests passed (idle reservation, UpCloud, different cgroup limits, host/disk, over budget, unlimited cgroup and locked systemd start), Bash syntax and Python compilation passed. All five queue services were active at the last measurement.

Material: [verification_summary.json](oom_slots_runtime/verification_summary.json), [verification.jsonl](oom_slots_runtime/verification.jsonl), [test_capacity.py](oom_slots_runtime/test_capacity.py), [deployed_units](oom_slots_runtime/deployed_units/), plus dmesg/ownership logs and backups in the same directory.
