# Quest: fix OOM-kills in `research.slice` (the BodyTwin SEED-500 run)

You are a lane_runner `build-lane-model` agent with full mandate from the coordinator (Anton). Complete the mission without asking; choose reasonable default choices and account for them. Network, cloud and writing in the workspace is allowed. No emails, pushes, publishing or credential changes.

## Symptoms (the new, newly introduced problem)
- `dmesg` shows OOM-kill in cgroup **`research.slice`** on both cloud hosts: **UpCloud 485 kills**, **OVH 37 kills**.
- Senaste exemplet:
  `oom-kill:constraint=CONSTRAINT_MEMCG oom_memcg=/research.slice task_memcg=/research.slice/agent-BT-FW48-SEED-072-1790714163.service task=opencode`
  `Memory cgroup out of memory: Killed process 518177 (opencode) anon-rss:660676kB`
- OVH: 35 agents started but `load average: 0.06` and ~29,8 GiB free → the agents are mostly waiting in admission (rate-limit) and are cgroup-killed.
- UpCloud: 12 agents on 7,9 GiB host; `research.slice` MemoryMax = 7 GiB → severe oversubscription.

## Sannolik rotorsak (verifiera!)
1. The slot calculation in the dispatcher counts **host memory** (`MemAvailable − BT_HEADROOM_MIB) / BT_SLOT_MIB`) and `BT_HOST_CAP`, but **ignores `research.slice` MemoryMax and that each agent has `MemoryMax=1500M`**.
   - OVH: `BT_SLOT_MIB=850` (default) but each agent gets up to 1500 M → oversubscription.
   - UpCloud: `BT_SLOT_MIB=700` but the agent is up to 1500 M → oversubscription → cgroup-OOM.
2. The coordinator just promoted OVH-dispatcherns `cap` 9→24 and `swarm_active_cap`/`swarm_worker_active_cap` to 8/8 → more concurrent agents → OOM.
3. Slice Limits: **OVH** `research.slice` MemoryMax=30064771072 (28 GiB); **UpCloud** MemoryMax=7516192768 (7 GiB). Host memory: OVH 31 GiB, UpCloud 7,9 GiB.

## Task
1. **Stop OOM now.** Ensure that the number of concurrent agents never exceeds what `research.slice` MemoryMax allows, given per-agent `MemoryMax=1500M` and a headroom left. Temporarily lower the `cap`/model ceiling if necessary.
2. **Fix the slot calculation** in the dispatcher so that it reads the actual cgroup limit:
   `ssh <host> "systemctl show research.slice -p MemoryMax -p MemoryCurrent"` and per-agent memory (1500 M), instead of using only host-`MemAvailable`/850.
   Formula should be approximately: `SLOTS = min(HOST_CAP − TOTAL, (slice_MemoryMax − slice_current − headroom)/1500M)`.
3. **Apply to both hosts** and via both service units (see below), so that the same logic applies to OVH and UpCloud.
4. **Keep the SEED-500 run**: no jobs must be dropped; no new OOM-kills in `research.slice` after the fix.
5. When the fix is ​​verified: restore a sustainable concurrency (preferably higher than 9 for OVH if slice/host allows, but never over the cgroup limit).

## Key files and facts
- **Dispatcher (shared by both hosts):** `tasks/lanes/ovh_agents/bt_queue_ovh.sh`
  - rad 12: `CAP=$(cat $D/cap 2>/dev/null || echo 30)`
  - rad 31–38: `SLOTS` from remote `MemAvailable`, `BT_HEADROOM_MIB`, `BT_SLOT_MIB` and `BT_HOST_CAP`.
  - rad 63: `sudo systemd-run --quiet --collect --slice=research.slice --unit=agent-$J-... --uid=ubuntu -p MemoryMax=1500M -p MemoryHigh=1250M -p CPUQuota=100% -p RuntimeMaxSec=5400 ...`
- **Service-units (lokala):**
  - `bt-queue-ovh.service`: default (BT_HOST_CAP=42, BT_HEADROOM_MIB=4096, BT_SLOT_MIB=850).
  - `bt-queue-upcloud.service`: `Environment=BT_CLOUD_HOST=root@212.147.226.179`, `BT_CLOUD_DIR=.../lanes/upcloud_agents`, `BT_HOST_CAP=10`, `BT_HEADROOM_MIB=1536`, `BT_SLOT_MIB=700`.
- **Local cap files:** `.../lanes/ovh_agents/cap` (now 24), `.../lanes/upcloud_agents/cap` (10), `cap_swarm_worker` .
- **Model policy:** `/opt/agents/FREE_MODEL_POLICY.json` on both hosts (`swarm_active_cap`, `swarm_worker_active_cap`, `desired_host_concurrency`, backoffs).
  - **IMPORTANT:** the files in `/opt/agents` MUST remain **ubuntu-owned** (0644/0664). Write them as **ubuntu**, never as root:
    - OVH: `ssh ubuntu@51.77.110.4 python3 -` (run as ubuntu).
    - UpCloud: `ssh root@212.147.226.179 'sudo -u ubuntu python3 -'`.
    - Protected with `flock` on `/opt/agents/rate_policy.lock`. Root-owned locks used to cause massive `exit=65` crashes — NOT reintroduces it.
- **SSH-nyckel:** `local_config_path/hunt_20260923`; `UserKnownHostsFile=external_research_path`; `IdentitiesOnly=yes`, `BatchMode=yes`.
- **Hostar:** OVH `ubuntu@51.77.110.4`; UpCloud `root@212.147.226.179`.
- **Logs:** `.../lanes/ovh_agents/queue_ovh.log` and `.../lanes/upcloud_agents/queue_ovh.log` (lines like `ovhstart`, `fetched ... exit=`, `start failed`).
- **Ramp:** `external_research_path` (controls `swarm_active_cap`; CEILING 42 OVH / 10 UpCloud; RESERVE 4096/1300 MiB; halves on rate-events, +3/+1 when productive and memory ≥ RESERVE).
- **Space guard:** `external_research_path` (+ `space-guard.timer`) keep `external_mount` ≥12 G.

## Context: The SEED-500 run
- 500 seed-derived mutations are first in queue: jobs `BT-FW48-SEED-001..500`, directories in `external_mount`.
- Gender: `.../bodytwin/tasks/lanes/bt_queue.txt`. Results are collected in `.../bodytwin/results/<J>/RESULTS.md`.
- Most non-ready jobs die from **provider-rate-limits** (`exit=75`: `explicit provider rate limit` , `admission wait expired` ). The free models (swarm/swarm_worker) are per-IP rate-limited; swarm_worker is **paused** (don't touch it). Egress is spread over OVH-direct + home IP via SSH tunnel (`EGRESS_POOL.json`).
- Done = `RESULTS.md` exists (hot or in archive via symlink).

## Important to NOT do
- Never change ownership of `/opt/agents` files to root.
- Don't touch swarm_worker-kvoten/konfigurationen.
- Delete no jobs or results; don't touch credentials.
- Do not add more agents than cgroup/host can handle.

## Definition of done
1. Slot calculation reads `research.slice` MemoryMax + per-agent memory and is respected on both hosts.
2. No new OOM-kills in `research.slice` for at least ~15 min while the queue is running (verify with `dmesg`/`journalctl -k` on both hosts).
3. Concurrency is reset to a sustainable value (higher than 9 for OVH if space exists).
4. SEED-jobben continues; new `RESULTS.md` are added and no jobs are lost.
5. Brief written summary of what was changed and measured effect, save in `tasks/lanes/LANE_RUNNER_FIX_OOM_SLOTS_RESULT.md`.

## Correction by coordinator (important — max, don't limit unnecessarily)
- **OVH has PLENTY of space.** `research.slice` MemoryMax=28 GiB, MemoryCurrent ~2 GiB, host 31 GiB, 16 cores, load ~0,06. OVH shall drive **many more** agents (~16–18 at 1500 M/agent), no less. `lanes/ovh_agents/cap` may be 24 or higher.
- **Only UpCloud is crowded:** slice 7 GiB, host 7,9 GiB → **max ~2–3 agents**; set its slots accordingly.
- **The local machine is a third execution machine and should be used:** 20 cores, ~20 GiB free. Driver is `tasks/lanes/bt_queue.sh` (not `bt_queue_ovh.sh` ), cap `tasks/lanes/bt_queue.max` (now 6 — raise reasonably, eg 10–12), gates `MemAvailable > 12 GB` and `mem-PSI full avg10 < 10` . Local jobs run as `nice 19 opencode` and are **not** in `research.slice` .
- Conclusion: do the slot calculation **per machine** (OVH large, UpCloud small, local own gate) so that all free capacity is used without any cgroup/OOM:ar.
