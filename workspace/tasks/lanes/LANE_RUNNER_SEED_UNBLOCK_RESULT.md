# LANE_RUNNER_SEED_UNBLOCK_RESULT

Completed 2026-09-30 with Anton’s full mandate. Both fixes are installed and verified. SEED jobs actually run via account D/swarm_worker; three such jobs have returned results and been retrieved. All 500 SEED IDs remain in the queue, including SEED-253. The previous 84 results remain; the latest result count is **93**.

## Marker inventory and cleanup

| Host | Inventoried markers | Live | Dead | Reused/inconsistent | Removed at installation | Latest markers / actual model processes |
|---|---:|---:|---:|---:|---:|---:|
| OVH | 233 | 15 | 218 | 0 | 219 | 15 / 15 |
| UpCloud | 1070 | 3 | 1067 | 0 | 1067 | 2 / 2 |

The inventory was performed before installation. OVH got one additional dead marker between inventory and cleanup; UpCloud had two live markers remaining at cleanup after a normal exit. Cleanup held `model_slots.lock`; code/policy installation also held `rate_policy.lock`. All remote writes and cleanup ran as **ubuntu**, UID 1000, also on UpCloud via `sudo -u ubuntu`. Backups were created as ubuntu before replacement.

The original counter already ignored PIDs that were actually missing. Therefore the 1286 removed dead files alone cannot be credited with today’s capacity blockage. PID reuse was a real code risk but **was not found** in the inventory. The measured improvement also includes fallback and SEED priority, not an isolated effect of file cleanup.

Latest check: 2026-09-30 11:16:00 CEST. All remaining markers correspond to the launcher’s actual PID, correct job/cgroup and valid start identity. **No orphaned OpenCode processes and no gaps between marker and model process** in the latest snapshot. Agent units can also contain wrappers/admission waiters that do not yet have a model slot; they do not count as ongoing model runs.

## Exact code and field changes

On **both** hosts:

- `/opt/agents/slot_markers.py`: new shared marker handling. Checks `/proc/<pid>/stat`, zombie status, cmdline (`run_profile_ovh.py`, correct job), agent cgroup, boot ID and the process’s start ticks. Legacy markers are rejected if the process’s start time is later than the marker’s reservation. Maximum age **5700 s = RuntimeMaxSec 5400 + 300 s margin**. A dead, reused, too old or inconsistent marker is removed under the slot lock. New markers are written via temporary file + atomic rename and include `start_ticks`/`boot_id`.
- `/opt/agents/route_runtime.py`: new `_active_counts`; `reserve_slot` uses the locked shared scan instead of only `os.kill(pid,0)`. Model selection accounts for egress availability before applying free preference/fill-target. Reserve can be chosen when no free option can be admitted, due to model cooldown, model cap or lack of healthy egress with available capacity. A single limited IP does not automatically send the job to swarm_worker if another healthy free route exists.
- `/opt/agents/egress_pool.py`: `_active_counts` uses the same locked cleanup. New `has_available` atomically reads egress state/health and existing cap values, without consuming a route allocation or changing egress records. Egress → slot is the only lock order; admission takes no egress lock under the slot lock.
- `/opt/agents/fallback_reserve.py`: `choose(..., free_unavailable=True)` replaces the requirement that both free providers’ global backoff must be active simultaneously. Reserve window, swarm_worker backoff and concurrency caps remain. Quota stop is also hard-limited to **at most 95 %**, even if a future policy field specifies a higher value.
- `/opt/agents/run_profile_ovh.py`: own marker is removed under the slot lock and own ephemeral runtime is cleaned up even upon setup error/SIGTERM/deadline. `run_guarded` has a `finally` that terminates only its own `start_new_session=True` process group. SIGTERM/INT/HUP give controlled exit; SIGKILL/OOM is cleaned up by the next locked scan. An internal deadline of **1575 s** gives a 25 s cleanup margin before the packet wrapper’s existing SIGKILL at 1600 s; deadline returns 124.

`FREE_MODEL_POLICY.json` was changed atomically under `rate_policy.lock`:

| Field | Before | After |
|---|---|---|
| `swarm_worker_value_only` | `true` | `false` |
| `allocation` | global both-cooldown | available free route first, D reserve when free admission is unavailable |
| `swarm_worker_reserve_authorization` | only high-value + both-cooldown | Anton’s SEED_UNBLOCK mandate, ordinary jobs and egress/cap admission included |
| `swarm_worker_reserve_account` | D | D |
| `swarm_worker_reserve_cap` OVH | 16 | **16** |
| `swarm_worker_reserve_cap` UpCloud | 3 | **3, existing stricter host cap preserved** |
| `swarm_worker_weekly_stop_pct` | 95 | **95** |
| `swarm_worker_weekly_used_pct` | 0.0, manually UNKNOWN | **0.0, manually UNKNOWN, untouched** |
| `swarm_worker_reserve_until` | existing 24h window | untouched |

The actual console value for D is unknown. No estimate or change to Anton’s manual quota value has been made. The fallback can now consume D quota faster; the 95% protection uses the manually stated value, not a new console measurement.

Locally, `tasks/lanes/ovh_agents/bt_queue_ovh.sh` was changed (shared by both drivers):

- The queue presentation uses `seed_unblock_evidence/seed_queue_order.py`: stable SEED priority, the same rows exactly once. The canonical queue file was not edited; no jobs or rows were removed. This was needed because most SEED jobs were after AUTO/PLAN jobs and therefore did not get the available host slots.
- The retry gate is checked before start but charged **after successful unit start**. Upon denied start, the own claim is released and the start pass breaks, so capacity and completed jobs are checked again. Previously a pass continued through the queue with an old slot snapshot and burned retries without real starts.
- Under `body_rate_retry.lock`, **one demonstrably unused retry each for 16 jobs** was restored in `BODY_RATE_RETRIES.json`. Only records whose `last_claim` temporally matched a denied start and lacked a subsequent successful start/active claim were changed. Executed attempt history and results were preserved; no general retry reset was done. Full before/after per job is in `failed_admission_repair.json`.

Both `bt-queue-ovh`/`bt-queue-upcloud` were cleanly stopped/started after dispatcher changes. Ongoing separate agent units were not stopped. Memory/OOM admission, host caps, slice limits and launch_reserved were not changed.

## Real fallback, actual provider run and results

For verification, free backoff was briefly set to at least `now + 120 s` under the policy lock on both hosts, twice while existing workers finished. Existing longer backoff was never shortened. The temporary additions have expired; current policy retains ordinary backoff history.

Actual automated SEED runs on OVH:

| Job | Requested model | Account / selected model | Verification |
|---|---|---|---|
| BT-FW48-SEED-493 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + retrieved result |
| BT-FW48-SEED-494 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + retrieved result |
| BT-FW48-SEED-495 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + retrieved result |

Example, actual `MODEL_ROUTE.json` for SEED-495:

```json
{"account":"D","requested":"swarm","selected":"swarm_worker","model":"opencode-go/swarm_worker-v4.1-flash","started":1790759089.4781682}
```

Process snapshot shows `opencode run --model opencode-go/swarm_worker-v4.1-flash` with these jobs and the marker’s launcher as parent. The snapshot for SEED-493/-494/-495 contains 7 completed tool calls each at measurement time. swarm_worker processes and tool calls were observed **also on UpCloud** in its ordinary AUTO jobs during the same fallback test; UpCloud also runs the existing SEED wave. Free routes resumed when verification backoff expired.

`CLOUD_COLLECTION.json` for completed jobs documents rsync checksum without differences, preserved source/inputs and matching result header. The results are produced reports, **not scientifically admitted or independently reviewed evidence**.

## Start rate and SEED integrity

Comparable 10-minute windows:

- Before: **2026-09-30 10:52:50 CEST → 2026-09-30 11:02:50 CEST**.
- After: **2026-09-30 11:02:50 CEST → 2026-09-30 11:12:50 CEST**. Snapshot covers 600 of 600 s.

| Host | All `ovhstart` before → after | SEED starts before → after | Retrieved results before → after | Retrieved exit=75 before → after |
|---|---:|---:|---:|---:|
| ovh | 10 → 14 | 0 → 14 | 9 → 13 | 0 → 0 |
| upcloud | 4 → 4 | 0 → 4 | 4 → 4 | 0 → 0 |

Total **14 → 18 starts/10 min**, 28.6% increase; SEED starts **0 → 18**. This is a short operational window, not a long-term estimate of provider capacity. Previous historical ~14 starts/results are not compared with an invented longer after period. The logs have only HH:MM:SS; the measurement reads the latest day’s tail backward to avoid double-counting the same clock time from older days.

**500/500** SEED IDs in the queue before and after, **500/500** task briefs remain, **SEED-253 remains**. All previous **84** result IDs remain. New retrieved SEED results: BT-FW48-SEED-484, BT-FW48-SEED-486, BT-FW48-SEED-493, BT-FW48-SEED-494, BT-FW48-SEED-495, BT-FW48-SEED-496, BT-FW48-SEED-497, BT-FW48-SEED-498, BT-FW48-SEED-499. No SEED results were deleted. Existing limits for real retries are preserved; this work does not promise that all 500 scientific tasks are already complete.

## Protections, command evidence and artifacts

**27 tests pass**: live/reused/boot/legacy/age/dead/inconsistent PID; permanent cleanup in both counters; egress backoff and healthy free alternative; non-high-value SEED → D/swarm_worker; 95% stop; cap 16; reserve window expiry; swarm_worker backoff; own marker cleanup; actual subprocess tests for SIGTERM/deadline, correct exit status and stopping only the own child process group while an unrelated process survives. Dispatcher also passed `bash -n`. The installed files’ SHA256 match the tested copies.

The following commands were used, with remote Python via stdin where the scripts contain the actual inventory/patch:

```sh
ssh ubuntu@51.77.110.4 'sudo -u ubuntu python3 -'
ssh root@212.147.226.179 'sudo -u ubuntu python3 -'
python3 tasks/lanes/seed_unblock_evidence/verify_patch.py
bash -n tasks/lanes/ovh_agents/bt_queue_ovh.sh
python3 tasks/lanes/seed_unblock_evidence/observe.py
systemctl --user stop bt-queue-ovh bt-queue-upcloud
systemctl --user start bt-queue-ovh bt-queue-upcloud
systemctl --user is-active bt-queue-ovh bt-queue-upcloud
```

Both dispatcher services are `active`. `/opt/agents/EGRESS_POOL.json` is byte-identical to before on both hosts; direct/home_via_ssh/warp records have not changed. Profile keys/config have not changed. No credentials, messages or external publications are included.

All evidence is in [seed_unblock_evidence](seed_unblock_evidence/):

- `ovh.inventory.before.json`, `upcloud.inventory.before.json`: marker inventory with PID/cmdline/cgroup.
- `ovh.deploy.json`, `upcloud.deploy.json`: UID, actual cleanup, policy before→after and unchanged egress hash.
- `ovh.patch.diff`, `upcloud.patch.diff`, `dispatcher.patch.diff`: exact diffs; `slot_markers.py` is the new file’s full source.
- `tests.json`, `verify_patch.py`: test outcomes and executable verification.
- `ovh.seed_fallback_proof.json`, `ovh.swarm_worker_processes.txt`, `upcloud.fallback_snapshot.json`, `upcloud.swarm_worker_processes.txt`: actual routes, markers, processes and tool calls.
- `BT-FW48-SEED-493.*`, `BT-FW48-SEED-494.*`, `BT-FW48-SEED-495.*`: retrieved MODEL_ROUTE, exit/attempts, results and collection receipts.
- `ovh.observation.json`, `upcloud.observation.json`, `metrics.json`, `seed.before.json`, `seed.after.json`: final check and measurement data.
- `failed_admission_repair.json`, `BODY_RATE_RETRIES.before_failed_admission_repair.json`: bounded restoration of unused retries.
