from pathlib import Path
import datetime,json
B=Path(__file__).resolve().parent
m=json.loads((B/'metrics.json').read_text());before=json.loads((B/'seed.before.json').read_text());after=json.loads((B/'seed.after.json').read_text());obs={n:json.loads((B/(n+'.observation.json')).read_text()) for n in ['ovh','upcloud']};tests=json.loads((B/'tests.json').read_text());repair=json.loads((B/'failed_admission_repair.json').read_text())
def dt(t):return datetime.datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S CEST')
t=m['measurement_start'];table=[]
for n in ['ovh','upcloud']:
 b=m['metrics'][n]['before_10min'];a=m['metrics'][n]['after_so_far'];table.append(f"| {n} | {b['all_starts']} → {a['all_starts']} | {b['seed_starts']} → {a['seed_starts']} | {b['fetched_results']} → {a['fetched_results']} | {b['fetched_exit75']} → {a['fetched_exit75']} |")
btot=sum(m['metrics'][n]['before_10min']['all_starts'] for n in obs);atot=sum(m['metrics'][n]['after_so_far']['all_starts'] for n in obs)
report=f'''# LANE_RUNNER_SEED_UNBLOCK_RESULT

Completed 2026-09-30 with Anton’s full mandate. Both fixes are installed and verified. SEED jobs actually run via account D/swarm_worker; three such jobs have returned results and been collected. All 500 SEED IDs remain in the queue, including SEED-253. The previous 84 results remain; the latest result count is **{len(after['done_ids'])}**.

## Marker inventory and cleanup

| Host | Markers inventoried | Live | Dead | Reused/inconsistent | Removed during installation | Latest markers / actual model processes |
|---|---:|---:|---:|---:|---:|---:|
| OVH | 233 | 15 | 218 | 0 | 219 | {obs['ovh']['markers']} / {len(obs['ovh']['runtimes'])} |
| UpCloud | 1070 | 3 | 1067 | 0 | 1067 | {obs['upcloud']['markers']} / {len(obs['upcloud']['runtimes'])} |

The inventory was performed before installation. OVH gained one more dead marker between inventory and cleanup; UpCloud had two live markers remaining at cleanup after a normal completion. Cleanup held `model_slots.lock`; code/policy installation also held `rate_policy.lock`. All remote writes and cleanup ran as **ubuntu**, UID 1000, also on UpCloud via `sudo -u ubuntu`. Backups were created as ubuntu before replacement.

The original counter already ignored PIDs that were actually missing. Therefore the 1286 removed dead files alone cannot be credited with today’s capacity blockage. PID reuse was a real code risk but **was not found** in the inventory. The measured improvement also includes fallback and SEED priority, not an isolated effect of file cleanup.

Latest check: {dt(min(o['time'] for o in obs.values()))}. All remaining markers correspond to the launcher’s actual PID, correct job/cgroup and valid start identity. **No orphaned OpenCode processes and no gaps between marker and model process** in the latest snapshot. Agent units can also contain wrappers/admission waiters that do not yet have a model slot; they do not count as ongoing model runs.

## Exact code and field changes

On **both** hosts:

- `/opt/agents/slot_markers.py`: new shared marker handling. Checks `/proc/<pid>/stat`, zombie status, cmdline (`run_profile_ovh.py`, correct job), agent cgroup, boot ID and process start ticks. Legacy markers are rejected if the process start time is later than the marker’s reservation. Maximum age **5700 s = RuntimeMaxSec 5400 + 300 s margin**. A dead, reused, too old or inconsistent marker is removed under the slot lock. New markers are written via temporary file + atomic rename and include `start_ticks`/`boot_id`.
- `/opt/agents/route_runtime.py`: new `_active_counts`; `reserve_slot` uses the locked shared scan instead of only `os.kill(pid,0)`. Model selection considers egress availability before applying free preference/fill target. Reserve can be chosen when no free option can be admitted because of model cooldown, model ceiling or lack of healthy egress with available capacity. A single limited IP does not automatically send the job to swarm_worker if another healthy free route exists.
- `/opt/agents/egress_pool.py`: `_active_counts` uses the same locked cleanup. New `has_available` atomically reads egress state/health and existing cap values, without consuming a route allocation or changing egress entries. Egress → slot is the only lock order; admission takes no egress lock under slot lock.
- `/opt/agents/fallback_reserve.py`: `choose(..., free_unavailable=True)` replaces the requirement that both free providers’ global backoff must be active simultaneously. Reserve window, swarm_worker backoff and concurrency ceilings remain. The quota stop is also strictly limited to **at most 95 %**, even if a future policy field were to specify a higher value.
- `/opt/agents/run_profile_ovh.py`: own marker is removed under slot lock and own ephemeral runtime is cleaned up also on setup error/SIGTERM/deadline. `run_guarded` has `finally` terminating only its own `start_new_session=True` process group. SIGTERM/INT/HUP give controlled termination; SIGKILL/OOM are cleaned up by the next locked scan. An internal deadline of **1575 s** gives 25 s cleanup margin before the packet wrapper’s existing SIGKILL at 1600 s; deadline returns 124.

`FREE_MODEL_POLICY.json` was changed atomically under `rate_policy.lock`:

| Field | Before | After |
|---|---|---|
| `swarm_worker_value_only` | `true` | `false` |
| `allocation` | global both-cooldown | available free route first, D reserve when free admission is missing |
| `swarm_worker_reserve_authorization` | only high-value + both-cooldown | Anton’s SEED_UNBLOCK mandate, ordinary jobs and egress/cap admission included |
| `swarm_worker_reserve_account` | D | D |
| `swarm_worker_reserve_cap` OVH | 16 | **16** |
| `swarm_worker_reserve_cap` UpCloud | 3 | **3, existing stricter host ceiling preserved** |
| `swarm_worker_weekly_stop_pct` | 95 | **95** |
| `swarm_worker_weekly_used_pct` | 0.0, manual UNKNOWN | **0.0, manual UNKNOWN, untouched** |
| `swarm_worker_reserve_until` | existing 24h window | untouched |

The actual console value for D is unknown. No estimate or change of Anton’s manual quota value has been made. Fallback can now consume D quota faster; the 95% protection uses the manually entered value, not a new console measurement.

Locally `tasks/lanes/ovh_agents/bt_queue_ovh.sh` was changed (shared by both drivers):

- Queue presentation uses `seed_unblock_evidence/seed_queue_order.py`: stable SEED priority, the same rows exactly once. The canonical queue file was not edited; no jobs or rows were removed. This was needed because most SEED jobs lay behind AUTO/PLAN jobs and therefore did not get available host slots.
- The retry gate is checked before start but charged **after successful unit start**. On denied start, the own claim is released and the start pass breaks, so capacity and completed jobs are checked again. Previously a pass continued through the queue with an old slot snapshot and burned retries without actual starts.
- Under `body_rate_retry.lock`, **one provably unused retry each for {len(repair)} jobs** was restored in `BODY_RATE_RETRIES.json`. Only entries whose `last_claim` matched a denied start in time and lacked a subsequent successful start/active claim were changed. Executed attempt history and results were preserved; no general retry reset was done. Full before/after per job is in `failed_admission_repair.json`.

Both `bt-queue-ovh`/`bt-queue-upcloud` were stopped/started cleanly after dispatcher changes. Ongoing separate agent units were not stopped. Memory/OOM admission, host ceilings, slice limits and launch_reserved were not changed.

## Real fallback, actual provider run and results

For verification, free backoff was briefly set to at least `now + 120 s` under the policy lock on both hosts, twice while existing workers finished. Existing longer backoff was never shortened. The temporary additions have expired; current policy retains the ordinary backoff history.

Actual automated SEED runs on OVH:

| Job | Requested model | Account / selected model | Verification |
|---|---|---|---|
| BT-FW48-SEED-493 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + collected result |
| BT-FW48-SEED-494 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + collected result |
| BT-FW48-SEED-495 | swarm | **D / swarm_worker** | process + completed tools + exit 0 + collected result |

Example, actual `MODEL_ROUTE.json` for SEED-495:

```json
{{"account":"D","requested":"swarm","selected":"swarm_worker","model":"opencode-go/swarm_worker-v4.1-flash","started":1790759089.4781682}}
```

Process snapshot shows `opencode run --model opencode-go/swarm_worker-v4.1-flash` with these jobs and the marker’s launcher as parent. The snapshot for SEED-493/-494/-495 contains 7 completed tool calls each at measurement time. swarm_worker processes and tool calls were observed **also on UpCloud** in its ordinary AUTO jobs during the same fallback test; UpCloud also runs the existing SEED wave. Free routes resumed when verification backoff expired.

`CLOUD_COLLECTION.json` for completed jobs documents rsync checksum without differences, preserved source/inputs and matching result header. The results are produced reports, **not scientifically admitted or independently reviewed evidence**.

## Start rate and SEED integrity

Comparable 10-minute windows:

- Before: **{dt(t-600)} → {dt(t)}**.
- After: **{dt(t)} → {dt(t+600)}**. Snapshot covers {min(m['elapsed_s'],600):.0f} of 600 s.

| Host | All `ovhstart` before → after | SEED starts before → after | Collected results before → after | Collected exit=75 before → after |
|---|---:|---:|---:|---:|
{chr(10).join(table)}

Total **{btot} → {atot} starts/10 min**, {((atot/btot)-1)*100:.1f}% increase; SEED starts **0 → {sum(m['metrics'][n]['after_so_far']['seed_starts'] for n in obs)}**. This is a short operating window, not a long-term estimate of provider capacity. Earlier historical ~14 starts/result are not compared against an invented longer after-period. The logs only have HH:MM:SS; the measurement reads the latest day’s tail backwards to avoid double-counting the same clock time from earlier days.

**500/500** SEED IDs in the queue before and after, **500/500** task briefs remain, **SEED-253 remains**. All previous **84** result IDs remain. Newly collected SEED results: {', '.join(after['new_done']) or 'none yet'}. No SEED results were deleted. Existing limits on real retries are preserved; this work does not promise that all 500 scientific assignments are already complete.

## Protection, command evidence and artefacts

**{tests['passed']} tests pass**: live/reused/boot/legacy/age/dead/inconsistent PID; permanent cleanup in both counters; egress backoff and healthy free alternative; non-high-value SEED → D/swarm_worker; 95% stop; ceiling 16; reserve window end; swarm_worker backoff; own marker cleanup; actual subprocess tests for SIGTERM/deadline, correct exit status and stopping only the own child process group while an unrelated process survives. Dispatcher has also passed `bash -n`. The installed files’ SHA256 match the tested copies.

The following commands were used, with remote Python via stdin where the scripts contain the inventory/patch itself:

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

Both dispatcher services are `active`. `/opt/agents/EGRESS_POOL.json` is byte-identical to before on both hosts; direct/home_via_ssh/warp entries have not been changed. Profile keys/config have not been changed. No credentials, mailings or external publishing are included.

All evidence is in [seed_unblock_evidence](seed_unblock_evidence/):

- `ovh.inventory.before.json`, `upcloud.inventory.before.json`: marker inventory with PID/cmdline/cgroup.
- `ovh.deploy.json`, `upcloud.deploy.json`: UID, actual cleanup, policy before→after and unchanged egress hash.
- `ovh.patch.diff`, `upcloud.patch.diff`, `dispatcher.patch.diff`: exact diffs; `slot_markers.py` is the new file’s full source.
- `tests.json`, `verify_patch.py`: test outcomes and executable verification.
- `ovh.seed_fallback_proof.json`, `ovh.swarm_worker_processes.txt`, `upcloud.fallback_snapshot.json`, `upcloud.swarm_worker_processes.txt`: actual routes, markers, processes and tool calls.
- `BT-FW48-SEED-493.*`, `BT-FW48-SEED-494.*`, `BT-FW48-SEED-495.*`: collected MODEL_ROUTE, exit/attempts, results and collection receipts.
- `ovh.observation.json`, `upcloud.observation.json`, `metrics.json`, `seed.before.json`, `seed.after.json`: final check and measured data.
- `failed_admission_repair.json`, `BODY_RATE_RETRIES.before_failed_admission_repair.json`: scoped restoration of unused retries.
'''
(B.parent/'LANE_RUNNER_SEED_UNBLOCK_RESULT.md').write_text(report)
print('written',B.parent/'LANE_RUNNER_SEED_UNBLOCK_RESULT.md')
