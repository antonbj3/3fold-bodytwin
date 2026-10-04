# Task: connect swarm_worker fallback on account D (+ clean up quota logic)

You are a lane_runner `build-lane model` agent with full authority from the coordinator (Anton). Execute without asking; report default choices. Priority: **correct, safe fallback** — not burning D's quota.

## The situation (verified by the coordinator)
- The free models (swarm/swarm_worker) are rate-limited per (model, egress IP) and are the bottleneck. **swarm_worker is fallback** and must be used when both free models are in cooldown.
- **Anton says: only account D has swarm_worker quota left right now. Accounts A and C are effectively exhausted.**
- But the guard logic points incorrectly:
  - OVH `ubuntu@51.77.110.4`: `swarm_worker_reserve_account = "A"`, `swarm_worker_weekly_used_pct = 65.87`
  - UpCloud `root@212.147.226.179`: `swarm_worker_reserve_account = "C"`, `swarm_worker_weekly_used_pct = 86.78`
  - Both have `swarm_worker_backoff_until ≈ 2790683370` (**year 2058** — an old hard-coded pause that blocks all swarm_worker), and `swarm_worker_reserve_until = 0`.
  - `swarm_worker_weekly_used_pct` is a **manual** field (the key: "manual: run swarm_worker-usage <pct> after checking the console") and measures nothing itself.

## Key facts
- **Profile D already exists** on both hosts in `/opt/agents/profiles/D/data/opencode/auth.json` (ubuntu-owned, mode 600), with providers `opencode` + `opencode-go`. Verified to work locally and on both clouds (`MODEL_ROUTE.json` → `account: D`). Runner `run_profile_ovh.py` has `choices=('A','B','C','D')`.
- **Fallback logic:** `/opt/agents/fallback_reserve.py` (`choose`), called from `route_runtime.reserve_slot` when no free models are available. swarm_worker allocation is governed by `FREE_MODEL_POLICY.json`: `swarm_worker_reserve_account`, `swarm_worker_reserve_cap`, `swarm_worker_reserve_until`, `swarm_worker_weekly_used_pct`, `swarm_worker_weekly_stop_pct`, `swarm_worker_backoff_until`, `swarm_worker_value_only`.
- **Helper:** `local_config_path/bin/swarm_worker-usage` (status/set/stop/resume/limit) — but it only sets the **same** pct on both hosts and does not redirect the account.
- **IMPORTANT about ownership:** write `/opt/agents` files as **ubuntu**, never root (OVH: run as ubuntu; UpCloud: `sudo -u ubuntu ...`). Protected by `flock` on `/opt/agents/rate_policy.lock`.
- **Fallback account D's key** exists in the profile; swarm_worker model `opencode-go/swarm_worker-v4.1-flash`.

## Task
1. **Verify first** that D's swarm_worker actually works from the hosts (run a minimal prompt through `run_profile_ovh.py D swarm_worker` in a temporary job directory). Record `MODEL_ROUTE.json`.
2. **Redirect fallback to D** on both hosts: set `swarm_worker_reserve_account = "D"` (and the corresponding go account if `fallback_reserve.choose` uses such a field — read the code first).
3. **Remove the old 2058 pause:** set `swarm_worker_backoff_until` to `0` (or an expired value) so swarm_worker is not globally blocked. **Justify in the results file.**
4. **Quota counter:** set a reasonable starting value for D. Anton's actual D quota is unknown — set **0 %** as a conservative starting value if an exact value cannot be read, and note that Anton updates the exact value tomorrow. Keep `swarm_worker_weekly_stop_pct = 95`. **Do not touch A/C's historical numbers unnecessarily** — they have their own accounts; update only what is needed for fallback to work.
5. **Keep D's quota protected:** swarm_worker must only be used as **fallback** (when the free models are in cooldown), with `swarm_worker_reserve_cap` reasonable (e.g. unchanged 16). Do not change to "always on".
6. **Verify:** a real run that forces fallback (e.g. temporarily putting the free models into backoff in a test directory, or running a HIGH_VALUE job) must land on `account: D`, `selected: swarm_worker`, and leave `MODEL_ROUTE.json`/`ADMISSION_WAIT` as evidence. Restore any test manipulation afterwards.
7. **Preserve the run:** the SEED jobs must not be disturbed; no new failed units; no critical single egress/account.

## Constraints
- Do not touch the OOM/slot logic (just fixed) or the egress pool (a separate Sol agent is running it).
- No emails/pushes/credentials. Do not change profile keys. No SEED results deleted.
- If D's swarm_worker does not work: leave the system unchanged (except for documentation) rather than pointing fallback at a dead account.

## Definition of done
1. Fallback points to **D** on both hosts; the 2058 pause gone; quota value set (0 % or measured).
2. Proven fallback run on `account: D` / `selected: swarm_worker`.
3. A/C otherwise untouched; free models still primary; swarm_worker only during cooldown.
4. Short summary in `external_research_path`: exact field changes before/after, command evidence, what was NOT changed.
