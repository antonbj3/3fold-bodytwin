# Task: two urgent fixes for the SEED-500 flow

You are a lane_runner `build-lane model` agent with full authority from the coordinator (Anton). Execute without asking; report default choices. Priority: **get more jobs completed**. Do not touch the OOM/slot-memory logic beyond what is stated here.

## Background (verified by the coordinator 2026-09-30)
- **498 of 500 SEED jobs have been started** (1150 start events, 462 jobs started >1 time, up to 6 times) but only **84 completed**. The churn is enormous: ~14 starts per completed job. Main cause: `exit=75` (provider rate limit / admission wait expired).
- **Slot accounting leaks:** OVH has **233** active slot markers in `/opt/agents/active_models/`, UpCloud **1070** — but **0 swarm_worker processes** and almost no real agents. Old markers for dead processes remain and block new starts (the queue thinks the slots are occupied). `route_runtime.reserve_slot`/`egress_pool._active_counts` only clear markers when `os.kill(pid,0)` fails — but markers from reused PIDs or unclean exits remain.
- **The swarm_worker fallback is not used at all:** `swarm_worker_value_only=true` AND the fallback requires **both** swarm AND swarm_worker to be in cooldown simultaneously → 0 swarm_worker runs. Account D has quota left (`swarm_worker_reserve_account=D`, `pct=0.0`, 24h window open).
- Anton wants: **swarm_worker may also be used for non-high-value SEED jobs** (it is the only capacity left). Continued stop at 95 %.

## Task 1: clear slot markers and stop the churn
1. **Inventory** `/opt/agents/active_models/*.json` on both hosts. Classify: live process (pip exists and corresponds to a real agent), dead (pid does not exist), or reused/inconsistent (the marker's `job`/`started` matches no running unit). Report the distribution.
2. **Clear dead markers** safely and atomically under `model_slots.lock` (and `rate_policy.lock` when writing policy). Feel free to add an **age limit** (e.g. marker older than `RuntimeMaxSec` 5400 s + margin ⇒ definitely dead) to make cleanup robust, not just pid-based.
3. **Make cleanup permanent:** patch `route_runtime._active_counts`/`reserve_slot` (and `egress_pool._active_counts`) so that markers whose process is dead OR older than a reasonable limit are ignored/cleaned, and completed jobs always remove their own marker (even on SIGTERM/timeout). Preserve `run_guarded` semantics and stopping only the own process group.
4. **Verify:** after the fix, the number of active markers on both hosts must correspond to actually running agents; the queue must fill free slots again. Measure start rate before/after (new `ovhstart` per 10 min).

## Task 2: make the swarm_worker fallback actually usable
1. **Lower the fallback threshold:** `swarm_worker_value_only` must **not** require high-value for SEED jobs. Set `swarm_worker_value_only=false` (thus allow the Go reserve for ordinary jobs too) or change the `_high_value` condition in `route_runtime.reserve_slot` so that SEED jobs (`BT-FW48-*`/`BT-DW48-*`) qualify. Keep the account stop at 95 % and `swarm_worker_reserve_cap` (16) as protection.
2. **Ensure the fallback triggers when needed** — not only when both free models are in *global* cooldown. If an egress IP is in backoff but the model is otherwise free, consider whether swarm_worker should still be able to take the job; document the choice. The goal: when the queue stalls due to rate limits, swarm_worker jobs must actually start.
3. **Verify** with a real run: force fallback (free models in backoff) and prove in `MODEL_ROUTE.json` that `account=D, selected=swarm_worker`, and that swarm_worker processes are actually visible on the hosts.
4. **Continued protection:** monitor D's quota (95 % stop). If D's real console value is unknown, leave `swarm_worker_weekly_used_pct` as Anton sets it, but note that fallback now consumes quota faster.

## Constraints
- **Write `/opt/agents` files as ubuntu** (OVH runs as ubuntu; UpCloud: `sudo -u ubuntu`). Never root. `rate_policy.lock`/`model_slots.lock` are held open during writing.
- Do not touch profile keys or the egress pool's `direct`/`home_via_ssh`/`warp` records (just verified).
- No emails/pushes/credentials. No SEED results deleted.
- Restart `bt-queue-ovh`/`bt-queue-upcloud` cleanly after a change if code has been reloaded.
- Verify that the SEED run does not lose jobs (499/500 must be in the queue; SEED-253 was once restored).

## Definition of done
1. Dead slot markers gone; active markers match reality on both hosts; permanent protection built in; measurable increase in start rate.
2. swarm_worker fallback tested and proved to be used for SEED jobs; `account=D/selected=swarm_worker` in at least one real run; quota protection (95 %) retained.
3. Exact file/field changes before→after, command evidence and measured effect in `tasks/lanes/LANE_RUNNER_SEED_UNBLOCK_RESULT.md`.
