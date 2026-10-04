# Egress/WARP — results 2026-09-30

**Working WARP is deployed in the egress pool on both OVH and UpCloud.** Endpoint `162.159.192.1:2408` shakes hands. Local direct traffic still gives `31.208.14.101`; WARP gives `104.28.234.225` and Cloudflare's `warp=on`. Over 24 minutes of observation after the first handshake; the scheduled check has **82/82 healthy cough probes**, with the same WARP-IP.

## Changes and operation

- `warp-min.conf` missing private key. The already generated wgcf key was restored, without re-registration or printing of the key. `Table=off`, table 51820 and mark `0xca6c` were retained; the rules got explicit priorities 100/101. `wg-quick@warp-min` is enabled at boot, with keepalive 25 s. The host's default route is unchanged.
- The proxy used the wrong socket option: **40 was replaced with `socket.SO_MARK` (36 on this Linux)**. It also binds the source address of the tunnel and terminates at mark-/bindfel. The service runs in system-scope as **anton**, with only `CAP_NET_ADMIN` and `Restart=always`; the old user service was archived. The port is still `127.0.0.1:18188`. [Linux socket documentation](https://www.man7.org/linux/man-pages/man7/socket.7.html) describes the right for SO_MARK.
- Two separate, auto-restarting SSH services expose this proxy on the same loopback port on the respective host. No additional WireGuard accounts or tunnels are needed on the cloud hosts.
- Both pools got `warp`, **weight=2, cap=2**. `direct` and `home_via_ssh` are exactly unchanged: weight 9/1, direct cap=0, home cap=8 on OVH and 6 on UpCloud. The Home tunnels remain. The low weight option/maxing home was not needed.
- `research-warp-health.timer` checks the entire respective host→SSH→WARP path every 30 second. A health record applies to 75 s. Missing, expired or negative health deselects WARP and sets its model backoff; other paths continue.
- `pick()` and `record()` now use the same model key (`swarm`/`swarm_worker`); old aliases are migrated without losing counters/backoff. Fractional weights work. A start failure backs off the affected path; the model's global pause only changes when all paths are unavailable.
- EThe cold OpenCode runtime could start with a native directory before the public directory fetch was complete. For WARP, **enable public model metadata** is therefore copied to the isolated runtime cache. It also works for admission waitlists who already had legacy routing code loaded. The catalog is updated hourly; the health check requires both free models to be present. Public directory/npm-bootstrap goes via normal route; the model request and IP-kontroller go via WARP.

## Parade riktiga probar

Same account **A**, model **swarm_worker**, task and input: read the factors 17/19, use Python and write `{"product":323}`. The final probes used the installed routing-/cachekoden from empty runtimes, without changing the production policy or slot markers.

| Host | Egress/IP | Outcome | Tid |
|---|---|---|---:|
| OVH | direct / 51.77.110.4 | Explicit provider rate-limit, no answer | 9,12 s |
| OVH | WARP / 104.28.234.225 | Correct answer, 3 completed utility calls | 21,01 s |
| UpCloud | WARP / 104.28.234.225 | Correct answer, 4 completed tool calls | 22,01 s |

The paired OVH outcome shows per-IP independence for this account/model/task. A previous check with preloaded directory gave the same rate-limit/pass difference. This is **a new shared IP bucket**, not two independent WARP buckets for the two hosts. Exact ratio and increase in jobs/minute have not been measured. Backoff state is per host.

**12 tests per host passed:** model alias, migration, fractional weight, health loss/expiry, healthy WARP when others IP:n backoff, live cap, public cache content, old routing calls, directory error→fallback, startup error and explicit rate-limit without model break. An unauthorized proxy process was verified aborting on SO_MARK error. Both hosts' added/changed egress files, directory and health records are ubuntu owned and written through ubuntu.

## SEED and delimitations

SEED jobs/results have not been manually modified, stopped, restarted or restarted by this lane. Before/after checks show continued SEED dispatch and waiting/running processes; the queues' regular collection and provider interruptions continued during the work. `run_profile_ovh.py` , `fallback_reserve.py` and `launch_reserved.py` have the same content hash as at launch. The admission/slot features' ASTs are unchanged. The routing code was changed only for egress cache and boot failure backoff. **No new SEED run was forced; true separate probes were used.** The baseline included jobs in progress, not a hash inventory of all already collected SEED results.

The probe uncovered a race in existing `research-runtime-cache.service`: a completed diagnostic runtime was lost between directory listing and `stat()`. After completing the tests, the existing service was successfully executed: `Result=success`. Its code and OOM/slot-logik were not changed. **Both cloud hosts have zero failed system units at final check; all egress services are healthy.**

Exception from a completely empty local failed list: `swarm_worker25f-20260928.service` was failed in the parallel swarm_worker task (00:34:02 CEST, exit 1) and was left untouched. The two previous local system errors remain. swarm_worker policy was simultaneously changed by the separate `LANE_RUNNER_reserve_worker_D` lane according to its policy note; this lane has neither written nor reset its policy/kvot, credentials or settings. No OOM-/slot-parametrar were changed. No emails, pushes, new accounts or credentials were created.

Source material: [stabilitet](WARP_EVIDENCE/stability.json), [parat OVH-results](WARP_EVIDENCE/51.77.110.4_paired_result.json), [UpCloud-results](WARP_EVIDENCE/212.147.226.179_paired_result.json), [OVH slutstatus](WARP_EVIDENCE/51.77.110.4_final.json), [UpCloud slutstatus](WARP_EVIDENCE/212.147.226.179_final.json), [OVH tester](WARP_EVIDENCE/51.77.110.4_pool_tests.txt), [UpCloud tester](WARP_EVIDENCE/212.147.226.179_pool_tests.txt). Backups are in `WARP_EVIDENCE`, next to the modified cloud files and i `/etc/wireguard/warp-min.conf.before-egress-20260930`.
