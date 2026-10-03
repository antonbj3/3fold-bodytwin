# Mission: stabilize egress/WARP for the free models' capacity

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton). Execute without asking; report default selection. Priority: **progression, not perfection**Never leave the system worse than you found it.

## Bakgrund (verifierat)
- The provider rate limit for the free models (swarm/swarm_worker) is **per (model, egress-IP)**. More IP ⇒ more capacity. This is the only remaining bottleneck; OOM, accounts and load are already resolved.
- Egress-poolen: `~/research/FREE_AUTONOMY_20260926/OVH_DIRECT_20260927/egress_pool.py`, config `/opt/agents/EGRESS_POOL.json` on both hosts. Current entries: `direct` (values' own IP) and `home_via_ssh` (proxy `http://127.0.0.1:18187`)Each entry has: `proxy`, `weight`, `cap`; per-post backoff hanteras automatiskt.
- `apply_route` uses `egress_pool.pick(model,directory)` and sets `HTTPS_PROXY` for the job. `record()` puts an entry in backoff at rate-limit so that a dead egress never blocks the model.

## What has already been done locally (build on, don't start over)
- Cloudflare WARP WireGuard-profil genererad med `wgcf` (`/tmp/opencode/wgcf/`), konto skapat, tunnel-IP `172.16.0.2`.
- `warp_proxy.py` (minimal HTTP CONNECT-proxy) i `~/research/FREE_AUTONOMY_20260926/OVH_DIRECT_20260927/warp_proxy.py`, systemd-user service `warp-proxy.service` (auto-restart), lyssnar `127.0.0.1:18188`.
- Customized tunnel-config `/etc/wireguard/warp-min.conf` (`Table = off`, `ip rule` fwmark 0xca6c → table 51820, plus `not fwmark 0xca6c → main`). `SO_MARK=40` in the proxy directs outgoing traffic into WARP-tabellen without touching the host's default-route.
- **Remaining bugs:** the WARP tunnel is not currently shaking hands (`latest-handshakes = 0`), so the proxy falls back to normal IP. the WARP endpoint `162.159.192.1:2408` responds to ping/UDP.

## Tasks (in order of priority)
1. **Get WARP to shake hands** (or change endpoint): try `162.159.192.1..8`, `engage.cloudflareclient.com`, port 2408/500/1701/4500; verify `wg show` handshake and that `curl -x http://127.0.0.1:18188 https://api.ipify.org` gives a WARP-IP (`104.28.x`/`2a09:`), while direct traffic keeps `31.208.x`.
2. **If WARP does not become stable within a reasonable time:** still put in the record with **very low weight** and verify that `record()`/`pick()` correctly backoff it so that `direct`/`home_via_ssh` takes over. Instead, focus on **maximizing `home_via_ssh`** (more simultaneous tunnels/proxy ports towards the home IP) or other existing egress.
3. **Always preserve fallback:** no single egress must be critical. Verify with paired probes (same account/model/task via two egress) that an egress-IP-rate-limit does not pause the entire model.
4. **Operating mode:** The proxy should run as a service (auto-restart). If the WARP tunnel is needed on the hosts, do the corresponding minimal setup there (OVH `ubuntu@51.77.110.4`, UpCloud `root@212.147.226.179`) and write `/opt/agents` files as **ubuntu** (never root). Verify that the SEED run is not interrupted.

## Constraints
- Do not touch the swarm_worker ratio (fallback, paused). Don't touch the OOM/slot logic (just fixed and verified).
- Do not change `EGRESS_POOL.json` to remove `direct`/`home_via_ssh`. Just add.
- No emails/pushers/credentials. Do not touch SEED jobs or results.
- After change: verify with a **real SEED run** (or paired probes) that new egress is indeed used and that the IP differs.

## Definition of done
1. Either a working, backoff protected WARP-egress **or** documented free-WARP not stable + a maxed `home_via_ssh`/`direct` path.
2. Parade probar som visar per-IP-oberoende.
3. SEED run unaffected; no new failed units; no critical individual egress.
4. Kort sammanfattning i `~/research/FREE_AUTONOMY_20260926/OVH_DIRECT_20260927/LANE_RUNNER_EGRESS_WARP_RESULT.md` (what was changed, IP:n, measured power, what is NOT gjordes).
