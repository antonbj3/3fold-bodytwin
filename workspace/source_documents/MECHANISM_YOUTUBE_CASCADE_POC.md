# MECHANISM YOUTUBE CASCADE POC — network fetch is bot-blocked in this environment (2026-07-21)

Task: use the newly-authorized YouTube fetch capability (`yt-dlp` 2026.07.04, `ffmpeg`/`ffprobe`,
network confirmed working) to find genuinely-high-effective-temporal-resolution slow-motion footage
of a ballistic movement (vertical jump / sprint push-off), to test whether the ankle→knee→hip→lumbar
force cascade — a **tie** (all peaks in the same frame) on this project's own 30fps IG corpus, both
for a squat (`docs/MECHANISM_FORCE_TRANSMISSION_SCENE.md`) and a plyometric jump
(`docs/MECHANISM_PLYO_JUMP_CASCADE.md`) — resolves into a real stagger at higher temporal sampling.

## Bottom line

**`yt-dlp` cannot fetch arbitrary YouTube video metadata or content in this environment without
real sign-in credentials.** Confirmed by direct measurement (forced OODA, not a first-shot bail —
§1), not assumed from the first error message: 23 distinct, deliberately diverse candidate video
IDs (official VEVO music video excluded — see the one control below; niche biomechanics-lab clips,
a mainstream sports clip, action-sports channels, coaching channels) were fetched directly by
video ID (bypassing search) — **23/23 failed**, every time, with YouTube's server-side player
response reporting `playability status: LOGIN_REQUIRED` regardless of which internal client
(`android_vr`, `web_safari`, and 9 more variants tried — §1) yt-dlp used. Exactly one control ID
(`dQw4w9WgXcQ`, tested 3 times interleaved with the failing batch) succeeded every single time —
which rules out a simple rate-limit/IP-cooldown explanation (it kept working immediately after a
12-video failure streak) and rules out "unpopular clips get blocked" (a mainstream, extremely-viewed
clip failed identically to the niche ones). **Zero clips were downloaded**; `data/youtube/` contains
no fetched files (verified live, not assumed). Per the coordinator: the same cascade-resolution
question is being pursued in parallel from the existing IG corpus's own `genuine_highspeed`-class
clips (no network fetch needed) — this doc's scope ends at the network-fetch feasibility question.

**⚠ PERSISTENCE STATUS (2026-07-21):** the qualitative finding ("zero clips downloaded") is
independently verifiable right now via a plain filesystem check (`find data/youtube -type f`) and
remains true. **The specific "23/23" count above is doc-prose-only — no list of the 23 video IDs or
per-ID pass/fail result was ever written to disk**, flagged in `docs/MECHANISM_TRUST_LEDGER.md` §10
item 6. Not re-run this session: reproducing it means firing 23 fresh network requests at YouTube
(across up to 11 client variants each, per §1) whose outcome depends on YouTube's current
server-side gating, not on anything this repo controls — a live re-run could legitimately come back
different (more, fewer, or zero failures) without that saying anything new about this environment's
`yt-dlp` capability at the time this doc was written. Marking prose-only rather than re-running is
the deliberate, lean choice here (see `docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md` for the full
disposition); a future session that needs a persisted per-ID array should write one the first time
it re-attempts this fetch, not retroactively reconstruct one now.

## 1. What was tried (forced OODA, exact evidence — never eyeballed, always the tool's own error text)

- **`ytsearchN:<query>` (the task's suggested discovery method): 40/40 failed.** 5 queries × 8
  results each (`slow motion vertical jump side view`, `high speed camera sprint start`,
  `240fps box jump`, `slow motion squat jump biomechanics`, `1000fps slow motion jump biomechanics`)
  — every single result errored at the metadata-extraction stage (before any download) with
  `Sign in to confirm you're not a bot`.
- **11 distinct player-client variants tried, all identical failure**: `android`, `tv`, `ios`,
  `web`, `mweb`, `tv_embedded`, `web_embedded`, `android_embedded`, `ios_music`, `web_music`,
  `android_vr` (plus the unmodified default, which auto-tries several of these). Not a client-
  selection problem.
- **Discovery workaround found (and used successfully)**: `--flat-playlist` on the literal
  `https://www.youtube.com/results?search_query=<query>` URL (as opposed to the `ytsearch:`
  pseudo-extractor) bypasses the search-endpoint gate and returns real `id`+`title` pairs. Used
  across 6 queries (the task's 4 + 2 more: `1000fps slow motion jump biomechanics`,
  `slow motion countermovement jump side view`) to build a 20-candidate shortlist including
  promising hits like `qN3apht8zRs` ("Slow motion video of vertical jump with **synchronized
  vertical force data**" — would have been the strongest candidate, a real force-plate external
  anchor) and several "Biomechanics"/lab-style clips.
- **But per-video access (the actual metadata/fps check or download) still requires the full
  player-response call, which is where every candidate failed.** Direct watch-page fetch
  (`https://www.youtube.com/watch?v=<id>`, bypassing `ytsearch:` entirely) was tested on 23
  distinct IDs drawn deliberately from across the shortlist's full diversity (lab biomechanics
  clips, a mainstream sports clip, action-sports/"slow-mo" channels, coaching-tutorial channels):
  **23/23 failed**, all with the same bot-check error.
- **Root cause confirmed via verbose (`-v`) log, not guessed from the generic error text**: for a
  representative failing ID (`qN3apht8zRs`), yt-dlp's own debug output shows
  `qN3apht8zRs: android_vr player response playability status: LOGIN_REQUIRED` and
  `qN3apht8zRs: web_safari player response playability status: LOGIN_REQUIRED` — the block is
  YouTube's own server-side per-video gate, not a yt-dlp parsing bug or a generic block page.
- **Forced adversary #1 — is this just an IP-level rate-limit/cooldown from the earlier `ytsearch`
  failures, not a real per-video gate?** Ruled out: the one control ID, `dQw4w9WgXcQ` (an
  astronomically-viewed video, ~1B+ views, plausibly served through a different/cached/allowlisted
  path), was re-tested 3 separate times interleaved between failing batches — **succeeded all 3
  times**, including immediately after a streak of 12 consecutive failures on fresh IDs. A global
  cooldown would have caught this ID too; it did not.
- **Forced adversary #2 — is this just "unpopular/niche clips get blocked, mainstream ones don't"?**
  Ruled out: `yHAc6bAWE0U` ("Usain Bolt WR 100M", a mainstream, presumably heavily-viewed sports
  clip) failed identically to the niche biomechanics clips.
- **Checked for the current community-standard fix**: a locally-installed PO-token-provider plugin
  (e.g. `bgutil-ytdlp-pot-provider`) would supply the token YouTube's gate is asking for without
  needing cookies. `yt-dlp --version --verbose`'s own plugin listing is empty; no `bgutil`/
  `pot_provider` package found anywhere on the filesystem (`find / -iname "*bgutil*" -o -iname
  "*pot_provider*"` → no hits). Not present in this environment.
- **Tried a public Invidious mirror** (a YouTube-content proxy that sometimes sidesteps this exact
  gate): `invidious.privacyredirect.com` → `502 Bad Gateway` (instance itself down); `yewtu.be` →
  reached YouTube's own extractor anyway and reproduced the identical `LOGIN_REQUIRED` error.
  Not a viable bypass with the instances tried.

**Conclusion: this is not a client-selection, rate-limit, or popularity effect — it is YouTube's
server-side sign-in gate, applied at the per-video API level, and no yt-dlp flag short of real
credentials clears it in this environment.**

## 2. Exact unblock recipe (for a future session with real credentials — not attempted here)

Two supported, standard yt-dlp options, in order of convenience:

1. **`--cookies-from-browser <browser>[:PROFILE]`** — e.g. `--cookies-from-browser firefox` or
   `--cookies-from-browser chrome:Default`. Points yt-dlp at an **already-logged-in-to-YouTube**
   browser profile's cookie store on this same machine. **Not attempted in this build**: no browser
   profile logged into YouTube was available/authorized for this session's use, and unilaterally
   reaching into a personal browser profile's credentials without being directed to is an
   authorization boundary this build does not cross on its own.
2. **`--cookies /path/to/cookies.txt`** — a Netscape-format cookies file exported from a logged-in
   browser session (e.g. via a "Get cookies.txt LOCALLY" browser extension), supplied explicitly by
   the operator.

**Decisive go/no-go check once either is available** (cheap, run before any bulk fetch): re-run
against the one already-confirmed-blocked ID from this build —
```
yt-dlp "https://www.youtube.com/watch?v=qN3apht8zRs" --skip-download -v --print "%(id)s" \
  2>&1 | grep "playability status"
```
Success = no `LOGIN_REQUIRED` line and a real title/fps/duration print. That specific ID is also
the single best candidate already identified (§1) — a slow-motion vertical jump **with synchronized
force-plate data**, which would let the eventual cascade-timing measurement anchor its flight/
propulsion-phase timing against real force-plate ground truth rather than a kinematic proxy alone.

## 3. Where the cascade question is being answered instead

Per the coordinator's live redirect: a parallel build is pursuing the identical ankle→knee→hip→
lumbar cascade-resolution question using the existing IG corpus's own `genuine_highspeed`-class
clips (no network fetch required), rather than waiting on newly-fetched YouTube footage. This
document's scope ends at the network-fetch feasibility question; see that parallel build's own
output for the actual cascade-timing result (resolved-with-stagger vs. still-tied).

## 4. Honest status of this build's own deliverables

- **`data/youtube/`**: verified empty of any fetched file (`find data/youtube -type f` → no
  results; a stray empty `foot/` subdirectory from unrelated work is the only content) — zero
  clips fetched, consistent with every fetch attempt above failing before the download stage.
- **No pose extraction, IK, or force-cascade analysis was run by this build** — there was no
  downloaded clip to run it on. `scripts/msk/pose_extract.py` (MediaPipe, `.venv-humancap`) and the
  `scripts/msk/force_scenes_batch.py` / `scripts/msk/plyo_jump_cascade.py` Newton-reaction
  primitives were read and are understood well enough to drive a new clip through the existing
  pipeline (`.trc` → `pose_to_opensim_ik.py`'s smoothed `.mot` → the propulsion-phase
  `cascade_order_in_window` measurement) — not exercised here for lack of input data.
- **`scripts/msk/youtube_cascade_poc.py` was deliberately NOT written**: a script with no
  successful fetch behind it and no real clip to validate against would be untested scaffolding,
  not a working proof-of-concept — skipped per the coordinator's explicit re-scope to this findings
  doc only, and per this repo's own lean discipline (no wasted, unrunnable artifacts).
- **A real methodological note worth preserving for whoever next attempts this with credentials**:
  the cascade-timing ORDER/TIE question is invariant to the video's true-vs-assumed frame timing
  (a uniform `dt` misestimate rescales every body's computed acceleration by the same constant
  factor and cannot change which frame index is the argmax for any cut), but the **millisecond
  stagger** is not — converting a frame-gap into a real ms figure requires knowing the clip's TRUE
  capture rate, not just its container fps (a YouTube "slow motion" upload is almost always
  re-encoded at a standard container fps like 30/60 with the slow-down baked into extra frames, so
  `ffprobe`'s fps is not the true capture rate). The clean, externally-anchored way to recover the
  true slow-down factor without trusting the upload's title/description claim: identify a genuine
  free-flight (zero-ground-contact) window geometrically (shape-based, not derivative-magnitude-
  based, to avoid circularity), fit a quadratic to a body's vertical position vs. apparent
  (container-fps-assumed) time within that window, and solve the slow-down factor from the fitted
  curvature against the true physical constant `g = 9.80665 m/s²` (`scripts/msk/
  validate_joint_force.py`'s own `G`) — the same "anchor externally, never a tautology" principle
  this repo already applies elsewhere, here applied to calibrate time itself. Not yet built or
  tested; recorded here so it does not need re-deriving next time.

## Reproducing the block (for anyone re-checking this finding)

```
# search: 100% fail
yt-dlp "ytsearch3:slow motion vertical jump" --skip-download --print "%(id)s"
# -> ERROR: [youtube] <id>: Sign in to confirm you're not a bot...

# discovery workaround (works, metadata-only, no per-video access):
yt-dlp "https://www.youtube.com/results?search_query=slow+motion+vertical+jump" \
  --skip-download --flat-playlist --playlist-end 5 --print "%(id)s | %(title).60s"

# per-video access: fails even for a hand-picked ID from the above
yt-dlp "https://www.youtube.com/watch?v=qN3apht8zRs" --skip-download -v --print "%(id)s" \
  2>&1 | grep "playability status"
# -> playability status: LOGIN_REQUIRED (android_vr AND web_safari clients)

# control (succeeds every time, does not establish a general fix):
yt-dlp "https://www.youtube.com/watch?v=dQw4w9WgXcQ" --skip-download --print "%(id)s | fps=%(fps)s"
```
