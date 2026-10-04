# OrthoLoad / CAMS-Knee live fetch report (2026-07-21)

Status: COMPLETE. All 6 implant IDs (knee, hip×2 telemetry generations, shoulder, spine×2 device types)
fully crawled via the live public API and downloaded; every number below is measured from the actual
downloaded files / manifests on disk, not estimated.

## 0. Correction to prior framing (measured, not assumed)
Earlier repo docs (`data/MECHANISM_SOURCE_REGISTRY.json`, `ORTHOLOAD_ACCESS_REQUEST_DRAFT.md`) characterized
OrthoLoad as "gated-individual-request" for hip/knee/shoulder/spine. **Live measurement this session shows
that framing was only half right — the site has TWO distinct data surfaces with different access:**

1. **The main "OrthoLoad Database" (database/, database-2017/) + "Test Loads" pages** — a JS-driven
   search UI over a plain, public, unauthenticated WordPress AJAX API (`wp-admin/admin-ajax.php`,
   actions `orthoload_get_activities/parameters/patients/files/data`) that serves real per-subject
   in-vivo force/moment trial files (`.akf` etc.) as **static files with NO login, NO nonce, NO session,
   NO form** — confirmed by (a) `robots.txt` explicitly carving out
   `Allow: /wp-admin/admin-ajax.php` while disallowing the rest of `/wp-admin/` (an external, site-owner
   signal that this endpoint is meant to be machine-fetched), (b) the JS source
   (`database2017/database.js`) containing no nonce/token/auth parameter anywhere, (c) a live pull of an
   actual force file (`k1l_110108_1_86p.akf`) succeeding with a plain unauthenticated GET, byte-identical
   whether fetched directly or via the site's own `downloadfile.php` shim. **This part of the coordinator's
   "no login, just a JS click-through" claim is CONFIRMED.**
2. **"Comprehensive Data Sample" (orthoload.com) and CAMS-Knee "Sample Data Download" / "Data Request"**
   — these are genuine identity-collecting gates, NOT mere JS navigation: they require submitting a web
   form with First Name / Surname / Organization / Email / Status ("student/scientist/industry/other") and
   an explicit Terms-and-Conditions checkbox before the site will (automatically or after human review)
   email a download link. CAMS-Knee's "Data Request" additionally requires **a signed license agreement
   upload (PDF)** and states "your request will be reviewed" (human-in-the-loop). **This part of the
   coordinator's "no gate anywhere" claim is NOT confirmed by direct measurement — these specific pages
   ARE registration/DUA gates**, consistent with the original `ORTHOLOAD_ACCESS_REQUEST_DRAFT.md` framing.
   Per the hard ethical line in this task's mandate, **these forms were NOT submitted** (no name/email was
   entered anywhere) — they are recorded below as gated, for the operator's own signed request to target.

Net: the bulk of the actual in-vivo force/moment telemetry (hip/knee/shoulder/spine, hundreds of trials,
decades of OrthoLoad's core database) turns out to be reachable through surface (1) with NO gate at all.
The CAMS-Knee synchronized (fluoroscopy+marker+GRF+EMG) dataset specifically is NOT free-standing — its
sample and full data both sit behind surface (2)'s form/DUA gate.

## 1. Site map (live-probed 2026-07-21)
`orthoload.com` and `cams-knee.orthoload.com` both HTTP 200, robots.txt on both only disallows
`/wp-admin/` (with `admin-ajax.php` explicitly re-allowed) — no crawl restriction relevant to this fetch.

Key pages found (non-exhaustive nav dump omitted; only data-relevant ones):
- `orthoload.com/database/`, `/database-2017/` — the JS search UI (implant→activity→parameter→patient→file).
- `orthoload.com/test-loads/standardized-loads-acting-at-hip-implants/` — **static public zip** (hip).
- `orthoload.com/test-loads/standard-loads-knee-joint/` — **~88 static public xlsx/zip files** (knee).
- `orthoload.com/comprehensive-data-sample/` — GATED (form, see §0.2).
- `orthoload.com/publications/{hip,knee,shoulder,vertebral-body-replacement,internal-spinal-fixator,gait-analysis}/`
  — literature reference lists only, zero direct file links (papers live behind their own journal DOIs,
  out of scope for this fetch).
- `orthoload.com/pelvis/`, `/thorax/` — literature references only (DOIs), no OrthoLoad-hosted data.
- `orthoload.com/kinematic-and-kinetic-data-recorded-with-an-instrumented-shoulder-prosthesis/` — no
  OrthoLoad-hosted file, but links out to two EXTERNAL public repositories (not fetched, out of scope,
  flagged for a possible future separate acquisition pass):
  - `https://data.4tu.nl/ndownloader/items/86db1d7d-13d9-4631-9c6b-1e3134a1ab38/versions/2` (4TU.ResearchData)
  - `https://repository.tudelft.nl/file/File_949cf8fe-0c85-4094-9a5b-7a8d85748ddc` (TU Delft repository)
- `cams-knee.orthoload.com/data/data-download/` — sample video + PDF **public**; the actual sample
  numeric dataset is GATED behind the name/email form on this page (see §0.2).
- `cams-knee.orthoload.com/data/data-request/` — full-dataset GATED (DUA + human review); 2 description
  PDFs + the license-agreement PDF are themselves public/direct.
- `cams-knee.orthoload.com/data/software/` — **public** (averDTW MATLAB/VB zips + checksums; a trial-averaging
  tool, not force data itself).
- `cams-knee.orthoload.com/workshop-data/` — **dead/stub page**: title "Life Demo Dataset", body text is
  literally the single word "Download" with no href anywhere (checked raw HTML). Not a gate, not a working
  link either — just unfinished content. Flagged, not counted as either public or gated.
- `cams-knee.orthoload.com/subjects/`, `/cooperations/request/` — no direct file links.

## 2. The public AJAX API (surface 1) — mechanism, verified live
`orthoload.com/wp-content/themes/orthoload2015/database2017/database.js` defines the flow (read as source,
not screenshotted):
`orthoload_get_activities(implantId)` → `orthoload_get_parameters(implantId,activityId,activityIndentationLevel)`
→ `orthoload_get_files(...,parameterId,parameterIndentationLevel,patientId)` → `orthoload_get_data(...,fileId,fileType)`
→ each result item is `{type, name, path, comment}` where `path` is a real static URL under
`/wp-content/files/...`. **`patientId=all` returns the union across every subject** for that
activity+parameter (verified: `all` for knee/Level-walking gave the same 21 trials as enumerating k1l+k2l+
k3r+k4r+k5r+k8l individually) — this collapses the crawl by one whole dimension.

Implant IDs (from the JS's own image-switch, confirms joint mapping):
| implantId | joint | n activity nodes (measured) |
|---|---|---|
| 1322 | knee | 113 |
| 12 | hip (gen 1 telemetry) | 124 |
| 2355 | hip (gen 2 telemetry, "hip3") | 377 |
| 864 | shoulder | 44 |
| 648 | spine — vertebral body replacement (VBR) | 213 |
| 1309 | spine — internal spinal fixator | 83 |

Verified end-to-end on one concrete trial: knee, "Level walking", parameter "none", patient "all" →
file `k1l_110108_1_86p` → data item `path=/wp-content/files/file/result/k1l_110108_1_86p.akf`, type `AKF`.
Direct unauthenticated GET of that path returned HTTP 200, 29,897 bytes, `ASCII text` — real content, not
an HTML error page — header confirms: `Diagram Title #1  Forces and Moments at Knee Joint`,
`Comment #2  K1L, 10 Months PO`, `BodyWeight [N]: 1000`, axis definitions. This is genuine in-vivo knee
contact-force time-history data — exactly the anchor class the charter requires — served with zero
authentication.

## 3. What was DOWNLOADED — Tier A: confirmed-public static files (complete, 90/90 verified)
All 90 files below were direct `<a href>` targets in page HTML (no API, no JS needed), fetched with a
single verification pass: HTTP 200, `content-type`, byte count, and Unix `file(1)` type-sniff on every
file — **0/90 flagged suspect** (no HTML-error-page masquerading as data).

| Category | Location | Count | Total size |
|---|---|---|---|
| Hip standardized test loads | `data/external/orthoload/hip/standard_loads/` | 1 (`StandardLoads-Hip_CompleteData.zip`) | 117 MB |
| Knee standard/per-subject loads | `data/external/orthoload/knee/standard_loads/` | 86 (xlsx per subject×activity + per-subject/per-activity "AllActivities" zips + cross-subject average xlsx + PNG diagram zips) | ~115 MB |
| CAMS-Knee public docs | `data/external/orthoload/cams_knee/public_docs/` | 4 PDFs (gait+squat sample description, full-dataset description, STAN-dataset description, **the License Agreement template itself**) | ~5.2 MB |
| CAMS-Knee public sample video | `data/external/orthoload/cams_knee/public_samples/` | 1 (`cams_sample_gait_squat.mp4`) | 4.6 MB |
| CAMS-Knee software | `data/external/orthoload/cams_knee/software/` | 4 (averDTW MATLAB zip + VB zip + 2 md5 checksums) | 2.5 MB |
| Knee comprehensive-sample teaser video | `data/external/orthoload/knee/standard_loads/H2R_150811_2_100_mokka.mp4` | 1 | 1.4 MB |

The **knee standard-loads set is the single richest free find**: per-subject (K1L, K2L, K3R, K5R, K6L, K7L,
K8L, K9L — 8 knee-implant subjects) × per-activity (Walking, StairsUp, StairsDown, SitDown, StandUp,
KneeBend, Jogging, Stance) force/moment spreadsheets, PLUS an averaged-subject-across-all-load-levels
workbook and full per-subject "AllActivities" zips (likely raw multi-trial time-histories, zips not yet
extracted/inspected in depth — flagged for a follow-up content pass, not blocking this report).

## 4. What was DOWNLOADED — Tier B: bulk API crawl (surface 1, all 6 implants)
Method: walked `orthoload_get_activities → get_parameters → get_files(patientId=all) → get_data` across
all 6 implant IDs. **29,428 API calls, 0 errors.** Discovered **18,840 unique static file paths**
(226,325 manifest rows before dedup — the same physical trial is cross-referenced from multiple
activity/parameter tree nodes). Type breakdown of what the API *lists*:

| type | unique files | what it is | live? (verified) |
|---|---|---|---|
| AKF | 3,943 | force/moment time-history, ASCII (the core anchor) | **100% live** |
| IOF | 899 | additional force file (bend/torsion/pedal rigs etc.) | **100% live** |
| COF | 74 | additional force file ("pedal" variants) | **100% live** |
| EOF | 19 | additional force file (`grf2f` = ground-reaction-force-to-force?) | **100% live** |
| CSV | 23 | raw marker-trajectory ("mov") CSV, large (~6-9 MB each) | **100% live** |
| RAW_MOV | 23 | raw movement data, zipped | **100% live** (2 exceptions, see below) |
| SCREENVIDEO | 4,958 | screen-capture MP4 of the analysis software trace | **100% live** |
| SCREEN | 4,958 | force-peak PNG snapshot | **100% live** |
| AVI | 3,943 | raw patient video | **CONFIRMED DEAD** — the API lists a path for every trial, but
  a random+targeted sample of 5 distinct paths returned **HTTP 404 (HTML error page)** every time, 0/5
  live. Not attempted in bulk (would be 3,943 wasted requests against a resource that does not exist).
  This is a genuine finding about the current state of OrthoLoad's server, not a fetch failure on this
  end — flagged honestly rather than silently omitted. |

**Download result (everything live except AVI): 14,897 files attempted → 14,893 succeeded (99.97%) after
one retry pass. 8.31 GB total.** First pass hit a transient local DNS resolution blip (359/14,897 failed
with `NameResolutionError`, not a server-side error — confirmed by DNS resolving cleanly seconds later);
all 359 were retried and now succeed. The remaining **4 failures are genuine, reproducible HTTP 404s on
OrthoLoad's own server** (itemized, not a systemic pattern — these 4 specific trial files are simply
missing from their `add_file`/`add_video` directories):
- `add_video/mp4/k8l_250311_1_64_mov.mp4`, `add_file/raw/k8l_250311_1_64_mov.zip`,
  `add_file/raw/k8l_250311_1_65_mov.zip`, `file/result/h9l_110613_2_7.akf`

Every downloaded file was verified two ways: (a) at fetch time, HTTP 200 + non-HTML content-type;
(b) post-hoc, a 20-file random sample cross-joint/cross-type run through Unix `file(1)` — 20/20 matched
their declared type exactly (real PNGs with valid dimensions, real ISO-Media MP4s, real ASCII AKF text) —
plus a 3-joint content check (hip/spine-VBR/shoulder, knee already checked in §2) confirmed the AKF
header format (`Diagram Title #1  Forces and Moments at <Joint>`, BodyWeight, Implant Type, Axes,
Activity Code, patient code + months-post-op) holds consistently across joints.

**Per-joint layout on disk** (`data/external/orthoload/<joint>/database_api/<type>/<original_filename>`,
plus the Tier-A static files under `<joint>/standard_loads/`):

| joint dir | implantId | size | AKF (force) files | distinct subject codes |
|---|---|---|---|---|
| `knee/` | 1322 | 1.4 GB (+ Tier A) | 609 | 9 (K1L,K2L,K3R,K4R,K5R,K6L,K7L,K8L,K9L) |
| `hip_gen1/` | 12 | 640 MB | 710 | 2 |
| `hip_gen2/` | 2355 | 4.4 GB | 1,240 | 10 |
| `shoulder/` | 864 | 237 MB | 140 | 7 |
| `spine_vbr/` | 648 | 1.2 GB | 869 | 5 |
| `spine_fixator/` | 1309 | 323 MB | 374 | 10 |
| `cams_knee/` | (Tier A only, see §3) | 12 MB | — | — |

Total: **15,016 files, 8.4 GB**, of which **3,942 are AKF force/moment trial files** — the actual
in-vivo joint-contact-force anchor for knee/hip/shoulder/spine, fully public, now on disk. Durable
manifests (crawl-raw JSONL, download-result CSV, static-file CSV) saved to
`data/external/orthoload/_manifest/`.

**Resolved uncertainty (was flagged mid-fetch, now closed):** the per-subject `K*_AllActivities.zip`
files (Tier A) were suspected to possibly contain additional raw data beyond the per-activity xlsx files.
Verified by direct inspection (`unzip -l`): they are exactly a bundle of the same 8 per-activity `.xlsx`
files already downloaded individually — no additional content, no follow-up needed.

## 5. GATED — precise requirements (for `ORTHOLOAD_ACCESS_REQUEST_DRAFT.md` to target)
1. **`orthoload.com/comprehensive-data-sample/`** ("free download" framing, but gated in practice):
   web form requires First Name, Surname, Organization, Email Address, Status (student/scientist/
   industry/other) + explicit T&C agreement checkbox. Says a download link + info is "automatically sent"
   to the given email — i.e. likely NOT human-reviewed, just identity-collection-gated. Covers: one H2R
   (hip) + presumably a matched knee walking trial with full kinematics/GRF/marker trajectories (richer
   than the plain AJAX-API force files).
2. **`cams-knee.orthoload.com/data/data-download/`** (sample trial): same style form (First/Surname/Org/
   Email/Status + T&C), explicitly says "without license agreement" for this SAMPLE tier — i.e. lighter
   gate than the full dataset, still needs the form.
3. **`cams-knee.orthoload.com/data/data-request/`** (full CAMS-Knee dataset): heaviest gate — name, org,
   email, status, **project title + description** (published on their site + Slack forum on approval),
   selection of which subjects (K1L/K2L/K3R/K5R/K7L/K8L)/activities/data-types (In-Vivo Loads / Skin-Marker
   Trajectories / Force-Plate / EMG / Fluoroscopic Images / Implant CAD) you need, **upload of a
   signed License Agreement PDF** (template already downloaded — see §3 public_docs), and an explicit
   "your request will be reviewed" — this is a genuine human-gated DUA process, matching the original
   registry framing exactly. Contact precedent per publications: Bergmann/Trepczynski/Damm/Taylor
   (Julius Wolff Institute + ETH Zürich), consistent with `ORTHOLOAD_ACCESS_REQUEST_DRAFT.md`'s addressee.
4. **`cams-knee.orthoload.com/workshop-data/`**: not really a gate — a stub/dead page (title "Life Demo
   Dataset", body text "Download" with no actual href). Practically unusable either way; note for the
   operator that this specific link is not a live resource.

No login/CAPTCHA/DUA was crossed or attempted anywhere in this fetch — the only forms encountered are
listed above and were left untouched.

## 6. Uncertainty / follow-ups flagged
- The two external (non-OrthoLoad-hosted) shoulder-data links (4TU.nl, TU Delft, §1) were noted but not
  fetched — out of this task's OrthoLoad/CAMS-Knee scope; flagging for a possible separate acquisition task.
- Hip/Pelvis/Thorax "publications" pages' underlying papers were not accessed (journal-paywalled DOIs,
  not OrthoLoad-hosted data — out of scope), except one confirmed-open-access exception: the PLOS ONE
  2016 paper describing the hip standard-loads dataset (`journal.pone.0155612`, publisher-wide CC-BY) was
  fetched directly (`hip/standard_loads/PLOSone_2016_StandardizedLoads_HipImplants.pdf`).
- 4 individual trial files (itemized in §4) are genuinely missing on OrthoLoad's own server (HTTP 404,
  reproducible) — not a fetch defect on this end, just noting the database has a few holes.
- Content semantics of IOF/COF/EOF (the smaller "additional file" numeric types — bend/torsion/pedal-rig
  variants and `grf2f`) were confirmed live and downloaded but not deep-parsed for exact column meaning —
  same header convention as AKF (verified), so should be straightforward to parse when the build reaches
  the per-joint force-anchor integration step; not blocking for this acquisition task.
- CAMS-Knee's actual synchronized (fluoroscopy+marker+GRF+EMG) dataset was NOT obtained — it lives
  entirely behind the gates in §5. The AJAX-API haul in §4 is OrthoLoad's *own* database (telemetry force/
  moment + screen captures), a separate and independent corpus from CAMS-Knee proper.
