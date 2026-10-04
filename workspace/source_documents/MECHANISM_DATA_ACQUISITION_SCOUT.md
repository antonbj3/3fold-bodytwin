# MECHANISM DATA ACQUISITION SCOUT — 5 targets, precise fetch status (2026-07-21)

Fetches only the highest-value in-vivo anchors + models the force-investigation docs
(`docs/MECHANISM_CONTACT_WAVEFORM.md`, `docs/MECHANISM_JAM_CONTACT_DECORR.md`,
`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`) flagged as open gaps. Report discipline per task: every claim below is
{fetched-and-verified | login-gated (exact step) | unavailable}, with byte-counts/row-counts/parseable-sample
proof — no claim of success without a file-layer check. This doc is a HYPOTHESIS for the coordinator to QC, not
a certified result.

## Summary table

| # | Target | Status | Path | Blocker / Next step |
|---|---|---|---|---|
| 1 | Grand Challenge (Fregly/D'Lima) in-vivo knee | **PARTIAL — fetched-and-verified for 1/6 subjects** | `data/external/opensim_jam_grand_challenge_DM/` (new, 104 MB), `data/external/simtk_kneeloads_DM_measured_force/` (copied+verified, 3.7 MB) | Other 5 subjects (Competitions 1-5): SimTK login-gated — free registration at `simtk.org/account/register.php`, join project `kneeloads` (group_id=413), then `frs/download_confirm.php` file IDs listed below become live |
| 2 | OrthoLoad hip contact-force waveform | **FETCHED-AND-VERIFIED — fixes a corrupted prior download** | `data/external/orthoload/hip/standard_loads_refetch/StandardLoads-Hip_CompleteData.zip` (205 MB, valid) | None — ready to consume; old corrupted 117 MB copy left untouched at the original path |
| 3 | Multi-segment foot OpenSim model | **login-gated (real project confirmed); OA literature anchor fetched** | No raw files (blocked) — but see manifest for the paper | SimTK project `footankle_model` (group_id=2092) — same free-registration unblock as #1; OR use the gold-OA paper directly (DOI 10.1080/10255842.2021.1968844, live at `research-repository.griffith.edu.au`) as a literature-only anchor now |
| 4 | Subject-specific/higher-fidelity deformable knee model | **FETCHED-AND-VERIFIED** (bundled with #1) | `data/external/opensim_jam_grand_challenge_DM/DM.osim` + implant `Geometry/` | None for this specific model — it is genuinely subject/implant-specific (real TKA CAD: Tibial Tray/Insert, Femoral Component, Patellar Button), not the generic lenhart2015 donor |
| 5 | Handsfield et al. 2014 MRI muscle-volume data | **FETCHED-AND-VERIFIED — primary source, upgrades the existing third-party transcription** | `data/external/handsfield2014/supplementary/mmc1.docx` (44,823 B) | None to fetch further (only 1 supplementary file exists, mmc2-5 confirmed 404); a clean cell-by-cell parse of Table 3 vs the existing CSV is a scoped follow-up, not a fetch gap |

**Isolation respected throughout:** all writes under `bodytwin/data/external/` (gitignored) or `bodytwin/docs/`;
no file created or modified outside this repo; the pre-existing corrupted OrthoLoad zip and the pre-existing
`data/external/handsfield2014/PROVENANCE.md` were read but **not** touched (only files this session created
were written to). No git add/commit/push. No SimTK registration/login performed by this agent — every
login-gated result below was confirmed gated by testing the actual mechanism (JS-redirect to
`/account/login.php`), never assumed from memory alone.

---

## 1. Grand Challenge Competition to Predict In Vivo Knee Loads (Fregly/D'Lima) — flagship anchor

**simtk.org/projects/kneeloads is reachable this session** (HTTP 200 — a change from a prior session's
byte-confirmed DNS/TCP-timeout outage recorded in `docs/MECHANISM_FOOT_FIDELITY_PLAN.md`). The project page
itself (description, file **names**) is openly browsable; the file **bytes** are not:

```
curl https://simtk.org/frs/download_confirm.php/file/3555/OpenSim%20Model.zip?group_id=413
-> HTTP 200, 195-byte body:
   <script>window.top.location='/account/login.php?triggered=1&return_to=...'</script>
```
This is a machine-confirmed login-redirect, not a bare 403 or a summarizer's guess — re-confirms the same gate
`bt_memory/knee-grand-challenge-unblock-recipe.md` documented in a prior session.

**The page's own file listing reveals the true scope**: 6 competitions (First-Sixth), each with its own
subject and up to 9 data-type zips (Raw CT, Raw Motion, Raw Strength, Synchronized Motion, Synchronized
Strength, Geometry, Raw X-Ray, OpenSim Model, Raw MR, Raw PreOp/PostOp CT for Competition 4 specifically).
Selected file IDs (group_id=413), for the operator to use once registered:

| Competition | File | file ID |
|---|---|---:|
| Sixth | Geometry Data.zip | 1925 (or 2795/3443, later versions) |
| Sixth | Raw CT Data.zip | 1926 (or 2796) |
| Sixth | OpenSim Model.zip | 3555 |
| Fourth | OpenSim Model.zip | (listed same page, "3. Data for Fourth Competition: OpenSim Model.zip") |
| First-Fifth | Geometry / Raw CT / Raw Motion / Synchronized Motion / Strength Data | all listed, all same login gate |

**Exact unblock recipe for the operator** (a human action this agent will not perform — no self-registration,
no ToS agreement, per this repo's own established doctrine): (1) `simtk.org/account/register.php` — free,
captcha-gated, no institutional affiliation required; (2) join/request project `kneeloads`; (3) re-fetch with
the session cookie; the exact file IDs above are ready to use, no further searching needed.

### What WAS obtained without any login (Apache-2.0, GitHub, `opensim-jam-org/jam-resources`)

The OpenSim-JAM team independently re-published subject **"DM"** (forensically confirmed = the **Sixth
Competition** subject, via the `.trc` files' own embedded path `...LMB_grand_challenge_database\...\
SixthCompetition\DM_ngait_og1_new.trc`) as an open, Apache-2.0-licensed model + experimental data bundle. Fetched
this session (see `data/external/opensim_jam_grand_challenge_DM/PROVENANCE.md` for full detail):
- Subject-specific implant model `DM.osim` (993,340 bytes, well-formed OpenSim 4.0 XML) + its own TKA
  CAD geometry (Tibial Tray/Insert, Femoral Component, Patellar Button — real prosthesis component meshes).
- **17 trials** of paired ground-reaction-force (`.sto`) + full-body marker motion (`.trc`, 120 Hz, 68 markers)
  — i.e. genuine "simultaneous gait" data for a real subject.
- 87 MB of motion data, verified parseable (header row/column counts match declared `nRows`/`nColumns` and
  frame counts exactly).

Combined with the **pre-existing** (not fetched by this session — copied+verified, see
`data/external/simtk_kneeloads_DM_measured_force/PROVENANCE.md`) **47-trial measured tibiofemoral contact-force
CSV set** (real `Fx,Fy,Fz,Tx,Ty,Tz,GON,GRFz` columns, MD5-verified byte-identical to its external-drive source),
**17 of those 47 trials cross-match by exact name** with the freshly-fetched open motion data — a real,
non-circular over-determination (two independently-obtained sources agree on trial identity).

**Net result for the flagship anchor: 1 of 6 subjects (DM/Sixth Competition) is now fully in-repo with
model + paired motion + measured force, zero login required. The other 5 subjects need the operator's own
SimTK registration** — this is the single highest-value remaining human action for the whole scouting task.

---

## 2. OrthoLoad hip contact-force waveform — corrupted download FIXED

The existing `data/external/orthoload/hip/standard_loads/StandardLoads-Hip_CompleteData.zip` (117,398,546
bytes) was already forensically confirmed truncated in `docs/MECHANISM_CONTACT_WAVEFORM.md` (valid ZIP magic,
no end-of-central-directory record). Re-fetched the **identical URL**
(`https://orthoload.com/wp-content/uploads/StandardLoads-Hip_CompleteData.zip`) fresh this session:

- **205,345,463 bytes** (1.75x the old file — confirms the old one was an interrupted transfer, not a
  server-side problem).
- Verified via Python `zipfile`: **504 entries**, `testzip()` returns `None` (zero CRC failures across every
  member — full-archive integrity, not a sample check).
- **No login/request-form needed** — straight anonymous HTTPS GET, HTTP 200.

New file: `data/external/orthoload/hip/standard_loads_refetch/StandardLoads-Hip_CompleteData.zip`. The old
corrupted file was left untouched (isolation: only touch files this agent creates) — coordinator decides
whether to supersede it.

---

## 3. Multi-segment foot OpenSim model

`docs/MECHANISM_FOOT_FIDELITY_PLAN.md` (prior session) flagged a specific project slug, `footankle_model`, as
**"UNVERIFIED — treat as FALSE until re-checked"** because simtk.org itself was down at the time (DNS
resolved but TCP timed out) and a WebFetch summarizer had plausibly hallucinated matching details.

**This session, simtk.org is reachable, and `footankle_model` is REAL** — HTTP 200, genuine project content
(not a generic error page): *"SimTK: Multi-segment Foot and ankle model validated using biplanar
videoradiography: Project Home"* — *"A multi-segment foot and ankle model consisting of the tibia, talus,
calcaneus, midfoot, forefoot and toes, with a total of 7 degrees of freedom"* — matching exactly the
"Maharaj JC model" description already in this repo's own foot-fidelity survey. Real downloadable file list
recovered (group_id=2092): `Maharaj2021.osim`, `Maharaj2021_BothLegs.osim`, `Maharaj2021_BothLegs_V2.osim`,
bilateral `FootModel_{calcn,talus,digits,forefoot,midfoot}_{l,r}.stl`, sample `static.c3d/.trc` +
`walk.c3d/.trc/_ik.mot`, `SampleData_SCALED.osim`.

**Login-gate confirmed by direct test** (not assumed from the kneeloads pattern):
```
curl https://simtk.org/frs/download_confirm.php/file/6470/Maharaj2021.osim?group_id=2092
-> same JS redirect to /account/login.php
```

**Citation recovered and its own paper is fully open** (a genuine bonus find, via NCBI eutils + Unpaywall,
zero login): **Maharaj et al., "Modelling the complexity of the foot and ankle during human locomotion: the
development and validation of a multi-segment foot model using biplanar videoradiography," *Computer Methods
in Biomechanics and Biomedical Engineering*, 2022. DOI: 10.1080/10255842.2021.1968844.** Unpaywall confirms
**gold open access**; live-verified (`hdl.handle.net/10072/415568` -> HTTP 200 ->
`research-repository.griffith.edu.au`) — the paper itself (methods, validation numbers) is freely readable
right now even though the raw `.osim`/mesh files remain SimTK-login-gated.

---

## 4. Subject-specific / higher-fidelity deformable-cartilage knee model

Resolved as part of #1: `DM.osim` (Grand Challenge Sixth-Competition subject) is a genuinely different,
implant-specific model from the generic 62 kg `lenhart2015.osim` already in use
(`docs/MECHANISM_JAM_CONTACT_DECORR.md`) — it carries real TKA prosthesis CAD geometry (Tibial Tray, Tibial
Insert, Femoral Component, Patellar Button — named exactly as such in the mesh files), not a native/healthy
knee. Apache-2.0, already fetched (see #1 for verification detail).

**Checked and ruled out as a further upgrade (an honest negative, not chased further):** `jam-resources` also
ships `models/knee_healthy/smith2019/smith2019.osim`, a candidate alternative to lenhart2015. Measured its
total body mass by summing every `<mass>` tag: **62.2 kg** — within 0.3% of lenhart2015's own implied 62.0 kg
body weight (`docs/MECHANISM_JAM_CONTACT_DECORR.md`: 608.01 N / 9.81). This is consistent with smith2019 being
the **same or a closely related donor/lineage** as lenhart2015 (same OpenSim-JAM research group), not an
independent subject — so it does **not** resolve the donor-diversity gap the way `DM.osim` does. Not fetched
(would not add real independence).

---

## 5. Handsfield et al. 2014 MRI muscle-volume data — primary source recovered

The existing `data/external/handsfield2014/PROVENANCE.md` (pre-existing, not created this session) already
established the main article (PMID 24368144, DOI 10.1016/j.jbiomech.2013.12.002) is paywalled (Unpaywall,
Semantic Scholar, PMC all checked negative) and that the in-repo CSV is a **third-party transcription**.

**This session forced that "unavailable" finding one step further** (per watertight discipline — a
leaning-negative claim needs its strongest fair channel tried) by re-confirming the paywall via 3 *different*
channels (Elsevier's own Content API, Europe PMC API, direct ScienceDirect fetch — all closed, consistent) —
**and then tried Elsevier's separate supplementary-content CDN directly**:
```
curl https://ars.els-cdn.com/content/image/1-s2.0-S0021929013006234-mmc1.docx
-> HTTP 200, application/word, 44,823 bytes
```
**This is the paper's own primary-source supplementary file, openly served despite the main article's
paywall.** Verified by extracting `word/document.xml` and reading the text (11,949 characters): it contains
**Supplemental Tables 1, 2, and 3** — Table 3 specifically is *"Coefficients for linear best fit for
individual muscle volumes against the metrics of total lower limb muscle volume, mass, and height*mass...
Volume = b1*metric + b2"* for all 35 muscles — **the exact table the existing third-party CSV states it was
derived from.** mmc2 through mmc5 (same URL pattern) all returned HTTP 404 — confirmed only one supplementary
file exists, not a fetch gap.

**Honest limitation, disclosed, not hidden:** a clean, cell-by-cell machine parse of Table 3 (needed to
numerically cross-validate every one of the third-party CSV's `b_slope`/`b_intercept` values against this
primary source) was attempted but not completed this session — the docx table did not cleanly separate under
a `</w:tc>`/`</w:tr>` boundary-marking approach. The content is verified to exist and be correct in substance
(human-readable in `mmc1_raw_flattened_extract.txt`); a rigorous numeric cross-check needs a proper
`python-docx`/XML-table parse, flagged as a concrete follow-up, not claimed as done.

---

## Files (all new, all under `bodytwin/data/external/` or `bodytwin/docs/`)

- `data/external/opensim_jam_grand_challenge_DM/` — DM.osim, markers, reserve actuators, selective implant
  Geometry (11 meshes), fluoroscopy data, 17 trials of motion/GRF data, `PROVENANCE.md`. 70 files, 104 MB.
- `data/external/simtk_kneeloads_DM_measured_force/` — 47 measured-force CSVs (copied+MD5-verified),
  `PROVENANCE.md`. 47 files, 3.7 MB.
- `data/external/orthoload/hip/standard_loads_refetch/StandardLoads-Hip_CompleteData.zip` +
  `PROVENANCE.md`. 205 MB, verified-valid ZIP (504 entries, 0 CRC failures).
- `data/external/handsfield2014/supplementary/mmc1.docx` + `mmc1_raw_flattened_extract.txt` +
  `PROVENANCE.md`. 44,823-byte primary source + 11,949-char verification extract.
- `docs/MECHANISM_DATA_ACQUISITION_SCOUT.md` — this file.
- `data/MECHANISM_DATA_ACQUISITION_SCOUT_manifest.json` — machine-readable version of the summary table.

## What this session deliberately did NOT do (scope discipline)

- Did not register a SimTK account or agree to any ToS on the operator's behalf (targets #1 remainder and #3).
- Did not attempt CAMS-Knee's full dataset (`cams-knee.orthoload.com`) beyond confirming its gate structure —
  already-fetched sample data + description PDFs pre-exist in this repo
  (`data/external/orthoload/cams_knee/`); full access needs a signed data-license request
  (`cams-knee.orthoload.com/data/data-request/`, contact form to Dr. Adam Trepczynski/Dr. Philippe Damm's
  group at Charité, `jwi.charite.de`) — a human-owned-identity action, not chased further as it would just
  re-confirm an already-known gate.
- Did not chase the "the reference model knee-grand-challenge-4th" mirror (subject "JW", proprietary license) already
  flagged in `bt_memory/SESSION_HANDOFF.md` — re-confirming an already-known gate would not add new value.
- Did not attempt a full cell-by-cell parse of Handsfield Table 3 (see honest limitation, §5) — out of scope
  for a fetch task; flagged as a follow-up.
- Did not run any consuming analysis script (e.g. re-running `contact_waveform_analysis.py` against the fixed
  hip zip, or building a new knee-force cert against `DM.osim`) — this task was scoped to fetch+verify, not
  to re-certify; those are natural next steps for whoever picks this up.
