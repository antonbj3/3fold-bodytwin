BT-DAT-Q014
# PUBLIC measured datasets constraining the Q014 model — inventory + one verified sample

`PREREG.md` (sha256 `0cb4b71948ef587a9394ae25f2b104a1e2e4d2f12f130c2b8e0779fb3d621b70`, frozen
before the first download) · `DATA_SOURCES.json` · `results.json` (all numbers) ·
`samples/` + `samples/LOAD_CHECK.json` · `analyse_samples.py`.
Run: 1 thread, `nice -n 19`, < 60 s, < 100 MB. No cloud run was needed.

## 0. Continuation, not restart

The previous interrupted session (read from `agent.log`, 7 lines, 2 steps) only did `ls` and
read `BRIEF.md`. **No artefacts existed** — no `PREREG.md`, no `samples/`,
no download. Everything below is new; nothing could be reused.

## 1. What is constrained (parameter audit of `inputs/Q014_model.py`)

Q014 has **zero anatomical parameters**. The entire parameter vector is
`{t_p, t_w, t_o, tau, a_A, a_B, G, sigma_w, n_seq, alpha}` — scalar or time-based
(`PARAMETER_TABLE`, lines 58–69). Therefore geometry/anatomy datasets by construction
cannot touch any of them (FC-D6), and therefore the choice of *dataset* is decisive.

| Model parameter | Owner-declared state | Which measured data can pin it | Found? |
|---|---|---|---|
| `tau` | **UNVERIFIED, free parameter, never measured** | within-subject decay time series, log-linear slope = `-1/tau` | measured in a *different* system, see §3 |
| `t_w` | SOURCE R1 Table 4 = 0 d | cohort/design metadata | **yes, 11 values** |
| `t_p` | DERIVED from "3 months" | design metadata | yes (e.g. 4 v + 1 v + 4 v) |
| `t_o` | ASSUMPTION = `t_p` | design metadata (outcome time point) | partial |
| `a_A` | IDENTIFIED from R1 Tab 2A | individual cell means in AB/BA cohort | no (arm level only) |
| `a_B` | ASSUMPTION = 0 | cell means (comparator) | no (arm level only) |
| `G` | CONVENTION = 1 | independent observable | **UNKNOWN, and searchable** |
| `sigma_w` | DERIVED from R1 SE | cell SD / within-subject | no (arm level only) |
| `n_seq` | SOURCE R1 (7, 9) | enrolled count per sequence | **yes** |

## 2. Source register (core kernels, IC1–IC6 in PREREG)

| # | Source | HTTP / size | Licence | Format | Subjects | Constrains |
|---|---|---|---|---|---|---|
| 1 | **ClinicalTrials.gov API v2**, PD + `crossover` | **200**, 822 770 B | NLM public data, no registration | JSON v2 | 50 hits, **46 CROSSOVER**, n = 8–94, median 24 | `t_w`, `t_p`, `n_seq` |
| 2 | **Remifentanil PK** (Pinheiro & Bates 2000, via R `nlme`) | **200**, 124 398 B (0.12 MB) | `nlme` GPL-2\|GPL-3; per-dataset licence **unspecified → UNKNOWN** | CSV 12 cols | **65 patients**, 2107 points | `tau`, the model’s single-`tau` assumption |
| 3 | Quinidine PK (via `nlme`) | **200**, 92 228 B | as above → UNKNOWN | CSV 14 cols | 136 patients | `tau` — **failed as N-B replicate** |
| 4 | Tetracycline1 (via `nlme`) | **200**, 881 B | as above → UNKNOWN | CSV 5 cols | 5 subjects | nothing (too few time points) |
| 5 | PhysioNet MIMIC-IV v2.2 | **200** (landing page 73 729 B) | Credentialed Health Data License + CITI/DUA | CSV/Postgres | ~65 000 patients (documented) | `tau` — **best candidate, not downloaded** |
| 6 | Zenodo 167808 (35 femur+tibia) | **200** (page 83 917 B) | per record, unverified → UNKNOWN | STL/PLY | 35 (documented) | **nothing — scoped out** (FC-D6) |
| 7 | PPMI | **DNS error** (`Could not resolve host`) | data-use agreement | CSV | UNKNOWN | `a_A`, `sigma_w` — **not admitted (IC1 failure)** |
| 8 | OpenCap | **DNS error** | UNKNOWN | time series | UNKNOWN | **scoped out** |
| 9 | R1 Schissler 2023 (`10.3389/fneur.2023.1197281`) | not retrieved again here | Frontiers OA, CC-BY | published table | 16 (7 + 9) | `a_A`, `sigma_w`, `n_seq`, `t_p`, `t_w` |

Full fields per row are in `DATA_SOURCES.json`. Rows 7–8 remain as **leads**,
not as sources, with HTTP status verbatim (`null` + reason).

## 3. Sample loaded and verified (IC7 / FC-D1)

`samples/Remifentanil.csv` — **124 398 B = 0.12 MB, far below the 50 MB cap.**

* **Form:** 2107 rader × 12 kolumner, index `rownames`; `ID, Subject, Time, conc, Rate,
  Amt, Age, Sex, Ht, Wt, BSA, LBM`. dtypes: 8 × float64, 2 × int64, 1 × object (`Sex`),
  1 × int64. Missing `conc` in 115/2107 rows (5.5 %, documented, not removed).
  `conc` 0.1–245.4 ng/ml. 20–54 obs/patient.
* **Units** (read from `rdrr.io/cran/nlme/man/Remifentanil.html`, HTTP 200, **not**
  guessed): `Time` = **min from start of infusion**, `conc` = **ng/ml**,
  `Rate` = µg/min, `Amt` = µg, `Age` = years, `Ht` = cm, `Wt` = kg.
* **Coordinate frame: `null`.** 1-D time series (time versus concentration) — a spatial frame is
  defined only for 1-D and **was not fabricated** (PREREG IC7). Same for the other three
  samples. The entire machine description in `samples/LOAD_CHECK.json` (shape, dtypes,
  head, sha256 per file).
* Other samples: `Quinidine.csv` (92 228 B), `Tetracycline1.csv` (881 B),
  `Nitrendipene.csv` (1 644 B, **discarded on inspection** — it is a tissue/enzyme table,
  `activity, NIF, Tissue, log.NIF`, not a PK series).

## 4. What the measurement gives — `tau` from decay (FC-D2, N-A)

The identification method is not mine: the owner’s own PREREG §4 requires "log-linear
slope = `-1/tau`". The segment is per patient from `argmax(conc)` to the last measurement,
`conc>0`, ≥ 5 punkter, `r² ≥ 0.5`.

| Quantity | Value | Source |
|---|---|---|
| Patients evaluated | **65 / 65** | `results.json:tau_remifentanil` |
| `tau` median (per patient) | **10.58 min**, IQR 9.70–12.09 | same |
| Median half-life | **7.33 min** | `tau*ln2` |
| `r²` median / min | **0.963 / 0.842** | same |
| Pooled slope (n = 1520) | **−0.0624 /min**, 95 % CI **[−0.0647, −0.0602]** | `results.json:pooled_remifentanil` |
| Pooled `tau` | 16.01 min, `r²` = 0.669 | same |

**Countertest N-A (placebo, labels shuffled within patient):** real median `r²`
**0.963** versus placebo mean **0.560** (Δ`r²` = **0.403**, 200 permutations) →
**PASS**. Identification thus reads signal, not noise. *Transparency:* placebo `r²`
reached 0.971 in one permutation — the mean is what distinguishes, individual rounds do
not do so unambiguously.

## 5. The most important finding — the port cannot be closed by a short memory time

Q014’s ports are **dimensionless**: `S = 1−e^{−t_p/tau}`, `R = e^{−(t_w+t_o)/tau}`,
`c = 1−e^{−t_o/tau}`. Thus a dataset in another system (ng/ml) can still constrain
the model — through the ratio `tau/t_p`, not through absolute days. Inserted into the model’s own
design block (`t_p = t_o = 90 d`, `t_w = 0`, from `Q014_model.py`):

| Measured `tau` | `tau/t_p` | `S` | `c` | `R` | Parenthesis `S(1+R)−c` |
|---|---|---|---|---|---|
| remifentanil median 10.58 min | **8.2e−5** | 1.0 | 1.0 | 0.0 | **0.0 exakt** |
| remifentanil max 27.19 min | 2.1e−4 | 1.0 | 1.0 | 0.0 | 0.0 exakt |
| kinidin median 4923 h | 3.8e−2 | 1.0 | 1.0 | 3.7e−12 | 3.7e−12 |

**Conclusion (derivation from the model’s own formula, not from measured data):** if memory time
is short — and it is short in every measured decay found, 7 minutes to 205 days
— then `S` and `c` both go to 1 simultaneously, the parenthesis collapses **exactly to 0**, and
`lam_hat = 0`. The order effect can exist only in the window `t_p ~ tau`. It is a
**negative result about Q014’s calibration point**: with a 90-day span and any
known adaptive memory time, port G2b is closed. In practice: the right dataset is not one
giving `tau` in days, but one where **the memory scale is comparable to the period**.

## 6. Hypotheses and verdicts (against N-B, N-C, N-D)

* **H1 — seek `tau` comparable to `t_p`: NOT MET.** No publicly measured dataset
  was found in the session with `tau` on the order of weeks–months. The closest is
  MIMIC-IV (row 5), which is *gated* and therefore could not be downloaded.
* **H2 — the design block against the background: MET, with numbers.** 46 registered
  AB/BA PD trials, `n` median 24 (the Q014 calibration has 16). Stated `t_w`:
  **48 h, 7 d, 14 d, 15±2 d, 1 v, 3 month**; NCT00305331: `t_p` = 4 v, `t_w` = 1 v, `t_p` = 4 v
  → **`t_w/t_p` = 0.25**, the only ratio the model actually needs. Only **11/46**
  mention any washout at all. **N-D: the calibration point `t_w = 0` is thus neither
  typical of the population nor even the most commonly stated value** — and "not mentioned" is
  not "sufficient", so 35/46 leave `t_w` as **UNKNOWN**, not as 0.
* **H3 / N-B — is a single `tau` sufficient: UNKNOWN according to PREREG (FC-D4 breached).**
  Remifentanil alone gives a clear positive dependence of measured `tau` on
  exposure duration: Spearman **ρ = 0.833, p = 7.9e−18, n = 65**; median `tau`
  9.94 min (short exposure) versus 12.10 min (long), OLS +0.108 min per min
  exposure, Mann–Whitney p = 2.7e−8. But that is **not** enough: quinidine gave only
  **2/136** usable patients (multi-dose, steady state), and nitrendipine revealed
  itself to be something else entirely. Without a second independent dataset, FC-D4 fails and H3 stands
  as **UNKNOWN** — I would rather report that than a half-supported falsification.
* **N-C — no `tau` from design metadata.** No washout value in §2 is used as
  memory time. The order `t_w` → `tau` does not occur in any row.

## 7. What failed

* **FC-D4 / H3** — lacks a second independent decay series. This is the single most important
  gap in the delivery.
* **H1** — negative result, see §5.
* **`G` remains UNKNOWN.** It is not a failure of the search: the latent state is seen only
  through the outcome, so `G` and `a` are identifiable **only as a product**. R1 gives
  arm-level means, not individual cells — therefore `sigma_w` had to be derived as `SE*sqrt(n)`.
  No public source was found that resolves `G`; claiming a value would be fabricated.
* **The licence for the `nlme` series is not established** (GPL on the package, no
  per-dataset licence). Reported as UNKNOWN, not as "free".
* PPMI and OpenCap: **not admitted** — DNS could not be resolved from the sandbox. They
  remain as leads, not as claimed sources.

## 8. Next steps (ordered by value, not by quantity)

1. **MIMIC-IV** (row 5, PhysioNet, gated) — the only candidate in the list whose time scale
   (hours–days) is on the same order as `t_p = 90 d`. Go through CITI/DUA.
2. **A second, independent decay dataset with the same design as remifentanil** (an
   infusion dataset, not multi-dose), to determine H3. Without it, FC-D4 remains
   UNKNOWN.
3. **Individual-level data for an AB/BA balance cohort** (PPMI if the wall opens) — the only
   way out of the `sigma_w` derivation and `G` insensitivity.
4. `t_o` is an **ASSUMPTION** in the model and it is a *design choice*; metadata for
   outcome time point (row 1) can test it without any new measurement.

---

*All numbers in §3–§6 are in `results.json` with a path; derivations (the gate parenthesis,
the ratio `t_w/t_p`) are labelled as derivation, not measurement. Source status per row: measured /
derived / assumption / UNKNOWN. No verdict words, no invented measured data.*
