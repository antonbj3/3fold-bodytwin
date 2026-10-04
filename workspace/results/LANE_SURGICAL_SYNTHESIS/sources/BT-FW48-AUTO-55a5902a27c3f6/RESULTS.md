# BT-FW48-AUTO-55a5902a27c3f6

Datum: 2026-09-30
Status: **PENDING_INDEPENDENT_REVIEW**. No physical validation claimed.

## 1. Decision addressed

Acquire the vessel map ONCE as a single measurement specification and stop. The
acquisition specification is delivered in §5 below, once. No vessel counts were
obtained and none are claimed: specimen access does not exist in this packet
(`JOB.json` `prereq`). The remaining budget was spent on a *different available
mechanism*: the weighting operator that converts a per-decade map into the flow
functionals the downstream SURG_* families consume.

## 2. Measured application gain (this is a code-correctness + operator result)

A dimension-algebra check (`unit_check()`, C1) **caught a real algebra error in
the PREREG text before results were reported**. PREREG §3.3 stated
`gamma_j = dP*R_j^3/(2 mu L_out)`, but `gamma = 4v/R` with
`v = dP R^2/(8 mu L)` gives `gamma = dP R/(2 mu L)`. The exponent was wrong by
two powers. This is recorded as **AMENDMENT 1** (§6); the corrected linear form
is used throughout. Had the check not been run, every shear-side number in this
job would have been wrong by `(R)^-2` — a factor of ~1e5 at 300 µm.

This is a genuine negative finding about the PREREG, preserved deliberately.

## 3. What was actually computed

**Declared baseline map** (ASSUMPTION, no provenance; PREREG §2.1):
`n = [1.5e4, 1.5e4, 3.0e2, 3.0e2]` per `1/m^2 per radius decade` for
`[3,10), [10,30), [30,100), [100,300) um`. Log-uniform within each of two
populations, amplitude ratio 50.

**Frozen relations** (PREREG §3), Poiseuille at a severed end:

    Q_vessel = pi*dP*R^4/(8*mu*L_out)
    q        = sum_j n_j * Q_vessel(R_j)            [m/s, per unit cut-face area]
    Q_total  = q * A_face                           [m^3/s]
    gamma_j  = dP*R_j/(2*mu*L_out)                  [1/s]   (AMENDED)

Ports: `L_out = 2.0e-3 m` (Q033 `Parameters.thickness_m`),
`A_face = pi*a*b = 4.698e-4 m^2` (Q033 `ellipse_geometry`, a=0.023, b=0.0065).
`mu = 3.5e-3 Pa s`, `dP = 1.0e4 Pa` are assumptions; both enter linearly and
do not change the *shares*, which are the quantity this job is about.

### Result 3.1 — Baseline flow (all four decades)

| quantity | value |
|---|---|
| `q` per unit cut-face area | **1.5375e-4 m/s** |
| `Q_total` over the Q033 ellipse | **7.2211e-8 m^3/s** |
| declared total count density | 3.06e4 hits/m^2 |

### Result 3.2 — THE load-bearing result: the two functionals weight the decades in OPPOSITE order

Flow share `s_j = Q_j/q` vs shear-weighted severed-end share `g_j` (the
quantity the platelet/hemostasis family needs for adhesion thresholds):

| decade (µm) | flow share `s_j` | shear-end share `g_j` |
|---|---|---|
| 3–10   | **0.0049 %** | **20.02 %** |
| 10–30  | 0.4926 %   | **63.31 %** |
| 30–100 | 0.9852 %   | 4.00 %  |
| 100–300| **98.517 %**| 12.66 % |

Leave-one-decade-out sensitivity on each functional (relative error in the
functional when the decade is removed):

| decade removed | error in `q` (flow) | error in shear-weighted count |
|---|---|---|
| 3–10 µm   | **−0.0049 %** | −20.02 % |
| 10–30 µm  | −0.4926 % | −63.31 % |
| 30–100 µm | −0.9852 % | −4.00 %  |
| 100–300 µm| **−98.52 %** | −12.66 % |

**Interpretation.** Flow is a `R^4`-weighted functional and is utterly
dominated by the 100–300 µm decade. The shear/platelet functional is `R`-weighted
and is dominated by the *two smallest* decades, which together carry 83 % of it.
So there is **no single "vessel map" that can be acquired to one accuracy**: the
map that fixes bleeding is nearly orthogonal to the map that fixes platelet
adduction. Any downstream port that consumes `q` alone and ignores `g` will be
systematically wrong about which vessels matter for clotting. This is a new
conjecture-level structural claim, not a measurement.

### Result 3.3 — Matched control (C5). The count-matched single log-uniform is WRONG by 25×

Strongest control available without a specimen: a **single** log-uniform
population over 3–300 µm, with its one free amplitude matched on the *total count
integral* `sum_j n_j` — the only scalar a section-plane histogram delivers
directly. (Matching the amplitude on `q` instead is circular: a single amplitude
reproduces any scalar `q` exactly, giving a trivially perfect ratio of 1.0. That
degenerate control was written first, and is recorded here as a rejected method,
§6 AMENDMENT 2.)

    control ratio  q_control / q_declared = 25.38

The count-matched single-population control **under-predicts flow by a factor of
25.4**. The JOB.json falsifier threshold is a factor of 2, so the falsifier is
**triggered**: the declared two-population structure IS load-bearing for flow,
and a single log-uniform fit is not an adequate substitute. Reported as an
observation, not as a pass/fail gate (PREREG C5).

### Result 3.4 — Box span and the bar the measurement must beat

The declared amplitude boxes (each a factor of 16 range) span a **factor of 16**
in `q` — i.e. the current assumption permits anything between 1/16 and 16× the
stated flow. The acquisition is worth doing only if it narrows this by more than
~2× per decade. That is the concrete bar, and it was not derivable without §3.2.

### Result 3.5 — Precision feasibility: one decade is unresolvable and must be bounded instead

Poisson counting, `rel SE = 1/sqrt(N)`, relative error contributed to `q` equal to
the decade's flow share, targeting 10 % on `q`:

| decade | counts needed | section area needed at declared `n` |
|---|---|---|
| 3–10 µm  | **4,121,306** | **274.8 m^2** — infeasible |
| 10–30 µm | 412 | 0.027 m^2 |
| 30–100 µm| 103 | 0.343 m^2 |
| 100–300 µm| <1 | negligible |

**This is the decisive practical finding.** The 3–10 µm decade would require
~275 m² of section plane — several whole-body-surface equivalents — to be
counted precisely enough to matter for flow. It is **not measurable** at the
precision anyone could achieve.

But §3.2 shows it contributes **0.0049 %** of `q`. So the correct engineering
answer is **not** to measure it, and **not** to assume a value for it either — it
is to carry it as an explicit **bounded interval** with a stated error budget, and
record the bound as a first-class result. The declared map's contribution from
this decade is far below any other uncertainty in the chain. This removes a
decade from the acquisition's burden and is the single most actionable output of
this job.

## 4. Capability delivered vs. unresolved

- **Useful capability (new)**: a closed-form map→functional weighting operator
  with two distinct functionals (`R^4` flow, `R` shear) whose decade
  sensitivities are computed exactly, giving per-decade precision targets and an
  explicit infeasibility verdict on the smallest decade.
- **Proved under assumptions** (PREREG §3, Poiseuille, declared map): flow
  dominance and shear dominance order as tabulated; the 25.4× control error; the
  274.8 m² infeasibility.
- **Unresolved question**: what the *actual* per-decade counts are. Untouched,
  by design.
- **No physical validation whatsoever.** Every number above is the output of the
  declared assumption map, and is only as good as that map.

## 5. THE ACQUISEMENT SPECIFICATION (delivered once; do not re-audit)

**Specimen** — resected human skin/subcutis specimen or fresh cadaver,
unfixed, or perfusion-fixed with the fixation stated.
**Section plane** — a single plane **perpendicular to the proposed incision
line** (i.e. a transverse cross-section of the intended cut). Report the plane
normal explicitly as a triad relative to (skin surface normal, Langer line
direction). This orientation is mandatory: a directional vessel map changes every
count, and it is also required to invert to volumetric density later.
**Layer partition** — epidermis excluded; **dermis**, **subcutis**, **fascia**
labelled separately per section, with layer boundaries marked on the image.
**Radius bins** — decades `[3,10)`, `[10,30)`, `[30,100)`, `[100,300)` µm, with
the 3–10 µm decade **bounded, not exhaustively counted** (§3.5).
**Reported unit** — counts per unit **cut-face area**, in absolute
`1/m^2 per radius decade`, per layer, per section plane. This is the frozen
definition; it is *not* volumetric density and must not be silently converted.
**Minimum replication** — ≥3 sections per specimen, ≥3 specimens, reported
individually before pooling.
**Orientation control (strongest control, from JOB.json)** — the same specimen
re-sectioned at a second plane orientation (e.g. 90° in-plane). Tests whether the
map is orientation-dependent. This is a required part of the acquisition, not
optional.
**Precision targets, per decade, for the flow functional** — from §3.2/§3.5:
`[3,10)` ≤ 0.01 % of flow (bound only, no target needed);
`[10,30)` ≤ 0.5 %; `[30,100)` ≤ 1 %; `[100,300)` ≤ **2 %** (dominant, and cheap
to measure). For the shear/platelet functional the targets invert:
`[3,10)` and `[10,30)` must be resolved to ~20 % and ~5 % respectively, or the
platelet-adhesion port stays uninterpretable.

Proposals 2, 8 and 9 become executable against this contract. The physical
acquisition itself is not executable in this packet.

## 6. Amendments and rejected methods (preserved)

- **AMENDMENT 1** — PREREG §3.3 shear exponent `R^3` → `R`, found by `unit_check()`
  C1 before reporting. PREREG text left as frozen; the correction lives in
  `vessel_map_flow.py:exit_shear` and here.
- **AMENDMENT 2** — the first control implementation matched the single-population
  amplitude on `q` itself, which is circular and returned a trivially perfect
  ratio of exactly 1.0. Replaced with matching on the total count integral. The
  degenerate version is recorded rather than deleted, because a 1.0 ratio from a
  control is exactly the kind of result that looks like success and is not.
- **REJECTED** — exhaustive counting of the 3–10 µm decade (275 m² required).
  Not run; infeasibility established analytically from `s_j`, not assumed.

## 7. Scope, lineage, limits

- Parents: `BT-HX-Q033` (`Q033_incision_model.py`, `Q033_PREREG.md`). Target
  context `BT-CTX-SURG-INCISION`; consumers `BT-CTX-SURG-HEMOSTASIS` and
  `BT-FW48-SURG-02`.
- Read: exactly 2 source files (Q033 model, Q033 PREREG) before code was written,
  per the bound.
- Found while reading: **Q033's transport model has no vessel term at all** —
  `transport_flow()` contains Darcy, evaporation and exudate only. So the vessel
  map is not merely imprecise in the current chain, it is **structurally absent**
  from `BT-HX-Q033`. The gap is larger than "an unmeasured parameter".
- Code unaudited. `mu` and `dP` are assumptions. All counts are declared
  assumptions with no provenance. No fitted constant was introduced after the
  PREREG freeze.

## 8. What remains incomplete

1. **Every physical count.** The entire map is unmeasured.
2. **Orientation dependence.** Unmeasured; §5 mandates the control but no data
   exists, so no volumetric inversion is attempted.
3. **The 3–10 µm decade bound.** Analytically negligible for flow but currently
   a declared value, not a bounded interval (§3.5 recommends the change; the
   bound itself is not yet specified numerically).
4. **Shear/platelet coupling.** `g_j` is computed but the downstream platelet
   adhesion threshold physics was not read or coupled in this job; the
   functional is offered as a port, not as a validated hemostasis model.
5. **No layer partitioning was computed at all** — the model is single-layer.
   The dermis/subcutis/fascia split is specified for acquisition (§5) but the
   weights cannot be computed without per-layer maps.