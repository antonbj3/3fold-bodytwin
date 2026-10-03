"""Run the preregistered checks under a wall-clock budget and write ../results.json.

Usage:  python3 run_all.py [budget_seconds]     (default 240)
Sections run in priority order; results.json is rewritten after every section,
and any section/case not reached before the deadline is listed in "not_run".
"""
import json, math, os, platform, random, sys, time
import cellcouple as c

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 240.0
T0 = time.time()
TOL_BAL, TOL_MESH = 1e-6, 0.02
MET = ["gradient_index", "effectiveness", "ATP_ADP_cons", "T_cons"]
STRESS = dict(c.BASE, Da_c=100.0, Da_p=1000.0)          # thin consumption layer
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results.json")
out = dict(task="CLOUD-W3-CELL-COUPLE", python=platform.python_version(), budget_s=BUDGET,
           baseline_params=c.BASE, stress_params=STRESS, not_run=[],
           note="All parameter values are assumed nondimensional groups; no empirical calibration.",
           data_access=dict(attempted_hosts=["bionumbers.hms.harvard.edu", "api.crossref.org", "doi.org",
               "europepmc.org", "www.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "www.ebi.ac.uk",
               "openorganelle.janelia.org", "www.allencell.org", "zenodo.org", "en.wikipedia.org", "arxiv.org",
               "www.biorxiv.org", "api.openalex.org"],
               outcome="all blocked by the session egress proxy; no external data or literature was read",
               matched_reference=None))


def left():
    return BUDGET - (time.time() - T0)


def save():
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)


def section(name, fn, need=5.0):
    if left() < need:
        out["not_run"].append(name); save(); return
    t = time.time()
    try:
        fn()
        out.setdefault("section_time_s", {})[name] = round(time.time() - t, 1)
    except Exception as e:  # record, never hide
        out.setdefault("errors", {})[name] = repr(e)
    save()


# 1. Conservation: adenylate total and ATP budget during transients
def s_conservation():
    runs = []
    for name, p in (("baseline", c.BASE), ("stress", STRESS)):
        for N, dt, ns in ((112, 0.01, 300), (448, 0.01, 100), (112, 0.1, 100), (112, 0.001, 300)):
            if left() < 3:
                out["not_run"].append("conservation:%s:N%d:dt%g" % (name, N, dt)); continue
            r = c.transient(N, p, dt=dt, nsteps=ns); r["case"] = name; runs.append(r)
    worst = max(max(r["max_rel_adenylate_drift"], r["max_rel_ATP_budget_residual"]) for r in runs)
    out["conservation"] = dict(runs=runs, worst=worst, pass_=worst < TOL_BAL)


# Counterexample A: non-telescoping fluxes (cell-centre instead of face areas)
def s_ce_nonconservative():
    ce = c.transient(112, c.BASE, dt=0.01, nsteps=300, wrong_area=True)
    ce["detected"] = max(ce["max_rel_adenylate_drift"], ce["max_rel_ATP_budget_residual"]) > TOL_BAL
    out["counterexample_nonconservative_scheme"] = ce


# Counterexample B: Newton + one huge first pseudo-time step -> spurious negative-ADP steady state
BAD = dict(c.BASE, Da_c=2.4631361076245835, Da_p=135.95554400119, kappa=0.0509693468405705, delta=1.0758512495615127)
def s_ce_spurious():
    out["counterexample_spurious_state"] = dict(params=BAD,
        newton_bigstep=c.run_metrics(112, BAD, dt0=1e2, newton=True), picard_ramped=c.run_metrics(112, BAD),
        note="both variants conserve to round-off; only a positivity check separates them")


# 2. Mesh convergence (grids aligned with region interfaces: h = 0.7/N)
def convergence(p, Ns, key):
    full, done = [], []
    for N in Ns:
        if left() < 3:
            out["not_run"].append("%s:N%d" % (key, N)); continue
        full.append(c.run_metrics(N, p)); done.append(N)
    res = {}
    for k in MET:
        f = [r[k] for r in full]
        d = dict(values=f)
        if len(f) >= 2:
            d["rel_change_finest"] = abs(f[-1] - f[-2]) / abs(f[-1])
        if len(f) >= 3 and f[-2] != f[-1] and f[-3] != f[-2]:
            d["observed_order"] = math.log(abs((f[-3] - f[-2]) / (f[-2] - f[-1])), 2)
        res[k] = d
    r = dict(grids=done, metrics=res, min_conc=min(x["min_conc"] for x in full) if full else None)
    if len(done) >= 2:
        r["worst_rel_change"] = max(v["rel_change_finest"] for v in res.values())
    return r


def s_mesh_baseline():
    out["mesh_convergence_baseline"] = convergence(c.BASE, [28, 56, 112, 224, 448], "mesh_baseline")


def s_mesh_stress():
    out["mesh_convergence_stress"] = convergence(STRESS, [28, 56, 112, 224, 448], "mesh_stress")


# 3. Local elasticities d ln(metric)/d ln(param), central +-1 %, N = 112
def s_sensitivity():
    out["baseline_metrics_N112"] = c.run_metrics(112, c.BASE)
    sens = {}
    for k in ["Da_c", "Da_p", "kappa", "delta", "m_out", "c_in", "r_nuc"]:
        if left() < 3:
            out["not_run"].append("elasticity:" + k); continue
        hi, lo = dict(c.BASE), dict(c.BASE); hi[k] *= 1.01; lo[k] /= 1.01
        mh, ml = c.run_metrics(112, hi), c.run_metrics(112, lo)
        sens[k] = {m: (math.log(mh[m]) - math.log(ml[m])) / (2 * math.log(1.01))
                   for m in ["gradient_index", "effectiveness", "ATP_ADP_cons"]}
    out["elasticities"] = sens


# 4. Diffusion-limitation sweep at fixed Da_p/Da_c = 10 (well-mixed state fixed), N = 112
def s_sweep():
    rows = []
    for Da in [0.01, 0.1, 1, 3, 10, 30, 100]:
        if left() < 3:
            out["not_run"].append("Da_sweep:%g" % Da); continue
        m = c.run_metrics(112, dict(c.BASE, Da_c=Da, Da_p=10 * Da))
        rows.append(dict(Da_c=Da, **{k: m[k] for k in ("gradient_index", "effectiveness", "T_wellmixed", "T_cons", "min_conc")}))
    out["Da_sweep"] = rows


# 5. Prior propagation (synthetic): groups log-uniform in [x/3, 3x], N = 56
def s_prior():
    rng = random.Random(20260924); gi, ef, mn = [], [], []
    for _ in range(200):
        if left() < 2:
            break
        p = dict(c.BASE)
        for k in ["Da_c", "Da_p", "kappa", "delta"]:
            p[k] = c.BASE[k] * math.exp(rng.uniform(-math.log(3), math.log(3)))
        m = c.run_metrics(56, p); gi.append(m["gradient_index"]); ef.append(m["effectiveness"]); mn.append(m["min_conc"])
    pct = lambda v, q: sorted(v)[min(len(v) - 1, int(q * len(v)))]
    out["prior_propagation"] = dict(n_completed=len(gi), n_planned=200, seed=20260924, factor=3,
        gradient_index={q: pct(gi, float(q)) for q in ["0.05", "0.5", "0.95"]} if gi else None,
        effectiveness={q: pct(ef, float(q)) for q in ["0.05", "0.5", "0.95"]} if ef else None,
        min_concentration_over_draws=min(mn) if mn else None)


for name, fn in [("conservation", s_conservation), ("ce_nonconservative", s_ce_nonconservative),
                 ("ce_spurious", s_ce_spurious), ("mesh_baseline", s_mesh_baseline),
                 ("mesh_stress", s_mesh_stress), ("sensitivity", s_sensitivity), ("Da_sweep", s_sweep),
                 ("prior_propagation", s_prior)]:
    section(name, fn)

cons = out.get("conservation", {}).get("pass_")
mesh = [out.get(k, {}).get("worst_rel_change") for k in ("mesh_convergence_baseline", "mesh_convergence_stress")]
mesh_ok = None if None in mesh else max(mesh) < TOL_MESH
v = lambda b: "UNKNOWN" if b is None else ("PASS" if b else "FAIL")
out["criteria"] = dict(balance_lt_1em6=v(cons), mesh_convergence_lt_2pct=v(mesh_ok),
                       numerical_criterion=v(None if cons is None or mesh_ok is None else (cons and mesh_ok)),
                       empirical_calibration="UNKNOWN")
save()
print(json.dumps(out["criteria"]), "not_run:", out["not_run"], "runtime", out["runtime_s"], "s")
