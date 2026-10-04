#!/usr/bin/env python3
"""CLOUD-A359: NHANES 1999-2002 knee-extensor strength vs DXA leg lean.

Usage:
  python3 a359.py fetch  DATA_DIR        # download CDC XPT files (needs cdc.gov access)
  python3 a359.py real   DATA_DIR OUT.json  # empirical analysis on downloaded files
  python3 a359.py synth  OUT.json        # SYNTHETIC pipeline self-test (not validation)

Design (frozen in RESULTS.md before any data was retrieved):
  train = 1999-2000 (DEMO, BMX, MSX, DXX), test = 2001-2002 (*_B), age >= 50,
  complete cases on every variable used by any model, same persons for all models.
  Primary DXA model: torque ~ legLean*height (+ intercept). Comparator: torque ~ mass.
  Pass if held-out RMSE improvement of primary DXA model vs mass-only >= 10%.
Requires numpy, pandas only. All randomness uses SEED.
"""
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

SEED = 20260924
N_BOOT = 2000
N_FOLDS = 5
IMPROVE_THRESHOLD = 0.10  # PREREG: >= 10 % RMSE improvement

BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes"
FILES = {  # name -> URL. Docs: same URL with .htm instead of .xpt
    "DEMO": f"{BASE}/Public/1999/DataFiles/DEMO.xpt",
    "BMX": f"{BASE}/Public/1999/DataFiles/BMX.xpt",
    "MSX": f"{BASE}/Public/1999/DataFiles/MSX.xpt",
    "DXX": f"{BASE}/Public/1999/DataFiles/DXX.xpt",
    "DEMO_B": f"{BASE}/Public/2001/DataFiles/DEMO_B.xpt",
    "BMX_B": f"{BASE}/Public/2001/DataFiles/BMX_B.xpt",
    "MSX_B": f"{BASE}/Public/2001/DataFiles/MSX_B.xpt",
    "DXX_B": f"{BASE}/Public/2001/DataFiles/DXX_B.xpt",
}
# Variable names as recalled from the CDC codebooks; the script asserts they exist
# and aborts otherwise (they MUST be checked against the .htm data dictionaries).
V = dict(
    id="SEQN", age="RIDAGEYR", sex="RIAGENDR", psu="SDMVPSU", strata="SDMVSTRA",
    wt="BMXWT", ht="BMXHT",               # kg, cm
    force="MSDPF",                         # knee-extensor isokinetic peak force, N
    imp="_MULT_",                          # DXA multiple-imputation index 1..5
    llean="DXDLLLE", rlean="DXDRLLE",      # left/right leg lean excl. BMC, g
)
# Lever-arm conversions force (N) -> torque (N m). No NHANES lever arm is published
# in the files; these are ASSUMPTIONS for sensitivity only.
LEVER = {
    "fixed_0.30m": lambda h_m: np.full_like(h_m, 0.30),
    "fixed_0.35m": lambda h_m: np.full_like(h_m, 0.35),
    "height_0.20H": lambda h_m: 0.20 * h_m,   # ~ 0.8 x shank length (0.246 H, Winter)
}
PRIMARY_LEVER = "fixed_0.30m"


def fetch(d):
    os.makedirs(d, exist_ok=True)
    for k, u in FILES.items():
        p = os.path.join(d, k + ".xpt")
        print("GET", u, flush=True)
        urllib.request.urlretrieve(u, p)


def load_cycle(d, sfx):
    rd = lambda k: pd.read_sas(os.path.join(d, k + sfx + ".xpt"), format="xport")
    demo, bmx, msx, dxx = rd("DEMO"), rd("BMX"), rd("MSX"), rd("DXX")
    for df, cols in [(demo, ["id", "age", "sex", "psu", "strata"]), (bmx, ["id", "wt", "ht"]),
                     (msx, ["id", "force"]), (dxx, ["id", "imp", "llean", "rlean"])]:
        miss = [V[c] for c in cols if V[c] not in df.columns]
        if miss:
            raise SystemExit(f"variables {miss} not found; check data dictionary")
    df = (demo[[V[c] for c in ["id", "age", "sex", "psu", "strata"]]]
          .merge(bmx[[V["id"], V["wt"], V["ht"]]], on=V["id"], how="left")
          .merge(msx[[V["id"], V["force"]]], on=V["id"], how="left"))
    dxx = dxx[[V["id"], V["imp"], V["llean"], V["rlean"]]].copy()
    return tidy(df, dxx)


def tidy(df, dxx):
    """Harmonise names; returns (person frame, long DXA frame with imputation index)."""
    df = df.rename(columns={V[k]: k for k in V})
    dxx = dxx.rename(columns={V[k]: k for k in V})
    dxx["leglean_kg"] = (dxx["llean"] + dxx["rlean"]) / 1000.0
    return df, dxx[["id", "imp", "leglean_kg"]]


def select(df, dxx, log):
    n0 = len(df); log["n_demo"] = n0
    df = df[df.age >= 50]; log["n_age50"] = len(df)
    log["missing_age50"] = {c: int(df[c].isna().sum()) for c in ["wt", "ht", "force"]}
    # same persons in every model and every imputation: leg lean present in all implicates
    ok = dxx.assign(v=dxx.leglean_kg.notna()).groupby("id").v.agg(["all", "size"])
    ids_dxa = set(ok.index[ok["all"] & (ok["size"] == dxx.imp.nunique())])
    log["missing_age50"]["leglean"] = int((~df.id.isin(ids_dxa)).sum())
    df = df.dropna(subset=["wt", "ht", "force", "sex", "psu", "strata"])
    df = df[(df.force > 0) & df.id.isin(ids_dxa)]
    log["n_complete"] = len(df)
    return df.reset_index(drop=True)


MODELS = {  # name -> design-matrix builder (intercept added automatically)
    "mass_only": lambda d: [d.wt],
    "mass_ht_sex_age": lambda d: [d.wt, d.ht_m, d.female, d.age],
    "leglean": lambda d: [d.leglean_kg],
    "leglean_x_ht": lambda d: [d.leglean_kg * d.ht_m],        # PRIMARY DXA model
    "leglean_x_ht_sex_age": lambda d: [d.leglean_kg * d.ht_m, d.female, d.age],
}
PRIMARY_DXA, COMPARATOR = "leglean_x_ht", "mass_only"


def X(name, d):
    cols = MODELS[name](d)
    return np.column_stack([np.ones(len(d))] + [np.asarray(c, float) for c in cols])


def fit_predict(name, tr, te, y):
    b, *_ = np.linalg.lstsq(X(name, tr), tr[y].to_numpy(), rcond=None)
    return X(name, te) @ b


def rmse(r):
    return float(np.sqrt(np.mean(np.square(r))))


def analyse(train_p, train_dxa, test_p, test_dxa, lever):
    """Per-imputation fit on train, predict test; pool predictions across imputations."""
    out = {}
    imps = sorted(set(train_dxa.imp) & set(test_dxa.imp))
    preds = {m: [] for m in MODELS}
    cv = {m: [] for m in MODELS}
    rng = np.random.default_rng(SEED)
    for k in imps:
        tr = prep(train_p, train_dxa[train_dxa.imp == k], lever)
        te = prep(test_p, test_dxa[test_dxa.imp == k], lever)
        # grouped CV inside training: groups = design PSU within stratum
        g = (tr.strata.astype(int) * 10 + tr.psu.astype(int)).to_numpy()
        ug = np.unique(g); fold_of = dict(zip(ug, rng.permutation(len(ug)) % N_FOLDS))
        f = np.array([fold_of[x] for x in g])
        for m in MODELS:
            res = np.empty(len(tr))
            for j in range(N_FOLDS):
                res[f == j] = tr.y[f == j] - fit_predict(m, tr[f != j], tr[f == j], "y")
            cv[m].append(rmse(res))
            preds[m].append(fit_predict(m, tr, te, "y"))
    te = prep(test_p, test_dxa[test_dxa.imp == imps[0]], lever)
    y = te.y.to_numpy()
    P = {m: np.mean(preds[m], axis=0) for m in MODELS}  # pooled over imputations
    boot_idx = np.random.default_rng(SEED + 1).integers(0, len(y), (N_BOOT, len(y)))
    for m in MODELS:
        r = y - P[m]
        br = np.sqrt(np.mean(r[boot_idx] ** 2, axis=1))
        out[m] = dict(train_groupcv_rmse=float(np.mean(cv[m])), test_rmse=rmse(r),
                      test_rmse_ci95=[float(np.percentile(br, 2.5)), float(np.percentile(br, 97.5))],
                      test_mean_resid_male=float(r[te.female == 0].mean()),
                      test_mean_resid_female=float(r[te.female == 1].mean()))
    rc = y - P[COMPARATOR]; rp = y - P[PRIMARY_DXA]
    imp_pt = 1 - rmse(rp) / rmse(rc)
    bimp = 1 - np.sqrt(np.mean(rp[boot_idx] ** 2, 1)) / np.sqrt(np.mean(rc[boot_idx] ** 2, 1))
    out["primary_improvement"] = dict(point=float(imp_pt),
                                      ci95=[float(np.percentile(bimp, 2.5)), float(np.percentile(bimp, 97.5))],
                                      n_test=int(len(y)), n_train=int(len(tr)), n_imputations=len(imps))
    return out


def prep(p, dxa, lever):
    d = p.merge(dxa, on="id", how="inner").dropna(subset=["leglean_kg"]).copy()
    d["ht_m"] = d.ht / 100.0
    d["female"] = (d.sex == 2).astype(float)
    d["y"] = d.force * LEVER[lever](d.ht_m.to_numpy())
    return d.sort_values("id").reset_index(drop=True)


def run(train, test):
    (trp, trd), (tep, ted) = train, test
    log = {"train": {}, "test": {}}
    trp = select(trp, trd, log["train"]); tep = select(tep, ted, log["test"])
    overlap = set(trp.id) & set(tep.id)
    assert not overlap, "leakage: SEQN overlap between cycles"
    log["seqn_overlap"] = 0
    log["results"] = {lv: analyse(trp, trd, tep, ted, lv) for lv in LEVER}
    pi = log["results"][PRIMARY_LEVER]["primary_improvement"]
    log["criterion_pass"] = bool(pi["point"] >= IMPROVE_THRESHOLD)
    return log


LEVER["unit"] = lambda h_m: np.ones_like(h_m)  # force in N, no lever-arm assumption


def synth_cycle(rng, n, id0):
    """SYNTHETIC data with an invented structure; used only to test the code path."""
    sex = rng.integers(1, 3, n); fem = sex == 2
    age = rng.integers(20, 86, n).astype(float)
    ht = np.where(fem, 161, 175) + rng.normal(0, 7, n)
    lean = np.where(fem, 13, 18) - 0.05 * (age - 50) + rng.normal(0, 2, n)
    fat = np.where(fem, 30, 22) + rng.normal(0, 8, n)
    wt = 2.9 * lean + fat
    torque_true = 5.0 * lean * ht / 100 + rng.normal(0, 18, n)
    force = torque_true / 0.30
    ids = np.arange(id0, id0 + n)
    p = pd.DataFrame(dict(id=ids, age=age, sex=sex, psu=rng.integers(1, 3, n),
                          strata=rng.integers(1, 14, n), wt=wt, ht=ht, force=force))
    p.loc[rng.random(n) < 0.15, "force"] = np.nan  # synthetic missingness
    dxa = pd.concat([pd.DataFrame(dict(id=ids, imp=k, leglean_kg=lean + rng.normal(0, .1, n)))
                     for k in range(1, 6)])
    dxa.loc[rng.random(len(dxa)) < 0.05, "leglean_kg"] = np.nan
    return p, dxa


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "synth"
    if cmd == "fetch":
        fetch(sys.argv[2]); return
    if cmd == "real":
        d = sys.argv[2]
        res = {"status": "EMPIRICAL", **run(load_cycle(d, ""), load_cycle(d, "_B"))}
        out = sys.argv[3]
    else:
        rng = np.random.default_rng(SEED)
        res = {"status": "SYNTHETIC_SELFTEST_NOT_VALIDATION",
               **run(synth_cycle(rng, 9000, 1), synth_cycle(rng, 11000, 20000))}
        res["criterion_pass"] = "NOT_APPLICABLE_SYNTHETIC"  # planted signal; not evidence
        out = sys.argv[2] if len(sys.argv) > 2 else "synthetic_selftest.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["results"][PRIMARY_LEVER]["primary_improvement"], indent=1))


if __name__ == "__main__":
    main()
