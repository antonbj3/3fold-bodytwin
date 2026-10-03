#!/usr/bin/env python3
"Contrasts (predictions) per scenario from pred_rows.json -> predictions.json.\nDefinitions (PREREG §Predictions). MBL in mm from insertion; '_fromload' = from loading (t_load 0,25 years)."
import json, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "pred_rows.json")))
# biology-only bone (no mechanics): MBL_bio(t) = L_bio(1 − e^(−t/τ_bio)); added as separate scenarios "<scen>|bio"
_extra = []
for r in rows:
    tb = 0.05 if r["scen"] == "tau_bio_0.05" else 0.3 if r["scen"] == "tau_bio_0.3" else 0.15
    q = dict(r, scen=r["scen"] + "|bio", mech_extra=0.0, failed=False)
    for t in ("05", "1", "5"):
        q["MBL_" + t] = r["L_bio"] * (1 - math.exp(-{"05": 0.5, "1": 1.0, "5": 5.0}[t] / tb))
    _extra.append(q)
rows = rows + _extra
THIN, THICK = [1.5, 2.0], [2.5, 3.0, 3.5, 4.0]


def sel(sc, geo, fo=1.0, tms=None):
    return [r for r in rows if r["scen"] == sc and r["geo"] == geo and r["f_over"] == fo and (tms is None or r["t_muc"] in tms)]


def mean(rs, k):
    return float(np.mean([r[k] for r in rs])) if rs else None


def fromload(rs, k, tau_bio=0.15, t_load=0.25):
    return float(np.mean([r[k] - r["L_bio"] * (1 - math.exp(-t_load / tau_bio)) for r in rs])) if rs else None


out = {}
for sc in sorted(set(r["scen"] for r in rows)):
    base = sc.split("|")[0]
    tb = 0.05 if base == "tau_bio_0.05" else 0.3 if base == "tau_bio_0.3" else 0.15
    P0 = sel(sc, "G0_matched")
    d = {}
    for t in ("1", "5"):
        k = "MBL_" + t
        d["P0_MBL_%sy" % t] = mean(P0, k)
        d["P0_MBL_%sy_fromload" % t] = fromload(P0, k, tb)
        d["P0_PSmix_MBL_%sy" % t] = mean(P0 + sel(sc, "G1_PS035"), k)
        pop = P0 + sel(sc, "G1_PS035")
        thin = [r for r in pop if r["t_muc"] in THIN]; thick = [r for r in pop if r["t_muc"] in THICK]
        d["Ca_thin_minus_thick_%sy" % t] = mean(thin, k) - mean(thick, k) if thin and thick else None
        g1 = sel(sc, "G1_PS035")
        d["Cb_PS_minus_matched_%sy" % t] = mean(g1, k) - mean(P0, k) if g1 else None
        g2, g3 = sel(sc, "G2_mach05"), sel(sc, "G3_mach10")
        mach = g2 + g3
        d["Cc_rough_minus_machined_%sy" % t] = mean(P0, k) - mean(mach, k) if mach else None
        f2 = sel(sc, "G0_matched", 2.0)
        d["Cd_load2x_minus_1x_%sy" % t] = mean(f2, k) - mean(P0, k) if f2 else None
        d["Cd_bruxism_%sy" % t] = 0.0   # location-specific reference: leg is adapted to individual load
        g4 = sel(sc, "G4_PS035_sub1")
        d["Ce_sub1_minus_equi_PS_%sy" % t] = mean(g4, k) - mean(g1, k) if g4 and g1 else None
    d["P0_increment_1to5"] = d["P0_MBL_5y"] - d["P0_MBL_1y"]
    d["failed_any_f1"] = any(r["failed"] for r in rows if r["scen"] == sc and r["f_over"] == 1.0)
    d["mech_extra_P0_mean"] = mean(P0, "mech_extra")
    d["L_bio_P0_mean"] = mean(P0, "L_bio")
    out[sc] = d
# envelope over sensitivity scenarios (excl. abs_* which are separate model variants)
keys = [k for k in out["primary"] if isinstance(out["primary"][k], float)]
env = {}
for k in keys:
    v = [out[s][k] for s in out if not s.startswith("abs_") and "|" not in s and out[s].get(k) is not None]
    env[k] = [min(v), max(v)] if v else None
json.dump({"by_scenario": out, "envelope_rel": env, "n_scen": len(out)}, open(os.path.join(HERE, "predictions.json"), "w"), indent=1)
p = out["primary"]
for k in keys:
    print("%-34s %8.3f  bio %8.3f  [%s]" % (k, p[k], out["primary|bio"][k], ", ".join("%.2f" % v for v in env[k]) if env[k] else ""))
