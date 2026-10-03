#!/usr/bin/env python3
"""Deterministically compute and audit the BT-DM raw-data packets."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tasks"))
import datamatrix as reader

RESULTS = ROOT / "results"
TABLES = RESULTS / "DATAMATRIX_TABLES"
LANE = RESULTS / "CX-DMCOMPUTE"
META = ["dataset", "unit", "person", "activity", "packet", "raw_sha256", "raw_source", "check_value", "missing_channels", "outlier_flags"]


def labels(spec):
    unit = spec["unit"]
    bits = unit.split("_", 1)
    person = bits[0].upper()
    activity = bits[1] if len(bits) > 1 else unit
    activity = re.sub(r"\d+$", "", activity).rstrip("_-.")
    return person, activity or unit


def missing_channels(proto, names):
    names = set(names)
    if proto.startswith("P-GRF"):
        ids = {n[2:] for n in names if re.fullmatch(r"fz\d+", n)}
        return [f"{axis}{i}" for i in sorted(ids) for axis in ("fx", "fy", "fz") if f"{axis}{i}" not in names]
    if proto.startswith("P-EMG"):
        return [x for x in ("semimem", "bifem", "vasmed", "medgas", "soleus") if x not in names]
    if proto == "P-CONTACT-COMP":
        a = ("pm", "am", "al", "pl") if "pm" in names else ("fx", "fy", "fz")
        return [x for x in a if x not in names]
    if proto == "P-MARKER-GEOM":
        bases = ("R.Asis", "L.Asis", "R.Knee.Lateral", "R.Knee.Medial", "L.Knee.Lateral", "L.Knee.Medial", "R.Ankle.Lateral", "R.Ankle.Medial", "L.Ankle.Lateral", "L.Ankle.Medial")
        return [b + a for b in bases for a in "xyz" if ("ranklemedial" + a if b == "R.Ankle.Medial" and "r.ankle.medialx" not in names else b.lower() + a) not in names]
    return [] if any(b + "x" in names for b in ("sacrum", "lumbar", "r.psis")) else ["pelvic_marker"]


def measured_values(proto, names, arr):
    """Follow the packet reader, treating all-zero marker triplets as absent.

    The reader converts eTibia/eKnee pound-force to N and takes GRF peaks
    from one plate. Keep these definitions here rather than rescaling tables
    after computation, so packet and aggregate outputs agree.
    """
    values, check = reader.measure(proto, names, arr)
    if proto not in ("P-MARKER-GEOM", "P-IMU-VIRTUAL"):
        return values, check
    lookup = {name: i for i, name in enumerate(names)}

    def point(base):
        key = base.lower()
        if key == "r.ankle.medial" and "r.ankle.medialx" not in lookup:
            key = "ranklemedial"
        p = np.stack([arr[:, lookup[key + axis]] for axis in "xyz"], axis=1).copy()
        p[np.all(p == 0, axis=1)] = np.nan
        return p

    if proto == "P-MARKER-GEOM":
        pairs = [("pelvis_width_mm", "R.Asis", "L.Asis"),
                 ("right_knee_width_mm", "R.Knee.Lateral", "R.Knee.Medial"),
                 ("left_knee_width_mm", "L.Knee.Lateral", "L.Knee.Medial"),
                 ("right_ankle_width_mm", "R.Ankle.Lateral", "R.Ankle.Medial"),
                 ("left_ankle_width_mm", "L.Ankle.Lateral", "L.Ankle.Medial")]
        for field, a, b in pairs:
            delta = point(a) - point(b)
            distance = np.linalg.norm(delta, axis=1)
            values[field] = float(np.nanmedian(distance)) if np.isfinite(distance).any() else math.nan
    else:
        base = next(b for b in ("Sacrum", "Lumbar", "R.Psis") if b.lower() + "x" in lookup)
        p = point(base)
        good = np.all(np.isfinite(p), axis=1)
        t = arr[good, lookup["time(sec)"]]
        p = p[good]
        if len(t) < 20 or len(np.unique(t)) != len(t):
            return {field: math.nan for field in values}, check
        v = np.gradient(p, t, axis=0) / 1000
        a = np.gradient(v, t, axis=0)
        speed = np.linalg.norm(v, axis=1)
        accel = np.linalg.norm(a, axis=1)
        values.update(sacrum_accel_rms_m_s2=float(np.sqrt(np.mean(accel * accel))),
                      sacrum_accel_peak_m_s2=float(np.max(accel)),
                      sacrum_speed_mean_m_s=float(np.mean(speed)),
                      sacrum_speed_peak_m_s=float(np.max(speed)),
                      sacrum_vertical_range_mm=float(np.ptp(p[:, 2])))
    return values, check


def compute():
    LANE.mkdir(exist_ok=True)
    TABLES.mkdir(exist_ok=True)
    packets = sorted(RESULTS.glob("BT-DM-*/inputs/PROTOCOL.json"))
    if len(packets) != 605:
        raise RuntimeError(f"Expected 605 packets, found {len(packets)}")
    if any("jw_lungef1" in p.as_posix().lower() for p in packets):
        raise RuntimeError("F-8 sealed force packet present")
    cache = {}
    by_protocol = defaultdict(list)
    stats = Counter()
    for i, path in enumerate(packets, 1):
        packet = path.parent.parent
        spec = json.loads(path.read_text())
        if "jw_lungef1" in spec["unit"].lower() or "jw_lungef1" in spec["raw_source"].lower():
            raise RuntimeError("F-8 sealed force source encountered")
        raw = path.parent / spec["raw_file"]
        data = raw.read_bytes()
        sha = hashlib.sha256(data).hexdigest()
        if sha != spec["raw_sha256"]:
            raise RuntimeError(f"SHA256 mismatch: {packet}")
        if sha not in cache:
            cache[sha] = reader.raw_table(data)
        names, arr = cache[sha]
        proto = spec["protocol"]
        missing = missing_channels(proto, names)
        if missing:
            values = {field: math.nan for field in spec["fields"]}
            check = math.nan
            stats["missing_packet"] += 1
        else:
            values, check = measured_values(proto, names, arr)
            if set(values) != set(spec["fields"]):
                raise RuntimeError(f"Field mismatch: {packet}")
            if not math.isclose(check, spec["check_value"], abs_tol=1e-10, rel_tol=1e-8):
                raise RuntimeError(f"Odd/even mismatch: {packet}")
        missing.extend(f"{field}:no_valid_samples" for field, value in values.items() if not math.isfinite(value))
        person, activity = labels(spec)
        row = dict(dataset=spec["dataset"], unit=spec["unit"], person=person,
                   activity=activity, packet=packet.name, raw_sha256=sha,
                   raw_source=spec["raw_source"], check_value=check,
                   missing_channels=";".join(missing), outlier_flags="", **values)
        by_protocol[proto].append(row)
        stats[proto] += 1
        if i % 100 == 0:
            print(f"computed {i}/{len(packets)}", flush=True)
    # Mark anomalies; never delete or modify the measurement itself.
    flags = []
    for proto, rows in sorted(by_protocol.items()):
        fields = reader.PROTOCOLS[proto][1]
        for field in fields:
            v = np.array([r[field] for r in rows], float)
            finite = v[np.isfinite(v)]
            med = np.median(finite) if len(finite) >= 10 else math.nan
            mad = np.median(np.abs(finite - med)) if len(finite) >= 10 else math.nan
            scale = 1.4826 * mad
            for row in rows:
                x = row[field]
                reasons = []
                if math.isfinite(x):
                    if scale > 0 and abs(x - med) > 6 * scale:
                        reasons.append("robust_z_gt_6")
                    if field.endswith("_fraction") and not 0 <= x <= 1:
                        reasons.append("fraction_range")
                    if (field.endswith("_rms") or field.endswith("_peak_N") or field.endswith("_duration_s") or field.endswith("_width_mm")) and x < 0:
                        reasons.append("negative_magnitude")
                    if field.endswith("_width_mm") and x > 500:
                        reasons.append("width_gt_500_mm")
                for reason in reasons:
                    tag = f"{field}:{reason}"
                    row["outlier_flags"] += (";" if row["outlier_flags"] else "") + tag
                    flags.append({"packet": row["packet"], "protocol": proto, "field": field,
                                  "value": float(x), "reason": reason})
        with (TABLES / (proto + ".csv")).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=META + fields)
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: (r["dataset"], r["unit"])))
    # Packet files share exactly the same audited row as the protocol table.
    for rows in by_protocol.values():
        for row in rows:
            packet = RESULTS / row["packet"]
            spec = json.loads((packet / "inputs/PROTOCOL.json").read_text())
            with (packet / "TABLE.csv").open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=META + spec["fields"])
                writer.writeheader()
                writer.writerow(row)
            (packet / "RESULTS.md").write_text(
                f"# {packet.name}\n\ncomputed deterministically by CX-DMCOMPUTE. "
                f"See `TABLE.csv` and `../DATAMATRIX_TABLES/{spec['protocol']}.csv`.\n\n"
                f"Raw SHA-256: `{row['raw_sha256']}`. Odd/even check: {row['check_value']:.10g}. "
                f"Missing channels: {row['missing_channels'] or 'none'}. "
                f"Outlier flags: {row['outlier_flags'] or 'none'}.\n"
            )
    (LANE / "OUTLIERS.json").write_text(json.dumps(flags, indent=2) + "\n")
    return by_protocol, stats, flags


def checks(by_protocol):
    summary = {}
    for proto, rows in sorted(by_protocol.items()):
        v = np.array([r["check_value"] for r in rows], float)
        summary[proto] = {"n": len(rows), "odd_even_median": float(np.nanmedian(v)),
                          "odd_even_p95": float(np.nanpercentile(v, 95)),
                          "odd_even_max": float(np.nanmax(v)),
                          "missing": sum(bool(r["missing_channels"]) for r in rows),
                          "flagged_rows": sum(bool(r["outlier_flags"]) for r in rows)}
    repeat = []
    for proto, rows in sorted(by_protocol.items()):
        field = reader.PROTOCOLS[proto][1][0]
        groups = defaultdict(list)
        for r in rows:
            if math.isfinite(r[field]):
                groups[(r["dataset"], r["person"], r["activity"])].append(r[field])
        cvs = [float(np.std(x, ddof=1) / abs(np.mean(x))) for x in groups.values()
               if len(x) >= 2 and abs(np.mean(x)) > 1e-12]
        repeat.append({"protocol": proto, "field": field, "groups": len(cvs),
                       "median_cv": float(np.median(cvs)) if cvs else None})
    bilateral = []
    for r in by_protocol["P-MARKER-GEOM"]:
        for part in ("knee", "ankle"):
            l, right = r[f"left_{part}_width_mm"], r[f"right_{part}_width_mm"]
            if math.isfinite(l) and math.isfinite(right) and 0 < l <= 500 and 0 < right <= 500:
                bilateral.append(abs(l - right) / ((l + right) / 2))
    peak = {(r["dataset"], r["unit"]): r["peak_resultant_N"] for r in by_protocol["P-GRF-PEAK"]}
    paired = [(peak[(r["dataset"], r["unit"])], r["contact_peak_N"])
              for r in by_protocol["P-CONTACT-COMP"] if (r["dataset"], r["unit"]) in peak]
    corr = float(np.corrcoef(np.array(paired).T)[0, 1]) if len(paired) >= 3 else None
    out = {"protocols": summary, "repeatability": repeat,
           "bilateral": {"n": len(bilateral), "median_relative_difference": float(np.median(bilateral)) if bilateral else None},
           "contact_vs_grf": {"n": len(paired), "pearson_peak": corr}}
    (LANE / "CHECKS.json").write_text(json.dumps(out, indent=2) + "\n")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--compute", action="store_true")
    args = p.parse_args()
    if args.compute:
        rows, _, _ = compute()
        print(json.dumps(checks(rows), indent=2))


if __name__ == "__main__":
    main()
