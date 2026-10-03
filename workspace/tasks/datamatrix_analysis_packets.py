#!/usr/bin/env python3
"""Package complete cross-trial tables for coordinator review, without queueing."""
from __future__ import annotations
import csv
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
LANE = OUT / "CX-DMCOMPUTE"
TABLES = OUT / "DATAMATRIX_TABLES"

SPECS = [
    ("BT-DA-001", "A", ["P-CONTACT-COMP", "P-GRF-PEAK"],
     "Fit k(activity) in N1g to measured implant peak from all paired activities. Which activities transfer across people? Exclude jw_lungef1/F-8.",
     "N1g/B24; CX-KNEEMERGE", "Leave out one person and report held-out MAE versus a single global k and |GRF| baseline."),
    ("BT-DA-002", "A", ["P-CONTACT-COMP", "P-GRF-PEAK", "P-GRF-IMPULSE", "P-GRF-BALANCE"],
     "Which GRF peak, impulse, and balance features explain implant peak residual beyond resultant GRF peak?",
     "N1g; Field U384", "Leave out one person; compare held-out MAE with |GRF| only and permuted feature columns."),
    ("BT-DA-003", "B", ["P-CONTACT-COMP", "P-GRF-PEAK", "P-EMG-TIMING", "P-EMG-AMPLITUDE"],
     "Does cross-muscle EMG activation add information on implant peak beyond GRF across matched trials? Raw EMG amplitude is not MVC normalized.",
     "A324 H3/A325; CX-KNEEMERGE", "Leave out one person; compare GRF-only, GRF+EMG, and within-person permuted EMG."),
    ("BT-DA-004", "B", ["P-CONTACT-COMP", "P-GRF-PEAK", "P-IMU-VIRTUAL"],
     "Can pelvis-derived motion summaries estimate measured implant peak across activities? Quantify the gap to GRF.",
     "A321 sparse sensors; N1g", "Leave out an activity and compare IMU-only with constant and GRF-only baselines; permute activity labels."),
    ("BT-DA-005", "B", ["P-MARKER-GEOM"],
     "Which person-level marker widths have enough valid coverage for a geometry prior? Assess pelvis separately from sparse knee and ankle widths; flag implausible static-trial values.",
     "N7c geometry manager; anatomical what-if priors", "Hold out an activity per person where coverage allows; compare predicted width to measured holdout and a pooled median. Report unavailable fields as UNKNOWN."),
    ("BT-DA-006", "A", ["P-CONTACT-COMP", "P-GRF-PEAK", "P-EMG-AMPLITUDE", "P-IMU-VIRTUAL", "P-MARKER-GEOM"],
     "Estimate trial repeatability as a measurement-noise floor by activity and person for each modality. Separate between-trial and between-person variation.",
     "N1g error budget; A321 sensor budget", "Bootstrap whole repeated trials; use a held-out repetition and compare to the within-group median."),
    ("BT-DA-007", "C", ["P-CONTACT-COMP", "P-GRF-BALANCE", "P-GRF-PEAK"],
     "Do active force-plate fraction and directional RMS improve implant peak estimates beyond GRF peak across activities?",
     "A321 sparse sensors; N1g", "Leave out one activity; compare peak-only and peak-plus-balance models with balance columns permuted within person."),
    ("BT-DA-008", "C", ["P-CONTACT-COMP", "P-GRF-PEAK", "P-GRF-IMPULSE"],
     "Does GRF-to-implant force calibration transfer between GC datasets after accounting for activity and individual?",
     "N1g external transfer; B24 null models", "Leave out one dataset; compare a pooled model with global-ratio and activity-permuted controls."),
]


def main():
    proposal = []
    manifest = []
    for ident, priority, protocols, question, consumers, check in SPECS:
        folder = OUT / ident
        inputs = folder / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        trial_rows = {}
        files = []
        for proto in protocols:
            src = TABLES / (proto + ".csv")
            if not src.exists():
                raise FileNotFoundError(src)
            dest = inputs / src.name
            shutil.copyfile(src, dest)
            files.append({"file": dest.name, "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()})
            with dest.open(newline="") as f:
                for row in csv.DictReader(f):
                    if "jw_lungef1" in row["unit"].lower():
                        raise RuntimeError("F-8 source in analysis packet")
                    key = (row["dataset"], row["unit"])
                    trial_rows.setdefault(key, {"dataset": row["dataset"], "unit": row["unit"],
                                                "person": row["person"], "activity": row["activity"]})[proto] = "1"
        with (inputs / "TRIAL_LABELS.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["dataset", "unit", "person", "activity"] + protocols)
            writer.writeheader()
            writer.writerows(sorted(trial_rows.values(), key=lambda r: (r["dataset"], r["unit"])))
        files.append({"file": "TRIAL_LABELS.csv", "sha256": hashlib.sha256((inputs / "TRIAL_LABELS.csv").read_bytes()).hexdigest()})
        complete = sum(all(proto in r for proto in protocols) for r in trial_rows.values())
        people = len({(r["dataset"], r["person"]) for r in trial_rows.values() if all(proto in r for proto in protocols)})
        brief = (f"# {ident}\n\n## Question\n\n{question}\n\n## Inputs\n\n"
                 f"Complete protocol tables: {', '.join(protocols)}. `inputs/TRIAL_LABELS.csv` gives all trial, activity, and person labels; "
                 f"a 1 in a protocol column means that trial is present. {complete} trials have every requested table; "
                 f"{people} dataset-person groups occur in that intersection. See `inputs/MANIFEST.json` for hashes.\n\n"
                 f"## Consumers\n\n{consumers}.\n\n## Independent check\n\n{check}\n\n"
                 "Report held-out sample counts and uncertainty. Preserve missing values and outlier flags; do not delete flagged rows silently. "
                 "Treat all model results as unvalidated until the stated check is passed.\n")
        (folder / "BRIEF.md").write_text(brief)
        (inputs / "MANIFEST.json").write_text(json.dumps(files, indent=2) + "\n")
        proposal.append(f"{priority} swarm {ident} | {question} | {consumers} | {check}")
        manifest.append({"packet": ident, "tables": protocols, "complete_trials": complete, "dataset_person_groups": people})
    (LANE / "QUEUE_PROPOSAL.txt").write_text("\n".join(proposal) + "\n")
    (LANE / "ANALYSIS_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
