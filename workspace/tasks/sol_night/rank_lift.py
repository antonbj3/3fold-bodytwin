#!/usr/bin/env python3
"Collect the observables the swarm says would LYFTA en rangdeficiens — projektets dyrbaraste lista.\n\nWhy the file exists. Measured 2026-10-02 over 5 984 swarm reports: 794 of them (13,3 %) make a\nexplicit claim of STRUKTURELL unidentifiability, while only 49 (0,8 %) complains of\nparameter uncertainty — a ratio of 16 to 1. The most common individual form is rank deficiency, 261\nreports. This means that the binding constraint of the project over subsystems is not calculation or\ncoverage without OBSERVABILITY: every seventh job says that the quantity it was asked cannot count\ndetermined from the observables it received.\n\nA rank deficiency is not a lack of data. It is a property of the observation set, and it\nsays which DIRECTION is unobserved. Therefore, a rank deficient job carries something of an ordinary\nnegative result does not carry: which SINGLE added observable would restore the rank. 157 off\nthe reports name it, and no one has collected them.\n\nAnd it reinterprets seed 1. The value of an external biology model is not that it is another model.\nIt's that it supplies an observable we don't have — so exactly what a rank deficiency needs.\nAn external model that only reproduces a quantity we already calculate does not raise any rank.\n\nStrict matching on purpose. A first loose count with the word \"structurally\" gave 38,1 % and\noverestimated by 2,9× against the strict 13,3 %. Same error class as when the substring \"gate\" happened\nmatch 309 of 400 short texts against 49 with word limits. Loose patterns do not count here.\n\nUse:\n  python3 tasks/build_night/rank_lift.py            # skriver notes/RANK_LIFT_OBSERVABLES.json + .md\n  python3 tasks/build_night/rank_lift.py --min-only  # bara de som anger ett ANTAL extra observabler\n"
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'results' / 'SWARM48_GRAPH_REPORTS'
OUT_JSON = ROOT / 'notes' / 'RANK_LIFT_OBSERVABLES.json'
OUT_MD = ROOT / 'notes' / 'RANK_LIFT_OBSERVABLES.md'

RANK_DEF = re.compile(r'\brank[- ]deficien\w+', re.I)
# How many extra observables are required, when the job says it outright.
MIN_EXTRA = re.compile(r'minimum\s+extra\s+observables?\s+(?:required|needed)\s*[:=]?\s*\**\s*(\d+)', re.I)
# A sentence like NAMNGER the lifting observable. Requires both a lift verb and a measure word, so that
# "Nothing here is a measurement" does not count as a suggestion.
LIFT_VERB = re.compile(r'\b(?:lift|lifts|lifting|restore[sd]?|break|breaks|breaking|resolve[sd]?|'
                       r'remove[sd]?|eliminat\w+)\b', re.I)
MEAS_WORD = re.compile(r'\b(?:observable|observables|measurement|readout|assay|pair|channel|'
                       r'profile|A-V|arterial-venous)\b', re.I)
# Negations that make a sentence a waiver rather than a proposition.
NEGATED = re.compile(r'\b(?:cannot|can not|does not|do not|never|no\s+\w+\s+(?:can|could)|'
                     r'nothing here|not a measurement)\b', re.I)
FAMILY = re.compile(r'\b(Q0\d\d|Q1\d\d|SURG_[A-Z]+|IMMUNITY|BIORESP|MITOSTRESS|SOLBENCH)\b')


def sentences(text: str) -> list[str]:
    return [s.strip().replace('\n', ' ') for s in re.split(r'(?<=[.!?])\s+', text)]


def harvest() -> dict:
    rows = []
    stats = Counter()
    for f in sorted(SRC.glob('*/*/RESULTS.md')):
        try:
            t = f.read_text(errors='replace')
        except OSError:
            continue
        stats['reports'] += 1
        if not RANK_DEF.search(t):
            continue
        stats['rank_deficient'] += 1
        job = f.parts[-3]

        m = MIN_EXTRA.search(t)
        n_extra = int(m.group(1)) if m else None

        proposals, refusals = [], []
        for s in sentences(t):
            if len(s) > 420 or not MEAS_WORD.search(s):
                continue
            if not LIFT_VERB.search(s):
                continue
            (refusals if NEGATED.search(s) else proposals).append(s[:300])

        if not proposals and n_extra is None and not refusals:
            stats['rank_deficient_silent'] += 1
            continue
        if proposals:
            stats['names_a_lift'] += 1
        if refusals and not proposals:
            stats['only_excludes'] += 1
        fams = sorted(set(FAMILY.findall(t)))
        rows.append({
            'job': job,
            'families': fams[:6],
            'minimum_extra_observables': n_extra,
            'proposed_lifting_observables': proposals[:3],
            'excluded_as_insufficient': refusals[:2],
            'path': str(f.relative_to(ROOT)),
        })
    # Rank: a job specifying ONE count and naming the observable is most actionable.
    rows.sort(key=lambda r: (r['minimum_extra_observables'] is None,
                             r['minimum_extra_observables'] or 99,
                             -len(r['proposed_lifting_observables'])))
    return {'_meta': dict(stats) | {
        'what': "Observables that the swarm's rank deficient jobs say would raise the rank. A rank deficiency is a property of the observation set, not a lack of data, so it tells which direction is unobserved.",
        'caveat': "The jobs are run by a weaker model. Each row is a CANDIDATE that must be verified against the report's own source before it is ordered. All PENDING_INDEPENDENT_REVIEW.",
        'strict_matching': "Loose patterns are not counted: the word \\\"structurally\\\" alone yielded 38,1% against the strict measurement's 13,3%.",
    }, 'rows': rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-only', action='store_true',
                    help='bara rader som anger ett antal extra observabler')
    a = ap.parse_args()
    data = harvest()
    rows = data['rows']
    if a.min_only:
        rows = [r for r in rows if r['minimum_extra_observables'] is not None]
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(data, ensure_ascii=False, indent=1))

    md = ["# Observables that would lift a rank deficiency — harvested from the swarm",
          '',
          "Generated by `tasks/build_night/rank_lift.py`. A rank deficiency tells which direction",
          "is unobserved, so a rankdeficient job carries something a regular negative does not carry: which",
          "SINGLE added observable that resets the rank.",
          '',
          f"Reports searched: {data['_meta'].get('reports', 0)}. Rangdeficienta: {data['_meta'].get('rank_deficient', 0)}. Namnger ett lyft: {data['_meta'].get('names_a_lift', 0)}Only exclude inadequate roads: {data['_meta'].get('only_excludes', 0)}.",
          '',
          "Each line is a CANDIDATE from a weaker model's run and should be verified against",
          "the report's own source before ordering.",
          '']
    for r in rows[:60]:
        n = r['minimum_extra_observables']
        md.append(f"## {r['job']}" + (f" — minsta antal extra observabler: {n}" if n is not None else ''))
        if r['families']:
            md.append(f"familjer: {', '.join(r['families'])}")
        for p in r['proposed_lifting_observables']:
            md.append(f"- LYFTER: {p}")
        for x in r['excluded_as_insufficient']:
            md.append(f'- UTESLUTET as insufficient: {x}')
        md.append(f"Source: '{r['path']}`")
        md.append('')
    OUT_MD.write_text('\n'.join(md))
    print(json.dumps(data['_meta'], ensure_ascii=False, indent=1))
    print(f"skrev {OUT_JSON.name} ({len(data['rows'])} rader) and {OUT_MD.name}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
