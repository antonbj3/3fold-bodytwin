#!/usr/bin/env python3
"""Does every number an edge cites still come out of the script that owns it?

Why this exists, and why now. Three decisions were built today and each one has an edge that quotes
its numbers in prose: T-E26 quotes the incision dwell chain, T-E27 the shaft inversion, T-E28 the
passage friction budget. Both sides were then edited several times -- the incision decision's refusal
was narrowed after the proof lane falsified its falsifier, the needle edge was corrected three times, the
meniscus edge was dissolved and rebuilt. Each edit is a chance for an edge to keep quoting a number
the code no longer produces, and a stale quoted number is worse than a missing one because it reads
as evidence.

The check is mechanical and deliberately dumb: run each owning script, collect every number it
writes into its own result JSON, and confirm that every number the edge's prose quotes appears
there. A number in the prose that the script cannot produce is reported with both values. Nothing
is repaired automatically, because the right repair differs -- sometimes the edge is stale, sometimes
the script regressed, and guessing which would hide the second case.

Numbers are compared as decimal strings at the precision the edge states, which is the endpoint
convention: an edge that says 0.137 matches 0.13712, and an edge that says 0.137125 does not match
0.1371. That is the rule PROOF_LANE_PRECISION_PROP established after a width-based rule failed all 1200
of its constructed boundary cases.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[2]
# Numbers an edge legitimately quotes from ANOTHER lane's measurement. They cannot come out of this
# edge's own script and reporting them as stale would train the gate to be ignored. Each one is
# listed with where it came from, so the list stays auditable rather than becoming a silencer.
CROSS_LANE = {
    '1273.2395447351628': 'LANE_TISSUE_ADHESION_CONTACT r1, default patch pressure',
    '0.7853981633974483': 'LANE_TISSUE_ADHESION_CONTACT r1, default patch area',
    '36.378273': 'LANE_TISSUE_ADHESION_CONTACT r1, ratio to the 35 kPa gel',
}

OWNED = {
    'T-E26-dwell_time-vs-collateral_damage_depth':
        ('tasks/assembly/incision_setting_decision.py', 'results/ASSEMBLY_INCISION_SETTING/DECISION_V1.json'),
    'T-E27-shaft_sensing-vs-tip_tangential_load':
        ('tasks/assembly/tip_load_decision.py', 'results/ASSEMBLY_TIP_LOAD/DECISION_V1.json'),
    'T-E28-pass_count-vs-forward_friction_budget':
        ('tasks/assembly/passage_friction_decision.py', 'results/ASSEMBLY_PASSAGE_FRICTION/DECISION_V1.json'),
}
# A number with at least three significant figures. Below that almost everything matches by accident,
# and counts like "two settings" or "4 N" are not claims about a computed quantity.
QUOTED = re.compile(r'(?<![\d.eE+-])(\d+\.\d{2,}(?:[eE][+-]?\d+)?|\d\.\d*[eE][+-]?\d+)')
# Years, PMIDs, DOIs and section numbers are not results.
SKIP_CONTEXT = re.compile(r'(PMID|doi|LNCS|pp\.|10\.\d{4})', re.I)


def numbers_in(obj) -> set[str]:
    """Every number the script's own output carries, as strings at full precision."""
    out: set[str] = set()
    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        elif isinstance(o, bool):
            return
        elif isinstance(o, (int, float)):
            out.add(repr(float(o)))
        elif isinstance(o, str):
            for m in QUOTED.finditer(o):
                out.add(repr(float(m.group(1))))
    walk(obj)
    return out


def matches(quoted: str, produced: set[str]) -> bool:
    """Endpoint convention: the quoted decimal must agree with some produced value to its own
    stated precision, and a quoted value with more digits than any produced value cannot match."""
    # The first version treated exponent notation as "no stated precision" and demanded agreement to
    # 1e-9, so an edge quoting 9.6333e-07 for a produced 9.633333e-07 was reported as stale when the
    # edge had simply rounded. Significant figures are the right measure for both notations.
    q = float(quoted)
    mant = quoted.lower().split('e')[0].replace('-', '').replace('.', '').lstrip('0')
    sig = max(1, len(mant.rstrip('0')) if '.' in quoted else len(mant))
    for p in produced:
        pv = float(p)
        if f'{pv:.{sig}g}' == f'{q:.{sig}g}':
            return True
    return False


def main() -> int:
    net = json.loads((W / 'CONSTRAINT_NETS.json').read_text())
    edges = {}
    def walk(o):
        if isinstance(o, dict):
            if o.get('id') in OWNED:
                edges[o['id']] = o
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(net)

    report = {'cross_lane_allowed': CROSS_LANE, 'checked': [], 'stale': [], 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    for eid, (script, out) in OWNED.items():
        e = edges.get(eid)
        if e is None:
            report['stale'].append({'edge': eid, 'problem': 'edge not found in the net'})
            continue
        r = subprocess.run([sys.executable, str(W / script)], capture_output=True, text=True, cwd=W)
        if r.returncode != 0:
            report['stale'].append({'edge': eid, 'problem': f'{script} exited {r.returncode}',
                                    'stderr': r.stderr[-400:]})
            continue
        produced = numbers_in(json.loads((W / out).read_text()))
        prose = ' '.join(str(v) for k, v in e.items()
                         if isinstance(v, str) and k not in ('id', 'evidence', 'external_reference'))
        missing = []
        for m in QUOTED.finditer(prose):
            tok = m.group(1)
            ctx = prose[max(0, m.start() - 40):m.start()]
            if SKIP_CONTEXT.search(ctx):
                continue
            if tok in CROSS_LANE:
                continue
            if not matches(tok, produced):
                missing.append(tok)
        report['checked'].append({'edge': eid, 'script': script,
                                  'numbers_produced': len(produced),
                                  'quoted_numbers_not_produced': missing})
        if missing:
            report['stale'].append({'edge': eid, 'quoted_but_not_produced': missing})

    d = W / 'results/ASSEMBLY_CONSISTENCY'
    d.mkdir(parents=True, exist_ok=True)
    (d / 'EDGE_VS_SCRIPT_V1.json').write_text(json.dumps(report, indent=1, ensure_ascii=False))
    for c in report['checked']:
        n = len(c['quoted_numbers_not_produced'])
        print(f"  {c['edge'][:46]:46} {c['numbers_produced']:4} numbers from the script, "
              f"{n} quoted numbers that the script does not produce")
        for t in c['quoted_numbers_not_produced']:
            print(f'      missing: {t}')
    print(f"  edges with problems: {len(report['stale'])} of {len(OWNED)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
