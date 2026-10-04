#!/usr/bin/env python3
"""Propose a net edge for every runnable decision that the net does not yet know about.

The measured gap, 2026-10-04. The net has 164 edges, of which 5 name a runnable decision, 4 carry an
external reference and 3 carry a falsifier. Separately, tasks/assembly holds 23 decision files and
ALL 23 run clean, spanning eye, knee, spine, needle and instrument, immune, mitochondrial, transport
and renal chains. So the scarce thing in this repo is not decisions and not edges -- it is the
junction between them. Eighteen working decisions are invisible to the net.

This proposes that junction rather than making it. For each decision it reads the result JSON the
decision writes, pulls the fields that carry a claim -- question, claim_type, control, falsifier,
external reference, and the decision's own headline number -- and emits an edge proposal with an
evidence pointer of the form the net already uses. Nothing is written into the net: proposals land in
results/ASSEMBLY_WIRING/PROPOSALS_V1.json for admission one at a time, because the 500-entry list and
the 69 templated jobs are what automatic admission produced last time.

What it deliberately does not do. It does not invent a falsifier for a decision that has none; a
missing falsifier is reported as missing, since writing one here would put the coordinator's guess
into the net as if the decision had earned it. It does not guess an external reference either: a
decision whose code carries no PMID or DOI is proposed as UNKNOWN with that stated.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path(__file__).resolve().parents[2]
OUT = W / 'results/ASSEMBLY_WIRING'
# Not decisions about the body: the two gates in this directory match *decision* by name, and the
# glass batch file belongs to cad-to-simulation, not BodyTwin.
NOT_A_BODY_DECISION = {'decision_edge_consistency.py', 'wire_decisions_to_net.py',
                       'glass_batch_energy_decision.py'}

LOCATOR = re.compile(r'(PMID\s*:?\s*\d{6,9}|10\.\d{4,9}/[^\s\'"`,)\]]+)')
# Fields a decision uses for its headline claim, in the order they are preferred.
CLAIM_KEYS = ('question', 'claim', 'headline', 'verdict', 'reading', 'problem')
CONTROL_KEYS = ('control', 'strongest_control', 'control_equal_friction_passes', 'control_coulomb')
NUM_KEYS = ('decision_value_N', 'decision_value_s', 'decision_value_passes', 'decision_value')


def first(obj, keys):
    """Depth-first search for the first of these keys carrying a non-empty value."""
    if isinstance(obj, dict):
        for k in keys:
            v = obj.get(k)
            if v not in (None, '', [], {}):
                return k, v
        for v in obj.values():
            got = first(v, keys)
            if got:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = first(v, keys)
            if got:
                return got
    return None


def main() -> int:
    net_ids = set()
    net_owners = set()
    def walk(o):
        if isinstance(o, dict):
            if o.get('id'):
                net_ids.add(o['id'])
            for k in ('owner', 'owner_lane'):
                if o.get(k):
                    net_owners.add(str(o[k]))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(json.loads((W / 'CONSTRAINT_NETS.json').read_text()))

    proposals, already, no_output = [], [], []
    for src in sorted((W / 'tasks/assembly').glob('*decision*.py')):
        if src.name in NOT_A_BODY_DECISION:
            continue
        rel = f'tasks/assembly/{src.name}'
        if rel in net_owners:
            already.append(rel)
            continue
        code = src.read_text(errors='ignore')
        # The first version looked only for OUT = W / 'results/...' and reported 13 of 23 decisions
        # as having no result, while every one of them had just run and printed the file it wrote.
        # Half of them name the directory as a plain local string instead. Any results/ASSEMBLY_*
        # literal in the file is the right thing to look for.
        cands = re.findall(r"['\"](results/ASSEMBLY_[A-Za-z0-9_]+)['\"]", code)
        if not cands:
            no_output.append({'decision': rel, 'problem': 'no results/ASSEMBLY_* path in the file'})
            continue
        sub = cands[0]
        d = W / sub
        files = sorted(d.glob('*.json')) if d.is_dir() else []
        if not files:
            no_output.append({'decision': rel, 'problem': f'{sub} holds no result JSON'})
            continue
        res = json.loads(files[0].read_text())
        claim = first(res, CLAIM_KEYS)
        ctrl = first(res, CONTROL_KEYS)
        fals = first(res, ('falsifier',))
        num = first(res, NUM_KEYS)
        locs = sorted({x.strip() for x in LOCATOR.findall(code)})
        proposals.append({
            'proposed_owner': rel,
            'evidence': f'{sub}/{files[0].name} :: '
                        + (f'{num[0]} = {num[1]}' if num else 'no single decision value emitted'),
            'claim': (claim[1] if claim else None),
            'control': (ctrl[1] if ctrl else None),
            'falsifier': (fals[1] if fals else None),
            'falsifier_missing': fals is None,
            'external_reference': (locs[0] if locs else None),
            'external_reference_status': 'FOUND_IN_CODE' if locs else 'UNKNOWN_NONE_IN_CODE',
            'all_locators_in_code': locs,
            'status_if_admitted': 'OPEN',
            'review_state': 'PENDING_INDEPENDENT_REVIEW',
        })

    report = {'net_edges_with_an_owner': len(net_owners),
              'decisions_already_wired': already,
              'decisions_without_a_result_json': no_output,
              'proposals': proposals,
              'proposals_missing_a_falsifier': sum(1 for p in proposals if p['falsifier_missing']),
              'proposals_without_any_locator': sum(1 for p in proposals
                                                   if p['external_reference'] is None),
              'admission_rule': 'one at a time, by hand, never in bulk',
              'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'PROPOSALS_V1.json').write_text(json.dumps(report, indent=1, ensure_ascii=False))
    print(f'  redan kopplade: {len(already)}   without result-JSON: {len(no_output)}   Proposal: {len(proposals)}')
    print(f"  of the proposals lack: {report['proposals_missing_a_falsifier']} falsifiers and {report['proposals_without_any_locator']} any locator at all")
    for p in proposals:
        print(f"   {'lok' if p['external_reference'] else '   '} "
              f"{'fals' if not p['falsifier_missing'] else '    '}  {p['proposed_owner'][16:]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
