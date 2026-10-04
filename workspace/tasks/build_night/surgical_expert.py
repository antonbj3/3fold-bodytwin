#!/usr/bin/env python3
"""The surgical dense expert as a SOURCE OF INFORMATION — normalized and questionable.

Why the file exists. Anton: the dense expert should be a source of information that the system asks,
not a document anyone reads. Four Sonnet agents harvested 155 records of ~40 specialty lenses,
but they are in four directories with FOUR OLIKA SCHEMAN — the same thing is called `threshold`,
`threshold_or_behaviour` and `expected_behaviour_or_threshold`, and verification status is called
`locator_status`, `locator_verified` and `locator_check`. As long as it's true, that's four reports
and not a source.

The verification mode differs DRASTISKT between the blocks and it must never be equalized:
  D 33 records — each locator retrieved as PubMed record with abstract. Strongest.
  C 17 entries — each locator confirmed to exist; a misremembered reference corrected.
  A 36 record — generated from memory, no web access, all locators unchecked.
  B 69 records — generated from memory, 14 with low reproduction confidence, returned for resolution.
B's 69 are therefore worth less than D's 33. The quantity lies about the quality, and whoever asks the source
must see it, so each record carries `trust` and a question can be filtered on it.

Edges are formed between STORHETER, not between anatomy (dental evidence: every link that gave a
useful result matched on the limiting quantity and its unit). The specialty is the lens
one collected with; the key is the greatness. Therefore, `--quantity` is the main query and `--lens` is a filter.

Use:
  python3 tasks/build_night/surgical_expert.py --build # normalize and write source
  python3 tasks/build_night/surgical_expert.py --quantity flow # query per quantity class
  python3 tasks/build_night/surgical_expert.py --quantity force --trust high
  python3 tasks/build_night/surgical_expert.py --timescale INTRAOPERATIVE --numeric
  python3 tasks/build_night/surgical_expert.py --stats
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

SRC = Path('external_research_path')
OUT = Path(__file__).resolve().parents[2] / 'notes' / 'SURGICAL_EXPERT_SOURCE.json'

# The block's verification rate, set by how the agent actually worked and not by its own self-image.
BLOCKS = {
    'A_neuro_spine':            ('A', 'recall_unverified', 'no web access; all locators unchecked'),
    'B_cardiothoracic_vascular':('B', 'recall_unverified', 'generated from memory; 14 low confidence; returned for resolution'),
    'C_abdominal_urologic':     ('C', 'locator_exists',    'each locator confirmed to exist; speech from search results, not full text'),
    'D_musculoskeletal_surface':('D', 'pubmed_retrieved',  'each locator retrieved as PubMed record with abstract'),
}
TRUST_ORDER = {'pubmed_retrieved': 'high', 'locator_exists': 'medium', 'recall_unverified': 'low'}

# Same thing under four names. First hit wins.
ALIAS = {
    'threshold': ('threshold', 'threshold_or_behaviour', 'expected_behaviour_or_threshold'),
    'cited_sources': ('cited_sources', 'cites_sources', 'cited_sources_for_root_trace'),
    'lens': ('lens', 'specialty_lenses', 'specialty'),
    'control_without_record': ('control_without_record', 'control_probe_question'),
    'quantity_class': ('limiting_quantity_class',),
    'locator_note': ('locator_status', 'locator_check', 'locator_verified', 'recall_confidence', 'confidence'),
}
# Infer magnitude class when the agent did not set it. The order matters: time before pressure, because
# "ischemia time" otherwise gets attached to "pressure" via "perfusion pressure" in the same sentence.
CLASS_RULES = (
    ('time',               r'\b(tid|time|duration|minut|minute|hour|ischem|dwell|clamp)\b'),
    ('temperature',        r'\b(temp|grad|celsius|heat|thermal|°c)\b'),
    ('flow',              r'\b(flow|flow|perfus|cardiac output|ml/min|l/min)\b'),
    ('pressure',          r'\b(pressure|tryck|mmhg|cmh2o|kpa|insuffl)\b'),
    ('force',             r'\b(force|kraft|newton|\bn\b|tension|torque|moment)\b'),
    ('clearance_distance',r'\b(distance|distance|margin|marginal|clearance|gap|thickness|tjocklek|mm\b|µm|micron|depth)\b'),
    ('concentration',     r'\b(concentration|koncentr|mol|mg/|µg|dos)\b'),
    ('current',           r'\b(current|current|milliamp|\bma\b|volt|stimulat)\b'),
    ('proportion',        r'\b(rate|percent|%|proportion|incidence|fraction|andel|sensitivity)\b'),
)


def pick(rec: dict, key: str):
    for name in ALIAS.get(key, (key,)):
        if rec.get(name) not in (None, ''):
            return rec[name]
    return None


def derive_class(rec: dict) -> str:
    explicit = pick(rec, 'quantity_class')
    if isinstance(explicit, str) and explicit.strip():
        return explicit.strip().lower()
    hay = ' '.join(str(rec.get(k, '')) for k in ('limiting_quantity', 'threshold', 'unit',
                                                 'threshold_or_behaviour',
                                                 'expected_behaviour_or_threshold')).lower()
    for name, pat in CLASS_RULES:
        if re.search(pat, hay):
            return name
    return 'unclassified'


def has_number(rec: dict) -> bool:
    """A threshold is only usable if a number is actually in it."""
    for field in ('numeric', 'threshold_has_number', 'numeric_threshold_present_in_source'):
        v = rec.get(field)
        if isinstance(v, bool):
            return v
    return bool(re.search(r'\d', str(pick(rec, 'threshold') or '')))


def build() -> dict:
    rows = []
    for folder, (block, verification, note) in BLOCKS.items():
        f = SRC / folder / 'RECORDS.json'
        if not f.exists():
            print(f'missing, skipping: {f}', file=sys.stderr)
            continue
        raw = json.loads(f.read_text())
        recs = raw if isinstance(raw, list) else (raw.get('records') or next(iter(raw.values())))
        for r in recs:
            if not isinstance(r, dict):
                continue
            # 2/10 22:35: the trust is PER POST where the agent did a resolution pass, not per block.
            # B was returned and 53 of 69 locators were verified, but 14 entries had a hard
            # errors: ten misattributed quotes and four figures that were not in the source. To let whole
            # the block be 'low' would hide the 29 whose number is confirmed, and to lift the whole
            # block would hide the 16 that still cannot be resolved.
            rec_trust = TRUST_ORDER[verification]
            loc_ok = r.get('locator_verified') is True or str(r.get('locator_status','')).lower() in ('verified','true','confirmed')
            num = str(r.get('number_status','')).lower()
            if loc_ok and ('confirm' in num or 'correct' in num):
                rec_trust = 'high'          # the source exists AND the speech is confirmed or corrected against it
            elif loc_ok:
                rec_trust = 'medium'        # source exists, number unchecked
            elif 'not_found' in num or 'absent' in num or 'contradict' in num:
                rec_trust = 'refuted'       # the number is not in the source or is contradicted by it
            rows.append({
                'id': f"{block}:{r.get('id')}",
                'block': block,
                'trust': rec_trust,
                'verification': verification,
                'lens': pick(r, 'lens'),
                'operative_step': r.get('operative_step'),
                'geometric_condition': r.get('geometric_condition'),
                'limiting_quantity': r.get('limiting_quantity'),
                'quantity_class': derive_class(r),
                'unit': r.get('unit'),
                'threshold': pick(r, 'threshold'),
                'has_number': has_number(r),
                'timescale': r.get('timescale'),
                'resolution_level': r.get('resolution_level'),
                'evidence_class': r.get('evidence_class'),
                'locator': r.get('locator'),
                'locator_note': pick(r, 'locator_note'),
                'cited_sources': pick(r, 'cited_sources'),
                'control_without_record': pick(r, 'control_without_record'),
            })
    meta = {
        'what': 'The surgical dense expert as a source of information. Normalized out of four'
                'persona blocks with four different schemes. Asked per STORHET, not per anatomy.',
        'records': len(rows),
        'blocks': {b: {'verification': v, 'note': n} for _, (b, v, n) in BLOCKS.items()},
        'trust_counts': dict(Counter(r['trust'] for r in rows)),
        'quantity_classes': dict(Counter(r['quantity_class'] for r in rows).most_common()),
        'timescale': dict(Counter(str(r['timescale']) for r in rows).most_common()),
        'with_number': sum(1 for r in rows if r['has_number']),
        'caveat': 'trust=low means the locator is uncontrolled rendering from memory.'
                  'Never use a low record as held data; it is a clue to verify.'
                  'Allt PENDING_INDEPENDENT_REVIEW.',
    }
    return {'_meta': meta, 'records': rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--quantity')
    ap.add_argument('--lens')
    ap.add_argument('--timescale')
    ap.add_argument('--trust', choices=('high', 'medium', 'low'))
    ap.add_argument('--numeric', action='store_true', help='only entries where the threshold carries a number')
    ap.add_argument('--stats', action='store_true')
    a = ap.parse_args()

    if a.build or not OUT.exists():
        data = build()
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(data, ensure_ascii=False, indent=1))
        print(f"skrev {OUT}: {data['_meta']['records']} poster")
        print(f"  tillit: {data['_meta']['trust_counts']}")
        print(f"  storhetsklasser: {data['_meta']['quantity_classes']}")
        print(f"  with numbers in the threshold: {data['_meta']['with_number']}")
        if a.build:
            return 0

    data = json.loads(OUT.read_text())
    rows = data['records']
    if a.stats:
        print(json.dumps(data['_meta'], ensure_ascii=False, indent=1))
        return 0
    sel = rows
    if a.quantity:
        q = a.quantity.lower()
        sel = [r for r in sel if q in r['quantity_class'] or q in str(r['limiting_quantity']).lower()]
    if a.lens:
        sel = [r for r in sel if a.lens.lower() in str(r['lens']).lower()]
    if a.timescale:
        sel = [r for r in sel if a.timescale.upper() in str(r['timescale']).upper()]
    if a.trust:
        sel = [r for r in sel if r['trust'] == a.trust]
    if a.numeric:
        sel = [r for r in sel if r['has_number']]
    print(f'{len(sel)} records by {len(rows)}')
    for r in sel[:40]:
        print(f"  [{r['trust']:6s} {str(r['timescale'])[:14]:14s}] {r['id']} {str(r['lens'])[:26]:26s} "
              f"{r['quantity_class']:18s} {str(r['threshold'])[:70]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
