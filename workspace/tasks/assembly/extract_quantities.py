"""Make a week of reports addressable: pull (quantity, value, unit) out of whatever shape they are in.

KNOWN DEFECT, measured on this index own output: the suffix rule reads a COMPOUND unit as its last
component. `mu_pa_s` is Pascal-seconds and was labelled seconds; `v_mol_s` is mol per second and was
labelled seconds; `darcy_flux_m_s` is metres per second, same. 7511 of 1165253 rows are certainly
wrong this way, 0.64 %, and that is a LOWER bound since it counts only the nine patterns I thought to
check. The errors concentrate in compound units, which is to say velocities, fluxes and viscosities --
exactly the quantities tissue transport is written in. Any dimensional filter built on this index is
weakest precisely where it would matter most.

The measurement that forced this, including my own false step. A mechanical test said 75 % of 13127
reports carried neither a value with a unit nor an external locator, which reads as a week of work with
doubtful content. The test was too narrow and I nearly reported it: a sample of the supposedly empty
group had 11 of 12 carrying more than three numbers and 12 of 12 carrying unit words in their text. The
numbers are there. Only 2.9 % of reports happen to put them in a `"value"` plus `"unit"` field pair,
which is the only shape anything here can query.

So the problem is not value, it is addressability. Thirteen thousand reports contain quantities that
nothing can ask for, which is why a quantity we needed could sit in the corpus while a job went out to
search the literature for it.

This reads any JSON shape and recovers triples: a numeric leaf, a unit taken from the key or a
neighbouring field, and the key path as the quantity name. It is deliberately conservative -- a number
with no recoverable unit is dropped rather than guessed, because a unitless number entering a chain is
the failure that produced a 3.33x error against a held-out measurement earlier in this project.
"""
from __future__ import annotations

import glob
import json
import re
from collections import Counter
from pathlib import Path

W = Path('')
OUT = W / 'results/ASSEMBLY_QUANTITY_INDEX'

# Unit suffixes as they appear in key names across this corpus, longest first so _mm_s beats _s.
# Compound units FIRST and matched longest-first, because the defect this fixes was reading only the
# last component: `mu_pa_s` is Pascal-seconds and came out as seconds, `v_mol_s` is mol per second,
# `darcy_flux_m_s` is metres per second. 7511 rows were wrong that way and they concentrated in
# velocities, fluxes and viscosities -- the quantities tissue transport is written in, so the
# dimensional filter built on this index was weakest exactly where it mattered most.
UNIT_SUFFIX = [
    # compound, longest first
    ('_w_m_k', 'W/(m*K)'), ('_j_kg_k', 'J/(kg*K)'), ('_mol_m2_s', 'mol/(m^2*s)'),
    ('_mol_m3_s', 'mol/(m^3*s)'), ('_percentage_points', 'percentage points'),
    ('_pmol_per_l', 'pmol/L'), ('_m_per_pa_s', 'm/(Pa*s)'), ('_pa_s', 'Pa*s'),
    ('_mol_s', 'mol/s'), ('_mol_m3', 'mol/m^3'), ('_kg_m3', 'kg/m^3'), ('_j_mol', 'J/mol'),
    ('_j_per_mol', 'J/mol'), ('_j_kg', 'J/kg'), ('_n_m', 'N*m'), ('_n_mm', 'N/mm'),
    ('_ml_min', 'mL/min'), ('_l_min', 'L/min'), ('_um_h', 'um/h'), ('_mm_h', 'mm/h'),
    ('_um_s', 'um/s'), ('_mm_s', 'mm/s'), ('_m_s', 'm/s'), ('_cm_s', 'cm/s'),
    ('_pa_1', '1/Pa'), ('_per_s', '1/s'), ('_per_h', '1/h'), ('_cm4', 'cm^4'),
    ('_mm2', 'mm^2'), ('_m2', 'm^2'), ('_m3', 'm^3'), ('_um2', 'um^2'),
    # single
    ('_mpa', 'MPa'), ('_kpa', 'kPa'), ('_gpa', 'GPa'), ('_pa', 'Pa'),
    ('_nm', 'nm'), ('_um', 'um'), ('_mm', 'mm'), ('_cm', 'cm'),
    ('_degc', 'degC'), ('_deg', 'deg'), ('_k', 'K'), ('_d', 'D'),
    ('_nn', 'nM'), ('_mol', 'mol'), ('_j', 'J'), ('_n', 'N'),
    ('_hz', 'Hz'), ('_s', 's'), ('_h', 'h'),
    ('_percent', 'percent'), ('_pct', 'percent'), ('_fraction', '1'), ('_ratio', '1'),
]
# Longest suffix wins regardless of list order, so a later single-unit entry cannot shadow a compound.
UNIT_SUFFIX.sort(key=lambda kv: -len(kv[0]))

DOI = re.compile(r'10\.\d{4,9}/[^\s",\)\]]+')
PMID = re.compile(r'pmid["\s:]*(\d{6,9})', re.I)
SKIP_KEY = re.compile(r'seed|version|epoch|elapsed|timestamp|index$|_id$|sha256|count$|^n_|_n$', re.I)


def unit_from_key(key: str) -> str | None:
    k = key.lower()
    for suf, unit in UNIT_SUFFIX:
        if k.endswith(suf):
            return unit
    return None


def walk(obj, prefix='', sibling_unit=None):
    """Yield (quantity_path, value, unit). A sibling `unit` field covers its whole object."""
    if isinstance(obj, dict):
        here = obj.get('unit') if isinstance(obj.get('unit'), str) else sibling_unit
        for k, v in obj.items():
            yield from walk(v, f'{prefix}.{k}' if prefix else k, here)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, f'{prefix}[{i}]', sibling_unit)
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        leaf = prefix.split('.')[-1]
        if SKIP_KEY.search(leaf):
            return
        unit = unit_from_key(leaf) or sibling_unit
        if unit:
            yield (prefix, float(obj), unit)


def main() -> int:
    rows, per_report, stats = [], Counter(), Counter()
    for p in glob.glob(str(W / 'results/*/results.json')):
        job = Path(p).parent.name
        try:
            raw = Path(p).read_text(errors='replace')
            d = json.loads(raw)
        except Exception:
            stats['unparsable'] += 1
            continue
        stats['reports'] += 1
        loc = DOI.search(raw)
        pmid = PMID.search(raw)
        got = 0
        for q, v, u in walk(d):
            rows.append({'job': job, 'quantity': q, 'value': v, 'unit': u,
                         'locator': (loc.group(0) if loc else
                                     (f'PMID:{pmid.group(1)}' if pmid else None))})
            got += 1
        per_report[job] = got
        stats['reports_with_at_least_one' if got else 'reports_with_none'] += 1
        stats['triples'] += got

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'QUANTITY_INDEX_V1.json').write_text(json.dumps(
        {'stats': dict(stats),
         'units_seen': dict(Counter(r['unit'] for r in rows).most_common(25)),
         'with_locator': sum(1 for r in rows if r['locator']),
         'method': ('unit taken from the key suffix or a sibling unit field; a number with no '
                    'recoverable unit is dropped rather than guessed'),
         'review_state': 'PENDING_INDEPENDENT_REVIEW',
         'rows': rows}, ensure_ascii=False))

    print(f"  reports read: {stats['reports']}  (unreadable {stats['unparsable']})")
    print(f"  med minst en storhet: {stats['reports_with_at_least_one']}  without: {stats['reports_with_none']}")
    print(f"  storheter utvunna: {stats['triples']}  "
          f"varav med DOI/PMID: {sum(1 for r in rows if r['locator'])}")
    print('  vanligaste enheter: ' + ', '.join(
        f'{u} {n}' for u, n in Counter(r['unit'] for r in rows).most_common(8)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
