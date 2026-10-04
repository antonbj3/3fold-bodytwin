"""Find numbers that are the SAME measurement in a different frame, and are being counted as two.

THE OBSERVATION THIS IS BUILT ON. Six collisions in this project on 2026-10-03/04 had one shape:
two numbers that are the same measurement expressed in different frames, readable as comparable
because the field name omits the frame. Each difference was an exact multiplicative factor:

  lens plane -> spectacle plane        x 0.6823 measured per eye, mean sensitivity 0.6916
  cylinder magnitude -> power vector   x 1.775 on the same 69 eyes
  stand-in axis 90 deg -> axis 0 deg   J0 changes sign, the arm's MAE 0.6209 -> 0.3636
  stated load -> area x peak pressure  x 0.132, 0.168, 0.1155, 0.133 on four edges
  SD -> SEM                            x sqrt(n)
  micromolar -> molecules              x 0.1037837 in the declared volume

THE DANGER IS NOT THE COLLISION. It is that a frame-shifted copy reads as an INDEPENDENT source. A
swarm report tonight claimed a value "confirmed by three independent sources for the same input";
if two of those three are the same measurement in two frames, the agreement is manufactured and the
confidence is false. That failure is invisible to a unit check, because both copies carry the right
unit.

WHAT THIS DOES. For every pair of numbers of the same dimension anywhere in the net, compute the
ratio and test it against the frame factors above plus plain decade and time-base slips. A hit means
the two are probably one measurement, and the pair is reported with the frame that would explain it.

WHAT A HIT IS NOT. It is not proof. Two genuinely different quantities can stand in a ratio of 1000
by accident, and the scan says so by reporting how many hits a factor produces: a factor that
explains everything explains nothing. The test of the instrument is whether it REDISCOVERS the six
collisions that were found by hand. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import itertools
import json
import math
import pathlib
import re

W = pathlib.Path('')
NET = W / 'CONSTRAINT_NETS.json'

UNITS = {
    'D': 'dioptre', 'N': 'force', 'mN': 'force', 'uN': 'force', 'nN': 'force',
    'Pa': 'pressure', 'kPa': 'pressure', 'MPa': 'pressure', 'mmHg': 'pressure',
    'mm': 'length', 'cm': 'length', 'um': 'length', 'µm': 'length', 'nm': 'length', 'pm': 'length',
    'mm2': 'area', 'mm²': 'area', 'um2': 'area', 'µm²': 'area',
    'uM': 'concentration', 'µM': 'concentration', 'mM': 'concentration', 'nM': 'concentration',
    's': 'time', 'ms': 'time', 'min': 'time', 'h': 'time', '%': 'fraction', 'deg': 'angle',
    '°': 'angle', 'K': 'temperature', 'N/mm': 'stiffness',
}
# factor -> what frame pair would explain it. Measured in this project unless marked generic.
FRAMES = {
    0.6823: 'lens plane -> spectacle plane, measured per eye in iol_decision.py',
    0.6916: 'lens plane -> spectacle plane at the mean sensitivity dR/dP',
    1.775: 'cylinder magnitude -> power-vector length, measured on the same 69 eyes',
    0.1037837: 'micromolar -> molecules in the declared 0.016 um^3 volume',
    0.132: 'stated 1000 N -> area x peak pressure, Rivarola FE row',
    1.5235: 'applied axial force -> nucleus pressure x disc area, Nachemson k',
    1000.0: 'decade slip of three (generic)',
    1e6: 'decade slip of six (generic)',
    60.0: 'per second -> per minute (generic)',
    3600.0: 'per second -> per hour (generic)',
}
for n in (3, 4, 5, 6, 7, 8, 10, 16, 20, 24, 30, 38, 69):
    FRAMES[round(math.sqrt(n), 6)] = f'SD -> SEM at n = {n} (generic)'
NUM = re.compile(r'(?<![\w.])(-?\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)\s*([A-Za-zµ°%/]+[0-9²³]?)')
REL_TOL = 0.015      # 1.5 percent: tight enough that a coincidence is unlikely to land on a factor


def numbers_in(text: str):
    for val, unit in NUM.findall(text):
        dim = UNITS.get(unit)
        if not dim:
            continue
        try:
            v = float(val.replace(',', '.'))
        except ValueError:
            continue
        if v != 0:
            yield abs(v), dim, f'{val} {unit}'


def count_hits(pool) -> tuple[list, dict]:
    hits, per_factor = [], {}
    for (v1, d1, l1, e1), (v2, d2, l2, e2) in itertools.combinations(pool, 2):
        if d1 != d2 or e1 == e2:
            continue
        ratio = max(v1, v2) / min(v1, v2)
        for f, why in FRAMES.items():
            probe = f if f > 1.0 else 1.0 / f
            if abs(ratio / probe - 1.0) <= REL_TOL:
                hits.append({'a': l1, 'a_edge': e1, 'b': l2, 'b_edge': e2, 'dimension': d1,
                             'ratio': round(ratio, 6), 'factor': f, 'frame': why})
                per_factor[f] = per_factor.get(f, 0) + 1
                break
    return hits, per_factor


def null_distribution(pool, draws: int = 40) -> dict:
    """How many hits does CHANCE give? Shuffle the values across records, keep the dimensions.

    Without this the scan is a sieve, not an instrument: with 207 numbers, 21 321 same-dimension
    pairs and 22 candidate factors at 1.5 percent tolerance, a few hundred hits are what randomness
    produces. A factor only means something if it beats its own null.
    """
    import random
    rng = random.Random(0)
    by_dim = {}
    for v, d, lit, eid in pool:
        by_dim.setdefault(d, []).append(v)
    counts = {}
    for _ in range(draws):
        shuffled = []
        pos = {d: 0 for d in by_dim}
        order = {d: rng.sample(vals, len(vals)) for d, vals in by_dim.items()}
        for v, d, lit, eid in pool:
            shuffled.append((order[d][pos[d]], d, lit, eid))
            pos[d] += 1
        _, pf = count_hits(shuffled)
        for f, n in pf.items():
            counts.setdefault(f, []).append(n)
    return {f: sum(ns) / len(ns) for f, ns in counts.items()}


def main() -> None:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    # The first pool was net edges only, 207 numbers, and nothing beat chance. That was the wrong
    # corpus: every collision found by hand sat between a LANE OUTPUT and the net, or between two
    # lanes -- the toric plane pair came from LENS_RESOLUTION r35 against iol_decision, the
    # magnitude/vector pair from CORNEA_SHAPE r24 against the assembly chain. Within the net each
    # measurement appears about once, so there is nothing for a frame shift to duplicate. Pool the
    # lane rounds and the assembly outputs alongside the net, which is where the copies live.
    pool = []
    for e in net['edges']:
        text = ' '.join(str(e.get(k, '')) for k in ('constraint', 'sensitivity_note'))
        for v, dim, lit in numbers_in(text):
            pool.append((v, dim, lit, 'net:' + e['id']))

    def flatten(obj, out):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    out.append((k, v))
                else:
                    flatten(v, out)
        elif isinstance(obj, list):
            for v in obj:
                flatten(v, out)

    import re as _re
    UNIT_IN_NAME = [('_D', 'dioptre'), ('_mm', 'length'), ('_um', 'length'), ('_N', 'force'),
                    ('_MPa', 'pressure'), ('_Pa', 'pressure'), ('_uM', 'concentration'),
                    ('_s', 'time'), ('_pct', 'fraction'), ('_percent', 'fraction'),
                    ('_molecules', 'count'), ('_K', 'temperature'), ('_pm', 'length')]
    sources = sorted(W.glob('results/LANE_*/night_rounds/r*.json')) + \
        sorted(W.glob('results/ASSEMBLY_*/*.json'))
    for f in sources:
        try:
            doc = json.loads(f.read_text())
        except Exception:
            continue
        pairs = []
        flatten(doc, pairs)
        tag = f.parent.name if f.parent.name.startswith('ASSEMBLY') else \
            f.parent.parent.name + '/' + f.stem
        for name, val in pairs:
            if val == 0 or abs(val) > 1e12:
                continue
            for suffix, dim in UNIT_IN_NAME:
                if name.endswith(suffix) or suffix.strip('_') in name.split('_')[-1:]:
                    pool.append((abs(float(val)), dim, f'{name}={val:.6g}', tag))
                    break
    hits, per_factor = count_hits(pool)
    null = null_distribution(pool)
    out = W / 'results/ASSEMBLY_FRAME_SHIFT'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'HITS.json').write_text(json.dumps(
        {'numbers_pooled': len(pool), 'pairs_tested': len(pool) * (len(pool) - 1) // 2,
         'hits': hits, 'hits_per_factor': {str(k): v for k, v in per_factor.items()},
         'null_expectation_per_factor': {str(k): v for k, v in null.items()},
         'rel_tol': REL_TOL, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=1,
        ensure_ascii=False))
    print(f'speech in the net with known dimension: {len(pool)}   couples tried: {len(pool) * (len(pool) - 1) // 2}   hits: {len(hits)}')
    print("\nhits against chance expectation (40 shuffles, dimensions preserved):")
    print(f"  {'faktor':>12} {'obs':>5} {'slump':>7} {'kvot':>6}  Declaration")
    rows = []
    for f, n in sorted(per_factor.items(), key=lambda kv: -kv[1]):
        exp = null.get(f, 0.0)
        lift = n / exp if exp > 0 else float('inf')
        rows.append((lift, f, n, exp))
    for lift, f, n, exp in sorted(rows, reverse=True):
        mark = "  <-- above chanceen" if lift > 2.0 and n >= 4 else ''
        print(f'  {f:>12} {n:>5} {exp:>7.1f} {lift:>6.2f}  {FRAMES[f][:52]}{mark}')
    print("\nthe ten most specific hitsna:")
    rare = sorted(hits, key=lambda h: per_factor[h['factor']])
    for h in rare[:10]:
        print(f"  {h['a']:>12} ({h['a_edge'][:18]:<18}) mot {h['b']:>12} "
              f"({h['b_edge'][:18]:<18}) kvot {h['ratio']:.4g} -> {h['frame'][:46]}")
    print('wrote', out / 'HITS.json')


if __name__ == '__main__':
    main()
