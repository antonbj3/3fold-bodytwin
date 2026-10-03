"""A decision for a robot that moves a spine segment: what range of motion keeps every ligament safe?

Why this is the seed rather than a side quest. The seed is tissue behaviour coupled to geometry as
support for medical robots, and a constitutive relation is only support for a robot if it answers the
question the robot has to answer: how far may I move this, and which structure stops me first. The
ligament table answers it and nobody has asked.

The material, read from `MECHANISM_SPINE_LIGAMENTS.md` and used as printed. Five ligaments, force in
newtons at every 20 degrees from -90 to +90:

  ALL  loads only in extension, peaking at 13.014 N near +50
  PLL  loads only in flexion, and barely, 0.525 N at its largest
  LF   loads in flexion, peaking 2.998 N near -10
  ISL  loads in flexion, peaking 0.869 N near -50
  SSL  loads in flexion, peaking 6.513 N near -50

The source also pre-registered its bands BEFORE sweeping, and states that the bands are inside the
model's permissive coordinate limit rather than a measured range of motion. That honesty note is why
this cell reports a MODEL range and never a clinical one.

The decision, and the part that needs the whole table rather than one ligament: for a tension limit,
find the widest angular interval containing the neutral position in which EVERY ligament stays under it,
and name the ligament that binds at each end. A limit that only consults the stiffest structure gets the
extension side right and the flexion side wrong, because the binding structure changes side.

The control is the single-ligament reading, which is what using ALL alone would give, and it is the
comparison that shows whether the coupling buys anything.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_LIGAMENT_ROM'
SRC = 'source_documents/MECHANISM_SPINE_LIGAMENTS.md'

# angle -> {ligament: newtons}, exactly as printed
TABLE = {
    -90: {'ALL': 0.000, 'PLL': 0.000, 'LF': 0.000, 'ISL': 0.000, 'SSL': 3.806},
    -70: {'ALL': 0.000, 'PLL': 0.000, 'LF': 0.003, 'ISL': 0.286, 'SSL': 6.422},
    -50: {'ALL': 0.000, 'PLL': 0.000, 'LF': 1.252, 'ISL': 0.869, 'SSL': 6.513},
    -30: {'ALL': 0.000, 'PLL': 0.018, 'LF': 2.862, 'ISL': 0.694, 'SSL': 4.075},
    -10: {'ALL': 0.000, 'PLL': 0.525, 'LF': 2.998, 'ISL': 0.060, 'SSL': 0.417},
     10: {'ALL': 4.811, 'PLL': 0.178, 'LF': 1.526, 'ISL': 0.000, 'SSL': 0.000},
     30: {'ALL': 12.388, 'PLL': 0.000, 'LF': 0.055, 'ISL': 0.000, 'SSL': 0.000},
     50: {'ALL': 13.014, 'PLL': 0.000, 'LF': 0.000, 'ISL': 0.000, 'SSL': 0.000},
     70: {'ALL': 6.141, 'PLL': 0.000, 'LF': 0.000, 'ISL': 0.000, 'SSL': 0.000},
     90: {'ALL': 0.022, 'PLL': 0.000, 'LF': 0.000, 'ISL': 0.000, 'SSL': 0.000},
}
LIGAMENTS = ('ALL', 'PLL', 'LF', 'ISL', 'SSL')


def widest_interval(limit: float, ligaments=LIGAMENTS) -> dict:
    """The widest run of sampled angles around 0 where every named ligament stays at or under limit."""
    angles = sorted(TABLE)
    ok = {a: all(TABLE[a][g] <= limit for g in ligaments) for a in angles}
    # walk outward from the two angles nearest neutral
    lo = hi = None
    for a in [x for x in angles if x > 0]:
        if ok[a]:
            hi = a
        else:
            break
    for a in [x for x in angles if x < 0][::-1]:
        if ok[a]:
            lo = a
        else:
            break
    binding_hi = binding_lo = None
    if hi is not None:
        nxt = [x for x in angles if x > hi]
        if nxt:
            binding_hi = max(ligaments, key=lambda g: TABLE[nxt[0]][g])
    if lo is not None:
        prv = [x for x in angles if x < lo]
        if prv:
            binding_lo = max(ligaments, key=lambda g: TABLE[prv[-1]][g])
    return {'low_deg': lo, 'high_deg': hi,
            'width_deg': (hi - lo) if (lo is not None and hi is not None) else None,
            'binding_in_extension': binding_hi, 'binding_in_flexion': binding_lo}


def main() -> int:
    rows = []
    for limit in (0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 13.0):
        allg = widest_interval(limit)
        only_all = widest_interval(limit, ('ALL',))
        rows.append({
            'tension_limit_N': limit,
            'all_five_ligaments': allg,
            'using_ALL_alone': only_all,
            'flexion_side_lost_by_using_one_ligament_deg': (
                (allg['low_deg'] - only_all['low_deg'])
                if (allg['low_deg'] is not None and only_all['low_deg'] is not None) else None),
        })

    peaks = {g: max(TABLE, key=lambda a: TABLE[a][g]) for g in LIGAMENTS}
    summary = {
        'question': ('what angular range keeps every ligament under a tension limit, and which '
                     'structure binds at each end'),
        'source': {'path': SRC, 'used_as_printed': True,
                   'source_honesty_note': ('the source states its bands are inside the model permissive '
                                           'coordinate limit and are NOT a measured range of motion, so '
                                           'every range here is a MODEL range')},
        'where_each_ligament_peaks_deg': peaks,
        'peak_force_N': {g: TABLE[peaks[g]][g] for g in LIGAMENTS},
        'the_binding_structure_changes_side': (
            'ALL carries load only in extension and peaks at 13.014 N, while SSL carries load only in '
            'flexion and peaks at 6.513 N. A limit set from the stiffest structure alone therefore '
            'governs one side of neutral and says nothing about the other'),
        'decision_rows': rows,
        'control': 'the same question answered from ALL alone, which is the single-ligament reading',
        'claim_type': 'information_link',
        'scope': ('arithmetic on a published table at 20-degree sampling; the interval endpoints are '
                  'sampled angles, not interpolated, so a true boundary lies within 20 degrees of each'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'LIGAMENT_ROM_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print('  topp per ligament (grader, N): ' + ', '.join(
        f'{g} {peaks[g]:+d}/{TABLE[peaks[g]][g]}' for g in LIGAMENTS))
    print(f"\n  {'limit N':>8s} {'alla fem':>18s} {'bindande flex/ext':>22s} {'bara ALL':>18s} {'tappat':>7s}")
    for r in rows:
        a, o = r['all_five_ligaments'], r['using_ALL_alone']
        rng = f"[{a['low_deg']}, {a['high_deg']}]"
        orng = f"[{o['low_deg']}, {o['high_deg']}]"
        bind = f"{a['binding_in_flexion']}/{a['binding_in_extension']}"
        print(f"  {r['tension_limit_N']:>8.1f} {rng:>18s} {bind:>22s} {orng:>18s} "
              f"{str(r['flexion_side_lost_by_using_one_ligament_deg']):>7s}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
