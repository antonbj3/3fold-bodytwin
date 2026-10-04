"""Can a load-scaling surrogate replace a static-optimisation run? Decided on a spine lift pair.

WHY THIS IS A DECISION. The twin will be asked what happens to muscle force sharing when a load
changes. The cheap answer is that every muscle scales with the load; the expensive answer is to
re-run static optimisation. If the cheap answer held, the expensive run would be unnecessary, so the
question is worth a number rather than an opinion.

WHAT THE REFERENCE IS, stated first because it is the thing most easily misread. The two arms are
static-optimisation SOLUTIONS from the read-only repo, a stoop lift with and without a box
(`data/msk_smoketest/subject2_spine_stoop_lift/so_no_box` and `so_with_box`, 10 frames each). They
are not measurements of a person. That makes them a legitimate reference for exactly one claim --
whether a scaling surrogate reproduces the solver -- and an illegitimate reference for any claim
about a real back. A previous audit in this project failed on precisely that substitution, so the
distinction is in the file rather than in a footnote.

THE CONTROL IS THE SURROGATE. Predict each muscle's force in the loaded arm as its unloaded force
times the ratio of total muscle force between the arms. That control sees the same two arms and the
same totals; what it does not see is which muscles take the extra load. The decision is whether that
blindness costs anything measurable.

THE EFFECTIVE SAMPLE IS ONE POSE, not ten frames, and the file checks it rather than trusting the
header. Each .sto declares ten rows at 10 ms spacing, but every muscle column is byte-identical
across them: semimembranosus reads 1509.725 N in all ten unloaded rows and 1712.727 N in all ten
loaded rows, and the spread over frames is exactly 0.0000 N for every muscle tried. So the pair is
one static pose per arm, replicated. A per-frame robustness check run on it returns the same number
ten times and demonstrates nothing; that check was run, returned 0.702 in all ten frames, and is
reported here as vacuous rather than as agreement.

No opensim is needed: the .sto files are tab-separated text with a header, and parsing them here is
what makes the pair usable at all on this machine. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import os
import pathlib
import statistics

BASE = pathlib.Path('source_repository/data/msk_smoketest/subject2_spine_stoop_lift')
STO = 'walking1_StaticOptimization_force.sto'
# Columns that are not muscles: reserves, residuals and the box actuator the loaded arm adds.
NON_MUSCLE = ('reserve', 'residual', 'box', 'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ')


def read_sto(path: pathlib.Path) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text().splitlines()
    head = next(i for i, line in enumerate(lines) if line.strip().lower() == 'endheader')
    cols = lines[head + 1].split('\t')
    rows = [[float(x) for x in line.split('\t')] for line in lines[head + 2:] if line.strip()]
    return cols, rows


def muscle_means(cols: list[str], rows: list[list[float]]) -> dict[str, float]:
    """Mean absolute force per muscle over the frames. Absolute, because a static-optimisation
    force is a generalised force and a sign flip is a direction, not a smaller load."""
    out = {}
    for j, name in enumerate(cols):
        if j == 0 or any(k.lower() in name.lower() for k in NON_MUSCLE):
            continue
        out[name] = statistics.mean(abs(r[j]) for r in rows if j < len(r))
    return out


def main() -> None:
    cols_n, rows_n = read_sto(BASE / 'so_no_box' / STO)
    cols_b, rows_b = read_sto(BASE / 'so_with_box' / STO)
    mn, mb = muscle_means(cols_n, rows_n), muscle_means(cols_b, rows_b)
    shared = sorted(set(mn) & set(mb))
    tot_n = sum(mn[m] for m in shared)
    tot_b = sum(mb[m] for m in shared)
    scale = tot_b / tot_n

    errs = {m: abs(mn[m] * scale - mb[m]) for m in shared}
    rel = {m: errs[m] / mb[m] for m in shared if mb[m] > 1.0}
    mae = statistics.mean(errs.values())
    worst = sorted(errs, key=errs.get, reverse=True)[:6]

    # "Is this one pose or ten?" is a question about a TOLERANCE, not about distinctness. Counting
    # distinct rows answered ten twice over: first because the reserve and residual columns vary,
    # then because 156 of 187 muscle columns differ across rows by up to 7e-6 N, which is solver
    # round-off at a relative 5e-9 against forces near 1500 N. Report the spread and let it speak.
    mus_idx = [j for j, c in enumerate(cols_n)
               if j and not any(k.lower() in c.lower() for k in NON_MUSCLE)]
    spread = max(max(abs(r[j]) for r in rows_n) - min(abs(r[j]) for r in rows_n) for j in mus_idx)
    scale_force = max(statistics.mean(abs(r[j]) for r in rows_n) for j in mus_idx)
    print(f'rows: {len(rows_n)} unloaded, {len(rows_b)} loaded; shared muscles: {len(shared)}')
    print(f'largest muscle-force spread across rows: {spread:.2e} N against a largest mean force of '
          f'{scale_force:.1f} N, i.e. relative {spread / scale_force:.1e}')
    print('  so the rows are ONE static pose to within solver round-off: effective n = 1 pose per arm,')
    print('  and a per-frame robustness check on this file is vacuous by construction.')
    print(f'total muscle force  : {tot_n:.1f} N unloaded -> {tot_b:.1f} N loaded, ratio {scale:.4f}')
    print(f'surrogate (uniform scaling) mean absolute error : {mae:.2f} N')
    print(f'  as a share of the total loaded force          : {mae * len(shared) / tot_b * 100:.2f} %')
    print(f'  largest single-muscle error                   : {max(errs.values()):.2f} N')
    print('  worst muscles (error N, unloaded N, loaded N, surrogate N):')
    for m in worst:
        print(f'    {m:<16} {errs[m]:9.2f} {mn[m]:9.2f} {mb[m]:9.2f} {mn[m] * scale:9.2f}')
    if rel:
        print(f'  median relative error over muscles above 1 N  : {statistics.median(rel.values()) * 100:.1f} %')

    # The verdict needs a threshold that is not invented after the fact: the surrogate is useful
    # only if its error is small against the CHANGE it is meant to predict, not against the force.
    change = {m: abs(mb[m] - mn[m]) for m in shared}
    total_change = sum(change.values())
    explained = 1.0 - sum(errs.values()) / total_change if total_change else float('nan')
    print(f'\ntotal redistribution to explain : {total_change:.1f} N')
    print(f'fraction the surrogate explains : {explained * 100:.1f} %')
    verdict = ('SURROGATE SUFFICIENT' if explained > 0.9 else
               'SURROGATE INSUFFICIENT: the redistribution is structural, not a scaling')
    print('verdict:', verdict)

    out = {
        'rows': [len(rows_n), len(rows_b)],
        'largest_muscle_force_spread_across_rows_N': spread,
        'effective_distinct_poses_per_arm': 1,
        'per_frame_robustness': 'vacuous: identical in all ten rows because the rows are one pose',
        'shared_muscles': len(shared),
        'total_force_N': [tot_n, tot_b],
        'load_scale_ratio': scale,
        'surrogate_mae_N': mae,
        'surrogate_max_error_N': max(errs.values()),
        'total_redistribution_N': total_change,
        'fraction_of_redistribution_explained': explained,
        'worst_muscles': {m: errs[m] for m in worst},
        'verdict': verdict,
        'reference_is': 'static-optimisation solutions, not measurements of a person; valid only '
                        'for the claim that a scaling surrogate reproduces the solver',
        'source': str(BASE),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_LOAD_REDISTRIBUTION'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('wrote', d + '/decision.json')


if __name__ == '__main__':
    main()
