#!/usr/bin/env python3
"""Admit the candidate tracks that are already consumable, and queue each against its consuming cell.

WEIGHTING, corrected 2026-10-03 after the operator put it right: the 500-entry program is a
MACHINE-GENERATED candidate list from a brainstorm session with a weaker model, not the operator's
own directives. The signal was in the data and I read past it -- the comparative improvement field, the status and
the time window are byte-identical across all 500, which a hand-written idea list is not. So these
carry no authority of their own: they are cheap, bound, falsifiable jobs and nothing more. Admitting
one is not recirculation against the operator's seeds, and the base rate for value should be set as
for generated candidates.

Why. The 500-entry candidate program has sat at `PROPOSED_NOT_ADMITTED` for every single seed since
2026-09-29, and dental confirmed on 2026-10-03 that no admission routine exists anywhere. A boundary
survey then measured which of them are consumable TODAY, by one strict criterion: the seed's observable
must carry a dimensioned unit that a named cell actually computes, matched against the 1057 unit strings
in the cells. Forty of 500 pass. The ten strongest have a named consuming cell.

The criterion earned its strictness. A first pass gave 86 matches and a second 42; the difference was
three traps that had to be fixed. Case-insensitive normalisation collapsed micromolar into micrometre
and produced eleven false pairs. Wall-clock seconds and minutes matched 46 of the 86 while being cost
rather than an observable. And a substring search for a phosphor material found 437 hits that were all
phosphorylation. The surviving forty are after all three corrections.

What admission means here: the seed becomes a job bound to the cell whose output IS its observable, with
the seed's own comparative improvement standard carried into the brief. It does not mean the seed is accepted as true.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('.')
SEEDS = W / 'tasks/free48/SEED_PROGRAM_500_IMPROVED_20260929.json'
QUEUE = W / 'tasks/lanes/bt_queue.txt'
RESULTS = W / 'results'

# Seed id -> (consuming cell or lane, the shared observable with its unit) from the boundary survey.
ADMIT = {
    'BT-FW48-SEED-043': ('Q154', 'balance residual, mol/s, over three mass pools'),
    'BT-FW48-SEED-210': ('Q168', 'balance residual, mol/s, consumed by an identifiability analysis'),
    'BT-FW48-SEED-345': ('Q156', 'balance residual, mol/s, with an explicit boundary flux'),
    'BT-FW48-SEED-083': ('Q146', 'organ-to-organ flux, mol/s, with a declared receiver'),
    'BT-FW48-SEED-211': ('Q146', 'organ-to-organ flux, mol/s, with a declared receiver'),
    'BT-FW48-SEED-228': ('SURG_INCISION', 'boundary force in N and deformation in mm'),
    'BT-FW48-SEED-127': ('Q100', 'frequency in Hz and amplitude in Pa, self-oscillation'),
    'BT-FW48-SEED-010': ('Q140', 'clearance in mL/min, already split into three pathways'),
    'BT-FW48-SEED-079': ('LANE_COUPLE_END_TO_END', 'clearance in L/h with preserved provenance'),
    'BT-FW48-SEED-085': ('LANE_EXTERNAL_SOLVER', 'pressure or modulus in Pa, independent unloading tangent'),
    'BT-FW48-SEED-126': ('Q052', 'pressure or modulus in Pa, stiffness and permeability from one specimen'),
    'BT-FW48-SEED-006': ('IMMUNITY', 'cytokine response in pg/mL, half-response concentration'),
    'BT-FW48-SEED-186': ('IMMUNITY', 'secondary response in pg/mL'),
}

BLOCK = re.compile(r'excluded_category|excluded_category|excluded_category|sinusoid|excluded_category|excluded_category|excluded_category', re.I)


def brief(seed: dict, cell: str, observable: str) -> str:
    def g(k, d=''):
        v = seed.get(k, d)
        return v if isinstance(v, str) else json.dumps(v, ensure_ascii=False)
    return f"""# {seed.get('job_id')} — admitted seed, bound to the cell that computes its observable

This candidate track has been sitting unadmitted since 2026-09-29 together with 499 others. It
comes from a generated program, not from the operator, so it carries no authority of its own. It is admitted now
because its observable is a quantity **{cell}** already computes: {observable}. That is the whole
reason, and it is the only kind of binding that does not require a new measurement.

**{cell} is shipped with this job, under `inputs/{cell}/`.** The first admitted seed stopped at step 1
because the cell was named but not shipped, and the sandbox correctly refused reads outside the job
directory. Read it there; do not look for it elsewhere on the host.

## The seed, as written
**Decision.** {g('decision')[:1200]}

**Observable.** {g('observable')[:400]}

**Regime.** {g('regime')[:300]}

**First principles.** {g('first_principles')[:600]}

**Twist.** {g('twist')[:400]}

## The standard the seed sets for itself
{g('outclass')[:400]}

Take that literally. The control is the **strongest eligible** supplied or domain control, never one
weakened to make a win visible. If you cannot name the strongest control, say so and stop — an
unnamed control makes the result unusable.

## What to do
1. **Read {cell} first** and report what it already computes in that unit, with the variable name. If
   the observable is not in fact its output, say so and stop: the admission was wrong and that is worth
   knowing.
2. **Run the seed's test** against the strongest control, and report one decisive number with its unit
   and what it is measured against.
3. **State the falsifier before computing**, and report it even when it fires. A fired falsifier is a
   result.
4. **If a summary statistic decides anything**, test its sufficiency: two states with the identical
   summary to machine precision, and the downstream difference. This failed in 83 per cent of a random
   sample of our own contracts, so it is not a formality.
5. **A missing measurement leaves as an orderable item** — quantity, unit, what it decides — in
   `ACQUISITION_TARGETS.json` in this directory, not as a sentence.

Status `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation. Write `RESULTS.md` starting
with the job id, plus `results.json` with your numbers.
"""


def main() -> int:
    seeds = {s.get('job_id'): s for s in json.loads(SEEDS.read_text())}
    made, missing, blocked = [], [], []
    for jid, (cell, observable) in ADMIT.items():
        s = seeds.get(jid)
        if s is None:
            missing.append(jid)
            continue
        text = brief(s, cell, observable)
        if BLOCK.search(text):
            blocked.append(jid)
            continue
        out = RESULTS / f'BT-SEED-ADMIT-{jid.replace("BT-FW48-SEED-", "")}'
        if (out / 'RESULTS.md').exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        (out / 'BRIEF.md').write_text(text)
        json.dump({'id': out.name, 'kind': 'admitted_seed', 'category': 'mechanism',
                   'seed_job_id': jid, 'consuming_cell': cell, 'observable': observable,
                   'seed_status_before': s.get('status'),
                   'breakthrough_priority': True,
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (out / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        # Ship the consuming cell with the job. Naming it is not enough: the sandbox denies reads
        # outside the job directory, which is what stopped BT-SEED-ADMIT-043 at step 1.
        src = W / 'tasks/free48/sources' / cell
        dst = out / 'inputs' / cell
        if src.is_dir() and not dst.exists():
            import shutil
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns('*.npz', '*.npy', '.bak*',
                                                                    '__pycache__'))
        elif not src.is_dir():
            lane = W / 'results' / cell
            if lane.is_dir():
                import shutil
                dst.mkdir(parents=True, exist_ok=True)
                for f in list(lane.glob('PORT*.json'))[:3] + list(lane.glob('RESULTS.md'))[:1]:
                    if f.stat().st_size < 2_000_000:
                        shutil.copy2(f, dst / f.name)
        made.append(out.name)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with QUEUE.open('a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    print(f'admitted {len(made)} of {len(ADMIT)}   seed id not found: {len(missing)}   '
          f'blocked by private filter: {len(blocked)}')
    for m in made:
        print('  ' + m)
    if missing:
        print('  not found: ' + ', '.join(missing))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
