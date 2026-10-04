"""More creative briefs, built from tonight's verified anomalies and from what each lane says blocks it.

The operator wants more of these because they produced the external anchors. The net-derived generator
is exhausted at 56 edges, and three attempts to mine the 13124-report corpus for more material all
failed for the same reason: the same key name means different things in different jobs. So the material
comes from two places that do carry meaning, because a person or a lane wrote it deliberately.

First, the anomalies measured and verified tonight. Each is a pair of numbers that cannot both be
comfortable, every one reproduced by me from the raw file rather than taken from a summary.

Second, each running lane's own stated obstacle. A lane that gates FAIL has written down precisely what
it could not get, in its own words, with its own number next to it. That is a better question than
anything I would invent, and nothing was consuming it.

Same form as the briefs that worked: the numbers bare, no interpretation of mine, no named fields to
search in, three framings before any search, and no required output shape.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('')
QUEUE = W / 'tasks/lanes/bt_queue.txt'
# The block list is not kept in the repo. tasks/assembly/excluded_terms.py says why, loads it from
# outside the tree, and matches everything if it cannot be read, so a missing list rejects rather
# than admits.
def _load_block_pattern(extra: str = '') -> 're.Pattern[str]':
    import importlib.util, pathlib
    for parent in pathlib.Path(__file__).resolve().parents:
        cand = parent / 'tasks' / 'assembly' / 'excluded_terms.py'
        if cand.exists():
            spec = importlib.util.spec_from_file_location('excluded_terms', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            pat = mod.block_pattern().pattern
            return re.compile(pat + ('|' + extra if extra else ''), re.I)
    return re.compile(r'(?s).*')  # loader gone: reject everything rather than pass everything


BLOCK = _load_block_pattern()

ASK = """

Before you look for anything: write down **three different ways of looking at this**, with
at least one that has nothing to do with tissue or the human body at all. Diverge first.

Then choose what you think goes furthest and pursue it. What do you find?

No particular answer format. Give DOI or PMID for what you build on. If you conclude
that the material above points somewhere other than the question, say so.

Status PENDING_INDEPENDENT_REVIEW.


## Diverge widely first

This is not a review job. I want ideas. Write three readings before searching, and let at least one
lie entirely outside tissue, body and medicine. Diverge, choose the one that goes furthest and follow
it all the way. What I want to see is something that is not in the material: a connection nobody has
posed, a quantity that would settle the question if someone measured it, or a way to make the question
decidable with something already measured elsewhere. DOI or PMID for what you build on.

If you see along the way that two numbers do not agree dimensionally,
say so in one line and move on. A bonus, not the task.

"""

# Verified tonight, each number reproduced from the raw file by me.
ANOMALIES = {
    'CORNEA-MODEL-OVERSTATES': """We put together the corneal astigmatism of a front side and a back side.
An apparatus measures the entire cornea in one step.

| | |
|---|---|
| our composite | 2,029 D |
| directly measured in one step | 1,831 D |
| mean difference per eye | 0,292 D |
| the error in the decision made from the model | 0,282 D |

The difference is greater than the error in the decision the model is used for. The decision still works.""",

    'ORTHOGONAL-AXES': """Two devices measured the back of the cornea on the same 69 eyes.

| | |
|---|---|
| axis difference, median | 86° |
| eyes within 15° | 0 of 69 |
| eyes between 75° and 105° | 65 of 69 |

On the front of the cornea, same eyes, same apparatus: 3,55° in median, 67 of 69 within 15°.""",

    'DEVICE-CHANGES-THE-LENS': """Same 83 eyes, same manufacturer catalog, two devices as input.

| | |
|---|---|
| eyes where the cylinder label becomes another | 50 of 83 |
| eyes where sphere+cylinder becomes another lens | 63 by 83 |
| worst-case failure, selection by device A | 0,559 D |
| worst-case-error, selection by device B | 0,544 D |
| worst-case failure, lens that holds under both | 0,372 D |

The residual versus outcome is almost equal for A and B: 0,282 versus 0,293 D.""",

    'AXIAL-LENGTH-SINGLE-SOURCE': """The length of the eye is the largest single term in the refraction calculation.
We have it from a device. No other value exists for any of the 89 eyes.

The quantities in which two devices exist differ systematically, not randomly:

| quantity | difference | dispersion |
|---|---|---|
| hornhinnans tjocklek | 9,90 µm | 6,26 |
| anterior chamber depth | −0,144 mm | 0,081 |
| medelhornhinnestyrka | +0,851 D | 0,361 |""",

    'CONFIDENCE-CARRIES-NOTHING': """Two decisions scored as probabilities against the cohort's own quota
as reference prognosis, 20 and 69 patients.

| decision | Brier | reference | skill | resolution |
|---|---|---|---|---|
| linsstyrka | 0,2559 | 0,2659 | +0,038 | 0,0817 |
| toriskt | 0,1384 | 0,1276 | −0,085 | 0,0034 |

The range of forecasts was 0,46 to 1,00. The decisions themselves beat the control clearly: 0,535 vs. 1,130 D
and 0,282 against 0,626 D.""",

    'ERROR-BELOW-ITS-OWN-SIGMA': """Three independent measured uncertainty terms in the input data give per eye one
composite standard deviation. The outcomes compared to it:

| | |
|---|---|
| eyes with errors UNDER own composite sigma | 50 of 89 |
| deras medelfel | 0,204 D |
| deras sammansatta sigma | 0,451 D |
| de 39 other, average error | 0,959 D |
| sammansatt sigma i medel | 0,456 D |
| faktiskt medelfel | 0,535 D |
| correlation sigma versus actual error, per eye | +0,127 |""",

    'PREDICTION-BEATS-MEASUREMENT': """The position of the lens after the operation determines the result: 1 mm becomes 1,349 D.

| | |
|---|---|
| decision with the MEASURED mode | 0,569 D|
| decision with a PREDICTED position from preoperative measurements | 0,535 D|
| own mean error of the prediction | 0,114 mm |
| the quantity's dispersion | 0,284 mm |

Guessing gave a better decision than measuring the same thing.""",

    'LABEL-WIDER-THAN-THE-STEP': """Lenses are manufactured in steps of 0,5 D and we choose a step.

| power range | permissible deviation according to standard |
|---|---|
| 0–15,00 D | ±0,30 D |
| 15,50–25,00 D | ±0,40 D |
| 25,50–30,00 D | ±0,50 D |
| over 30 D | ±1,00 D |

68 of our 89 lenses are in the ±0,40 band. Measured deviations in the literature: 0,18 ± 0,12 D vs. exact
label, 0,24 ± 0,11 D at 20 D and above. We have never measured a lens.""",
}


def main() -> int:
    made = []
    for tag, body in ANOMALIES.items():
        if BLOCK.search(body):
            continue
        jid = f'BT-ANOM-{tag}'
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / 'BRIEF.md').write_text(body.strip() + ASK)
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'creative_from_verified_anomaly', 'anomaly': tag,
                   'brief_style': 'raw_material_no_interpretation_no_named_analogies_diverge_first',
                   'required_output_form': None, 'category': 'exploratory',
                   'claim_type': 'information_link',
                   'numbers_reproduced_by_coordinator_from_raw_files': True,
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    # Each running lane's own stated obstacle, in its own words, with its own number beside it.
    for state in sorted((W / 'tasks/build_night/state').glob('*.round')):
        lane = state.stem
        rnd = int(state.read_text().strip() or 0)
        rep = None
        for r in range(rnd, max(rnd - 3, 0), -1):
            p = W / 'results' / lane / 'night_rounds' / f'r{r}.json'
            if p.exists():
                rep = p
                break
        if rep is None:
            continue
        try:
            dd = json.loads(rep.read_text())
        except Exception:
            continue
        obstacle = str(dd.get('obstacle') or '').strip()
        if len(obstacle) < 40 or BLOCK.search(obstacle):
            continue
        jid = 'BT-OBST-' + re.sub(r'^LANE_', '', lane).replace('_', '-')
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / 'BRIEF.md').write_text(
            f"# {re.sub(r'^LANE_', '', lane).replace('_', ' ').lower()}\n\n"
            f"One of our own runs stopped and wrote down what it couldn't get. Literally: \n\n"
            f"> {obstacle}\n\n"
            f"What it tried to do, also verbatim:\n\n"
            f"> {str(dd.get('operation') or '-')[:400]}\n\n"
            f"What it compared against:\n\n"
            f"> {str(dd.get('control') or '-')[:300]}\n"
            + ASK)
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'creative_from_lane_obstacle', 'lane': lane, 'round': rnd,
                   'brief_style': 'raw_material_no_interpretation_no_named_analogies_diverge_first',
                   'required_output_form': None, 'category': 'exploratory',
                   'claim_type': 'information_link',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    if made:
        slots = ('A', 'B', 'C', 'D')
        lines = [f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)]
        existing = QUEUE.read_text().splitlines() if QUEUE.exists() else []
        QUEUE.write_text('\n'.join(lines + existing) + '\n')
    print(f'{len(made)} new creative jobs, first in line')
    for m in made:
        print('  ' + m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
