"""Briefs written as a creative practitioner would write them, with the graph supplying the tension.

Three corrections from the operator, in order, and each one landed.

1. Stop steering. My briefs ran to thirty lines of required form and the jobs filled the form.
2. The graph should be part of it. It was not, and generating briefs from it immediately exposed two
   defects: a node named in a constraint but not an endpoint of it, invisible to any walk; and an edge
   still handing out a number I had withdrawn hours earlier.
3. Prompt like a human creativity expert, not an auditor. This is the one that reaches furthest. I was
   handing out status labels and gap figures and asking which edge a finding moves. That is an audit
   request wearing the clothes of a question, and it can be satisfied by lookup.

What a practitioner does instead: put a concrete contradiction in front of someone, told as a fact
rather than a status code; ask something that cannot be answered by looking it up; leave the output
shape alone; and invite the move that actually worked in the last batch -- importing a practice from a
field that solved this long ago. That lens, "what does another field do routinely that we do not do at
all", produced nine of the thirty-one anchored proposals and the only one worth building.

So the graph's job here is narrow and useful: find the pairs of numbers that cannot both be comfortable,
and state them plainly. The question stays open.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')

# Each: the anomaly as a plain fact, then the open question. The numbers come from the net and from
# tonight's cells; the framing does not.
BRIEFS = {
    'AXIAL-LENGTH': (
        'Eye length',
        """The eye's length is the largest single term in the entire refraction calculation. We have it from a
single device, and we have never seen a second value for the same eye.

The two devices we *can* compare do not disagree randomly. They differ systematically:
9,90 mikrometer in corneal thickness, 0,144 millimeter in anterior chamber depth. Every time. This is
not noise, these are two devices that agree with themselves and disagree with each other.

If the same applies to length, there is an error larger than anything we measured tonight — and we cannot
see it, because we only have one measurement.

What do you do about a quantity you can only measure in one way? Other fields have been in that situation:
metrology, surveying, astronomy, manufacturing control. What do they do?""",
    ),
    'POSTERIOR-CORNEAL-CYLINDER': (
        'Two devices, orthogonal answers',
        """Two devices measured the same 69 eyes. On the back of the cornea they report axes that lie
**86 grader apart** at the median — 65 of 69 eyes between 75 and 105 grader, not a single one within 15. They are
therefore systematically perpendicular to each other.

When we took the axis as given our result was 0,758 dioptrier, worse than guessing a
population average. When we rotated 90 grader it became 0,293 — almost exactly like the other device.
So it was a convention, not a disagreement.

But we discovered it by chance, because the number became unreasonable. In the same eyes the axes agree for
the *front* of the cornea, so there was no rule to carry over.

How do you know that two measurements of the same thing are expressed in the same convention, when no one has written down
the convention? This is not an eye problem.""",
    ),
    'CORNEAL-THICKNESS': (
        'The model overestimates',
        """We build the cornea from a front plus a back and calculate its astigmatism. CASIA
measures the whole cornea in one step. We get 2,029 dioptrier, the direct measurement gives 1,831. The difference is
0,292 — **larger than the error in the decisions we make from the model**.

So the model is wrong by more than what it is meant to decide, and it still works: the decision
lands at 0,28 against practice at 0,63.

That is uncomfortable. Either two errors cancel each other out, or the two things do not measure the same quantity
despite having the same name.

Which? And how could you determine that?""",
    ),
    'LENS-POSITION': (
        'Where the lens ends up',
        """Where the artificial lens ends up in the eye determines the final result: one millimeter becomes 1,35
dioptrier. But you only know afterwards, when it is in place.

We predict it from measurements taken before surgery and get within 0,114 millimeter, against the quantity's own
spread of 0,284. That is enough — 0,114 millimeter becomes 0,15 dioptrier, below the clinical threshold.

The odd thing: predicting the position gave a **better** decision than entering the measured position.
0,535 against 0,569 dioptrier. So the prediction is more useful than the measurement of the same thing.

Why would a guess beat a measurement? And where is the limit to that?""",
    ),
    'IOL-POWER-LABEL': (
        'Is the label correct',
        """The lenses come in steps of 0,5 dioptrier and we choose a step. We treat the number on
the lens as exact.

Our decision moves to another step in 44 of 89 eyes. So almost half the cases are decided by
half a dioptri.

If the manufacturing tolerance approaches half a step, there is a floor that no model can go
below — and we would not notice it, because we have never measured a lens.

How do other industries find out whether the component matches its label, without measuring every specimen?""",
    ),
    'STROMAL-INDEX': (
        'A number we calculated ourselves',
        """The cornea's refractive index enters every optical calculation we make. The number comes from our own
hydration model — we have never compared it with an independently measured value.

Around that point there are four relations we have not managed to determine. The nearest says that thickness,
index, hydration and surface distribution *together* leave 0,306 dioptrier undetermined, and that each individual
pair is exactly identical in four out of four cases. So the quantities cannot be distinguished with what we measure
today.

Measuring one of them better does not help if they cannot be separated.

How do you separate quantities that only appear together? That problem has been solved in several fields.""",
    ),
}

TAIL = """
You do not need to answer in any particular form. Say what you conclude, and give DOI or PMID for what
you build on so we can follow it. If you conclude that the question is wrongly posed, say that instead.

Status PENDING_INDEPENDENT_REVIEW.
"""


def main() -> int:
    for tag, (title, body) in BRIEFS.items():
        d = W / 'results' / f'BT-2ND-{tag}'
        d.mkdir(parents=True, exist_ok=True)
        (d / 'BRIEF.md').write_text(f'# {title}\n\n{body.strip()}\n{TAIL}')
        j = d / 'JOB.json'
        meta = json.loads(j.read_text()) if j.exists() else {'id': f'BT-2ND-{tag}'}
        meta.update({'brief_style': 'creative_practitioner_graph_supplies_the_tension',
                     'required_output_form': None,
                     'review_state': 'PENDING_INDEPENDENT_REVIEW'})
        json.dump(meta, j.open('w'), ensure_ascii=False, indent=1)
        print(f"  {tag:28s} {len((d / 'BRIEF.md').read_text()):>5d} characters")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
