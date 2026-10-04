"""A brief per cell asking what it could DECIDE, not what is wrong with it.

The operator's correction, and it lands. Every creative brief so far was built from a contradiction --
two of our numbers that disagree -- and a free-form brief built on a contradiction returns a resolution
of the contradiction. That is an audit in free form, and the audits have been good, but I selected for
them by choosing discrepancies as the raw material. 30 of 37 completed briefs mention a refutation of
one of our own numbers, which is the shape of the question coming back.

Generative material is not "these two numbers disagree". It is "here is what this cell computes, at the
resolution it actually computes it; what could it decide that nobody has asked it to".

Why that is the right axis. Every result that moved the project tonight came from inverting a chain into
a decision and scoring it against something outside us: the implant power, the toric cylinder, the lens
that holds under both instruments, the laser fluence. None of those needed new data. They needed someone
to ask what the forward model could be run backwards to choose. There are around 45 cells and that
question has been asked of four of them.

The material handed over is the cell's own highest-resolution content: what it computes, the quantities
with units it carries, and the external referents it already names. Nothing of mine is added -- no
framing of what matters, no suggested answer, no list of fields to look in.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('')
CELLS = W / 'tasks/free48/sources'
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
UNITY = re.compile(r'["\']([a-z0-9_]*_(?:um|mm|mm2|nm|pa|mpa|kpa|j|n|k|deg|degc|s|h|d|mol|nm_conc|'
                   r'ml_min|percent|ratio|fraction|per_s|m_s|mol_s|pa_s|kg_m3|j_mol))["\']', re.I)
DOI = re.compile(r'10\.\d{4,9}/[^\s"\',\)\]]+')

ASK = """

What can this cell **decide** that nobody has asked it to decide? A forward model can be run
backwards: instead of calculating what happens, choose what produces a desired outcome.

Write three different ways of looking at what it could decide, with at least one that
has nothing to do with tissue or the human body at all. Diverge first, then choose.

For what you choose: what is the decision, what is the goal, and what is outside
us to score it against? A decision without an external reference is an opinion.

No particular answer format. Give DOI or PMID. If the cell cannot
decide anything, say what is missing for it to do so.

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


def main() -> int:
    made, skipped = [], 0
    for cell in sorted(p for p in CELLS.iterdir() if p.is_dir()):
        text = ''
        for f in sorted(cell.rglob('*')):
            if f.is_file() and f.suffix in ('.py', '.md', '.json') and f.stat().st_size < 400_000:
                try:
                    text += f.read_text(errors='replace')
                except Exception:
                    pass
        if not text or BLOCK.search(text):
            skipped += 1
            continue

        quantities = sorted({m.group(1) for m in UNITY.finditer(text)})[:28]
        dois = sorted(set(DOI.findall(text)))[:6]
        # the cell's own stated purpose: first docstring line of its main module
        purpose = ''
        for f in sorted(cell.glob('*.py')):
            m = re.search(r'"""(.{20,400}?)(?:\n\n|""")', f.read_text(errors='replace'), re.S)
            if m:
                purpose = ' '.join(m.group(1).split())[:360]
                break
        if not quantities and not purpose:
            skipped += 1
            continue

        jid = f'BT-CAP-{cell.name}'
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)

        body = [f'# {cell.name}', '']
        if purpose:
            body += ['The cell\'s own description of what it does:', '', f'> {purpose}', '']
        if quantities:
            body += [f'Quantities it carries, with unit in the name ({len(quantities)} of them):', '']
            body += ['  ' + ', '.join(f'`{q}`' for q in quantities), '']
        if dois:
            body += ['External sources it already names:', '']
            body += ['  ' + ', '.join(dois), '']

        (d / 'BRIEF.md').write_text('\n'.join(body).rstrip() + ASK)
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'creative_capability_from_cell', 'cell': cell.name,
                   'quantities_handed_over': len(quantities), 'external_sources_named': len(dois),
                   'brief_style': 'capability_framing_no_contradiction_no_interpretation',
                   'required_output_form': None, 'category': 'exploratory',
                   'claim_type': 'capability',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    if made:
        slots = ('A', 'B', 'C', 'D')
        lines = [f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)]
        existing = [l for l in (QUEUE.read_text().splitlines() if QUEUE.exists() else [])
                    if l.split() and l.split()[-1] not in set(made)]
        QUEUE.write_text('\n'.join(lines + existing) + '\n')
    print(f'{len(made)} capability briefs, {skipped} cells excluded (empty or filtered)')
    for m in made[:12]:
        print('  ' + m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
