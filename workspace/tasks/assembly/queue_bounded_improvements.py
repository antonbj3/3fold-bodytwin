#!/usr/bin/env python3
"""Queue bounded improvement jobs, where each outcome is either a number with a unit or a stated absence.

Why bounded and not open brainstorming. The 500-entry candidate program was open brainstorming by a
weaker model, and the measured base rate was poor: 40 of 500 passed a strict consumability test, and 4
of the first 10 admitted turned out to rest on a false premise. All 500 carried byte-identical template
fields. The same capacity does produce real results when the question is bounded and carries an external
anchor — that is where the night's refutations came from.

So these jobs each start from something already measured and ask one question whose answer is checkable:

1. **Open edges of the tissue constraint net.** Eight edges carry a gap as a number. For each: name the
   measurement that would close it, and establish whether it exists.
2. **Deferred cells.** Five cells failed a sufficiency test and the quantity they need is one they do not
   carry. For each: is that quantity published anywhere?
3. **Unconsumed references.** Our own lanes cited 135 sources and nothing harvested 131 of them.

An outcome of "verified absence, here is what I searched" is as useful as a value: it says to change the
readout rather than wait for data. An improvement proposal with no anchor is not.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('')
QUEUE = W / 'tasks/lanes/bt_queue.txt'
RESULTS = W / 'results'
NET = W / 'CONSTRAINT_NETS.json'
REFS = W / 'notes/UNCONSUMED_REFERENCES.json'

# The five cells whose needed quantity they do not carry, from the repair analysis.
DEFERRED = {
    'Q052': ('pore diameter', 'the co-factor the permeability readout needs'),
    'Q054': ('contact area', 'the co-factor between prescribed load and local contact pressure'),
    'Q077': ('graph distance', 'the locator the transport-time readout needs, a topology not a volume'),
    'Q080': ('peak density', 'the weight that turns a population mean into a spatial peak'),
    'Q044': ('two angles', 'the directional coordinates a magnitude cannot carry'),
}

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
TABLE = (
    '| quantity | value | unit | where in the source | validity | n |\n'
    '|---|---|---|---|---|---|\n'
    '| what was measured | the number alone | the unit alone | table or figure | '
    'species/tissue/temp/method | sample size |\n'
)
RULES = (
    '\n## Rules\n'
    '1. **A stated absence is a result.** If no measurement exists, answer `VERIFIED_ABSENCE_OF_EVIDENCE`'
    ' and list what you searched. That tells us to change the readout instead of waiting.\n'
    '2. **Every value carries its validity range** — species, tissue, temperature, method, sample size. '
    'A value used outside its range is the commonest reason a record turns out unusable.\n'
    '3. **Say whether the number is measured or derived.** A value the authors computed from a model is '
    'not a measurement and must not arrive labelled as one.\n'
    '4. **Do not propose an improvement without an anchor.** One number with a locator beats five ideas.\n'
    '\nStatus `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation. Write `RESULTS.md` '
    'starting with the job id, plus `results.json`.\n'
)


def write_job(jid: str, text: str, meta: dict, web: bool = True) -> bool:
    if BLOCK.search(text):
        return False
    d = RESULTS / jid
    if (d / 'RESULTS.md').exists():
        return False
    d.mkdir(parents=True, exist_ok=True)
    (d / 'BRIEF.md').write_text(text)
    if web:
        (d / 'ALLOW_WEB').write_text('1\n')
    json.dump({'id': jid, 'review_state': 'PENDING_INDEPENDENT_REVIEW', **meta},
              (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
    return True


def main() -> int:
    made = []

    # 1. Open edges
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    for e in net['edges']:
        if e['status'] not in ('OPEN', 'UNKNOWN'):
            continue
        jid = 'BT-EDGE-' + re.sub(r'^T-E\d+-', '', e['id'])[:44].upper().replace('_', '-')
        text = (
            f"# Close one open edge of the tissue net\n\n"
            f"**Edge.** `{e['id']}` — status **{e['status']}**\n\n"
            f"**What is known.** {e['constraint']}\n\n"
            f"**Evidence behind that.** {e['evidence']}\n\n"
            f"The gap is already quantified, so this job is not to re-measure it. It is to answer one "
            f"question: **what single measurement would move this edge to TIGHT, and does it exist?**\n\n"
            f"Give the quantity, its unit, and what the measurement would decide. Then search for it and "
            f"report with the table below.\n\n" + TABLE + RULES)
        if write_job(jid, text, {'kind': 'edge_closure', 'edge': e['id'], 'status': e['status'],
                                 'category': 'calibration', 'claim_type': 'information_link'}):
            made.append(jid)

    # 2. Deferred cells
    for cell, (quantity, why) in DEFERRED.items():
        jid = f'BT-NEED-{cell}-{quantity.split()[0].upper()}'
        text = (
            f"# Find one quantity a cell needs and does not have\n\n"
            f"**Cell.** {cell}\n**Quantity.** {quantity}\n**Why it is needed.** {why}\n\n"
            f"This cell fails a sufficiency test: two states identical in the summary its readout uses "
            f"give different results. The missing coordinate has been named, and unlike five other cells "
            f"this one does not already carry it. So the question is whether it is published.\n\n"
            f"Search for it and report with the table below. If it is not published for a relevant "
            f"tissue, say so with what you searched — that is the result.\n\n" + TABLE + RULES)
        if write_job(jid, text, {'kind': 'needed_quantity', 'cell': cell, 'quantity': quantity,
                                 'category': 'calibration', 'claim_type': 'information_link'}):
            made.append(jid)

    # 3. The remaining unconsumed references, beyond the first 15 already queued
    refs = json.loads(REFS.read_text())['not_consumed'][15:55]
    for it in refs:
        doi = it['doi']
        slug = re.sub(r'[^A-Z0-9]+', '-', doi.upper())[:46].strip('-')
        jid = f'BT-FETCH-{slug}'
        text = (
            f"# Extract the values from one named reference\n\n"
            f"**Reference.** `{doi}`\n\n"
            f"Cited **{it['mentions']} times** in our own lane files while never harvested. The source is "
            f"given, so this is extraction and not a search.\n\n"
            f"Report with the table below. Also say whether this is a paper or a data archive: if raw "
            f"data is downloadable, give the access route, licence and size, because a dataset is worth "
            f"more than a value. And if the reference does not support what our files cite it for, say so "
            f"with the quote — that is more valuable than a table.\n\n" + TABLE + RULES)
        if write_job(jid, text, {'kind': 'reference_fetch', 'doi': doi,
                                 'mentions_in_our_files': it['mentions'],
                                 'category': 'calibration', 'claim_type': 'information_link'}):
            made.append(jid)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with QUEUE.open('a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    kinds = {}
    for m in made:
        k = m.split('-')[1]
        kinds[k] = kinds.get(k, 0) + 1
    print(f'queued {len(made)} bounded jobs: {kinds}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
