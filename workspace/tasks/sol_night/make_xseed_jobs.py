"Make Anton's X bookmarksseeds (BodyTwin-delen, ~/research/X_BOOKMARKS_20260930) to self-supporting swarming jobs.\n\nOne job per selected seed: results/BT-XSEED-<batch><B>/ with BRIEF.md, inputs/SEED.md (seed's section verbatim), JOB.json, ALLOW_WEB.\nThe runs are placed FIRST in tasks/lanes/bt_queue.txt (Anton's seeds before the backlog), under the queue lock. Idempotent."
import fcntl, hashlib, json, re, time
from pathlib import Path

def _private_input_pattern(public_pattern):
    """Combine public exclusions with required operator-maintained private terms."""
    import json
    from pathlib import Path
    import re
    config = Path.home() / ".bodytwin" / "private_source_terms.json"
    terms = json.loads(config.read_text())
    if not isinstance(terms, list) or not terms or any(not isinstance(t, str) or not t.strip() for t in terms):
        raise ValueError("A nonempty private-source exclusion list is required")
    return "(?:" + public_pattern + ")|(?:" + "|".join(re.escape(t) for t in terms) + ")"


X = Path('~/research/X_BOOKMARKS_20260930'); B = Path('.')
Q = B / 'tasks/lanes/bt_queue.txt'
import os
PICK = [tuple((int(x[0]), x[1:])) for x in os.environ['XSEED_PICK'].split(',')] if os.environ.get('XSEED_PICK') else [(4, 'B1'), (4, 'B2'), (3, 'B1'), (2, 'B3'), (2, 'B4'), (3, 'B4'), (3, 'B5'), (3, 'B6'), (3, 'B7'),
        (4, 'B6'), (4, 'B7'), (2, 'B5'), (4, 'B4')]
PROFILES = 'ABCD'

RULES = "\n## How you work (swarm job, The swarm/swarm_worker)\n\nThis is a seed from Anton's bookmarks, expanded and fact-checked by a stronger model. Seeds are ideas with a high base frequency of value, not ordered conclusions. Break it down to its smallest components and make the smallest crucial calculation first.\n\n1. Write WORK_STATUS.md directly with the next concrete calculation.\n2. PREREG.md before the first run: the decisive number, threshold, strongest equally informed control and what rejects the hypothesis. Hash it (PREREG.sha256).\n3. Build executable code (numpy/scipy) for the seed's core mechanism from first principles, with unit checking and at least one analytic boundary case in a test.\n4. Public data may be retrieved if it is small (< 200 MB) and open; cite DOI/URL and table/figure. Never make up measurement data; mark VERIFIED if a number cannot be looked up.\n5. RESULTS.md starting with the job's ID. results.json with the fields decision, actual_gain (against the control, with numbers), matched_accuracy, validity, missing_prerequisite, next_test, negative_result (boolean).\n6. FOLLOWUPS.json: 1–4 concrete next experiment if the result justifies it (decision, changed_input_or_operator, discriminating_test, strongest_comparator) — and EACH post should have its OWN external_referent (same scheme as above), neither inherited nor omitted. Measured 1/10: 45 of 63 new jobs in two hours declared nothing, so the child jobs lost the reference parent had. Declare synthetic_only or our_own_fixture honestly when the follow-up lacks an external reference.\n\n**Form (Anton 1/10, mandatory):** aim for A) a guarantee an equally informed control cannot give, B) a sharp statement where the budget is spent on bringing it down (report tried instances and counterexamples), or C) a measurement against a reference that does NOT come from us (public dataset, independent measurement, someone else's code, closed expression from a piece of paper). Write in the results.json field external_referent: {kind, locator (arXiv/DOI/PMID/URL), compared_quantity (exakt tal/utsaga), refutes_us}. If the seed has no external reference: declare our_own_fixture honestly. Beating our own control on our own lineup is not a goal in itself.\n\nLimits: work only in this directory; a thread; < 1 GB RAM; no internal data; no clinical recommendations. Everything is PENDING_INDEPENDENT_REVIEW. An honest negative finding that locates the obstacle is progress.\n"

def section(batch, tag):
    lines = (X / f'batch{batch}/SEEDS_BODYTWIN.md').read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if re.match(rf'## {tag}[. (]', l))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith('## ')), len(lines))
    return lines[start].lstrip('# ').strip(), '\n'.join(lines[start:end]).strip() + '\n'

def main():
    made, lines = [], []
    for n, (batch, tag) in enumerate(PICK):
        title, text = section(batch, tag)
        if re.search(_private_input_pattern('\\b(the collaborator|Grand Challenge|AnyBody)\\b'), text):
            print("skips (internal data mentioned):", batch, tag); continue
        jid = f'BT-XSEED-{batch}{tag}'; d = B / 'results' / jid
        if (d / 'RESULTS.md').exists(): continue
        (d / 'inputs').mkdir(parents=True, exist_ok=True)
        (d / 'inputs/SEED.md').write_text(text)
        (d / 'BRIEF.md').write_text(f'# {jid} — {title}\n\nRead inputs/SEED.md (Seed in full: ability, degradation, surgery, control and frozen targets). Perform the first executable step of the seed and go as far as the budget is sufficient.\n' + RULES)
        (d / 'ALLOW_WEB').write_text('public literature and small open datasets only\n')
        (d / 'JOB.json').write_text(json.dumps({'id': jid, 'wave': 'BT-XSEED', 'owner': 'BodyTwin coordinator anton-5f',
            'source': f'X_BOOKMARKS_20260930/batch{batch}/SEEDS_BODYTWIN.md {tag}', 'decision': title,
            'seed_sha256': hashlib.sha256(text.encode()).hexdigest(), 'created': time.time()}, ensure_ascii=False, indent=1))
        made.append(jid); lines.append(f'{PROFILES[n % 4]} swarm {jid}')
    with (B / 'tasks/lanes/bt_queue.lock').open('a') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        old = Q.read_text().splitlines()
        have = {l.split()[-1] for l in old if l.strip() and not l.startswith('#')}
        new = [l for l in lines if l.split()[-1] not in have]
        head = [l for l in old[:5] if l.startswith('#')]
        rest = old[len(head):]
        tmp = Q.with_name(Q.name + '.tmp'); tmp.write_text('\n'.join(head + new + rest) + '\n'); tmp.replace(Q)
    print(json.dumps({'made': made, 'queued_first': len(new)}))

if __name__ == '__main__':
    main()
