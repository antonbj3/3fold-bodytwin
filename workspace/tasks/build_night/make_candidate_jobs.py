"""Turn CANDIDATES.json (Sonnet brainstorm, reviewed by the coordinator) into swarm jobs results/<PREFIX>-<key>/.

Usage: python3 make_candidate_jobs.py <CANDIDATES.json> <PREFIX> [key,key,...]
Each candidate: key, title, form, decision, brief, external_referent, builds_on. Idempotent; the run lines are added to
tasks/lanes/bt_queue.txt under the queue's lock (the order is set by research_value.order, not by the file position).
"""
import fcntl, hashlib, json, re, sys, time
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


B = Path(''); Q = B / 'tasks/lanes/bt_queue.txt'
sys.path.insert(0, str(Path(__file__).parent))
from make_xseed_jobs import RULES  # same work rules and A/B/C form as XSEED

HIDDEN = """
**Hidden variables: vibrate them, deliver a cert (Anton 1/10).** A biological anchor is a projection of a state we don't see, and its SPRIDNING is the fingerprint of the hidden variables, not measurement noise. Hitting the mean can be wrong about each individual. Therefore, do not describe the hidden variables in text — find them with the graph engine tool and show if the decision holds when they are vibrated.

The tools are there, don't build them yourself (path: local_path):
- `hidden_variable.from_mechanism(model, x_range, params)` — MEKANISMSIDAN, before each report. Sweep each parameter over its published range and measure how much of (x, θ) TECKNET at dy/dx differs from the mean. A parameter that can reverse the relation within its own ranges is a hidden axis that the anchor's validity box must carry. Also returns the threshold at which the reversal first occurs.
- `hidden_variable.candidates(reports, min_side=2, min_known=0.7, n_perm=200)` — DATASIDAN. **Two thresholds that silence away (the graph 1/10):** it returns EMPTY LIST when n < 2*min_side, so **four** reports are required, not two, and the attribute must be present on at least 70% of them. Empty list is therefore NOT absence of hidden variable — check n and attribute coverage before interpreting it. Input: [{"margin": float, "sigma": float, "attributes": {...}}] or "sign": ±1 instead of margin. If two reports on the same thing disagree even though their validity boxes overlap, read that as an undeclared variable, not a contradiction. Rank candidate splits on the chi2 case, with permutation tests. One variable that BOTH sides point out is the throw.
- `claim_federation` is the ENTRY for pitting conflicting studies against each other (it traces each source to its root sources and counts N_eff, so four articles citing the same original count as ONE, and it distinguishes REGIME-BOUNDARY from CONTRADICTION). `disagreement_field` is downstream of hidden_variable, not as an alternative. `decision_cert.certify(decision, state)` — deliver the decision as a cert with a margin, not as a number. `flip_attribution` tells what quantity would reverse the decision, and `margin_net` gives the margin and alerts when it is negative.

Delivery requirements:
1. Run the mechanism side of your own model and report the parameters that can reverse sign within their ranges, with the threshold. Zero turns is also a result: then the anchor's box is sufficient for the sign question.
2. Try against the anchor's SPRIDNING, not just the average. If the model can produce the reported the dispersion by vibrating the designated the axes within published intervals, write that as its own result — it's stronger than the average hit. If it only passes the average: mark POPULATION_FITTED and say resolution not reached.
3. Deliver the decision as cert with a margin and with the quantity that would reverse it. Two mechanisms giving the same mean are not separated by the armature; name the one that would give the same number for other reasons and what part of the dispersion separates them.

An anchor without declared dispersion is weaker, not stronger: then write that the dispersion is unknown and that the average hit cannot be graded.
"""


TIME_WARNING = '''
**Time limit: ** the run can be interrupted after about 90 minutes (hard wall clock). Three of the first jobs in this wave died on 124 without having time to write RESULTS.md. Therefore, work interrupt-proof: write WORK_STATUS.md first, freeze PREREG.md early, and write a PARTIELL RESULTS.md with what you have as soon as the first digit is found — then update it. An aborted run that leaves code, frozen sources, and a partial report can be resumed; one who leaves only one log must start over.
'''

PATH_WARNING = '''
**Paths:** the job can run on a cloud host where local_path NOT exists. Everything you need must be in this job directory (inputs/) or obtained from public sources. Pointing the brief to a local path that doesn't exist: book it as missing_prerequisite and build from first principles instead. Never copy internal data here.
'''

DELIVERY = '''
**Delivery of reference (applies from 1/10): ** in results.json external_referent must contain kind, locator, compared_quantity (the number from the source with unit, species, regime), our_value (the number we calculated), command or excerpt (how it is recalculated, or the quote line with the value) and refutes_us. Reference is a hold-out data anchor: it must never be input to the model. List the lower level inputs and confirm that the anchor value is not among them. Do you not find reference: found: false with searched: [search phrases].
'''

def main():
    src, prefix = Path(sys.argv[1]), sys.argv[2]
    pick = set(sys.argv[3].split(',')) if len(sys.argv) > 3 else None
    made, lines = [], []
    for n, c in enumerate(json.loads(src.read_text())):
        if pick and c['key'] not in pick: continue
        text = json.dumps(c, ensure_ascii=False)
        if re.search(_private_input_pattern('\\b(the collaborator|Grand Challenge|AnyBody)\\b'), text):
            print('skips (internal data mentioned):', c['key']); continue
        jid = f"{prefix}-{c['key']}"; d = B / 'results' / jid
        if (d / 'RESULTS.md').exists(): continue
        (d / 'inputs').mkdir(parents=True, exist_ok=True)
        (d / 'inputs/CANDIDATE.json').write_text(json.dumps(c, ensure_ascii=False, indent=1))
        # 1/10 (anton-5f): a job that references a local path can be run on a cloud host where
        # local_path does not exist, and then loses its whole basis. The candidate can therefore list
        # 'copy_inputs': the files are copied into inputs/ so the job always carries them along.
        import shutil
        # 2/10: the loop variable was named src and shadowed the CANDIDATES.json path bound in
        # main(), so every JOB.json written after the first copy recorded a copied input file as its
        # provenance instead of the candidates file. Renamed to extra_input.
        for extra_input in (c.get('copy_inputs') or []):
            sp = Path(extra_input)
            if sp.is_file():
                shutil.copy2(sp, d / 'inputs' / sp.name)
            else:
                print('copy_inputs missing, skipping:', extra_input)
        (d / 'BRIEF.md').write_text(f"# {jid} — {c['title']}\n\nForm {c['form']}. {c['decision']}\n\n{c['brief']}\n\nBased on: {c.get('builds_on','')}\n" + HIDDEN + TIME_WARNING + PATH_WARNING + DELIVERY + RULES)
        (d / 'ALLOW_WEB').write_text('public literature and small open datasets only\n')
        # 2/10 17:15 (anton-5f): target_id, source_keys and category are now passed through from the
        # candidate instead of being hardcoded. Measured cause: impact_value.py scores a job by the
        # completed-result count of its target_id and applies the PRIORITY_TARGETS.json bonus by the
        # same key, so every job this script wrote was invisible to both. They all received the median
        # node load and a zero bonus, which means the BodyTwin priorities I had just configured could
        # never reach my own hand-written jobs. category defaults to mechanism as before; a candidate
        # that declares itself adversarial now lands in the challenge pool where it belongs.
        record = {'id': jid, 'wave': prefix, 'owner': 'BodyTwin coordinator anton-5f',
            'source': str(src), 'decision': c['decision'],
            'category': c.get('category', 'mechanism'),
            'external_referent': c.get('external_referent'),
            'candidate_sha256': hashlib.sha256(text.encode()).hexdigest(), 'created': time.time()}
        if c.get('target_id'):
            record['target_id'] = c['target_id']
        if c.get('source_keys'):
            record['source_keys'] = list(c['source_keys'])[:3]
        # 2/10 17:25 (anton-5f): research_value, and ONLY when the candidate really supplies all seven
        # fields. Measured reason: research_value.assess derives a job's family by keyword-matching the
        # decision text, and 192 of 300 pending jobs match nothing and land in 'other'. The ordering
        # then applies a per-family diversity penalty, so that one bucket holds 64 % of the queue and
        # buries jobs regardless of score - BT-FW48-FIXX-OCT2-KCAT-ADJUDICATE scored 60.77 against
        # another job's 55.98 and still sat at position 1814. A job may declare mechanism_family and
        # escape the bucket, but assess() honours that only when the research_value block is complete,
        # which also earns +12. Nothing is padded to reach the 15-character minimum: a candidate that
        # does not state all seven gets no block, because inventing the fields to clear a threshold
        # would be gaming the scorer rather than informing it.
        rv_fields = ('capability', 'obstacle', 'changed_operation', 'consumer',
                     'metric', 'strongest_control', 'falsifier')
        rv = {k: c[k] for k in rv_fields if isinstance(c.get(k), str) and len(c[k].strip()) >= 15}
        if len(rv) == len(rv_fields):
            if isinstance(c.get('mechanism_family'), str):
                rv['mechanism_family'] = c['mechanism_family']
            record['research_value'] = rv
        elif rv:
            print(f'{jid}: research_value omitted, {len(rv)} of {len(rv_fields)} fields specified '
                  f'(saknas: {sorted(set(rv_fields)-set(rv))})')
        (d / 'JOB.json').write_text(json.dumps(record, ensure_ascii=False, indent=1))
        made.append(jid); lines.append(f"{'ABCD'[n % 4]} swarm {jid}")
    with (B / 'tasks/lanes/bt_queue.lock').open('a') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        old = Q.read_text().splitlines()
        have = {l.split()[-1] for l in old if l.strip() and not l.startswith('#')}
        new = [l for l in lines if l.split()[-1] not in have]
        head = [l for l in old[:5] if l.startswith('#')]
        tmp = Q.with_name(Q.name + '.tmp'); tmp.write_text('\n'.join(head + new + old[len(head):]) + '\n'); tmp.replace(Q)
    print(json.dumps({'made': made, 'queued': len(new)}))

if __name__ == '__main__':
    main()
