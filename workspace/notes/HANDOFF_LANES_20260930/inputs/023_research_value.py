"""Shared research scheduling advice; never scientific evidence admission.

Legacy classification uses the actual decision, not inherited BRIEF context or
category alone. Every queue row survives. Scores are uncalibrated heuristics.
"""
from collections import Counter
from pathlib import Path
import argparse
import json
import re
import time

ROOT = Path('.')
VERSION = 'research-value-20260930-v1'
VALUE_FIELDS = ('capability', 'obstacle', 'changed_operation', 'consumer',
                'metric', 'strongest_control', 'falsifier')
DIRECTIVE = '''Anton2026-09-30: drive ambitious research autonomously. Build the strongest
available construction toward a useful capability beyond current feasibility.
Connect supplied mechanisms through actual equations and compatible ports.
Execute a discriminating calculation early. At a tie or failure, diagnose the
load-bearing obstacle and change a consequential operation within this task's
scope and budget; carry the original goal forward. Compare the strongest equally
informed control at matched accuracy, including setup, update and checking cost.
A large gain is a search aim, never permission to invent a gain or evidence.
Prioritize coupled physical/design capabilities, representation or solver
changes, and calibration that unlocks a named computation. Routine verdict,
hash, preregistration, AST and reporting audits have a small supporting role.
Do not turn one missing input into successive renamed audits. Give its exact
measurement specification once and attack another available mechanism.
FOLLOWUPS/SWARM_QUEUE_ADD proposals must include research_value with capability,
obstacle, changed_operation, consumer, metric, strongest_control and falsifier
(concrete strings), plus mechanism_family and original source/parent IDs.
Prefer substantive cross-result constructions to more cataloguing. Preserve
negative results, graph lineage and scope. No scientific admission is implied.
'''
PLANNER_DIRECTIVE = '''
Planner allocation target: 80% constructive capability experiments/calibration,
15% decisive physical or mathematical challenges, 5% blocking evidence checks.
Routine report/label/hash audits must not refill themselves indefinitely.
Rank by capability unlocked, consequential uncertainty resolved, reusable
operator and complete cost. These are hypotheses, not calibrated expected values.
Read VALUE_BACKLOG.json and select diverse mechanism/consumer families. Reuse
existing equivalent packets before proposing successors. A missing prerequisite
can yield a conditional model test or one acquisition specification; never
claim physical validation. Explain why the next changed operation could beat
the strongest same-information control; if it cannot, change the construction.
'''

META = ('verdict string', 'reported string', 'vocabulary', 'bookkeeping',
        'sidecar', 'mtime', 'ctime', 'provenance gate', 'gate-independence',
        'static operator', 'dead declared constants', 'dead model constants',
        'corpus-level statistic', 'co-located hash', 'freeze-before',
        'k1=none', 'missing-regret', 'abstention channel', 'subsample identity')
FAMILIES = (
    ('history-coupling', ('history', 'binding', 'transport', 'perfusion', 'receptor', 'immune', 'mitochond', 'clearance', 'recirculat')),
    ('contact-interface', ('compliance', 'interface', 'contact', 'retention', 'mucosa')),
    ('geometry-stress', ('shell', 'crown', 'eigenstrain', 'stress', 'modulus', 'free surface', 'geometry', 'finite-element')),
    ('process-material', ('fatigue', 'surface-finishing', 'material', 'manufactur', 'sinter', 'thermal', 'heat', 'composition')),
    ('measurement-decision', ('sensor', 'observation', 'measurement', 'ambiguity', 'minimax', 'regret', 'identified set')),
    ('graph-computation', ('operator', 'coupling', 'inverse', 'graph', 'solver', 'certificate', 'representation')),
)


def assess(job, parent=None):
    decision = str(job.get('decision', ''))
    text = decision.lower()
    value = job.get('research_value', {})
    if not isinstance(value, dict): value = {}
    structured = all(isinstance(value.get(k), str) and len(value[k].strip()) >= 15 for k in VALUE_FIELDS)
    # Consumer words in copied context cannot turn a reporting audit into physics.
    meta = any(word in text for word in META) or bool(re.search(r'\bprereg(?:istration)?(?:\.md|\b)',text))
    physical = any(word in text for _, words in FAMILIES[:4] for word in words)
    category = job.get('category', '')
    if '-PLAN-' in str(job.get('id', '')): role = 'planner'
    elif meta: role = 'background'
    elif category == 'audit':
        role = 'review' if parent and assess(parent)['role'] in ('constructive', 'challenge') else 'background'
    elif category == 'adversarial': role = 'challenge'
    elif physical or structured or category in ('mechanism', 'calibration', 'integration'): role = 'constructive'
    else: role = 'challenge'
    family_text = str(parent.get('decision', '')).lower() if parent and category == 'audit' else text
    family = next((name for name, words in FAMILIES if any(word in family_text for word in words)), 'other')
    if meta: family = 'reporting-support'
    declared_family = job.get('mechanism_family', value.get('mechanism_family'))
    if structured and isinstance(declared_family, str): family = declared_family[:100]
    base = {'planner': 100, 'constructive': 60, 'challenge': 45, 'review': 35, 'background': 5}[role]
    keys = job.get('source_keys', [])
    score = base + 12*structured + 6*(isinstance(keys, list) and len(keys) >= 2)
    score += 5*bool(job.get('graph_candidate_ids')) + 4*bool(job.get('after'))
    if role == 'constructive': score += 5*physical
    return {'version': VERSION, 'role': role, 'score': score, 'family': family,
            'structured_value': structured, 'target_id': job.get('target_id'),
            'parent': job.get('after'), 'scheduling_advice_only': True,
            'reason': 'primary decision concerns reporting' if meta else 'capability, sources and mechanism diversity',
            'decision': decision[:600]}


def enrich(job, parent=None, planner=False):
    result = dict(job)
    result['value_advice'] = assess(result, parent)
    result['body'] = result.get('body', '') + '\n\n' + DIRECTIVE + (PLANNER_DIRECTIVE if planner else '')
    return result


def read_job(directory):
    try:
        p = directory/'JOB.json'
        if p.stat().st_size > 1024**2: return {}
        result = json.loads(p.read_text())
        return result if isinstance(result, dict) else {}
    except (OSError, ValueError): return {}


def order(lines, results=ROOT/'results'):
    """80/15/5 portfolio; diversified families/projects; preserve all rows."""
    pools = {r: [] for r in ('planner', 'constructive', 'challenge', 'review', 'background')}
    tail = []; advice = {}; cache = {}
    for index, line in enumerate(lines):
        parts = line.split()
        if len(parts) != 3 or not re.fullmatch(r'BT-[A-Za-z0-9_-]+', parts[2]): tail.append(line); continue
        ident = parts[2]; directory = results/ident
        if not directory.is_dir() or (directory/'RESULTS.md').exists() or (directory/'.ovh_claim').exists():
            tail.append(line); continue
        if (directory/'.local_claim').exists():
            try:
                pid = int((directory/'.local_claim').read_text()); live = (Path('/proc')/str(pid)).exists()
            except (OSError, ValueError): live = False
            if live: tail.append(line); continue
        job = cache.setdefault(ident, read_job(directory))
        job = dict(job, id=ident)
        parent_id = job.get('after')
        parent = None
        if isinstance(parent_id, str) and re.fullmatch(r'BT-[A-Za-z0-9_-]+', parent_id):
            if parent_id not in cache: cache[parent_id] = read_job(results/parent_id)
            parent = cache[parent_id]
        a = assess(job, parent); advice[ident] = a
        plan_stamp=ident.rsplit('-',1)[-1]
        created=int(plan_stamp)/1e9 if plan_stamp.isdigit() and int(plan_stamp)>10**12 else (int(plan_stamp) if plan_stamp.isdigit() else directory.stat().st_mtime)
        if a['role'] == 'planner' and time.time() - created > 5400:
            tail.append(line); continue
        project = 'bodytwin' if ident.startswith('BT-FW') else 'dental'
        seed = '-SEED-' in ident
        pools[a['role']].append((a['score']+3*seed, index, line, project, a['family'], parent_id))
    selected = []; projects = Counter(); families = Counter(); parents = Counter()
    def take(role):
        pool = pools[role]
        if not pool: return False
        best = max(range(len(pool)), key=lambda i: (pool[i][0]-8*families[(pool[i][3],pool[i][4])]
                     -6*projects[pool[i][3]]-12*parents[pool[i][5]]*(pool[i][5] is not None), -pool[i][1]))
        _, _, line, project, family, parent = pool.pop(best)
        selected.append(line); projects[project] += 1; families[(project,family)] += 1
        if parent: parents[parent] += 1
        return True
    for _ in range(2): take('planner')
    # Sixteen construction slots, three decisive challenges/reviews, one support
    # slot. Empty roles fall back to available work; nothing is dropped/starved.
    pattern = ['constructive']*4+['challenge']+['constructive']*4+['review']+['constructive']*4+['challenge']+['constructive']*4+['background']
    while any(pools.values()):
        for role in pattern:
            if not take(role):
                for alternative in ('constructive','challenge','review','planner','background'):
                    if take(alternative): break
    return selected+tail, advice


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--queue',type=Path,default=ROOT/'tasks/lanes/bt_queue.txt')
    parser.add_argument('--inspect',action='store_true');args=parser.parse_args()
    lines=args.queue.read_text().splitlines(keepends=True); ranked,advice=order(lines)
    if args.inspect:
        print(json.dumps({'version':VERSION,'time':time.time(),'rows_preserved':Counter(lines)==Counter(ranked),
              'role_counts':dict(Counter(a['role'] for a in advice.values())),
              'next':[{ 'job':l.split()[2],**advice[l.split()[2]]} for l in ranked if len(l.split())==3 and l.split()[2] in advice][:20]},indent=2))
    else:
        for line in ranked: print(line,end='')
