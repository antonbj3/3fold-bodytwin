"""Reuse the tested graph/queue controller for BodyTwin's own result-driven planning."""
from pathlib import Path
import importlib.util, json, sys, time
source=Path('local_path')
sys.path.insert(0,str(source.parent))
spec=importlib.util.spec_from_file_location('shared_controller',source);c=importlib.util.module_from_spec(spec);sys.modules['shared_controller']=c;spec.loader.exec_module(c)
c.D=Path('');c.B=c.D;c.A=c.D/'tasks/free48';c.RESULTS=c.D/'results'
c.STORE=Path('external_mount');c.CAT=json.loads((c.A/'CATALOG.json').read_text());c.STATE=c.A/'STATE.json';c.PREFIX='BT-FW48'
_accepted=c.accepted
def accepted_legacy_target(j):
    # Graph lane 2026-09-30: jobs briefed before the source-family target switch still name the old
    # catch-all BT-C4-ERROR-BUDGET; bind them to their first source family's catalog target instead.
    if isinstance(j,dict) and j.get('target_id')=='BT-C4-ERROR-BUDGET':
        keys=[k for k in j.get('source_keys') or [] if k in c.CAT]
        if keys and 'BT-C4-ERROR-BUDGET' not in {c.CAT[k]['target_id'] for k in keys}:
            j['legacy_target_id']=j['target_id'];j['target_id']=c.CAT[keys[0]]['target_id']
    return _accepted(j)
c.accepted=accepted_legacy_target
def planner():
    return {'id':'BT-FW48-PLAN-'+str(int(time.time())),'source_keys':list(c.CAT),'target_id':c.CAT['MITOSTRESS']['target_id'],'category':'exploratory',
    'decision':'Choose the next BodyTwin mechanisms and cross-scale couplings that change a defined observable at matched error and cost.',
    'baseline':'Existing supplied model families, completed followups and preserved negative controls.',
    'novelty_delta':'Expand actual completed findings into testable new operators, reusable interfaces and independent checks.',
    'body':'''BEFORE you suggest something: write three different ways of looking at what is worth doing next, where at least one does not involve tissue, cells or the human body at all. Diverge first, then choose. Write the three in RESULTS.md so they can be read afterwards, and say which one you chose and why. This applies CHOICE of direction, not the form of the answer: the requirements below remain unchanged for each proposal you submit.\n\nMeasured 2026-10-03 by anton-5f, therefore this step exists: 69 bound jobs with template form gave 1 DOI in 175 table rows and zero quantities any cell consumes, while six briefs without forced answer form and no given interpretation gave 8-24 DOI with PMID each, and the first to be consumed set a floor throughout the decision chain. The difference was that no one had decided what the answer would be about. At the same time, the dental lane 12 of 30 measured COUNTING jobs that failed on stringency, wrong tail and a universal quantifier out of a grid, so stringency requirements are not removed anywhere.\n\nRead inputs/SOURCE_CATALOG.json, RECENT_RESULTS.json and EXISTING_DECISIONS.json first. Read at most six further source files. Write SWARM_QUEUE_ADD.json early as a list of 6-12 ranked distinct proposals, each with source_keys (1-3 catalog keys), target_id, category, decision, baseline, novelty_delta, body. Category: mechanism/adversarial/audit/integration/calibration/exploratory. Aim 80% constructive mechanisms, implementations, integrations or capability-unlocking calibration, 15% decisive physical or mathematical challenges, 5% blocking evidence checks. Expand from first principles: state/control volume -> conservation and constitutive law -> observable -> representation -> assumption -> decisive falsifier. Connect cells, tissue geometry, material response, barrier/transport, organism exposure and measurement only through dimensionally defined ports. At least one third should connect two supplied mechanisms. Preserve negative evidence. No patient or licensed anatomy inputs are available. No broad literature-only jobs or duplicate identities. A hypothesis is useful only with a new discriminating test and a downstream model decision. These are definition/review tasks until physical inputs and native-code bindings are known; conditional synthetic experiments are permitted but are not validated. Keep full private-code coverage and graph gaps explicit. Return fewer proposals if novelty or sources are missing. RESULTS.md explains what each family could unlock and what is not yet identifiable. FORM (Anton 2026-10-01, mandatory): every proposal targets one of A) a guarantee the strongest equally informed control cannot give at any budget, B) a sharp claim with the budget spent trying to refute it (report instances tried and counterexamples found), or C) a measurement against a reference that does NOT come from this project. 'Beat an equally informed control on a setup we wrote ourselves' is forbidden as the only goal. Each proposal MUST includes an external_referent object: {\"kind\": independent_measurement|published_dataset|published_code|closed_form|external_review|our_own_fixture|synthetic_only, \"locator\": arXiv id, DOI, PMID, URL or existing path, \"compared_quantity\": the exact number or statement taken from the referent, \"refutes_us\": true if the referent could show our current answer is wrong}. Declare our_own_fixture honestly when there is no outside reference; such proposals rank lower. Prefer C with independent measurements. The referent is a held-out data anchor: the model must compute that quantity from lower-level mechanisms (molecules, cells, tissue physics), never take it or a textbook constant as an input; proposals that fit or insert the referent value are invalid.'''}
c.planner_job=planner

# 30/9 23:40 (anton-5f): FOLLOWUPS.json from workers is often rich but has another schema (decision, changed_input_or_operator,
# discriminating_test, strongest_comparator ...) without source_keys/target_id/category, and older proposals point at
# BT-C4-ERROR-BUDGET, which the families no longer have as a target. Both fell out as "invalid metadata" (128 of 164 in the last hour).
# Normalize in place before accepted(): inherit keys from the parent job, map the target to the key's current target_id.
_sn_orig_accepted, _sn_orig_ingest, _origin = c.accepted, c.ingest, {}
# 2/10 15:50 (anton-5f): the requirement text is a module constant because it now has a SECOND
# consumer. _normalize only reaches jobs accepted after 05:11; measured on 400 queued BT-FW jobs,
# 16/16 accepted after that carry the requirement and 33/384 accepted before it do. The backlog is
# therefore ~96 % jobs that cannot produce an external-reference result, and
# tasks/build_night/backfill_referent.py appends this same block to the queued, unstarted ones.
# One literal, two callers: if it lived in both places they would drift and the backfilled jobs
# would be asked for a different field than the live ones.
_SN_REFERENT_BLOCK = ('\n\nEXTERNAL REFERENCE (required in results.json): external_referent ='
    '{kind: independent_measurement|published_dataset|published_code|closed_form|external_review|'
    'our_own_fixture|synthetic_only, locator: DOI/PMID/arXiv/URL/path, compared_quantity: det exakta '
    'the number or statement from the call with unit and regime, our_value: the number we calculated, command or '
    'excerpt: how it is shaved, refutes_us: true if the call can prove us wrong}. The parent declared '
    'nothing, you must declare your own; don\'t inherit it. Missing external reference: declare our_own_fixture or '
    'synthetic_only ally. The anchor is the hall data and must never be input to the model.')
def _normalize(j):
    try:
        if not isinstance(j, dict): return
        parent = c.STATE_CACHE.get(_origin.get('id'), {}) if hasattr(c, 'STATE_CACHE') else {}
        if not j.get('source_keys') and parent.get('source_keys'):
            j['source_keys'] = [k for k in parent['source_keys'] if k in c.CAT][:3]
        keys = [k for k in (j.get('source_keys') or []) if k in c.CAT]
        if keys and j.get('target_id') not in {c.CAT[k]['target_id'] for k in keys}:
            if j.get('target_id'): j['legacy_target_id'] = j['target_id']
            j['target_id'] = c.CAT[keys[0]]['target_id']
        if j.get('category') not in ['mechanism','adversarial','audit','integration','calibration','exploratory']:
            d = str(j.get('decision', '')).lower()
            j['category'] = 'audit' if 'audit' in d else 'adversarial' if any(w in d for w in ('falsif', 'attack', 'refute')) else 'mechanism'
        if not (isinstance(j.get('baseline'), str) and len(j['baseline']) >= 15):
            j['baseline'] = str(j.get('strongest_comparator') or 'Parent result and its strongest equally informed control')[:6000]
        if not (isinstance(j.get('novelty_delta'), str) and len(j['novelty_delta']) >= 15):
            j['novelty_delta'] = str(j.get('changed_input_or_operator') or j.get('why_it_matters') or j.get('unlocks') or '')[:6000]
        if not (isinstance(j.get('body'), str) and len(j['body']) >= 15):
            parts = [f'{k}: {j[k]}' for k in ('why_it_matters','changed_input_or_operator','discriminating_test','strongest_comparator',
                     'outcome_if_negative','cost','unlocks','addresses','parent_artifact') if j.get(k)]
            j['body'] = ('Followup proposed by the parent job (unaudited). Freeze PREREG before running; strongest equally informed '
                         'control; preserve negatives.\n' + '\n'.join(parts))[:6000]
        # 2/10 05:15 (anton-5f): measured that 42 of 69 new jobs lacked external_referent. The cause is that
        # research_value.DIRECTIVE only reaches PLANNER JOBS; ordinary followup jobs get the dental controller's
        # RULES, which do not mention the field. The requirement is therefore put in the job's OWN body here, in
        # BodyTwin's code path, and only when it is missing. No followup is rejected: an empty field already costs
        # zero in the scorer.
        if not isinstance(j.get('external_referent'), dict) and isinstance(j.get('body'), str):
            if 'external_referent' not in j['body']:
                j['body'] = (j['body'] + _SN_REFERENT_BLOCK)[:9000]
    except Exception as e:
        c.log(f'normalize error {type(e).__name__}: {e}')
_SN_FAMILY_CAP, _SN_FAMILY_WINDOW_S = 10, 7200
def _sn_family_saturated(j):
    # 1/10 07:00 (anton-5f): followups of one family bred more of the same (78/78 new jobs on Q080 in 2 h).
    # Hold proposals whose source key already got >= _SN_FAMILY_CAP jobs in the last _SN_FAMILY_WINDOW_S; planners are unaffected.
    try:
        now=time.time(); keys=set(j.get('source_keys') or [])
        if not keys or len(keys)>3: return False
        jobs=getattr(c,'STATE_CACHE',{}) or {}
        recent=[r for r in jobs.values() if not r.get('planner') and now-r.get('created',0)<_SN_FAMILY_WINDOW_S]
        return any(sum(k in (r.get('source_keys') or []) for r in recent)>=_SN_FAMILY_CAP for k in keys)
    except Exception as e:
        c.log(f'family cap error {type(e).__name__}: {e}'); return False
def _sn_accepted(j):
    _normalize(j)
    if isinstance(j,dict) and _sn_family_saturated(j): return False
    return _sn_orig_accepted(j)
def _sn_ingest(path, state):
    _origin['id'] = Path(path).parent.name; c.STATE_CACHE = state.get('jobs', {})
    return _sn_orig_ingest(path, state)
c.accepted, c.ingest = _sn_accepted, _sn_ingest

# 1/10 09:45 (anton-5f): main() schedules audits directly via enqueue() with the job's old target_id
# (BT-C4-ERROR-BUDGET); after the target switch no key matches -> StopIteration in enqueue -> every run
# since 03:44 was aborted (exit 1 / timeout). Translate old targets here too; a remaining error must not stop the run.
_sn_orig_enqueue = c.enqueue
def _sn_enqueue(j, state, planner=False):
    try:
        if isinstance(j, dict):
            keys = [k for k in (j.get('source_keys') or []) if k in c.CAT]
            if keys and j.get('target_id') not in {c.CAT[k]['target_id'] for k in keys}:
                j['legacy_target_id'] = j.get('target_id'); j['target_id'] = c.CAT[keys[0]]['target_id']
        return _sn_orig_enqueue(j, state, planner=planner)
    except StopIteration:
        c.log(f"enqueue skipped {j.get('id') if isinstance(j, dict) else j}: no catalog key matches target")
        return False
c.enqueue = _sn_enqueue

# 1/10 09:55 (anton-5f): chain_backlog() reads the whole queue_ovh.log and stats ~6000 RESULTS.md per call and
# is called once per audit candidate in main(); once the crash above was gone, the run hit the 240 s timeout.
# Cache the value for 60 s within a run (enqueue only adds a few jobs, so the error is at most a handful).
_sn_orig_chain_backlog, _sn_backlog_cache = c.chain_backlog, {}
def _sn_chain_backlog(state):
    now = time.time(); hit = _sn_backlog_cache.get('v')
    if hit and now - hit[0] < 60: return hit[1]
    v = _sn_orig_chain_backlog(state); _sn_backlog_cache['v'] = (now, v); return v
c.chain_backlog = _sn_chain_backlog

if __name__=='__main__':
    if sys.argv[1:]==['--value-refresh']:
        import fcntl
        with (c.A/'controller.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if time.time()<c.UNTIL and not (c.A/'STOP').exists():
                state=json.loads(c.STATE.read_text())
                pending=[j for j,r in state['jobs'].items() if r.get('value_planner') and not (c.RESULTS/j/'RESULTS.md').exists() and time.time()-r['created']<5400]
                if not pending:
                    j=planner();j['id']=c.PREFIX+'-PLAN-'+str(time.time_ns())
                    if c.enqueue(j,state,planner=True):state['jobs'][j['id']]['value_planner']=True;state['last_plan']=time.time();c.save(c.STATE,state)
    else:c.main()
