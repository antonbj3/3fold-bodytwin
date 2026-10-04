"""Queues targeted The swarm planners for source families that have not yet received jobs (the controller’s own planners start only when the queue is nearly empty).
Uses the free_controller module’s planner()/enqueue() under controller.lock. Usage: python3 seed_new_families.py [--dry]"""
import collections, fcntl, importlib.util, json, sys, time
from pathlib import Path

B = Path('')
spec = importlib.util.spec_from_file_location('fc', B / 'tasks/free_controller.py'); argv = sys.argv; sys.argv = ['x']
fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc); sys.argv = argv; c = fc.c
DRY = '--dry' in sys.argv

with (c.A / 'controller.lock').open('a') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    c.CAT = json.loads((c.A / 'CATALOG.json').read_text())
    state = json.loads(c.STATE.read_text())
    used = collections.Counter(k for r in state['jobs'].values() for k in (r.get('source_keys') or []) if not r.get('planner') and time.time()-r.get('created',0)<6*3600)
    unused = [k for k in c.CAT if used[k] == 0]
    rest = [k for k in unused if not k.startswith('SURG_')]
    groups = [[k for k in unused if k.startswith('SURG_')]] + [rest[i::3] for i in range(3)]
    out = []
    for g in groups:
        if not g: continue
        j = fc.planner(); j['id'] = c.PREFIX + '-PLAN-' + str(time.time_ns())
        j['body'] += ('\nPriority source partition (families with no jobs in the last 6 hours): ' + ', '.join(g) +
                      '. Propose runnable jobs mainly for these, 1-3 source_keys each, target_id = the catalog target of one of the keys. '
                      'At least one third should connect two of these families through a dimensionally defined port. '
                      + ('Surgical chain (Anton 30/9): increase resolution where a scalpel opens tissue - materials, cells, tissues, bindings; '
                         'SURG_INCISION / SURG_COLLAGEN / SURG_HEMOSTASIS / SURG_HEALING share ports: cut geometry, damage zone, vessel map, time. '
                         if g is groups[0] else '') +
                      'Check EXISTING_DECISIONS to avoid duplicates.')
        out.append({'id': j['id'], 'keys': g})
        if not DRY and c.enqueue(j, state, planner=True):
            state['jobs'][j['id']]['value_planner'] = True
    if not DRY:
        state['last_plan'] = time.time(); c.save(c.STATE, state)
print(json.dumps({'dry': DRY, 'unused_keys': unused, 'planners': out}, ensure_ascii=False))
