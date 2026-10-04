from common import *
from contact import compare, gates, Q
import csv, gzip, collections, resource

def run():
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R1.json')['metrics']
    rows = list(csv.DictReader(gzip.open(V4 / 'raw/QUALITY_ROWS.csv.gz', 'rt')))
    old = {(r['task_id'], r['participant']): r for r in rows if r['split'] == 'test' and r['participant'] in FRONTIER}
    cases = sorted({r['case_key'] for r in old.values()})
    out = []
    checks = []
    for key in cases:
        tasks = read(V4 / 'payload/public/tasks' / f'{key}.json')
        ref = npz(V4 / 'payload/private/references' / f'{key}.npz')
        pred = {m: npz(V4 / 'payload/predictions' / m / f'{key}.npz') for m in FRONTIER}
        cache = {}
        for (i, t) in enumerate(tasks):
            for m in FRONTIER:
                r = old[t['task_id'], m]
                row = dict(case_key=key, task_id=t['task_id'], dataset=r['dataset'], family=t['family'], level=t['level'], participant=m, L1=r['L1'], resolution='PER_TOOTH')
                if f'outer_{i}' not in pred[m]:
                    row.update(status='NO_PREDICTION', L2='FAIL' if r['L1'] == 'FAIL' else 'UNKNOWN')
                else:
                    if t['family'] not in cache:
                        cache[t['family']] = npz(V4 / 'payload/public' / t['geometry_file'])
                    a = cache[t['family']]
                    pg = a['ceiling'] - pred[m][f'outer_{i}']
                    rg = a['ceiling'] - ref[t['family']]
                    z = compare(a['xy'], a['faces'], pg, rg)
                    gg = gates(z, pr)
                    row.update(status=z['status'], contact=z, gates=gg, L2='PASS' if r['L1'] == 'PASS' and all(gg.values()) else 'UNKNOWN' if z['status'] != 'SCORED' else 'FAIL')
                    if i == 0 and m == 'field_generator':
                        control = Q.contact_map(a['xy'], a['faces'], pg, rg)
                        checks.append(dict(case=key, error_mm2=abs(control['contact_symdiff_mm2'] - z['symdiff_mm2']) if control['contact_symdiff_mm2'] is not None else None))
                out.append(row)
        if len(out) % 576 == 0:
            print('R1', len(out), 'rows', flush=True)
    by = collections.defaultdict(list)
    for r in out:
        by[r['task_id']].append(r)
    items = []
    for (key, rr) in by.items():
        items.append(dict(task_id=key, family=rr[0]['family'], case_key=rr[0]['case_key'], dataset=rr[0]['dataset'], L1_discriminates=len({r['L1'] == 'PASS' for r in rr}) > 1, L2_discriminates=len({r['L2'] == 'PASS' for r in rr}) > 1, L2_statuses={r['participant']: r['L2'] for r in rr}, all_four_scored=all((r['status'] == 'SCORED' for r in rr))))
    summary = dict(requested_items=len(items), discriminating_L1=sum((x['L1_discriminates'] for x in items)), discriminating_L2=sum((x['L2_discriminates'] for x in items)), all_four_scored=sum((x['all_four_scored'] for x in items)), participants={m: dict(collections.Counter((r['L2'] for r in out if r['participant'] == m))) for m in FRONTIER}, families={f: dict(items=sum((r['family'] == f for r in items)), discriminating=sum((r['family'] == f and r['L2_discriminates'] for r in items))) for f in sorted({r['family'] for r in items})}, dropout=dict(requested_rows=len(out), not_scored=sum((r['status'] != 'SCORED' for r in out)), reasons=dict(collections.Counter((r['status'] for r in out if r['status'] != 'SCORED')))), seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/R1_ROWS.json', out)
    dump(ROOT / 'raw/R1_ITEMS.json', items)
    dump(ROOT / 'raw/R1_CONTROLS.json', checks)
    dump(ROOT / 'rounds/R1.json', summary)
    state('R1_DECIDED', f"{summary['discriminating_L2']}/{len(items)} spatially discriminating items", 'Replace nominal-only spatial answer by continuous pose sets on complete crowns')
    (ROOT / 'HANDOFF.md').write_text('R1 complete. Frozen nominal spatial-contact gates distinguish ' + str(summary['discriminating_L2']) + '/' + str(len(items)) + ' paired items, versus ' + str(summary['discriminating_L1']) + ' level-1 differences. These are represented roofs, not loaded complete crowns. Next construction: complete crown surfaces, same spatial referent, continuous pose uncertainty; then consume contact incidence in TANDLAST.\n')
    return summary
if __name__ == '__main__':
    print(clean(run()))
