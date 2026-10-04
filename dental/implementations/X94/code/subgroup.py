import json, time, resource
from fractions import Fraction as Q
import numpy as np
from geometry import R, write, sha, load, FDI
from validate import metrics
REPORT = R.parent / 'LANE_X7_BITE2TEXT/raw/ALL_STRUCTURED_REPORTS.json'
FIELDS = {'molar_right': 'I', 'molar_left': 'I', 'canine_right': 'I', 'canine_left': 'I', 'overbite': 'normal', 'overjet': 'normal', 'crossbite': 'absent'}

def eligible(reports, strict=True):
    if not reports:
        return (False, ['no_report'])
    fields = FIELDS if strict else {k: v for (k, v) in FIELDS.items() if not k.startswith('canine')}
    bad = []
    for (i, r) in enumerate(reports):
        labels = r['labels']
        for (k, v) in fields.items():
            if labels.get(k) != v:
                bad.append(k + ':' + ('UNKNOWN' if labels.get(k) is None else str(labels[k])))
        if labels.get('missing_teeth') is None:
            bad.append('missing_teeth:UNKNOWN')
        elif labels['missing_teeth']:
            bad.append('missing_teeth:present')
    return (not bad, sorted(set(bad)))

def main():
    start = time.perf_counter()
    assert sha(R / 'PREREG_R4.json') == (R / 'PREREG_R4.json.sha256').read_text().split()[0]
    (manifest, cases, sources) = load()
    reports = json.load(open(REPORT))
    partition = manifest['inherited_partition']
    profiles = {}
    selection = {}
    rejected = {}
    for c in cases:
        e = c['edges']
        w = [c['area'][edge] for edge in e]
        total = sum(w, Q(0))
        p = {j: {t: Q(0) for t in FDI[j]} for j in FDI}
        if total:
            for (edge, weight) in zip(e, w):
                p['upper'][edge[0]] += weight / total
                p['lower'][edge[1]] += weight / total
            profiles[c['case']] = {j: {str(t): float(v) for (t, v) in pj.items()} for (j, pj) in p.items()}
        selection[c['case']] = {}
        for mode in ['strict', 'moderate']:
            (ok, bad) = eligible(reports.get(c['case'], []), mode == 'strict')
            selection[c['case']][mode] = ok
            if not ok:
                rejected.setdefault(mode, {}).setdefault(c['case'], bad)
    frozen = dict(prereg_sha256=sha(R / 'PREREG_R4.json'), report_source_sha256=sha(REPORT), profiles=profiles, selection=selection, source_status='Full16 proximity-area predictions and eligibility frozen before subgroup quantile computation', frozen_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat())
    p = R / 'FROZEN_PREDICTIONS_R4.json'
    if p.exists():
        old = json.load(open(p))
        assert old['profiles'] == profiles and old['selection'] == selection
    else:
        write('FROZEN_PREDICTIONS_R4.json', frozen)
        (R / 'FROZEN_PREDICTIONS_R4.json.sha256').write_text(sha(p) + '  FROZEN_PREDICTIONS_R4.json\n')
    out = {}
    for (group, ids) in [('test', partition['test']), ('all', manifest['cases'])]:
        out[group] = {}
        for mode in ['all', 'strict', 'moderate']:
            selected = [c for c in ids if mode == 'all' or selection[c][mode]]
            used = [c for c in selected if c in profiles]
            observ = {}
            for j in ['upper', 'lower']:
                for tag in ['R', 'posterior']:
                    values = np.array([sum((v for (t, v) in profiles[c][j].items() if (int(t) // 10 in (1, 4) if tag == 'R' else int(t) % 10 >= 4))) for c in used])
                    observ[j + '_' + tag] = metrics(values, tag) if len(values) else None
            allowed = len(used) >= 30
            shape = allowed and all((observ['upper_' + tag]['sample_quantile_gate'] == 'PASS' and observ['upper_' + tag]['gain_gate'] == 'PASS' for tag in ['R', 'posterior']))
            out[group][mode] = dict(n_source=len(ids), n_selected=len(selected), n_usable=len(used), excluded_eligibility=len(ids) - len(selected), excluded_zero_full16_area=len(selected) - len(used), adequacy_gate='PASS' if allowed else 'UNKNOWN_N_LT30', shape_gate=('PASS' if shape else 'FAIL') if allowed else 'NOT_ADJUDICATED', metrics=observ)
    normal = {'labels': {**FIELDS, 'missing_teeth': []}}
    bad = {'labels': {**FIELDS, 'missing_teeth': [], 'overjet': None}}
    assert eligible([normal])[0] and (not eligible([normal, bad])[0])
    result = dict(claim_type='information_link', groups=out, physical_population_validation='UNKNOWN even if eligibility passes: reporter truth/age/dentition≥24/calibrated loaded pose not measured', eligibility_rejections=rejected, source_hash=sha(REPORT), mutated_missing_field_rejected=True, full16_zero_area_count=len(cases) - len(profiles), wall_s=time.perf_counter() - start, max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write('raw/SUBGROUP_DENOMINATOR.json', result)
    print(json.dumps({g: {m: {k: v for (k, v) in d.items() if k != 'metrics'} for (m, d) in ms.items()} for (g, ms) in out.items()}))
    write('CURRENT_WORK_STATE.json', dict(lane='X94-force-validation', status='REPORT_SELECTION_AND_DENOMINATOR_ROUND_COMPLETE', latest_gate={g: {m: d['shape_gate'] for (m, d) in ms.items()} for (g, ms) in out.items()}, next_operation='Final exact-integral verification, source ledger, figure, review feedback and same-patient measurement handoff'))
if __name__ == '__main__':
    main()
