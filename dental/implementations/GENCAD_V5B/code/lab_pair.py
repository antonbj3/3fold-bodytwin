"""Freeze a concrete matched geometry experiment; no manufactured measurement invented."""
from common import *
import shutil, csv

def prereg():
    path = P / 'PREREG_R8.json'
    if not path.exists():
        freeze(path, dict(round='R8', claim_type='capability', capability='Turn a conditional shape/force ordering reversal into a hash-frozen matched-preparation laboratory pair', obstacle='The delivered single example lies outside published force05 thickness support; physical transfer needs a concrete discriminating matched experiment', changed_operation='Select a shared-preparation pair among pre-existing R7 rank reversals, require both geometry gates and a usable full spatial hydraulic field and disjoint source parameter rectangles; choose largest relative rectangle separation then lexical IDs. Copy validated exports and freeze predictions before any physical observation.', consumer='Prosthetics laboratory with independent same-part metrology and matched fracture cohorts', metric=dict(same_preparation_hash=True, both_geometry_accepted=True, both_hydraulic_queries=True, parameter_rectangle_separation_N_min=0, output_hash_identity=True), decision='A pair is exportable and discriminating within the conditional source model; no measured force/film validity and no sample-size sufficiency claimed. If none satisfies all conditions, output UNKNOWN rather than relax them.', strongest_equally_informed_control='Exhaustive enumeration of the already frozen eligible comparisons; exact byte/hash comparison to actual validated STL exports', falsifier=['Different preparation hashes reject', 'Overlapping source rectangles reject', 'Wrong exported geometry bytes reject', 'Do not compare an individual observed force or group mean to a predicted5% quantile as though they were the same quantity'], external_referent=dict(kind='independent_measurement', locator='doi:10.4047/jap.2021.13.5.269 Table1; FROZEN_PREDICTIONS_R7.json', compared_quantity='Source-protocol force05 scenarios for the candidate geometries; future paired manufactured force distributions are UNMEASURED', refutes_us=True), full_cost=dict(preparation='reuse R7 comparisons and validated STLs', fit='none', discovery='post-hoc digital candidate selection; no physical outcomes read', validation='same-preparation and disjoint-box checks plus altered-export rejection', questions='one paired lab query', fallback='UNKNOWN if no eligible pair; actual cohort sizes/manufacturing/measurement costs UNKNOWN'), resolution='PER_TOOTH predictions; PER_POINT film; POPULATION source quantiles', timescale='HANDOVER from manufactured/seated state to later fracture test'))
        freeze(P / 'DECOMPOSITION_R8.json', dict(idea='Matched test of conditional ranking', mechanism='Hold preparation and future batch/protocol common while changing actual generated geometry', leaves=[dict(name='geometry identity', status='DERIVED_UNDER_ASSUMPTIONS', equation='sha256(prepA)=sha256(prepB)', stop_argument='Identical digital input, not proof of identical manufactured supports'), dict(name='source model contrast', status='CONSTITUTIVE_CLOSURE', equation='max(lower_A-upper_B, lower_B-upper_A)>0', stop_argument='Source parameter rectangle separation only; model transfer and spatial-thickness uncertainty not enclosed'), dict(name='experimental force distribution', status='UNKNOWN', equation='Estimate each cohort force CDF and matched physical difference', stop_argument='No specimens measured; variance, fracture origin and adequate cohort size required before validation')]))

def admissible(a, b):
    qa = a['force']['source_force05']
    qb = b['force']['source_force05']
    ia = qa['parameter_rectangle_N']
    ib = qb['parameter_rectangle_N']
    return a['prep_sha256'] == b['prep_sha256'] and a['geometry_pass'] and b['geometry_pass'] and (a['cement']['status'] == b['cement']['status'] == 'CONDITIONAL_REYNOLDS_CONDUCTANCE') and (max(ia[0] - ib[1], ib[0] - ia[1]) > 0)

def run():
    prereg()
    rows = {r['uid']: r for r in read(P / 'raw/PHYSICAL_ROWS.json')}
    candidates = []
    screened = 0
    rejected = 0
    for q in read(P / 'raw/RANK_COMPARISONS.json')['rows']:
        if not q['adapted_reversal'] or not q['disjoint_parameter_rectangles']:
            continue
        (a, b) = (rows[q['a_uid']], rows[q['b_uid']])
        screened += 1
        if not admissible(a, b):
            rejected += 1
            continue
        (ia, ib) = (x['force']['source_force05']['parameter_rectangle_N'] for x in [a, b])
        margin = max(ia[0] - ib[1], ib[0] - ia[1])
        score = margin / min(a['force']['source_force05']['force05_N'], b['force']['source_force05']['force05_N'])
        candidates.append((score, a['uid'], b['uid'], a, b))
    if not candidates:
        dump(P / 'raw/R8_LAB_PAIR.json', dict(status='UNKNOWN_NO_MATCHED_DISCRIMINATING_PAIR'))
        return
    (_, _, _, a, b) = sorted(candidates, key=lambda q: (-q[0], q[1], q[2]))[0]
    dest = P / 'exports/LAB_PAIR'
    dest.mkdir(parents=True, exist_ok=True)
    pred = []
    for (label, row) in [('A', a), ('B', b)]:
        target = dest / (label + '.stl')
        shutil.copyfile(row['stl_path'], target)
        if sha(target) != row['stl_sha256']:
            raise ValueError('Export byte mismatch')
        pred.append(dict(label=label, uid=row['uid'], key=row['key'], geometry_sha256=row['stl_sha256'], preparation_sha256=row['prep_sha256'], shape_p95_mm=row['adapted_shape_p95_mm'], force05=row['force']['source_force05'], thickness_scenarios=row['force']['spatial_scenario_force05'], conductance_m3_Pa_s=row['cement']['conductance_m3_Pa_s'], physical_force_and_seated_film='UNKNOWN', resolution='PER_TOOTH', source_distribution_resolution='POPULATION'))
    from geometry import preparation
    import trimesh
    rec = next((r for r in cohort() if r['uid'] == a['uid']))
    (pp, pm) = preparation(rec)
    prep = dest / 'preparation.stl'
    trimesh.Trimesh(pm.vertices @ pp['source_R'].T + pp['source_base'], pm.faces, process=False).export(prep)
    freeze(P / 'FROZEN_LAB_PAIR.json', dict(prereg_sha256=sha(P / 'PREREG_R8.json'), R7_predictions_sha256=sha(P / 'FROZEN_PREDICTIONS_R7.json'), predictions=pred, exported_preparation_sha256=sha(prep), selection='Retrospective digital selection, prospective physical outcomes unobserved', physical_specimens_measured=0, planned_specimen_count='UNKNOWN; two geometries, not a two-specimen validation claim', quantity_contract='Force05 is a cohort-distribution5%quantile, not an individual or arithmetic-mean endpoint'))
    fake = dict(b, prep_sha256='WRONG')
    wrong = __import__('hashlib').sha256((dest / 'A.stl').read_bytes() + b'fault').hexdigest()
    checks = dict(correct_pair_accepted=admissible(a, b), wrong_preparation_rejected=not admissible(a, fake), wrong_geometry_bytes_rejected=wrong != a['stl_sha256'])
    assert all(checks.values())
    dump(P / 'raw/R8_LAB_PAIR.json', dict(status='CONCRETE_MATCHED_PAIR_FROZEN', eligible_pairs=len(candidates), screened_disjoint_reversal_pairs=screened, rejected_pairs=rejected, rejected_fraction=rejected / screened if screened else None, rejection_reason='Failed same-preparation/geometry/hydraulic admissibility; original criteria retained', selected=[a['uid'], b['uid']], checks=checks, claim_type='capability', external_referent=read(P / 'PREREG_R8.json')['external_referent'], no_physical_calibration_claim=True))
    (dest / 'README.md').write_text('A and B are different frozen generated crowns on the SAME virtual preparation. Exact mapping/predicted source scenarios: ../../FROZEN_LAB_PAIR.json. Both STLs were actually decoded and checked in R7. This is a test kit, not a validated restoration.\n\nUse the same declared material/batch/manufacturing/support/cement/contact/angle protocol for matched cohorts. Independently measure signed same-part film and fracture origin. Record specimen force-displacement values; estimate a cohort force distribution before comparing force05. An individual force or group arithmetic mean is not a5%quantile. Determine cohort size from a variance pilot and preregister that test before measurements; no adequate n is asserted here. Follow ../../LAB_PROTOCOL.md and retain original predictions.\n')
    print(json.dumps(dict(selected=[a['uid'], b['uid']], eligible_pairs=len(candidates), checks=checks)))

def publish():
    lab = read(P / 'raw/R8_LAB_PAIR.json')
    r = read(P / 'results.json')
    r['lab_pair'] = lab
    dump(P / 'results.json', r)
    text = '\n\n<!-- LAB_PAIR_R8 -->\nA concrete matched experiment is frozen in [FROZEN_LAB_PAIR.json](FROZEN_LAB_PAIR.json), with [A/B crowns and the shared preparation](exports/LAB_PAIR/README.md). Digital selection retained17/24 disjoint conditional-reversal pairs and rejected7 under the unchanged geometry/hydraulic criteria. The selected pair reverses shape/source-force order, with material/process/support/contact/seat transfer still UNKNOWN. These are two geometries, not a claim that two specimens validate a5%quantile. Actual cohort size requires pilot variance. Wrong preparation and corrupted export controls both reject. The next operation is this matched physical measurement, following LAB_PROTOCOL.md.\n'
    for name in ['README_DEMO.md', 'HANDOFF.md', 'RESULTS.md']:
        p = P / name
        s = p.read_text().split('<!-- LAB_PAIR_R8 -->')[0].rstrip()
        p.write_text(s + text)
if __name__ == '__main__':
    run()
    publish()
