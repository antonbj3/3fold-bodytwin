import os, sys, resource, time, json, hashlib, datetime
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def verify():
    lock = json.loads((ROOT / 'RELEASE_LOCK.json').read_text())
    for r in lock['files']:
        q = ROOT / r['path'] if r['location'] == 'executed_copy_relative' else Path(r['path'])
        if not q.is_file() or sha(q) != r['sha256']:
            raise RuntimeError('RELEASE_HASH_MISMATCH ' + str(q))
    for pred in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        q = pred.with_suffix('.sha256')
        if q.exists() and q.read_text().split()[0] != sha(pred):
            raise RuntimeError('PREDICTION_HASH_MISMATCH ' + pred.name)
    print('release hashes checked against executed copy and read-only external inputs', flush=True)

def state(phase, gate, nextop):
    (ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X95-layer-thickness-nesting', 'status': phase, 'latest_gate': gate, 'next_operation': nextop, 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=2) + '\n')

def main():
    t0 = time.perf_counter()
    verify()
    if '--verify-only' in sys.argv:
        return
    import measure_geometry, pose_join, virtual_prep, validate_round2, height_core, coupled_height_core, exported_point_join, final_validation, report
    stages = [('R0', measure_geometry.run), ('R1', pose_join.run), ('R2', virtual_prep.run), ('R2_sufficiency', validate_round2.run), ('R3', height_core.run), ('R4', coupled_height_core.run), ('R5', exported_point_join.run), ('source_and_export_validation', final_validation.run)]
    phases = {}
    for (name, fn) in stages:
        state(name + '_RUNNING', 'Frozen gates unchanged', 'Complete ' + name + ' under fixed inputs')
        start = time.perf_counter()
        fn()
        phases[name] = time.perf_counter() - start
    cost = {'replay_phase_wall_seconds': phases, 'replay_wall_seconds_before_report': time.perf_counter() - t0, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_max': 4, 'GPU': False, 'fit_seconds': 0, 'source_discovery_and_manual_image_binding_cost': 'included in bounded4h session; separately measured human/agent active time UNKNOWN', 'historical_data_collection_cost': 'UNKNOWN', 'physical_lab_cost': 'UNKNOWN_NOT_RUN', 'commercial_time_or_remake_benefit': 'UNKNOWN', 'all_library_threads_enforced': True, 'address_space_limit_MiB': 3500}
    tr = time.perf_counter()
    out = report.run(cost)
    cost['report_wall_seconds'] = time.perf_counter() - tr
    cost['replay_total_wall_seconds'] = time.perf_counter() - t0
    cost['peak_rss_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    cost['own_file_bytes'] = sum((q.stat().st_size for q in ROOT.rglob('*') if q.is_file()))
    assert cost['own_file_bytes'] < 3000000000
    out['cost'] = cost
    pose_join.dump(ROOT / 'results.json', out)
    pose_join.dump(ROOT / 'COST.json', cost)
    feedback = {'target_id': 'DENT-MFG-PROCESS-MODEL', 'graph_binding_status': 'UNBOUND_DISPATCH_REJECTED_MISSING_EXACT_COVERAGE', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'result_file': str(ROOT / 'results.json'), 'sha256': sha(ROOT / 'results.json'), 'measured_quantity': 'YML half-disc wall slack at exterior points and exported virtual heightpreparation, not empirical manufacture/gap accuracy', 'units': 'mm and mm3', 'resolution': 'PER_POINT; consumer volume aggregation PER_TOOTH', 'uncertainty': 'source and actual preparation UNKNOWN; rigorous float/source enclosure MISSING', 'population_regime': 'three frozen R4 virtual scenarios; no physical measurements', 'preregistered_gate': 'PREREG_B05 + PREREG_R1..R5; originalfixedintaglio allposesfail, sourcecoupledvirtualcore3/3 only', 'baseline': 'inherited engineering0.5mm ablation; separately matched independent pose/triangle/column controls', 'outcome': out['outcome'], 'negative_result': True, 'negative_scope': 'R1 originalallpose failure; R2 nonmanifold export; R3 no positive anterior fixed-footprint core; no physicalqualification', 'dispatch_failure': 'numerical_dispatch_not_allowed and unresolved_prerequisites; no bypass', 'coverage_proposal': {'proposed_id': 'DENT-DESIGN-YML-LAYER-WALL-POSE', 'status': 'PROPOSED_NOT_NATIVE', 'description': 'sourcebound local minimumwall as function of crown indication, localcoordinate, discpose and half-disc region', 'inputs': ['case/specimen', 'externalmesh/intaglio hash', 'point/triangleID', 'sourceRuleID', 'frame', 'pose', 'units'], 'outputs': ['local requiredwall', 'slack', 'rejectedpose witness', 'conditional virtual preparation geometry'], 'consumer': ['CAD', 'lab', 'material'], 'proposed_edges': [{'from': 'YML source rule and locked crown/intaglio', 'to': 'local CAD/nesting constraint', 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'status': 'PENDING_INDEPENDENT_REVIEW'}, {'from': 'local constraint', 'to': 'DENT-MFG-PROCESS-MODEL consumer', 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'status': 'PROPOSED_CONDITIONAL_NOT_ADMITTED'}]}}
    pose_join.dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    state('COMPLETE_SCOPED_DIGITAL_DEMO', 'R1 negative preserved; R4/R5 virtualgeometry capability, physical qualificationUNKNOWN', 'Independent review; next replace virtualcore by matched actualprep and batch IFU; close arithmetic/source enclosure')
    print('demo complete:', out['outcome'], flush=True)
if __name__ == '__main__':
    main()
