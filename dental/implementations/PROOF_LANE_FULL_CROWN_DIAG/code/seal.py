from common import *
import subprocess

def run():
    replay = read(ROOT / 'raw/REPLAY_PARITY.json')
    assert replay['mode'] == 'full_remeasurement' and replay['all_shape_and_contact_gates_identical']
    validation = read(ROOT / 'raw/VALIDATION.json')
    assert validation['status'] == 'PASS'
    res = read(ROOT / 'results.json')
    res['reproducibility'] = replay
    demo = read(ROOT / 'raw/DEMO.json')
    res['cost']['one_command_replay_seconds'] = demo['seconds']
    initial = {n: read(ROOT / 'history/initial_measurements' / n)['seconds'] for n in ['LOCAL_PAIRS.json', 'ANNOTATED_PAIRS.json', 'RESCORE.json']}
    res['cost']['initial_experiment_seconds'] = initial
    res['cost']['total_recorded_numerical_seconds_including_replay'] = sum(initial.values()) + demo['seconds']
    res['cost']['source_curation_and_coding_walltime'] = 'Not separately timed; UNKNOWN, not zero'
    res['original_threshold_measurement'] = dict(path='history/initial_measurements/ANNOTATED_PAIRS.json', sha256=sha(ROOT / 'history/initial_measurements/ANNOTATED_PAIRS.json'))
    res['data_licences'] = dict(Teeth3DS=dict(licence='CC BY-NC-ND4.0', locator='https://osf.io/download/9dutn/', local_file='raw/sources/teeth3ds_license.txt', sha256=sha(ROOT / 'raw/sources/teeth3ds_license.txt')), Bits2Bites='CC BY-NC-SA according to task brief', Bite2Text='UNKNOWN for local copy')
    res['data_artifacts'] = {str(p): dict(sha256=sha(p), bytes=p.stat().st_size) for p in DATA.glob('*.npz')}
    ownbytes = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    databytes = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert ownbytes + databytes < 3000000000
    res['cost']['own_disk_bytes'] = ownbytes
    res['cost']['data_disk_bytes'] = databytes
    dump(ROOT / 'results.json', res)
    import rescore
    dependencies = {}
    for module in list(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and '/3fold-workspaces/dental/' in path and Path(path).is_file():
            dependencies[str(Path(path).resolve())] = sha(path)
    for p in (ROOT / 'code').glob('*.py'):
        dependencies[str(p)] = sha(p)
    dependencies[str(ROOT / 'run_all.sh')] = sha(ROOT / 'run_all.sh')
    dump(ROOT / 'SOURCE_CODE_MANIFEST.json', dict(files=dependencies, source_access='read-only predecessor imports; no source edits', runtime_file='raw/RUNTIME.json'))
    snapshot = dict(results=res, code_manifest=read(ROOT / 'SOURCE_CODE_MANIFEST.json'), result_sha256=sha(ROOT / 'results.json'), created_utc=now(), frozen_protocol_hashes={p.name: sha(p) for p in ROOT.glob('PREREG_*.json')})
    freeze(ROOT / 'REVIEW_SNAPSHOT.json', snapshot)
    feedback = dict(target_id='DENT-VAL-BASELINE-COMPARISON', source_generation='gen-dca5b6abb1a2c582', result_file='results/PROOF_LANE_FULL_CROWN_DIAG/REVIEW_SNAPSHOT.json', sha256=sha(ROOT / 'REVIEW_SNAPSHOT.json'), review_state='PENDING_INDEPENDENT_REVIEW', measured_quantity='Bilateral natural-crown point/triangle p95; rigid shape and static contact diagnostic of72 frozen main constructions; exact-summary insufficiency witnesses', units='mm;mm2;counts', uncertainty='Clinical accepted-pair floor UNKNOWN.18 annotated pairs in6 upper arches; lower R3 transfer unvalidated. Finite probes/four-start ICP; no rigorous numeric or scanner enclosure.3/18 versus4/18 bilateral threshold exceedances across probe densities.', population_regime='Local digital geometry only; Teeth3DS annotated natural bilateral surfaces, Bite2Text/Bits2Bites inferred labels and virtual preparation. No clinical recommendation.', preregistered_gate='PREREG_LOCAL;PREREG_ANNOTATED;PREREG_RESCORE. Frozen bilateral maxima0.4371454294/0.7253699179/0.2398965819mm. Original0.35 and function gates retained.', baseline='Original single-reference0.35p95; ordinary rigidICP; equally informed R3 scalar_control actually rescored', outcome='0/72 main designs pass the descriptive bilateral shape envelope;0/72 shape+nominal function. Whole-surfacep95 insufficient for contact: identity error0,contact difference1mm2 and separate location shift3mm.Clinical floor remains UNKNOWN.', negative_result=True)
    dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    dump(ROOT / 'GRAPH_COVERAGE_PROPOSAL.json', dict(status='PROPOSAL_ONLY', parent_target='DENT-VAL-BASELINE-COMPARISON', missing_coverage='Reference identity and same-preparation accepted-design distribution separated from geometry reconstruction; actual finish-line and loaded-pose prerequisites', resolution='PER_POINT', time_scale='SIMULTANEOUS', proposed_edge='accepted same-preparation registered surfaces -> task-specific crown evaluation', no_native_mutation=True))
    dump(ROOT / 'ATTEMPTS.json', dict(attempts=[dict(id='LOCAL', outcome='18 paired measurements; rejected as clinical calibration', obstacle='inferred labels/fragments and no accepted same-prep pairs', next='ANNOTATED'), dict(id='ANNOTATED', outcome='18 original-label pairs;4/18 over.35 with8192probes,3/18 with2048', obstacle='natural asymmetry is not clinical crown variability; only6clusters', next='RESCORE'), dict(id='RESCORE', outcome='0/72shape and0/72shape+function; equally informed scalar0/18shape', obstacle='residual form/target support; actual prep margin/loaded contact absent', next='NEXT_MEASUREMENT.json')], external_referent=res['external_referent']))
    handoff = (ROOT / 'HANDOFF.md').read_text().replace("First measuring cycle complete. The full one command run is ongoing; its status is in raw/run_all.log. Clinical calibration is UNKNOWN, not a passed result.", "Both measuring cycles and the complete single command run are complete. Measurement values and form/contact decisions reproduced according to raw/REPLAY_PARITY.json. Clinical calibration is UNKNOWN, not a passed result.")
    handoff = handoff.replace("Update the final graph binding only after successful full rerun ; review_state should be PENDING_INDEPENDENT_REVIEW .", "GRAPH_FEEDBACK.json points to the immutable REVIEW_SNAPSHOT.json. Register/inspect receipt in raw/GRAPH_FEEDBACK_RECEIPT.json ; review_state is PENDING_INDEPENDENT_REVIEW . No scientific admission has taken place.")
    handoff += "\nDriving: `./run_all.sh`. Final running times, peak memory, disk budget and source hashes are available in results.json. Checks: all injected errors are dropped; component selections added no surface i36/36 annotated teeth. Full literature numbers and locators are available in LITERATURE.md and raw/PUBLISHED_MEASUREMENTS.json.\n"
    (ROOT / 'HANDOFF.md').write_text(handoff)
    state('COMPLETE_DIAGNOSTIC_WITH_OPEN_CLINICAL_CALIBRATION', res['main_counts'], 'Independent review; acquire actual same-preparation accepted designs via NEXT_MEASUREMENT.json before another generator')
    print(json.dumps(dict(result_sha256=sha(ROOT / 'results.json'), snapshot_sha256=sha(ROOT / 'REVIEW_SNAPSHOT.json'), counts=res['main_counts'], replay=replay, peak_MiB=res['cost']['max_observed_peak_rss_MiB'], disk_bytes=ownbytes + databytes)))
if __name__ == '__main__':
    run()
