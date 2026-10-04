from common import *
import shutil
import integrity

def run():
    integrity.verify()
    r = read(P / 'results.json')
    attempts = []
    for k in range(1, 8):
        p = read(P / f'PREREG_R{k}.json')
        g = read(P / f'raw/R{k}_GEOMETRY.json')
        rows = g['rows']
        executed = g.get('rescue_requested') if k == 5 else g.get('executed', sum((x['status'] != 'SOURCE_UNSCORED' for x in rows)))
        attempts.append(dict(round=f'R{k}', prereg_sha256=sha(P / f'PREREG_R{k}.json'), parent=p.get('parent_prereg_sha256'), operation=p['changed_operation'], requested=189, constructed=sum((x['status'] == 'ADAPTED' for x in rows)), geometry_accepted=sum((x.get('geometry_pass', False) for x in rows)), executed=executed, not_run=g.get('not_run', 0), outcome='R7_final' if k == 7 else 'Preserved prior construction', seconds=g['seconds'], external_physical_measurements=0))
    dump(P / 'ATTEMPTS.json', attempts)
    names = ['README_DEMO.md', 'LEADERBOARD.md', 'RESULTS.md', 'HANDOFF.md', 'FACIT.csv', 'raw/PHYSICAL_ROWS.json', 'raw/QUALITY_VECTOR.csv', 'raw/QUALITY_VECTOR_FULL.json.gz', 'raw/RANK_COMPARISONS.json', 'raw/FAULT_INJECTIONS.json', 'raw/SUFFICIENCY.json', 'raw/SOURCE_VALIDATION.json', 'raw/R7_GEOMETRY.json', 'raw/FRESH_GEOMETRY_REPLAY.json', 'raw/RUN_ACCOUNTING_REBUILD.json', 'raw/RUN_ACCOUNTING.json', 'INPUT_LOCK.json', 'CODE_SOURCE_LOCK.json', 'ARTIFACT_MANIFEST.json', 'FROZEN_PREDICTIONS_R7.json', 'FROZEN_PREDICTIONS_R4.json', 'ATTEMPTS.json', 'figures/physical_axes.png', 'figures/crown_field.png', 'exports/EXAMPLE.json', 'raw/EXPORT_CLI_CHECK.json', 'COMMANDS.md', 'LAB_PROTOCOL.md', 'RUNTIME_LOCK.json', 'LICENSES.json', 'SOURCE_REUSE.md'] + [f'PREREG_R{k}.json' for k in range(1, 8)] + [f'DECOMPOSITION_R{k}.json' for k in range(1, 8)]
    artifacts = []
    for name in names:
        src = P / name
        dest = P / 'review_artifacts' / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if sha(dest) != sha(src):
                raise ValueError('Review copy already exists with different bytes: ' + name)
        else:
            shutil.copyfile(src, dest)
        artifacts.append(dict(path=str(dest.relative_to(P)), source=name, sha256=sha(dest), bytes=dest.stat().st_size))
    snapshot = dict(claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', result=r, artifacts=artifacts, scientific_admission=False)
    freeze(P / 'REVIEW_SNAPSHOT.json', snapshot)
    packet = read(P / 'raw/GRAPH_PACKET.json')
    f = dict(target_id='DENT-VAL-BENCH-CROWN-FIT', result_file='results/PROOF_LANE_GENCAD_V5B/REVIEW_SNAPSHOT.json', sha256=sha(P / 'REVIEW_SNAPSHOT.json'), review_state='PENDING_INDEPENDENT_REVIEW', outcome='GEOMETRY_EXPORT_AND_CONDITIONAL_PHYSICAL_QUERIES; MATCHED_PHYSICAL_CALIBRATION_UNKNOWN', measured_quantity='Virtual preparation/crown gap and normal-ray thickness, decoded export topology/volume, source-conditioned force05 and Reynolds conductance, paired shape/force orderings', units='mm; um; source-conditioned N; m3/(Pa s); fractions', uncertainty='Discrete geometry without continuous enclosure; source parameter rectangles not jointCI; X1B provisional intervals and failed transfer; actual seated film/material/support/contact calibration UNKNOWN', population_regime='183 existing complete crowns,189 requests,18 sites/6case clusters/13participants; source-protocol curves; no matched manufactured measurements', preregistered_gate=[f'PREREG_R{k}.json' for k in range(1, 8)], baseline='Actual read-only X55/X1B/X13/BTE1 operations, independent numeric algebra and closed-form controls; source gates preserved. No algorithm-superiority claim.', negative_result=True, claim_type='capability', dispatch_scope='Review only: numerical experiment not admitted under this existing benchmark target', resolution=dict(geometry='PER_POINT', cement_hydraulics='PER_SURFACE_REGION', force_vector='PER_TOOTH', source_groups_and_counts='POPULATION', unmeasured_physical_transfer='PHENOMENOLOGICAL'), edge_contracts=r['edges'], missing_coverage_proposal='DENT-GENCAD-PREP-FILM-FRACTURE in HANDOFF.md; no native/source graph mutation', external_referent=r['external_referent'], scientific_admission=False, source_generation=packet.get('result_binding_template', {}).get('source_generation'))
    dump(P / 'GRAPH_FEEDBACK.json', f)
    state('REVIEW_PACKAGE_SEALED', dict(geometry=r['complete_panel'], ranking=r['ranking'], faults=r['controls']['count']), 'Register review-only feedback, independent scientific review and matched lab measurement')
    print(json.dumps(dict(snapshot_sha256=sha(P / 'REVIEW_SNAPSHOT.json'), review_artifacts=len(artifacts))))
if __name__ == '__main__':
    run()
