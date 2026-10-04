"""Add explicit final reporting contracts; numeric round files stay unchanged."""
from dental_release.paths import expand as _release_expand
import datetime as dt
import hashlib
import importlib.metadata as im
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X49_design_gate'))

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda : f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

def dump(p, d):
    p.write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    r = json.loads((HERE / 'results.json').read_text())
    raw = Path(r['run_directory'])
    refs = r['external_referent']
    for ref in refs:
        ref['can_refute'] = True
        ref['refutes_us'] = False
    r['external_referents'] = refs
    r['external_referent'] = dict(kind='published_code', locator='https://trimesh.org/trimesh.proximity.html', compared_quantity='Witness-to-all-target-facet distance and exact oriented geometry round-trip through published trimesh 3MF importer', refutes_us=False, can_refute=True, limitation='Independent code verification only; no physical same-object laboratory measurement')
    r['data_source'] = dict(kind='published_dataset', locator='https://doi.org/10.5281/zenodo.10597292', name='STS-Tooth3D, Wang et al.', license='CC BY 4.0', local_consumer='Five derived X1b crown and preparation exports', compared_quantity='Original exported oriented triangle geometry and numerical topology; not an independent fabricated Crown reference')
    r['coverage_accounting'] = dict(synthetic_pre_export=dict(rule_slots=190, unknown=84, fraction=84 / 190, all_injected_cases_included=True), X1b_pre_export=dict(rule_slots=50, unknown=40, fraction=0.8), X1b_after_crown_transport=dict(rule_slots=50, unknown=35, fraction=0.7, pass_count=13, fail_count=2), physical_validation=dict(cases=0, commercial_CAD_triples=0, registered_X1b_antagonists=0), rejected_X1b_preparations=dict(numerator=2, denominator=5, fraction=0.4, reason='M1/M2 have tiny facets below frozen cross-product-norm cutoff 1e-12 mm2. Numerical policy, not proven clinical/manufacturing defect.'), rejected_X1b_crown_transports=dict(numerator=0, denominator=5), independent_facet_queries=dict(rejected_records=0, reason='All selected labelled synthetic faces have locator and mm unit; completeness/anatomy remain asserted'))
    r['resolution_contract'] = dict(witness_coordinates='PER_POINT', surface_intervals='PER_SURFACE_REGION', mean_wall_summary='PER_TOOTH', case_rule_slot_fractions='descriptive aggregate over PER_SURFACE_REGION rule decisions; no population extrapolation', IFU_threshold='PER_SURFACE_REGION product specification', lab_bands='PHENOMENOLOGICAL', sinter_factor='PHENOMENOLOGICAL metadata; measured batch length ratio required to replace it', full_pose_uncertainty='UNKNOWN', timescale=dict(geometry_and_measurement='SIMULTANEOUS', design_to_future_manufacture='HANDOVER'))
    r['cost_accounting'] = dict(total_validation_wall_s=r['wall_s'], total_validation_cpu_s=r['cpu_s'], peak_RSS_MiB=r['peak_RSS_MiB'], fit_s=0, manual_discovery_cost='UNMEASURED', inherited_development_and_data_acquisition_cost='UNMEASURED; not omitted from claimed speedup because no speedup claimed', fallback_lab_measurement_cost='UNKNOWN', GPU_used=False, thread_limits='OMP/OPENBLAS/MKL/NUMEXPR=1; unchanged parent HiGHS explicitly sets 4, user hard limit <=4', prereg_deviation='planned threads=1 did not account for unchanged parent HiGHS setting=4; all executions remained within user <=4')
    r['limitations'] = ['No fabricated part or complete commercial CAD triple', 'Unsigned Hausdorff and fixed-XY axial analytic bounds rigorous under supplied assumptions; floating-point enclosure lacks outward-rounding proof', 'Physical geometry/registration error envelopes not measured', 'General insertion, rotational sweeps, complete machine/holder paths and self-intersection proof UNKNOWN', 'Physical cement film and fit UNKNOWN; group averages cannot replace local field', '3MF XSD and commercial CAM/import validation UNKNOWN', 'Local replay depends on pinned predecessor code paths; no portable lab-machine install claimed']
    r['manifest'] = 'ARTIFACT_MANIFEST.json'
    dump(HERE / 'results.json', r)
    dump(HERE / 'RUNTIME_LOCK.json', dict(python=sys.version, executable=sys.executable, packages={k: im.version(k) for k in ('numpy', 'scipy', 'trimesh', 'rtree', 'shapely')}, plotting='PYTHONNOUSERSITE=1 system numpy/matplotlib; no shared environment changes'))
    dump(HERE / 'CURRENT_WORK_STATE.json', dict(lane='X49-design-gate', status='ROUND_COMPLETE_REVIEW_PENDING', round='R3', latest_gate=f"{r['passed_tests']}/{r['total_tests']} frozen verification controls; physical lab validation UNKNOWN", ongoing='Final reader package and graph binding/coverage handoff', next_operation='Acquire one complete commercial CAD triple and measured pose/scan envelope; freeze 3D swept-solid/actual CAM predictions before bench observations', result_sha256=sha(HERE / 'results.json'), updated_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    files = []
    for root in (HERE, DATA):
        for p in sorted(root.rglob('*')):
            if not p.is_file() or p.name == 'ARTIFACT_MANIFEST.json' or '__pycache__' in p.parts:
                continue
            files.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)))
    dump(HERE / 'ARTIFACT_MANIFEST.json', dict(files=files, large_array_policy='No individual array >50MB; own data <3GB', local_bytes=sum((x['bytes'] for x in files if x['path'].startswith(str(HERE)))), external_bytes=sum((x['bytes'] for x in files if x['path'].startswith(str(DATA)))), review_state='PENDING_INDEPENDENT_REVIEW'))
    print('Final results contract and manifest saved; original round JSON unchanged')
if __name__ == '__main__':
    main()
