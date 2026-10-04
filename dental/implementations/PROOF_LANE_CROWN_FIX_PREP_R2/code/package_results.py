from common import *
from collections import Counter
from text_stl import write as text_stl

def ref(p):
    p = Path(p)
    return dict(path=str(p), sha256=sha(p))

def run():
    selected = {}
    allrows = {}
    for tag in ['C', 'E']:
        rows = read(R / f'RESULTS_{tag}.json')['rows']
        allrows[tag] = rows
        paths = {(x['key'], x['material']): x['certificate'] for x in read(R / ('RESULTS_PATH.json' if tag == 'C' else 'RESULTS_PATH_E.json'))['rows']}
        films = {(x['key'], x['material']): x for x in read(R / f'RESULTS_FILM_{tag}.json')['rows']}
        dies = {(x['key'], x['material']): x['preparation'] for x in read(R / f'RESULTS_PREPARATION_EXPORT_{tag}.json')['rows']}
        for row in rows:
            key = (row['key'], row['material'])
            if not row.get('restored_assembly_geometric_conjunction'):
                continue
            path = paths[key]
            film = films[key]
            die = dies[key]
            assert path['isolated_path_pass'] and path['subset_exact'] and film['continuous_pass']
            assert die['watertight'] and die['self_intersections'] == 0
            selected[key] = (row, path, film, die)
    exports = []
    for (row, path, film, die) in selected.values():
        m = load(row['mesh_path'])
        stem = Path(row['mesh_path']).stem
        stl = text_stl(D / (stem + '_EXACT_TEXT_RESEARCH.stl'), m['vertices'][m['faces']])
        assert stl['exact_triangle_roundtrip']
        part = dict(key=row['key'], fdi=row['fdi'], family=row['family'], material=row['material'], round=row['round'], resolution='PER_TOOTH', status='FROZEN_DIGITAL_GEOMETRY_RESEARCH_ONLY', full_crown_qualified=False, mesh=ref(row['mesh_path']), ascii_stl=stl, three_mf=ref(row['exports']['three_mf']), three_mf_exact_roundtrip=row['exports']['exact_roundtrip'], geometry_certificate=ref(R / 'exports' / (stem + '_CERTIFICATE.json')), continuous_film=ref(R / 'raw' / (stem + '_FILM.json')), path_certificate=ref(path['path']), full_preparation_formula=path['full_preparation_formula'], research_die=die, assembly_p95_mm=row['assembly_shape']['p95_mm'], assembly_p95_limit_mm=read(R / 'PREREG_B.json')['metrics']['shape_p95_mm'][row['family']], wall_lower_bound_mm=row['material_contract']['wall_mm'], active_nominal_film_interval_mm=film['interval_mm'], scope='Digital source-surface restored assembly. Original crown-only form still fails; material/biological/CAM/physical qualification UNKNOWN.')
        cp = R / 'exports' / (stem + '_ADMITTED_CERTIFICATE.json')
        dump(cp, part)
        part['admitted_certificate'] = ref(cp)
        for p in [stl['path'], row['exports']['three_mf'], die['path']]:
            dest = R / 'exports' / Path(p).name
            if not dest.exists():
                dest.symlink_to(p)
        exports.append(part)
    out = dict(frozen_utc=now(), claim_type='capability', selection='For each tooth/material, last C/E construction passing all assembly gates, exact path, continuous film and closed nonintersecting research die. No threshold retuned.', variants=len(exports), unique_teeth=len(set((x['key'] for x in exports))), rows=exports)
    dump(R / 'EXPORTS.json', out)
    families = {}
    for family in ['anterior', 'premolar', 'molar']:
        keys = {x['key'] for x in exports if x['family'] == family}
        families[family] = dict(n_teeth=6, assembly_geometry_candidates=len(keys), complete_digital_crowns=0, complete_chain=0, resolution='PER_TOOTH')
    cdie = next((x['preparation'] for x in read(R / 'RESULTS_PREPARATION_EXPORT_C.json')['rows'] if x['key'] == '06c45662b50389bd_premolar'))
    edie = next((x['research_die'] for x in exports if x['key'] == '06c45662b50389bd_premolar' and x['material'] == 'KATANA_ML'))
    (cv, ev, tv) = (cdie['retained_preparation_volume_mm3'], edie['retained_preparation_volume_mm3'], edie['native_closed_model_volume_mm3'])
    retention = dict(key='06c45662b50389bd_premolar', material='KATANA_ML', resolution='PER_TOOTH', native_model_volume_mm3=tv, flat_C_full_preparation_mm3=cv, adaptive_E_full_preparation_mm3=ev, additional_retained_mm3=ev - cv, removed_volume_reduction_percent=100 * (ev - cv) / (tv - cv), scope='Signed tetrahedral volume of emitted full-preparation approximations in same digital native tooth at same plane; source/physical uncertainty UNKNOWN. Not a population effect. Exact subset belongs to implicit CSG P.')
    rounds = {}
    for (tag, filename) in [('A', 'RESULTS_A_PILOT.json'), ('B', 'FROZEN_PREDICTIONS_B_PILOT.json'), ('C', 'RESULTS_C.json'), ('D', 'RESULTS_D.json'), ('E', 'RESULTS_E.json')]:
        rows = read(R / filename)['rows']
        counts = Counter((x['status'] for x in rows))
        reasons = Counter()
        for x in rows:
            if x['status'] == 'GENERATED':
                continue
            z = str(x.get('reason', x.get('error', '')))
            reason = 'EXACT_NATIVE_CLEARANCE' if x['status'] == 'REJECTED_EXACT' else 'HEIGHT_BELOW_PREREG_MINIMUM' if 'FIELD_PEAK' in z else 'SECTION_OR_CAP_TOPOLOGY' if 'LOOP' in z or 'CAP' in z or 'RADIAL' in z else z
            reasons[reason] += 1
        generated = sum((x['status'] == 'GENERATED' for x in rows))
        accepted = sum((x.get('restored_assembly_geometric_conjunction', False) for x in rows))
        rounds[tag] = dict(requested=len(rows), status_counts=dict(counts), rejected_before_scoring=len(rows) - generated, rejected_before_scoring_fraction=(len(rows) - generated) / len(rows), reasons=dict(reasons), assembly_geometry_accepted=accepted, complete_chain=0, evidence=ref(R / filename))
    per_tooth = []
    for rec in inputs():
        admitted = [dict(round=x['round'], material=x['material'], certificate=x['admitted_certificate']) for x in exports if x['key'] == rec['key']]
        per_tooth.append(dict(key=rec['key'], family=rec['family'], fdi=rec['fdi'], resolution='PER_TOOTH', geometry_candidates=admitted, complete_digital_crown=False, complete_chain=False, original_crown_only_form='FAIL_OR_NO_CANDIDATE'))
    costs = {p.stem: read(p) for p in (R / 'raw').glob('*RESOURCES.json') if p.name != 'DEMO_RESOURCES.json'}
    for tag in ['C', 'D', 'E']:
        rows = read(R / f'RESULTS_{tag}.json')['rows']
        costs[tag + '_recorded_seconds'] = dict(generation=sum((x.get('generation_seconds', 0) for x in rows)), scoring=sum((x.get('score_seconds', 0) for x in rows)), scope='Retained per-row timing; inherited E failures and pilot repetitions are not an additive total. Missing component-rescore timing UNKNOWN.')
    result = dict(lane='PROOF_LANE-crown-fix-prep-r2', updated_utc=now(), claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', outcome='Native-contained digital preparation and restored-assembly geometry for4/18 teeth; original complete crown chain0/18.', external_referent=read(R / 'PREREG_A.json')['external_referent'], referent_scope='Original independently acquired Teeth3DS IOS geometry, used during construction and checking; not a held-out preparation or physical-fit measurement. Native closure/CEJ/pulp not independently measured.', baseline=dict(lane='../PROOF_LANE_CROWN_FIX_PREP', complete_chain=0, teeth=18, exact_outside_witnesses=28, scope='Previous local rule preparations. New assembly metric must not be substituted for its original crown-only form metric.'), by_family=families, teeth=per_tooth, exports=ref(R / 'EXPORTS.json'), rounds=rounds, retained_substance=retention, sufficiency=read(R / 'raw/CONTROLS.json')['sufficiency'], validation=dict(fresh_recompute=ref(R / 'raw/FRESH_RECOMPUTE.json'), paired_valid_and_fault_checks=read(R / 'raw/CONTROLS.json')['count'] + read(R / 'raw/REAL_CONTROLS.json')['count'], controls=[ref(R / 'raw/CONTROLS.json'), ref(R / 'raw/REAL_CONTROLS.json')], certificate_replay='PENDING_FINAL_REPLAY'), cost=dict(known_receipts=costs, compile_peak_RSS_KiB=1116188, threads=1, GPU=False, memory_ceiling_MiB=3500, full_historical_peak='UNKNOWN; monitored process-tree and compiler receipts reported separately', preparation='Reuse closed native/outer supports and pinned parent code; original preprocessing not retimed', fit='No fitted clinical model', discovery='C/D/E fresh replay289.4s including all108 requests; A/B pilots separate retained logs', validation='Exact geometry certificates and mesh/shape/path/gap/export checks; receipts retained', questions=0, fallback='Rejected requests preserved; no heuristic acceptance fallback', full_token_cost='UNKNOWN'), limitations=['No complete original digital or physical crown qualification; old full-native crown-only p95 fails all emitted C/E candidates.', 'Restored-assembly p95 measures preserved lower native exterior plus crown exterior; area-stratified sampled statistic with no rigorous full-surface shape enclosure.', 'Exact P subset T concerns implicit CSG in completed digital support; exported research die is welded/rounded approximation.', 'Source finish curve measured from IOS model, but above-CEJ placement and pulp preservation UNKNOWN.', 'KATANA wall minima published; inherited0.5mm shoulder-width scenario not numerically supported by cited guide. Lithium disilicate1mm shoulder documented, but rounding/anatomical support/regular thickness unresolved.', 'Full-arch insertion, certified CAM, cement/seating forces, fatigue/fracture and physical fit UNKNOWN.', 'Film40–60um is continuous digital nearest-surface distance on active intaglio; shoulder0–50um separate CAD ramp. Neither is measured cement.', 'C079 anterior fails form; E098 premolar has98 crown self-intersections; neither is admitted.', 'No molar passes construction/certification; no small topology component discarded.'], phenomenological_debts=[dict(quantity='CAD nominal cement film', resolution='PHENOMENOLOGICAL', value_mm=0.05, replacing_measurement='Registered seated-gap field from independently measured crown/die or micro-CT', target_resolution='PER_POINT'), dict(quantity='Biological finish-line placement', resolution='PER_SURFACE_REGION', value='UNKNOWN', replacing_measurement='Independent CEJ and pulp annotations in same registered tooth frame')], edges=[dict(source='Native IOS facets and closed-support certificate', consumer='Preparation/CAD generator', resolution='PER_POINT', timescale='SIMULTANEOUS', quantity='Spatial containment and finish-curve registration'), dict(source='Frozen crown/die geometry', consumer='Manufacture then independent scan registration', resolution='PER_POINT', timescale='HANDOVER', quantity='Nominal surface geometry; as-built state unmeasured')], next_construction='Registered nonplanar source finish curve with independently labelled CEJ/pulp, then exact clearance-preserving shoulder fillet and collision-free collar; retain multi-loop molar failures as counterexamples.')
    dump(R / 'results.json', result)
    print(json.dumps(dict(variants=len(exports), teeth=out['unique_teeth'], by_family=families, retention=retention), indent=2))
if __name__ == '__main__':
    run()
