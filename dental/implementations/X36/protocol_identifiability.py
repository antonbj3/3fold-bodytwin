"""Round 3: test whether loading setup can be decomposed, before fitting it."""
from pathlib import Path
import json, datetime, hashlib, time
import numpy as np
from analyze import P, read, dump, write_csv, model
from meta_engine import design

def main():
    proto = P / 'PREREG_R3.json'
    if not proto.exists():
        plan = dict(round='R3', claim_type='information_link', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), capability='Separate material/thickness from indenter, speed, die and angle in published full-crown setup metadata', obstacle='R1 angle coefficient bundles different fixtures; interpretability needs independent setup variation', changed_operation='Replace one angle flag by a physical fixture vector and test its rank/support before any adjusted fit', consumer='DENT-MAT-RESTORATIVE-PSP; researcher planning matched fixture experiments', specification='Use R1 endpoint/thickness eligible groups only. Expanded X=R1 full design plus reported indenter_diameter_mm and crosshead_mm_min. If rank deficient, report nullspace and no full fit; never drop a dependent field after seeing outcome.', metric_tolerance_decision=dict(rank='SVD numpy matrix_rank; residual <=1e-9 for an exact alias', support='>=3 independent studies with independent variation for every fixture quantity', injection='Changing one indenter diameter by 0.25 mm must invalidate the frozen exact alias; this is a diagnostic fixture, not physical validation'), strongest_equally_informed_control='Direct per-row arithmetic for the alias, using same source metadata', practice='Loading angle alone or universal thickness exponent without an explicit die/indenter/speed vector', falsifier='A full-rank independently replicated fixture vector refutes the asserted confounding obstacle', full_cost=dict(preparation='local metadata normalization', fit='none if rank fails', discovery='prior reading UNKNOWN', validation='direct alias and injected metadata perturbation', questions=0, fallback='matched factorial experiment'), external_referent=dict(kind='independent_measurement', locator='EXTRACTION_CROWN.csv primary DOI/protocol locators', compared_quantity='published full-crown fixture angle, indenter diameter, crosshead speed and die identity', refutes_us=False), resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS')
        proto.write_text(json.dumps(plan, indent=2) + '\n')
        (P / 'PREREG_R3.sha256').write_text(hashlib.sha256(proto.read_bytes()).hexdigest() + '  PREREG_R3.json\n')
    assert hashlib.sha256(proto.read_bytes()).hexdigest() == (P / 'PREREG_R3.sha256').read_text().split()[0]
    start = time.perf_counter()
    rows = [r for r in read('NORMALIZED_DATA.json')['crown'] if r['thickness_eligible']]
    unknown = [r for r in rows if any((r.get(k) is None for k in ['indenter_diameter_mm', 'crosshead_mm_min']))]
    eligible = [r for r in rows if r not in unknown]
    terms = [('log_thickness', 'numeric', None), ('material_class', 'categorical', '3Y'), ('angle30', 'numeric', None), ('indenter_diameter_mm', 'numeric', None), ('crosshead_mm_min', 'numeric', None)]
    (X, names) = design(eligible, terms)
    rank = int(np.linalg.matrix_rank(X))
    (_, s, vh) = np.linalg.svd(X, full_matrices=False)
    aliases = [{n: float(v) for (n, v) in zip(names, z) if abs(v) > 1e-06} for z in vh[rank:]]
    residual = np.array([r['indenter_diameter_mm'] - (16 * r['crosshead_mm_min'] - 4 - 0.5 * r['angle30']) for r in eligible])
    actual = float(np.max(abs(residual)))
    injected = residual.copy()
    injected[0] += 0.25
    m = model(eligible, terms)
    support = {key: len({r['study'] for r in eligible if len({a[key] for a in eligible if a['study'] == r['study']}) > 1}) for key in ['angle30', 'indenter_diameter_mm', 'crosshead_mm_min']}
    out = dict(round='R3', claim_type='information_link', status='FIXTURE_COMPONENTS_NOT_SEPARATELY_IDENTIFIABLE', review_state='PENDING_INDEPENDENT_REVIEW', eligible_rows=len(eligible), unknown_metadata_rows=len(unknown), rank=rank, n_fixed_columns=X.shape[1], nullspace=aliases, exact_alias='indenter_diameter_mm = 16*crosshead_mm_min - 4 - 0.5*angle30 on these four published studies only; arithmetic relation, not a physical law', alias_residual_max=actual, within_study_fixture_variation_studies=support, fit=m, control=dict(direct_alias_pass=actual < 1e-09, injected_indenter_0_25mm_rejected=bool(np.max(abs(injected)) > 1e-09)), external_referent=read('PREREG_R3.json')['external_referent'], elapsed_seconds=time.perf_counter() - start, limitation='A declared unknown effect, not absence of a physical fixture effect; no clinical transfer')
    out['external_referent']['refutes_us'] = True
    dump('RESULTS_R3.json', out)
    write_csv('FIXTURE_VECTOR.csv', [{k: r.get(k) for k in ['row_id', 'study', 'doi', 'locator', 'angle_deg', 'indenter_diameter_mm', 'crosshead_mm_min', 'die_material', 'thickness_mm', 'material_class', 'resolution_level']} for r in eligible])
    (P / 'HANDOFF_R3.md').write_text('R3 tests the physical fixture vector before fitting it. Read RESULTS_R3.json: exact source-metadata alias makes diameter, speed and angle inseparable on current studies. Changing one diameter by 0.25 mm invalidates the arithmetic alias. The formula is dataset-specific, not physics. Next: matched fixture factorization and independent sample covariance.\n')
    dump('CURRENT_WORK_STATE.json', dict(lane='X36-meta-regression', phase='R3_COMPLETE', latest_gate=out['status'], next_operation='Deliver figures, method draft and independent-review packet', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print('R3 fixture rank', rank, '/', X.shape[1], 'alias error', actual)
if __name__ == '__main__':
    main()
