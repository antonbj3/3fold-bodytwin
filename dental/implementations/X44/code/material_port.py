"""K13: typed region -> constitutive operator. Unknown laws remain unknown.
This is a local research material registry, not a patch to BodyTwin's source registry.
"""
from dental_release.paths import expand as _release_expand
import numpy as np
from common import *

class InvalidMaterial(ValueError):
    pass

def isotropic(E, nu):
    if not (E > 0 and -1 < nu < 0.5):
        raise InvalidMaterial('non-positive elastic energy or incompressibility needs mixed formulation')
    la = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = la
    np.fill_diagonal(C[:3, :3], la + 2 * mu)
    C[3:, 3:] = np.eye(3) * mu
    return C

def resolve(entry, region, context):
    if region != entry['region']:
        raise InvalidMaterial('wrong anatomical phase')
    if entry['parameter_status'] == 'UNKNOWN':
        raise InvalidMaterial('No measured conditional law for this source body')
    unit = entry['E']['unit']
    value = entry['E']['value']
    if unit not in ['MPa', 'GPa']:
        raise InvalidMaterial('unrecognized stress unit')
    E = value * (1000 if unit == 'GPa' else 1)
    v = entry.get('validity', {})
    for key in ['jaw', 'tooth_class', 'phase', 'rate_protocol']:
        if key in v and context.get(key) != v[key]:
            raise InvalidMaterial('outside ' + key + ' validity')
    if 'force_N' in v and (not v['force_N'][0] <= context.get('force_N', float('-inf')) <= v['force_N'][1]):
        raise InvalidMaterial('outside force validity')
    if 'E_MPa' in v and (not v['E_MPa'][0] <= E <= v['E_MPa'][1]):
        raise InvalidMaterial('parameter/unit not in evidence-supported interval')
    if 'nu' in v and context.get('nu') != v['nu']:
        raise InvalidMaterial('E and nu are not separately identified')
    return isotropic(E, entry['nu'])

def run():
    p = verify('PREREG_K13.json')
    st = time.perf_counter()
    g = json.loads((L / 'raw/K09.json').read_text())
    data = np.load(D / 'region_field.npz')
    tag = data['tag']
    body = data['body']
    source = ROOT / 'results/F1_pdl_constitutive/runs.jsonl'
    rr = [json.loads(x) for x in source.read_text().splitlines()]
    curves = {}
    for (case, unn) in [('P12_max', 8), ('P12_max', 9), ('P14_max', 8), ('P14_max', 9)]:
        r = [x for x in rr if x.get('tag') == 'lin' and x['case'] == case and (x['unn'] == unn) and (x['etype'] == 'C3D10') and (x['nu_pdl'] == 0.45) and x['ok']]
        r = sorted(r, key=lambda x: x['E_pdl'])
        E = np.array([x['E_pdl'] for x in r])
        F = 0.15 / np.array([x['u_crowncentre_along_F_mm'] / x['F_N'] for x in r])
        bounds = np.exp(np.interp(np.log([7, 16.2]), np.log(F), np.log(E)))
        curves[f'{case}:{unn}'] = dict(E_MPa=E.tolist(), F_at_0p15_N=F.tolist(), conditional_E_interval_MPa=bounds.tolist(), resolution='PER_TOOTH')
    lower = max((x['conditional_E_interval_MPa'][0] for x in curves.values()))
    upper = min((x['conditional_E_interval_MPa'][1] for x in curves.values()))
    E0 = float(np.sqrt(lower * upper))
    assert lower < upper
    base = dict(model='linear_isotropic_elastic', parameter_evidence_resolution='PHENOMENOLOGICAL', resolution='PER_POINT', timescale='SIMULTANEOUS', parameter_status='CONSTITUTIVE_CLOSURE')
    entries = {'tooth': dict(base, region=1, E=dict(value=18600.0, unit='MPa'), nu=0.31, locator='https://doi.org/10.1371/journal.pone.0069990', validity={}, debt='Homogeneous whole tooth closure; no enamel/dentin split or patient E measurement'), 'bone': dict(base, region=3, E=dict(value=11500.0, unit='MPa'), nu=0.3, locator='https://doi.org/10.1371/journal.pone.0069990', validity={}, debt='Homogeneous source bone scenario, no cortical/cancellous calibrated map'), 'pdl_mandible_unknown': dict(region=2, parameter_status='UNKNOWN', resolution='PER_POINT', timescale='SIMULTANEOUS', debt='Same-tooth mandibular axial and lateral displacement with paired geometry/rate'), 'pdl_upper_incisor_conditional': dict(base, region=2, parameter_status='CONDITIONAL_POPULATION_CALIBRATION', E=dict(value=E0, unit='MPa'), nu=0.45, locator=p['external_referent']['locator'], validity=dict(jaw='maxilla', tooth_class='central_incisor', phase='loading', rate_protocol='F1_static_Keilig_0p1s_population_anchor', force_N=[5, 16], nu=0.45, E_MPa=[lower, upper]), debt='Geometry mismatch, constructed OFJ PDL and nu remain; effective population transfer only')}
    code = np.where(tag == 1, 1, np.where(tag == 3, 3, 2)).astype(np.int8)
    np.savez_compressed(D / 'material_binding.npz', body=body, region=tag, law_code=code, known_closure=tag != 2)
    ctx = dict(jaw='maxilla', tooth_class='central_incisor', phase='loading', rate_protocol='F1_static_Keilig_0p1s_population_anchor', force_N=10, nu=0.45)
    C = resolve(entries['pdl_upper_incisor_conditional'], 2, ctx)
    Cs = [resolve(entries['tooth'], 1, {}), resolve(entries['bone'], 3, {}), C]
    positive = all((np.linalg.eigvalsh(c).min() > 0 for c in Cs))
    rng = np.random.default_rng(44)
    energy_errors = []
    for entry in [entries['tooth'], entries['bone'], entries['pdl_upper_incisor_conditional']]:
        E = entry['E']['value']
        nu = entry['nu']
        la = E * nu / ((1 + nu) * (1 - 2 * nu))
        mu = E / (2 * (1 + nu))
        A = rng.normal(size=(30, 3, 3))
        eps = (A + A.transpose(0, 2, 1)) / 2
        strain = np.column_stack([eps[:, 0, 0], eps[:, 1, 1], eps[:, 2, 2], 2 * eps[:, 0, 1], 2 * eps[:, 1, 2], 2 * eps[:, 0, 2]])
        numerical = 0.5 * np.einsum('bi,ij,bj->b', strain, isotropic(E, nu), strain)
        reference = 0.5 * la * np.trace(eps, axis1=1, axis2=2) ** 2 + mu * np.einsum('bij,bij->b', eps, eps)
        energy_errors.append(float(np.max(np.abs(numerical - reference) / np.maximum(reference, 1e-12))))
    bad = []
    entry = entries['pdl_upper_incisor_conditional']
    for (name, ent, region, context) in [('wrong_region', entry, 3, ctx), ('MPa_mislabeled_GPa', dict(entry, E=dict(value=E0, unit='GPa')), 2, ctx), ('wrong_jaw', entry, 2, dict(ctx, jaw='mandible')), ('wrong_force', entry, 2, dict(ctx, force_N=100)), ('wrong_nu', entry, 2, dict(ctx, nu=0.49)), ('wrong_phase', entry, 2, dict(ctx, phase='unloading')), ('unknown_mandible', entries['pdl_mandible_unknown'], 2, {}), ('invalid_elasticity', dict(entry, nu=0.5), 2, ctx)]:
        try:
            resolve(ent, region, context)
            bad.append(dict(name=name, rejected=False))
        except InvalidMaterial as ex:
            bad.append(dict(name=name, rejected=True, reason=str(ex)))
    forces = {k: float(np.exp(np.interp(np.log(E0), np.log(v['E_MPa']), np.log(v['F_at_0p15_N'])))) for (k, v) in curves.items()}
    defaults = {k: {str(e): float(np.exp(np.interp(np.log(e), np.log(v['E_MPa']), np.log(v['F_at_0p15_N'])))) for e in [0.0689, 68.9]} for (k, v) in curves.items()}
    gates = dict(all_native_elements_typed=len(code) == len(tag) and set(np.unique(code)) == {1, 2, 3}, positive_energy=positive, independent_energy=max(energy_errors) < 1e-12, invalid_calls_rejected=all((x['rejected'] for x in bad)), conditional_anchor=all((7 <= x <= 16.2 for x in forces.values())), source_default_rejected=all((not 7 <= x <= 16.2 for v in defaults.values() for x in v.values())))
    report = dict(chain='K13', claim_type='capability', outcome='PASS_SCOPED_COMPILER' if all(gates.values()) else 'FAIL', gates=gates, external_referent=p['external_referent'], resolution='PER_POINT', timescale='SIMULTANEOUS', material_registry=entries, conditional_E_interval_MPa=[float(lower), float(upper)], conditional_E_mid_MPa=E0, conditioned_on='nu=0.45; F1 geometry and loading; not patient posterior', uncertainty='Conditional source-force interpolation of six old FE samples. Interpolation error between samples is UNKNOWN; parameter interval is not a formal clinical/material certificate', calibrated_F_at_0p15_N=forces, source_default_F_at_0p15_N=defaults, source_FE_curves=curves, independent_energy_max_relative_errors=energy_errors, invalid_probes=bad, coverage=dict(elements=len(tag), typed_elements=len(tag), closure_supported_elements=int((tag != 2).sum()), unknown_PDL_elements=int((tag == 2).sum()), unknown_fraction=float((tag == 2).mean()), body_count=int(len(np.unique(body))), resolution='PER_POINT'), rejection=dict(parameter_candidates=3, retained=1, rejected=2, reason='Both inherited PDL defaults fail the published force band; mandatory explicit unsupported mandibular law'), source_FE=dict(path=str(source), sha256=sha(source), new_FE_runs=0), cost=cost(st), physical_chain_status=_release_expand('PARTIAL: K09/K13 body-law binding runs, but @DENTAL_CASE_ID@ mandibular PDL cannot be calibrated from upper-incisor population measurements; global validated physical solve refused'))
    write(L / 'MATERIAL_REGISTRY.json', entries)
    write(L / 'raw/K13.json', report)
    state('K13_COMPLETE', gates, 'K50 conditional fit/freeze and external model-form alarm')
    print(json.dumps(dict(chain='K13', gates=gates, E_interval=report['conditional_E_interval_MPa'], unknown_fraction=report['coverage']['unknown_fraction'])))
if __name__ == '__main__':
    run()
