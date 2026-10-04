"""Geometry preport for PROOF_LANE, carrying explicit refusal of unsupported inputs."""
from dental_release.paths import expand as _release_expand
import json, zipfile, hashlib, importlib.util
from fractions import Fraction as Q
import numpy as np
from occlusion_operator import H, DATA, ZIP, read_stl, dump, REGIONS
PROOF_LANE = __import__('pathlib').Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/engines/3fold-motion-engine/_private/romi_collab/build/PROOF_LANE_OCKLUSION_20261001/code/occlusion.py'))

def enc(x):
    return str(Q(float(x)))

def V(x):
    return [enc(a) for a in x]

def exact_point(xy, tri):
    (a, b, c) = [[Q(float(x)) for x in v] for v in tri]
    (x, y) = map(lambda v: Q(float(v)), xy)
    e = [b[i] - a[i] for i in range(3)]
    f = [c[i] - a[i] for i in range(3)]
    det = e[0] * f[1] - e[1] * f[0]
    s = ((x - a[0]) * f[1] - (y - a[1]) * f[0]) / det
    t = (e[0] * (y - a[1]) - e[1] * (x - a[0])) / det
    return [str(x), str(y), str(a[2] + s * e[2] + t * f[2])]

def main(case=1):
    r = json.loads((H / 'raw/cases_3d' / f'{case:03d}.json').read_text())
    z = zipfile.ZipFile(ZIP)
    (U, _) = read_stl(z, r['upper_member'])
    (L, _) = read_stl(z, r['upper_member'].replace('upper.stl', 'lower.stl'))
    with np.load(DATA / f'{case:03d}_map.npz') as m:
        uf = m['upper_face']
        xy = m['xy']
        ur = m['upper_region']
    with np.load(DATA / f'{case:03d}_clearance.npz') as a:
        dist = a['distance_mm']
        lf = a['nearest_lower_face']
        lp = a['lower_point_mm']
        up = a['upper_point_mm']
    from occlusion_operator import predict_regions
    models = json.loads((H / 'raw/segmentation_model.json').read_text())
    (lr, _, _) = predict_regions(L.mean(1), models['lower'], True)
    ok = np.flatnonzero(dist <= 0.1)
    witnesses = []
    seen = set()
    for i in ok:
        pair = (int(ur[i]), int(lr[lf[i]]))
        if pair in seen:
            continue
        seen.add(pair)
        ut = U[uf[i]]
        lt = L[lf[i]]
        un = np.cross(ut[1] - ut[0], ut[2] - ut[0])
        ln = np.cross(lt[1] - lt[0], lt[2] - lt[0])
        witnesses.append(dict(upper_region=REGIONS[pair[0]], lower_region=REGIONS[pair[1]], upper_source_face=int(uf[i]), lower_source_face=int(lf[i]), upper_triangle_mm=[V(v) for v in ut], lower_triangle_mm=[V(v) for v in lt], upper_source_point_mm=exact_point(xy[i], ut), lower_nearest_point_mm=V(lp[i]), unsigned_clearance_mm=enc(dist[i]), upper_normal_unscaled=V(un), lower_normal_unscaled=V(ln), lower_point_status='IEEE64 closest-point, rational encoding of computed value; exact source-plane incidence not certified'))
    blockers = ['No validated individual FDI/instance segmentation in Bits2Bites', 'No root effective area, rotational support or anatomically reduced PDL matrix in IOS', 'No acquisition force/load direction or jaw moment', 'Contacts are separated source surfaces, not paired coincident central contacts; normals differ', 'PROOF_LANE v1 accepts one central contact/tooth and a specified closure; our source patches do not satisfy that family', 'Proximity field samples at0.2mm; four-region validation and coarse/fine force gates failed']
    port = dict(schema='dental-occlusion-geometric-contact/v1', consumer='PROOF_LANE_OCKLUSION2', status='UNCERTAIN', scope='OBSERVED_PAIRED_GEOMETRY_AND_UNSIGNED_SAMPLED_CLEARANCE_ONLY', case=case, frame='Bits2Bites registered RAS; no jaw transform', length_unit='mm', force_unit='N', modulus_unit='MPa=N/mm2', numeric_encoding='Dimensioned fields rational strings from measured float32 STL and computedIEEE64 values; no new measurement precision claimed', fdi_status='UNKNOWN_NOT_VALIDATED', load=None, pdl=None, contact_band_mm='1/10', source=dict(zip=str(ZIP), upper_member=r['upper_member'], upper_sha256=r['upper_sha256'], lower_sha256=r['lower_sha256']), regions=REGIONS, witnesses=witnesses, conditional_model=dict(support_per_tooth_N_per_mm=['500', '1130'], support_status='Incisor-derived hypothetical prior; posterior UNKNOWN;35scenario envelopes in round3', shares_pp=[enc(v) for v in r['metrics'][1]['mechanics']['shares_pp']], actual_patient_force_N=None), unresolved=blockers, port_proposal='Accept source-pair geometry with explicit region IDs, separate contact points, unsigned clearance and missing physical leaves; never promote to positive PROOF_LANE mechanical certificate')
    spec = importlib.util.spec_from_file_location('proof_lane_existing_occlusion', PROOF_LANE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        module.validate(port)
        accepted = True
        reason = 'ACCEPTED'
    except Exception as e:
        accepted = False
        reason = str(e)
    port['existing_proof_lane_v1_validation'] = dict(accepted=accepted, reason=reason, validator_path=str(PROOF_LANE), validator_sha256=hashlib.sha256(PROOF_LANE.read_bytes()).hexdigest())
    dump(H / 'PORT.json', port)
    dump(H / 'PROOF_LANE_INPUT_BLOCKERS.json', dict(status='UNCERTAIN', blockers=blockers, exact_existing_validator_accepted=accepted, reason=reason))
    print('PROOF_LANE preport exported; existing v1 accepted=', accepted, 'reason=', reason, 'witnesses=', len(witnesses))
if __name__ == '__main__':
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
