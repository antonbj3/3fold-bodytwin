from dental_release.paths import expand as _release_expand
import sys, time, zipfile, resource, gc
import numpy as np
from scipy import ndimage as ndi
from shared import *
if __name__ == '__main__':
    run = Path(sys.argv[1])
    case = sys.argv[2]
    pid = 'P' + str(int(case.rsplit('_', 1)[1]))
    out = run / case
    out.mkdir(exist_ok=True)
    d = DATA / run.name / case
    d.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    x12 = RESULTS / _release_expand('X12')
    sys.path.insert(0, str(x12 / 'code'))
    import archive_io as ar, pair_atlas_r6 as pa, calibrated_atlas as ca, tf2_io as io, consumer_ports as cp
    lookup = {r['name']: r for r in ar.index()}
    n = 'Pulpy3D/' + pid + '/data.nii.gz'
    (img, pulpctsha) = ar.nifti(lookup[n])
    pulpct = np.asanyarray(img.dataobj)
    with zipfile.ZipFile(io.ZIP) as z:
        member = io.ROOT + '/imagesTr/' + case + '_0000.mha'
        raw = z.read(member)
        (ct, sp, h) = io.read_mha_bytes(raw)
        ctsha = __import__('hashlib').sha256(raw).hexdigest()
        matches = pa.exact_match(pulpct, ct)
        assert len(matches) == 1, 'No unique exact patient pairing'
        bad = ct.copy()
        bad.flat[len(bad) // 2] += 1
        fault = not bool(pa.exact_match(pulpct, bad))
        del bad, pulpct, img, ct, raw
        lm = io.ROOT + '/labelsTr/' + case + '.mha'
        raw = z.read(lm)
        (lab, sp, lh) = io.read_mha_bytes(raw)
        labsha = __import__('hashlib').sha256(raw).hexdigest()
    assert all((h.get(k) == lh.get(k) for k in ['TransformMatrix', 'Offset', 'ElementSpacing', 'DimSize']))
    (img, psha) = ar.nifti(lookup['Pulpy3D/' + pid + '/gt_instance.nii.gz'])
    p = pa.transform(np.asanyarray(img.dataobj), *matches[0])
    (sem, cal) = ca.semantic_map(p, lab)
    assert sem is not None
    lut = np.arange(65536, dtype=np.uint16)
    if sem == 'left_right_swap':
        for fdi in pa.FDI:
            lut[fdi] = fdi + 10 if fdi < 40 else fdi - 10
    mapped = lut[p]
    del p, img
    np.savez_compressed(d / 'labels.npz', labels=lab, spacing=sp)
    rows = []
    excluded = []
    controls = []
    crops = {}
    boxes = ndi.find_objects(lab, max_label=48)
    for fdi in sorted(pa.FDI):
        if boxes[fdi - 1] is None:
            excluded.append(dict(fdi=fdi, reason='MISSING_TOOTH'))
            continue
        (ans, err) = pa.tooth_measure(pid, fdi, lab, mapped, sp, lh, True, box=boxes[fdi - 1])
        if err:
            excluded.append(err)
            continue
        (row, ctrl, crop) = ans
        row['semantic_map'] = sem
        row['source_pulpy_fdi'] = fdi if sem == 'identity' else fdi + 10 if fdi < 40 else fdi - 10
        rows.append(row)
        controls.append(ctrl)
        if fdi in [36, 46, 34, 44]:
            file = d / (str(fdi) + '_paired.npz')
            np.savez_compressed(file, **crop)
            crops[fdi] = dict(file=artifact(file), shape=list(crop['tooth'].shape), origin_zyx=crop['crop_origin_zyx'], crown_axis_sign=row['crown_axis_sign'])
    dump(out / 'X12_MEASUREMENTS.json', dict(patient_id=case, paired_pulpy_id=pid, rows=rows, excluded=excluded, controls=controls, crops=crops))
    sp = np.array(sp)
    Dmat = np.array([float(v) for v in lh.get('TransformMatrix', '1 0 0 0 1 0 0 0 1').split()]).reshape(3, 3)
    offset = np.array([float(v) for v in lh.get('Offset', '0 0 0').split()])
    index = np.array([[1, 2, 3], [17, 20, 14]], float)
    world = index * sp[::-1] @ Dmat.T + offset
    back = (world - offset) @ np.linalg.inv(Dmat).T / sp[::-1]
    rt = float(np.max(np.abs(back - index)) * sp.max())
    frame = dict(patient_id=case, unit='mm', array_axes='zyx', canonical_axes='xyz', array_to_world='world_xyz=Offset+TransformMatrix@(index_zyx[::-1]*spacing_zyx[::-1])', offset_xyz_mm=offset, direction=Dmat, spacing_zyx_mm=sp, shape_zyx=lab.shape, roundtrip_error_mm=rt, paired_CT_transform=matches[0], semantic_map=sem, crops=crops, IOS_registration=None, clinical_frame_accuracy='UNKNOWN')
    dump(out / 'FRAME.json', frame)
    pr = load(R / 'PREREG_R1.json')
    ss = pr['numerical_scenarios']
    nvpath = RESULTS / 'NV1_canals/per_case' / f'{case}.json'
    nv = load(nvpath)
    fg = module('full_geometry', RESULTS / 'LANE_X8_GUIDE_NERVE_RISK/full_geometry.py')
    risk = module('risk_operator', RESULTS / 'LANE_X8_GUIDE_NERVE_RISK/risk_operator.py')
    dec = module('decidability', RESULTS / 'LANE_X26_DECIDABILITY/decidability.py')
    profiles = load(RESULTS / 'LANE_X31_CLINICAL_ANSWERS/inputs/guide_profiles.json')['profiles']
    fdi = next((f for f in ss['CBCT_tooth_priority'] if str(f) in nv['teeth']))
    tr = nv['teeth'][str(fdi)]
    entry = np.array(tr['crest_entry'])
    axis = -np.array(tr['axis'])
    axis /= np.linalg.norm(axis)
    centers = np.argwhere(np.isin(lab, [3, 4])) * sp
    geo = fg.voxel_cylinder_bracket(centers, sp, entry, axis, 6.0, 2.0)
    guide = []
    for profile in profiles:
        margin = risk.scenario_margin(profile)
        bound = risk.bound_margin(profile)
        prob = float(risk.scenario_probability(geo['lower_mm'], profile))
        bnd = float(risk.body_bound(geo['lower_mm'], profile))
        scalar = risk.scalar_control(geo['lower_mm'], profile)
        guide.append(dict(guide_type=profile['id'], source=profile['source'], population_resolution='POPULATION', patient_geometry_resolution='PER_SURFACE_REGION', scenario_1pct_margin_mm=margin, moment_bound_1pct_margin_mm=bound, scenario_probability=prob, moment_bound_probability=bnd, scalar_control_error=abs(bnd - scalar), scope='CONDITIONAL_CLOSURE; not patient injury risk', decidability=dec.ScalarPort(margin=geo['lower_mm'] - margin, unit='mm', bias_bound=0.3, coefficient=0, physical_model_validated=False).evaluate(0)))
    remaining = []
    for row in rows:
        for delta in [0, 0.15, 0.3]:
            val = cp.uniform_shell_decision(row['horn_boundary_min_mm'], row['digital_two_surface_radius_mm'], 0.5, 0.5, delta)
            remaining.append(dict(fdi=row['fdi'], resolution='PER_SURFACE_REGION', delta_each_mm=delta, theta_ref='annotation_boundary_mm', **val))
    freeze(out / 'FROZEN_PREDICTIONS.json', dict(patient_id=case, pulpy_id=pid, source=dict(CT_sha256=ctsha, pulpy_CT_sha256=pulpctsha, label_sha256=labsha, pulp_sha256=psha, TF2_archive=io.ZIP, Pulpy_archive=str(ar.ARCHIVE), TF2_label_member=lm), frame=frame, virtual_implant=dict(fdi=fdi, entry_zyx_mm=entry, axis_zyx=axis, length_mm=6, radius_mm=2, geometry=geo), guide=guide, remaining_hard_tissue=remaining, physical_observation='NOT_RUN', labels_used_as_published_geometric_reference=True, clinical_predictions='ABSTAIN'))
    x15 = module('measure_r1', RESULTS / 'LANE_X15_BONE_DEHISCENCE/measure_r1.py')
    x15.DATA = d / 'bone'
    x15.DATA.mkdir(exist_ok=True)
    (bone, bprov, bfiles) = x15.process(case)
    dump(out / 'X15_BONE.json', dict(rows=bone, provenance=bprov, files=bfiles))
    x20 = module('measure', RESULTS / 'LANE_X20_SHORT_IMPLANT_SINUS/measure.py')
    eligible = []
    rejs = []
    sinus = []
    present = set((int(v) for v in np.unique(lab)))
    for f in x20.POSTERIOR:
        if f in present:
            rejs.append(dict(fdi=f, reason='PRESENT_TOOTH_NOT_EDENTULOUS'))
            continue
        if not all((str(q) in nv['teeth'] for q in [f - 1, f + 1])):
            rejs.append(dict(fdi=f, reason='MISSING_NEIGHBOR_FRAME'))
            continue
        eligible.append(f)
        (anchor, a, normal, tan, nh) = x20.frame(case, f)
        local = x20.Local(lab, sp, anchor, a)
        res = x20.trace(local, anchor, a, normal)
        sinus.append(dict(fdi=f, geometry=res, decision=x20.scalar_decide(res), physical='UNKNOWN'))
    if not eligible:
        f = next((f for f in x20.POSTERIOR if str(f) in nv['teeth'] and 'crest_entry' in nv['teeth'][str(f)]), None)
        if f:
            tr = nv['teeth'][str(f)]
            e = np.array(tr['crest_entry'])
            a = -np.array(tr['axis'])
            a /= np.linalg.norm(a)
            t = np.cross(a, [1, 0, 0] if abs(a[0]) < 0.9 else [0, 1, 0])
            t /= np.linalg.norm(t)
            local = x20.Local(lab, sp, e, a)
            probe = x20.trace(local, e, a, t)
        else:
            probe = None
    else:
        probe = None
    dump(out / 'X20_SINUS.json', dict(eligible_sites=eligible, rejected=rejs, rows=sinus, occupied_site_diagnostic=probe, decision='ABSTAIN_NO_MATCHED_EDENTULOUS_SITE' if not eligible else 'CONDITIONAL_GEOMETRY', maxillary_bone_label_present=2 in present, sinus_labels=sorted(present & {5, 6})))
    chosen = next((f for f in ss['CBCT_tooth_priority'] if f in crops), None)
    tooth_export = None
    if chosen:
        q = np.load(crops[chosen]['file']['path'])
        mask = q['tooth']
        localtri = fg.voxel_triangles(mask, sp) - sp[::-1] / 2
        localzyx = localtri[:, :, ::-1] / sp
        if crops[chosen]['crown_axis_sign'] < 0:
            localzyx[:, :, 0] = mask.shape[0] - 1 - localzyx[:, :, 0]
        idx = localzyx + np.array(crops[chosen]['origin_zyx'])
        globaltri = idx[:, :, ::-1] * sp[::-1] @ Dmat.T + offset
        if crops[chosen]['crown_axis_sign'] < 0:
            globaltri = globaltri[:, [0, 2, 1]]
        fg.write_stl(d / 'paired_tooth_world.stl', globaltri)
        tooth_export = artifact(d / 'paired_tooth_world.stl')
        np.savez_compressed(d / 'figure_geometry.npz', tooth_xyz=(np.argwhere(lab == chosen)[::12] * sp)[:, ::-1], pulp_xyz=(np.argwhere(mapped == chosen)[::2] * sp)[:, ::-1], canal_xyz=centers[::40, ::-1], sinus_xyz=(np.argwhere(np.isin(lab, [5, 6]))[::40] * sp)[:, ::-1])
    valid_bone = [v for row in bone for k in ['buccal', 'lingual'] for v in [row.get(k, {})] if v.get('status') == 'MEASURED_ANNOTATION_SURROGATE']
    result = dict(patient_id=case, pulpy_id=pid, claim_type='capability', identity=dict(exact_CT_pair=True, transform=matches[0], mutated_voxel_rejected=fault, voxel_count=int(np.prod(lab.shape)), source_CT_hashes=[ctsha, pulpctsha], semantic_calibration=cal), frame_roundtrip_mm=rt, pulp_teeth=len(rows), pulp_excluded=excluded, selected_fdi=fdi, paired_tooth_export=tooth_export, canal_clearance_mm=[geo['lower_mm'], geo['upper_mm']], canal_bracket_width_mm=geo['gap_mm'], guides=guide, remaining_hard_tissue=remaining, bone=dict(sections=len(bone), measured_sides=len(valid_bone), clearance_range_mm=[min((v['clearance_mm'] for v in valid_bone)), max((v['clearance_mm'] for v in valid_bone))] if valid_bone else None, physical='UNDETERMINED_NO_CEJ_OR_OUTCOME'), sinus=load(out / 'X20_SINUS.json'), thermal=dict(X16='ABSTAIN: population study temperatures not same patient/protocol', X33='ABSTAIN: no heat input, spray boundary, dentin/material calibration for selected patient; prospective simulation in next round'), controls=dict(CT_corruption_rejected=fault, EDT_KDTree_max_error_mm=max((c['max_abs_error_mm'] for c in controls)), EDT_injected_1mm_rejected=all((c['mutated_distance_rejected'] for c in controls)), guide_scalar_max_error=max((r['scalar_control_error'] for r in guide))), external_referent=dict(kind='published_dataset', locator=str(ar.ARCHIVE) + '::Pulpy3D/' + pid + '; ' + str(io.ZIP) + '::' + lm, compared_quantity='Exact CT identity and pulp/outer-tooth/nerve/sinus annotation geometry', refutes_us=False), runtime_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(out / 'CBCT_RESULTS.json', result)
    print('CBCT', case, 'teeth', len(rows), 'clearance', result['canal_clearance_mm'], 'bone sides', len(valid_bone), 'sinus sites', len(eligible), flush=True)
