from dental_release.paths import expand as _release_expand
import copy, hashlib, sys, time, zipfile
from scipy import ndimage as ndi
from common import *
from region_field import *

def run(out, data_dir):
    t0 = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    freeze(out / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', prereg=source(ROOT / 'PREREG_C2_REGION.json'), source_crop_predictions=[dict(source=source(PARENT_DATA / pid / '36_paired.npz'), expected_raw_tooth_mismatches=0, expected_raw_pulp_mismatches=0) for pid in ['ToothFairy2P_048', 'ToothFairy2P_242']], predicted_query_invariants=dict(max_frame_error_mm=1e-09, max_source_height_error_mm=1e-06, owner_mismatches=0), physical_measurement='NOT_RUN', prior_exposure='Parent source-derived geometries known; descriptive reference replay, not blinded prediction'))
    descriptors = []
    fields = []
    checks = []
    rows = []
    members = []
    rejected = []
    source_sufficiency = []
    manifest = RESULTS / 'LANE_X7_BITE2TEXT/raw/DATA_MANIFEST.json'
    source(manifest)
    archive = load(manifest)['zip']
    x21 = RESULTS / _release_expand('X21')
    expansion = []
    with zipfile.ZipFile(archive) as zfile:
        available = set(zfile.namelist())
        for p in sorted((x21 / 'raw/cases').glob('*.json')):
            case = p.stem
            if case == _release_expand('@DENTAL_CASE_ID@'):
                continue
            m = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bite2Text')) / case / 'labels+landmarks.json'
            good = m.is_file() and all((case + '/ios/ios_' + j + '.stl' in available for j in ['upper', 'lower']))
            if good:
                expansion.append(case)
            else:
                rejected.append(dict(patient_id=case, reason='Missing local labels or raw paired arches'))
            if len(expansion) == 2:
                break
        for case in [_release_expand('@DENTAL_CASE_ID@')] + expansion:
            pid = 'Bite2Text_' + case
            if case == _release_expand('@DENTAL_CASE_ID@'):
                metapath = PARENT_DATA / pid / 'labels+landmarks.json'
                framepath = PARENT_RUN / pid / 'FRAME.json'
                fr = load(framepath)
                R = np.array(fr['local_to_native_rotation'])
                center = np.array(fr['native_origin_xyz_mm'])
            else:
                metapath = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bite2Text')) / case / 'labels+landmarks.json'
                framepath = x21 / 'raw/cases' / (case + '.json')
                fr = load(framepath)['frame']
                R = np.array(fr['R_columns'])
                center = np.array(fr['center_xyz_mm'])
            source(metapath)
            source(framepath)
            meta = load(metapath)
            A = np.eye(4)
            A[:3, :3] = R.T
            A[:3, 3] = -R.T @ center
            frameid = pid + ':X11_common_mm'
            for jaw in ['upper', 'lower']:
                member = case + '/ios/ios_' + jaw + '.stl'
                blob = zfile.read(member)
                memberhash = hashlib.sha256(blob).hexdigest()
                m = meta['arches'][jaw]
                if memberhash != m['input_sha256']:
                    raise ValueError('IOS_SOURCE_MEMBER_HASH')
                labelpath = Path(m['labels']['path'])
                source(labelpath)
                if sha(labelpath) != m['labels']['sha256']:
                    raise ValueError('IOS_LABEL_HASH')
                tri = read_stl(blob)
                (v, inv) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
                lab = np.load(labelpath)['labels']
                if len(v) != len(lab):
                    raise ValueError('IOS_LABEL_INDEX')
                vertex_labels = lab[inv.reshape(-1, 3)]
                pure = np.all(vertex_labels == vertex_labels[:, 0, None], axis=1)
                memberrow = dict(archive=archive, member=member, sha256=memberhash, bytes=len(blob))
                members.append(memberrow)
                mfield = dict(field_id=pid + '_' + jaw, patient_id=pid, frame_id=frameid, modality='IOS', backend='IOS_SURFACE', jaw=jaw, native_to_patient=A, source_sha256=memberhash, source=memberrow, label_source=source(labelpath), region_semantics='Predicted FDI only for unanimous facet vertices; 0 unassigned soft-tissue/gingiva label; mixed facets UNKNOWN', distance_semantics='SIGNED_VERTICAL_GAP_NOT_SDF', anatomical_accuracy='UNKNOWN_FDI', absent_modalities=['CBCT'], quantity_debt='Loaded pose/region forces, scan error, independent FDI, preparation')
                (field, path) = save_field(out, data_dir, mfield, dict(triangles=tri, facet_vertex_labels=vertex_labels))
                fields.append(field)
                descriptors.append(str(path))
                pts = field.tri[np.linspace(0, len(tri) - 1, 129, dtype=int)].mean(1)
                q = field.query(pts, patient_id=pid, frame_id=frameid)
                sys.path.insert(0, str(RESULTS / 'LANE_X18_CROWN_ANTAGONIST/code'))
                import geometry as g
                (h, face) = g.query_height(field.tri, pts[:, :2], jaw == 'upper')
                known = np.isfinite(h)
                err = float(np.max(np.abs(h[known] - q['height_mm'][known])))
                original_labels = vertex_labels[q['source_address'][q['known_geometry']]].astype(np.int64)
                own = np.where(np.all(original_labels == original_labels[:, 0, None], axis=1), original_labels[:, 0], -1)
                mismatch = int(np.count_nonzero(own != q['region_id'][q['known_geometry']]))
                rt = float(np.max(np.abs(transform(transform(v[:100], A), np.linalg.inv(A)) - v[:100])))
                checks += [check(pid + '_' + jaw + '_source_height', err < 1e-06, err + 0.01 > 1e-06, dict(max_error_mm=err, injection='Add 0.01 mm to predicted heights')), check(pid + '_' + jaw + '_source_owner', mismatch == 0, bool(np.any(own != own + 1)), dict(mismatches=mismatch, injection='Increment returned source owner')), check(pid + '_' + jaw + '_frame', rt < 1e-09, rt + 0.01 > 1e-09, dict(roundtrip_mm=rt, injection='Add 0.01 mm to inverse transform'))]
                for (key, value) in [('patient_id', 'WRONG_PATIENT'), ('frame_id', 'WRONG_FRAME'), ('unit', 'um')]:
                    kwargs = dict(patient_id=pid, frame_id=frameid)
                    kwargs[key] = value
                    checks.append(check(pid + '_' + jaw + '_' + key, True, rejects(field.query, pts, **kwargs)))
                rows.append(dict(patient_id=pid, field_id=field.meta['field_id'], modality='IOS', source_faces=len(tri), mixed_faces=int((~pure).sum()), unassigned_pure_faces=int((pure & (vertex_labels[:, 0] == 0)).sum()), owner_known_faces=int((pure & (vertex_labels[:, 0] > 0)).sum()), rejected_fraction=float(np.mean(~pure | (vertex_labels[:, 0] == 0))), exclusion_reason='No unanimous positive source FDI; geometry retained with UNKNOWN owner', tested_points=len(pts), maximum_source_height_error_mm=err, roundtrip_error_mm=rt))
    io = imported('r2_tf2_io', RESULTS.parent / 'cells/geometry/tf2_io.py')
    ar = imported('r2_archive_io', RESULTS / 'LANE_X12_PULPY3D/code/archive_io.py')
    lookup = {r['name']: r for r in ar.index()}
    for pid in ['ToothFairy2P_048', 'ToothFairy2P_242']:
        framepath = PARENT_RUN / pid / 'FRAME.json'
        source(framepath)
        fr = load(framepath)
        cropfile = PARENT_DATA / pid / '36_paired.npz'
        source(cropfile)
        crop = np.load(cropfile)
        labmember = io.ROOT + '/labelsTr/' + pid + '.mha'
        with zipfile.ZipFile(io.ZIP) as zfile:
            blob = zfile.read(labmember)
        (rawlab, sp, hdr) = io.read_mha_bytes(blob)
        lh = hashlib.sha256(blob).hexdigest()
        o = crop['crop_origin_zyx']
        shape = np.array(crop['tooth'].shape)
        sl = tuple((slice(int(a), int(a + b)) for (a, b) in zip(o, shape)))
        tooth = rawlab[sl] == 36
        pulpname = 'Pulpy3D/P' + str(int(pid.rsplit('_', 1)[1])) + '/gt_instance.nii.gz'
        (img, ph) = ar.nifti(lookup[pulpname])
        pul = np.asanyarray(img.dataobj)
        (perm, flips) = fr['paired_CT_transform']
        pul = np.transpose(pul, perm)
        for (axis, flip) in enumerate(flips):
            if flip:
                pul = np.flip(pul, axis)
        original_pulp_fdi = 46 if fr['semantic_map'] == 'left_right_swap' else 36
        pulp = pul[sl] == original_pulp_fdi
        crown_sign = int(fr['crops']['36']['crown_axis_sign'])
        if crown_sign < 0:
            tooth = tooth[::-1].copy()
            pulp = pulp[::-1].copy()
        tm = int(np.count_nonzero(tooth != crop['tooth']))
        pm = int(np.count_nonzero(pulp != crop['pulp']))
        mutated = tooth.copy()
        mutated.flat[0] = ~mutated.flat[0]
        checks.append(check(pid + '_raw_source_masks', tm == pm == 0, not np.array_equal(mutated, crop['tooth']), dict(tooth_mismatches=tm, pulp_mismatches=pm, injection='Flip one crop owner voxel', injected_mismatches=int(np.count_nonzero(mutated != crop['tooth'])))))
        D = np.array(fr['direction'])
        sp = np.array(sp)
        offset = np.array(fr['offset_xyz_mm'])
        M = D @ np.diag(sp[::-1]) @ np.eye(3)[::-1]
        J = np.diag([crown_sign, 1, 1])
        shift = np.array([shape[0] - 1 if crown_sign < 0 else 0, 0, 0])
        A = np.eye(4)
        A[:3, :3] = M @ J
        A[:3, 3] = offset + M @ (o + shift)
        owner = np.where(pulp, 2, np.where(tooth, 1, 0)).astype(np.uint8)
        edt = (ndi.distance_transform_edt(~tooth, sampling=sp) - ndi.distance_transform_edt(tooth, sampling=sp)).astype(np.float32)
        src = dict(archive=io.ZIP, member=labmember, sha256=lh, bytes=len(blob))
        members.append(src)
        members.append(dict(archive=str(ar.ARCHIVE), member=pulpname, sha256=ph, bytes=lookup[pulpname]['uncompressed']))
        meta = dict(field_id=pid + '_tooth36', patient_id=pid, frame_id=pid + ':CBCT_world_mm', modality='CBCT', backend='CBCT_VOXELS', native_to_patient=A, source_sha256=digest([lh, ph]), source=[members[-2], members[-1]], crop_origin_zyx=o, spacing_zyx_mm=sp, crown_axis_sign=crown_sign, index_to_full_zyx=dict(matrix=J, offset=o + shift), region_semantics={'0': 'outside selected tooth in crop; not certified anatomical air', '1': 'total_hard_tissue_36', '2': 'pulp_36'}, distance_semantics='NEAREST_GRID_CENTER_EDT_NOT_CONTINUOUS_SDF', absent_modalities=['IOS'], quantity_debt='DEJ, actual preparation, image-to-anatomy bound and registered IOS')
        (field, path) = save_field(out, data_dir, meta, dict(owner=owner, tooth_center_edt_mm=edt))
        fields.append(field)
        descriptors.append(str(path))
        index = np.indices(owner.shape).reshape(3, -1).T
        points = transform(index, A)
        query = field.query(points, patient_id=pid, frame_id=meta['frame_id'])
        mismatches = int(np.count_nonzero(query['region_id'] != owner.ravel()))
        unknown = int(np.count_nonzero(~query['known_region']))
        rt = float(np.max(np.abs(transform(points, np.linalg.inv(A)) - index)) * np.max(sp))
        mutated = query['region_id'].copy()
        mutated[0] += 1
        checks.append(check(pid + '_all_voxel_centres', mismatches == unknown == 0, not np.array_equal(mutated, owner.ravel()), dict(mismatches=mismatches, unknown=unknown, injection='Change one owner returned by query', injected_mismatches=int(np.count_nonzero(mutated != owner.ravel())))))
        checks.append(check(pid + '_coordinate_roundtrip', rt < 1e-09, rt + 0.3 > 1e-09, dict(error_mm=rt, injection='One-voxel crop-origin error')))
        wrong_meta = copy.deepcopy(field.meta)
        wrong_A = np.eye(4)
        wrong_A[:3, :3] = M
        wrong_A[:3, 3] = offset + M @ o
        wrong_meta['native_to_patient'] = wrong_A
        wrong = RegionField(wrong_meta, field.arrays).query(points, patient_id=pid, frame_id=meta['frame_id'])
        difference = int(np.count_nonzero(query['region_id'] != wrong['region_id']))
        source_sufficiency.append(dict(patient_id=pid, summary='Same owner array/histogram, voxel-volume sum and sorted EDT; only coordinate map differs', summary_identity_error=0, downstream_world_region_mismatches=difference, queries=len(index), fraction_changed=difference / len(index), minimum_extension='Signed local-z affine with crown_axis_sign and crop extent; crop origin alone is insufficient', resolution='PER_POINT', external_referent=dict(kind='published_dataset', locator=io.ZIP + '::' + labmember, compared_quantity='Source-coordinate region label after local-crop z reversal', refutes_us=difference == 0)))
        query_owner = query['region_id'].reshape(owner.shape)
        from_carrier_tooth = query_owner > 0
        from_carrier_pulp = query_owner == 2
        dist = ndi.distance_transform_edt(from_carrier_tooth, sampling=sp)
        retained = from_carrier_tooth & (dist > 0.65)
        orig_dist = ndi.distance_transform_edt(crop['tooth'], sampling=sp)
        ref_retained = crop['tooth'] & (orig_dist > 0.65)
        exportdata = dict(tooth=from_carrier_tooth, pulp=from_carrier_pulp, retained_after_virtual_erosion=retained, spacing=sp, crop_origin_zyx=o)
        maskpath = data_dir / (pid + '_X33_inputs.npz')
        np.savez_compressed(maskpath, **exportdata)
        checks.append(check(pid + '_X33_mask_handover', np.array_equal(retained, ref_retained) and np.array_equal(from_carrier_pulp, crop['pulp']), not np.array_equal(~retained, ref_retained), dict(injection='Complement retained mask', time_scale='HANDOVER', source=source(maskpath))))
        rows.append(dict(patient_id=pid, field_id=meta['field_id'], modality='CBCT', tested_points=len(index), mismatches=mismatches, crown_axis_sign=crown_sign, rejected_queries=unknown, maximum_roundtrip_error_mm=rt, raw_mask_mismatches=tm + pm, X33_input=source(maskpath), erosion_mm=0.65, erosion_status='PHENOMENOLOGICAL; inherited virtual preparation, no physical cut'))
        del rawlab, pul, img, points, index, query, blob
    base = next((f for f in fields if f.meta['backend'] == 'CBCT_VOXELS'))
    angle = 0.23
    (c, s) = (np.cos(angle), np.sin(angle))
    T = np.array([[c, -s, 0, 7], [s, c, 0, -3], [0, 0, 1, 2], [0, 0, 0, 1.0]])
    newmeta = copy.deepcopy(base.meta)
    newmeta['frame_id'] += ':synthetic_reexpression'
    newmeta['native_to_patient'] = T @ base.A
    moved = RegionField(newmeta, base.arrays)
    idx = np.array([[2, 3, 4], [10, 3, 4], [2, 12, 4], [2, 3, 15], [20, 15, 10]], float)
    p = transform(idx, base.A)
    q = transform(p, T)
    bridge = dict(patient_id=base.meta['patient_id'], unit='mm', source_frame_id=base.meta['frame_id'], target_frame_id=moved.meta['frame_id'], source_source_sha256=base.meta['source_sha256'], target_source_sha256=moved.meta['source_sha256'], source_to_target=T, identity_evidence=dict(locator='same source field, synthetic rigid re-expression; not an IOS-CBCT pair', sha256=base.meta['source_sha256']), fit_landmark_ids=['fit_a', 'fit_b', 'fit_c'], heldout_landmarks=dict(ids=['v1', 'v2', 'v3', 'v4', 'v5'], source_xyz_mm=p, target_xyz_mm=q), heldout_acceptance_mm=1e-09, uniform_registration_error_bound_mm=None, time_scale='SIMULTANEOUS', acquisition_state='Identical source state synthetically re-expressed; future real pairing must record acquisition times and anatomy changes', evidence_kind='our_own_fixture')
    query_index = np.concatenate([idx, np.array([np.argwhere(base.arrays['owner'] == label)[0] for label in [0, 1, 2]])])
    query_points = transform(query_index, base.A)
    joined = registered_query(bridge, base, moved, query_points)
    reference = base.query(query_points, patient_id=base.meta['patient_id'], frame_id=base.meta['frame_id'])
    bridge_checks = []
    for (name, key, value) in [('patient', 'patient_id', 'wrong'), ('unit', 'unit', 'um'), ('hash', 'source_source_sha256', '0' * 64), ('frame', 'target_frame_id', 'wrong')]:
        bad = copy.deepcopy(bridge)
        bad[key] = value
        bridge_checks.append(check('bridge_' + name, True, rejects(validate_bridge, bad, base, moved)))
    for (name, matrix) in [('reflection', np.diag([-1.0, 1, 1, 1])), ('scale', np.diag([1000.0, 1000, 1000, 1])), ('translation', T + np.array([[0, 0, 0, 0.1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]))]:
        bad = copy.deepcopy(bridge)
        bad['source_to_target'] = matrix
        bridge_checks.append(check('bridge_' + name, True, rejects(validate_bridge, bad, base, moved)))
    bad = copy.deepcopy(bridge)
    bad['fit_landmark_ids'] = ['v1']
    bridge_checks.append(check('bridge_independent_landmarks', True, rejects(validate_bridge, bad, base, moved)))
    bad = copy.deepcopy(bridge)
    bad['evidence_kind'] = 'published_dataset'
    bridge_checks.append(check('bridge_real_identity_artifact', True, rejects(validate_bridge, bad, base, moved)))
    different_patient = next((ff for ff in fields if ff.meta['patient_id'] != base.meta['patient_id']))
    bridge_checks.append(check('bridge_actual_cross_patient', True, rejects(validate_bridge, bridge, base, different_patient)))
    checks += bridge_checks
    checks.append(check('bridge_query_invariance', np.array_equal(joined['region_id'], reference['region_id']) and set(reference['region_id']) == {0, 1, 2}, not np.array_equal(joined['region_id'] + 1, reference['region_id']), dict(injection='Increment joined owner', evidence_kind='our_own_fixture', region_ids=reference['region_id'])))
    dump(out / 'REGISTRATION_FIXTURE.json', bridge)
    dump(out / 'REGISTRATION_FIXTURE_TARGET_FIELD.json', moved.meta)
    dump(out / 'REGISTRATION_FIXTURE_POINTS.json', query_points)
    template = {k: None for k in bridge}
    template.update(schema='Patient360-IOS-CBCT-registration-v2', dataset_locator='https://doi.org/10.6084/m9.figshare.26965903.v3', dataset_status='METADATA_ONLY_NOT_DOWNLOADED', evidence_kind='published_dataset', unit='mm', uniform_registration_error_bound_mm=None)
    dump(out / 'FUTURE_IOS_CBCT_PAIR_TEMPLATE.json', template)
    aa = np.array([1, 2], np.uint8)
    bb = aa[::-1].copy()
    sdfa = np.array([-1.0, 1.0])
    sdfb = sdfa[::-1].copy()
    suff = dict(summary_histogram_identity_error=int(np.max(np.abs(np.bincount(aa) - np.bincount(bb)))), sorted_distance_summary_identity_error_mm=float(np.max(np.abs(np.sort(sdfa) - np.sort(sdfb)))), downstream_fixed_point_region_ids=[int(aa[0]), int(bb[0])], downstream_difference=1, minimum_extension='Retain spatial source address and frame with every region/distance value.', evidence_kind='our_own_fixture', resolution='PER_POINT')
    for field in fields:
        missing = 'CBCT' if field.meta['modality'] == 'IOS' else 'IOS'
        ans = query_modality(fields, missing, np.zeros((1, 3)), field.meta['patient_id'], field.meta['frame_id'])
        checks.append(check(field.meta['field_id'] + '_missing_modality', ans['status'] == 'ABSTAIN_MISSING_MODALITY', ans['region_id'] != [0], dict(injection='Relabel missing modality as background 0')))
    dump(out / 'FIELD_MANIFEST.json', dict(descriptors=[source(p) for p in descriptors]))
    result = dict(claim_type='capability', descriptors=descriptors, rows=rows, patients=sorted({r['patient_id'] for r in rows}), expansion_cases=expansion, selection_rejected=rejected, source_members=members, checks=checks, registration=joined['registration'], actual_cross_modal_pairs=0, sufficient_summary_test=suff, source_sufficiency=source_sufficiency, external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bite2text/; https://ditto.ing.unimore.it/toothfairy2/; exact local members in source_members', compared_quantity='Original source facet heights/owners and all tooth36/pulp crop voxels', refutes_us=not all((x['pass_'] for x in checks))), runtime_seconds=time.perf_counter() - t0)
    dump(out / 'REGION_RESULTS.json', result)
    return (result, fields)
