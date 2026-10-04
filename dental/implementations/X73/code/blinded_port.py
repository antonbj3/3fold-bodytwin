"""Frozen native-grid ROIs for future independent contours. No human data invented."""
from common import *
import copy, resource
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

def validate_reader(r, roi):
    assert r['roi_id'] == roi['roi_id'] and r['roi_sha256'] == roi['roi_sha256']
    assert r['frame_id'] == roi['frame_id'] and r['unit'] == 'mm' and (r['spacing_mm'] == [0.3, 0.3])
    assert r['blinded'] is True and r['preconsensus'] is True and r['reader_id']
    assert Path(r['mask_path']).exists() and sha(r['mask_path']) == r['mask_sha256']
    a = np.load(r['mask_path'], allow_pickle=False)
    assert a.shape == (72, 72) and a.dtype == np.dtype(bool) and a.any()
    return a

def compare_readers(records, rois):
    if not records:
        return dict(status='UNKNOWN_NO_INDEPENDENT_READER_DATA', measured_rois=0, independent_spread_mm=None)
    out = []
    for rid in sorted(set((r['roi_id'] for r in records))):
        rr = [r for r in records if r['roi_id'] == rid]
        assert len(rr) >= 2 and len({r['reader_id'] for r in rr}) == len(rr)
        masks = [validate_reader(r, rois[rid]) for r in rr]
        values = []
        for (i, a) in enumerate(masks):
            pa = np.argwhere(a & ~ndi.binary_erosion(a)) * 0.3
            for b in masks[i + 1:]:
                pb = np.argwhere(b & ~ndi.binary_erosion(b)) * 0.3
                values.extend(cKDTree(pb).query(pa, workers=1)[0])
                values.extend(cKDTree(pa).query(pb, workers=1)[0])
        out.append(dict(roi_id=rid, reader_count=len(rr), resolution='PER_POINT', time_scale='SIMULTANEOUS', spread=summary(values), limits='2D contour repeat spread only; declared blinding/reader independence require external review; not 3D nerve anatomy or physical guide validation'))
    return dict(status='LOCAL_CONTOURS_MEASURED_PENDING_EXTERNAL_REVIEW', measured_rois=len(out), rois=out)

def main():
    check_frozen()
    start = time.perf_counter()
    sites = json.loads((ROOT / 'raw/LINEAGE_SITES.json').read_text())
    witnesses = [r for r in sites if r['classes']['tf2'] == 'ABOVE' and r['union']['class_2mm'] == 'BELOW']
    bycase = {}
    for r in witnesses:
        bycase.setdefault(r['case'], []).append(r)
    base = DATA / 'blinded_rois'
    base.mkdir(parents=True, exist_ok=True)
    rois = []
    predictions = []
    with zipfile.ZipFile(TF1) as z:
        for (case, rr) in sorted(bycase.items()):
            patient = rr[0]['patient']
            (im, ih) = npy(z, f'ToothFairy_Dataset/Dataset/{patient}/data.npy')
            for r in sorted(rr, key=lambda r: r['fdi']):
                rid = f'ROI_{len(rois) + 1:03d}'
                pt = np.array(r['distances']['tf2']['witness']['witness_zyx_mm'])
                voxel = np.rint(pt / 0.3).astype(int)
                zi = int(voxel[0])
                low = voxel[1:] - 36
                high = low + 72
                assert 0 <= zi < im.shape[0] and np.all(low >= 0) and np.all(high <= im.shape[1:]), 'ROI lies outside scan'
                image = im[zi, low[0]:high[0], low[1]:high[1]].copy()
                path = base / (rid + '.npy')
                np.save(path, image, allow_pickle=False)
                frame = dict(axis='native axial z-index fixed', z_index=zi, origin_yx_indices=low.tolist(), spacing_mm=[0.3, 0.3], voxel_center_world_zyx_mm=[zi * 0.3, float(low[0] * 0.3), float(low[1] * 0.3)], shape=[72, 72])
                frame_id = hashlib.sha256(json.dumps(frame, sort_keys=True).encode()).hexdigest()
                rois.append(dict(roi_id=rid, roi_path=str(path), roi_sha256=sha(path), frame_id=frame_id, frame=frame, reader_input_has_reference_masks=False, unit='mm', resolution='PER_POINT', time_scale='SIMULTANEOUS'))
                masks = {}
                for sb in r['source_bindings']:
                    with np.load(r['point_artifact']['path']) as d:
                        pp = d[sb['dataset'] + '_voxels_zyx']
                        sel = pp[(pp[:, 0] == zi) & np.all(pp[:, 1:] >= low, 1) & np.all(pp[:, 1:] < high, 1)][:, 1:] - low
                        m = np.zeros((72, 72), bool)
                        m[tuple(sel.T)] = True
                        masks[sb['dataset']] = dict(contour_occupancy_sha256=arrsha(m), foreground_pixels=int(m.sum()))
                predictions.append(dict(roi_id=rid, patient_dataset_alias=patient, image_group=r['image_group'], case=case, fdi=r['fdi'], segment=r['segment'], source_image_npy_sha256=ih, source_mask_bindings=r['source_bindings'], alias_source_bindings=r['alias_source_bindings'], roi_reference_mask_predictions=masks, source_fixed_cylinder_clearance_mm=r['distances']['tf2'], observed_union=r['union'], loss_upper_mm=r['loss_upper_mm'], interpretation='3D query prediction frozen separately. A single axial contour measures a local 2D repeat quantity and cannot validate the global 3D cylinder gap.', roi_sha256=sha(path), frame_id=frame_id))
            del im
    dump(ROOT / 'raw/BLINDED_ROI_INPUT.json', dict(roi_count=len(rois), rois=rois, reader_protocol='Two different readers, blinded to each other and all predicted masks; retain each original contour before consensus; use opaque IDs, no personal data.', selection='Source-exposed selection on actual reversals, not a random patient population sample'))
    predfile = ROOT / 'FROZEN_PREDICTIONS_R5_LOCAL_CONTOURS.json'
    payload = dict(prereg_sha256=sha(ROOT / 'PREREG_R5_BLINDED_PORT.json'), predictions=predictions, physical_measurements_performed=False)
    if predfile.exists():
        old = json.loads(predfile.read_text())
        assert old['prereg_sha256'] == payload['prereg_sha256'] and old['predictions'] == predictions, 'Frozen predictions drift; new prereg/version required'
    else:
        dump(predfile, payload | dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        predfile.with_suffix('.sha256').write_text(sha(predfile) + '  ' + predfile.name + '\n')
    template = ROOT / 'READER_SUBMISSION_TEMPLATE.json'
    if not template.exists():
        dump(template, dict(schema='X73-local-contour-readers-v1', frozen_predictions_sha256=sha(predfile), records=[], required_record_fields=['roi_id', 'roi_sha256', 'frame_id', 'reader_id', 'blinded', 'preconsensus', 'unit', 'spacing_mm', 'mask_path', 'mask_sha256'], note='mask .npy boolean72x72, independent reader pairs per ROI; no synthetic fixtures submitted as measurements'))
    assert rois
    fixture = base / 'CONTROL_ONLY_MASK.npy'
    np.save(fixture, np.eye(72, dtype=bool), allow_pickle=False)
    row = dict(roi_id=rois[0]['roi_id'], roi_sha256=rois[0]['roi_sha256'], frame_id=rois[0]['frame_id'], reader_id='FIXTURE_A', blinded=True, preconsensus=True, unit='mm', spacing_mm=[0.3, 0.3], mask_path=str(fixture), mask_sha256=sha(fixture))
    valid = validate_reader(row, rois[0]).shape == (72, 72)
    probes = {}
    for (key, val) in [('roi_sha256', '0' * 64), ('frame_id', '0' * 64), ('unit', 'pixels'), ('spacing_mm', [3.0, 0.3]), ('blinded', False), ('preconsensus', False)]:
        bad = copy.deepcopy(row)
        bad[key] = val
        try:
            validate_reader(bad, rois[0])
            probes[key + '_poison_rejected'] = False
        except AssertionError:
            probes[key + '_poison_rejected'] = True
    try:
        compare_readers([row, row], {r['roi_id']: r for r in rois})
        probes['duplicate_reader_rejected'] = False
    except AssertionError:
        probes['duplicate_reader_rejected'] = True
    measured = compare_readers([], {})
    assert valid and all(probes.values()) and (measured['status'] == 'UNKNOWN_NO_INDEPENDENT_READER_DATA')
    artifact = dict(path=str(predfile), sha256=sha(predfile))
    result = dict(claim_type='capability', status='BLINDED_PORT_READY; INDEPENDENT_SPREAD_UNKNOWN', rois=len(rois), frozen_predictions=artifact, reader_result=measured, contract_controls=probes, positive_fixture_accepted=True, control_fixture_is_measurement=False, cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, scratch_bytes=diskcheck()), external_referent=dict(kind='published_dataset', locator=str(TF1) + ' and raw/LINEAGE_SITES.json exact source hashes', compared_quantity='Published image ROI and digital reference prediction only; blinded wall repeat has no external data yet', refutes_us=False))
    dump(ROOT / 'raw/R5_RESULT.json', result)
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
