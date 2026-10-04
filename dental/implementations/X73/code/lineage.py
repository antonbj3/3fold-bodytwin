from common import *
import copy, resource
sys.path.insert(0, str(X8))
from full_geometry import voxel_cylinder_bracket

def main():
    check_frozen()
    start = time.perf_counter()
    patients = json.loads((ROOT / 'raw/PATIENTS.json').read_text())
    sites = json.loads((ROOT / 'raw/SITES.json').read_text())
    byid = {p['patient']: p for p in patients}
    groups = {}
    for p in patients:
        groups.setdefault(p['image_identity']['tf1']['sha256'], []).append(p['patient'])
    canonical = {p: min(ps, key=lambda k: int(k[1:])) for ps in groups.values() for p in ps}
    duplicates = []
    new = []
    for ps in groups.values():
        if len(ps) > 1:
            duplicates.append(dict(canonical=canonical[ps[0]], aliases=sorted(ps), image_sha256=byid[ps[0]]['image_identity']['tf1']['sha256']))
    for source in sites:
        r = copy.deepcopy(source)
        p = byid[r['patient']]
        alias_ids = groups[p['image_identity']['tf1']['sha256']]
        r['image_group'] = canonical[r['patient']]
        r['within_id_union'] = copy.deepcopy(r['union'])
        r['alias_source_bindings'] = []
        for alias in alias_ids:
            if alias == r['patient']:
                continue
            ap = byid[alias]
            assert ap['shape'] == p['shape']
            assert ap['image_identity']['tf1']['sha256'] == p['image_identity']['tf1']['sha256']
            r['alias_source_bindings'].append(dict(patient=alias, source_bindings=ap['source_bindings'], point_artifact=ap['point_artifact'], image_sha256=ap['image_identity']['tf1']['sha256']))
            assert sha(ap['point_artifact']['path']) == ap['point_artifact']['sha256']
            with np.load(ap['point_artifact']['path']) as d:
                for sb in ap['source_bindings']:
                    key = alias + '::' + sb['dataset']
                    points = d[sb['dataset'] + '_voxels_zyx'] * SP
                    pose = r['pose']
                    x = voxel_cylinder_bracket(points, np.full(3, SP), np.array(pose['entry_zyx_mm']), np.array(pose['axis_zyx']), pose['length_mm'], pose['radius_mm'])
                    br = {k: x[k] for k in ['lower_mm', 'upper_mm', 'gap_mm', 'active_voxel_boxes', 'witness']}
                    assert br['gap_mm'] <= 0.0001
                    r['distances'][key] = br
                    r['classes'][key] = classify(br['lower_mm'], br['upper_mm'])
        lo = min((x['lower_mm'] for x in r['distances'].values()))
        hi = min((x['upper_mm'] for x in r['distances'].values()))
        r['union'] = dict(lower_mm=lo, upper_mm=hi, class_2mm=classify(lo, hi))
        ref = r['distances']['tf2']
        r['loss_lower_mm'] = max(0.0, ref['lower_mm'] - hi)
        r['loss_upper_mm'] = max(0.0, ref['upper_mm'] - lo)
        new.append(r)
    dump(ROOT / 'raw/LINEAGE_SITES.json', new)
    dump(ROOT / 'raw/IMAGE_GROUPS.json', dict(canonical=canonical, groups=groups, duplicate_groups=duplicates, physical_image_groups=len(groups), dataset_patient_ids=len(patients), duplicate_group_ids=sum((len(p['aliases']) for p in duplicates)), cross_alias_query_sites=sum((bool(r['alias_source_bindings']) for r in new)), cross_alias_extra_2mm_reversals=sum((r['within_id_union']['class_2mm'] == 'ABOVE' and r['union']['class_2mm'] == 'BELOW' for r in new)), cross_alias_distance_reductions=sum((r['within_id_union']['upper_mm'] > r['union']['upper_mm'] + 1e-10 for r in new)), cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    print('Image groups', len(groups), 'duplicate groups', duplicates, flush=True)
if __name__ == '__main__':
    main()
