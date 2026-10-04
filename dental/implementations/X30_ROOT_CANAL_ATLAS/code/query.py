"""One tooth query with branch-level source frame and uncertainty preserved."""
from dental_release.paths import expand as _release_expand
import argparse, csv, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X30'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='P1')
    ap.add_argument('--fdi', type=int, default=36)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--point-clearance', type=Path)
    args = ap.parse_args()

    def match(r):
        return r['case'] == args.case and int(r['fdi']) == args.fdi
    teeth = [r for r in csv.DictReader((ROOT / 'raw/R1_teeth.csv').open()) if match(r)]
    if len(teeth) != 1:
        raise SystemExit('Tooth absent/rejected; inspect upstream exclusions. No substituted tooth.')
    topo = [r for r in csv.DictReader((ROOT / 'raw/R2_topology.csv').open()) if match(r)]
    budgets = [r for r in csv.DictReader((ROOT / 'raw/R3_measurement_budget.csv').open()) if match(r)]
    paths = [r for r in map(json.loads, (DATA / 'R1_paths.jsonl').open()) if match(r)]
    for r in paths:
        r['rootward_graph_endpoint_mm'] = r.pop('apical_annotation_endpoint_mm')
        r['graph_endpoint_IAN_boundary_mm'] = r.pop('annotation_endpoint_IAN_boundary_mm', None)
        r['anatomical_apex_IAN_mm'] = None
    output = {'case': args.case, 'canonical_FDI': args.fdi, 'claim_type': 'capability', 'tooth': teeth[0], 'topology_sensitivity': topo, 'branch_window_metrology_budget': budgets, 'partial_annotation_branches': paths, 'coordinate_contract': 'Points in TF2 full-volume voxel-index frame(zyx)*0.3mm; add MHA origin/direction to obtain world frame. Raw graph endpoint is10%rootward-region cut,not apex.', 'status': 'ANNOTATION_GEOMETRY_ONLY', 'clinical_working_length_mm': None, 'anatomical_apex_distances_mm': None, 'Schneider_Pruett_landmarks': None, 'fracture_probability': None, 'external_reference': 'https://doi.org/10.12659/PJR.901840;Tables1-4; population plausibility gate failed'}
    text = json.dumps(output, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')
    if args.point_clearance:
        import numpy as np
        from operators import measure, H
        manifest = json.loads((ROOT / 'raw/R1_manifest.json').read_text())
        item = next((r for r in manifest if match(r)))
        from atlas import sha
        if sha(item['path']) != item['sha256']:
            raise SystemExit('Local mask hash drift')
        with np.load(item['path']) as a:
            (p, t, origin) = (a['pulp'], a['tooth'], a['origin_zyx'])
        (_, _, points, _, distance) = measure(p, t, origin)
        with args.point_clearance.open('w') as f:
            writer = csv.writer(f)
            writer.writerow(['case', 'fdi', 'resolution_level', 'z_index', 'y_index', 'x_index', 'inside_tooth', 'unsigned_boundary_center_distance_mm', 'digital_voxel_union_lower_mm', 'digital_voxel_union_upper_mm', 'anatomical_segmentation_error_mm', 'clinical_perforation_risk'])
            for (q, d) in zip(points, distance):
                idx = q + origin
                writer.writerow([args.case, args.fdi, 'PER_POINT', *idx.tolist(), bool(t[tuple(q)]), float(d), max(0, float(d) - np.sqrt(3) * H), float(d) + np.sqrt(3) * H, 'UNKNOWN', 'UNKNOWN'])
if __name__ == '__main__':
    main()
