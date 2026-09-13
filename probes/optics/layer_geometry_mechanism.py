"""Measure adjacent box geometry before constructing a region-labelled mesh."""
import json
from pathlib import Path
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2]


def measure():
    rows=[];vertices=[];interface=[]
    for region,lo,hi in ((1,0.,10.),(2,10.,60.)):
        m=trimesh.creation.box(extents=[60,60,hi-lo]);m.apply_translation([30,30,(lo+hi)/2])
        mask=np.all(m.triangles[:,:,2]==10.,axis=1)
        rows.append({'region':region,'volume':float(m.volume),'vertices':len(m.vertices),'faces':len(m.faces),
            'closed':bool(m.is_watertight),'winding_consistent':bool(m.is_winding_consistent),
            'interface_triangles':int(mask.sum()),'interface_area':float(m.area_faces[mask].sum()),
            'interface_normals':m.face_normals[mask].tolist()})
        vertices.append(m.vertices);interface.extend([tuple(sorted(map(tuple,t))) for t in m.triangles[mask]])
    all_vertices=np.concatenate(vertices)
    return {'regions':rows,'input_vertices':len(all_vertices),'unique_coordinate_vertices':len(np.unique(all_vertices,axis=0)),
        'input_interface_triangles':len(interface),'unique_interface_triangles':len(set(interface))}


def main():
    legs=[measure(),measure()]
    gates={'two_tables_exact':legs[0]==legs[1],
        'closed_oriented_regions':all(r['closed'] and r['winding_consistent'] for r in legs[0]['regions']),
        'analytic_volumes_exact':[r['volume'] for r in legs[0]['regions']]==[36000.,180000.]}
    report={'legs':legs,'gates':gates,'scope':'synthetic source geometry only; no assembled transport mesh yet'}
    (ROOT/'reports/layer_geometry_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
