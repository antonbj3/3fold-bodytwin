#!/usr/bin/env python3
"""Small region-payload hand-off example: a region-bound geometry/material/mass/heat record.

Synthetic fixture only (the repo's own layered box); it is an interface example, not a
physiological claim. Every field carries units, provenance and any missing value, so the result
envelope (``bodytwin.framework``) can seal it without the geometry package owning that envelope.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / 'src'))

from bodytwin.geometry.region_mass_v1 import (  # noqa: E402
    MaterialProperty, Region, area_mm2, integrate_region, interface_area_by_pair_mm2,
    lookup_density, mass_moments_surface, region_boundary_volume_mm3, total_quantity)
from bodytwin.geometry.layered_box_v1 import layered_box  # noqa: E402


def outward_faces(box, label):
    return np.concatenate([box['faces'][box['back'] == label],
                           box['faces'][box['front'] == label][:, ::-1]])


def quantity(q):
    return {'name': q.name, 'value': q.value, 'units': q.units,
            'provenance': q.provenance, 'missing': list(q.missing)}


def material(m):
    if m is None:
        return None
    return {'name': m.name, 'value': m.value, 'units': m.units, 'provenance': m.provenance,
            'source': m.source, 'source_sha256': m.source_sha256}


def build():
    box = layered_box()
    interface = area_mm2(box, 'material_interface_area')
    cp = MaterialProperty('muscle_cp_assumption', 3490.0, 'J/(kg*K)', 'assumption',
                          'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K')
    states = []
    for label in (1, 2):
        region = Region(f'SYNTH-BOX-{label}', region_boundary_volume_mm3(box, label), 'mm3',
                        'outward boundary integral of a synthetic closed surface',
                        area_kind='material_interface_area', area_mm2=interface)
        states.append(integrate_region(region, density_key='MSK-MUSCLE', specific_heat=cp))
    moments = mass_moments_surface(box['vertices'], outward_faces(box, 2),
                                   density=lookup_density('MSK-MUSCLE'),
                                   region_id='SYNTH-BOX-2', basis='region_volume')
    return {
        'schema': 'i2_region_example_v1',
        'synthetic': True,
        'source_sha256': {'layered_box_v1.py':
                          '7a8e9ba717c4a6b4b75b88911f08d97616bb978611b1d793e6bf46344f5cd6e5',
                          'mesh_ingest_v1.py':
                          'a05396a93879319e9004b985cefd550b1d43476eb513a47e9c7c709c3c7277cd'},
        'area': {'kind': 'material_interface_area', 'value': interface, 'units': 'mm2',
                 'method': 'shared-face count once, back->front',
                 'by_pair': {f'{a}-{b}': v for (a, b), v in sorted(interface_area_by_pair_mm2(box).items())}},
        'regions': [
            {'region_id': s.region_id, 'material_region_id': 'MSK-MUSCLE',
             'volume_mm3': s.geometry.volume_mm3,
             'volume_method': s.geometry.volume_method, 'mass_basis': s.basis,
             'density': material(s.density), 'mass': quantity(s.mass_kg),
             'specific_heat': material(s.specific_heat),
             'heat_capacity': quantity(s.heat_capacity_j_per_k),
             'missing': list(s.missing)}
            for s in states],
        'totals': {'mass': quantity(total_quantity(states, 'mass_kg')),
                   'heat_capacity': quantity(total_quantity(states, 'heat_capacity_j_per_k'))},
        'mass_moments_region_2': {
            'region_id': moments.region_id, 'basis': moments.basis,
            'mass_kg': moments.mass_kg,
            'centroid_mm': None if moments.centroid_mm is None else moments.centroid_mm.tolist(),
            'inertia_kg_mm2': (None if moments.inertia_kg_mm2 is None
                               else moments.inertia_kg_mm2.tolist()),
            'units': {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'},
            'inertia_SI_units': moments.inertia_SI_units,
            'missing': list(moments.missing)},
    }


if __name__ == '__main__':
    out = build()
    text = json.dumps(out, indent=2, sort_keys=True)
    (HERE / 'region_example.json').write_text(text + '\n')
    print(text)
