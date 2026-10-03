#!/usr/bin/env python3
"""Seed 1, the interoperability test, done on one quantity with a published anchor.

The seed is that external biology models should be connectable so more questions become askable. Its
blocker was never a modelling question: no finite-element solver existed on this machine at all, so
interoperability had never been tested. scikit-fem 12.0.2 is now installed, which is an independent
implementation rather than a second copy of our own assumptions.

The quantity. The literature collection carries 62 records in the unit *ratio (bending/stretch at
matched peak tensile strain)*, with a computed floor of **0.25** at a peak strain of 0.1 for pure
bending with the neutral axis at mid-wall, for a response that is pointwise linear in tensile strain and
inert in compression.

Where 0.25 comes from, and why that matters here. Under uniform stretch the whole thickness carries the
peak strain, so the mean tensile response is e0. Under pure bending with the neutral axis at mid-wall
the strain varies linearly from -e0 to +e0, so only half the thickness is in tension and its mean
tensile strain is e0/2, giving a mean response of e0/4 — hence 0.25. The number is therefore ANALYTIC,
and it rests on the strain varying linearly through the thickness.

So the test is sharper than a comparison: an independent solver computes the strain field instead of
assuming it, which means it can say where the analytic 0.25 stops holding. A thick wall, a short span or
a clamped end all break the linear-through-thickness assumption.

What this delivers: three numbers for one quantity — ours (analytic), the external solver's (computed),
and the published floor — plus the itemised list of what had to be supplied to make an external solver
answer our question at all. That list is the actual interoperability result.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from skfem import (Basis, ElementVector, ElementTetP1, MeshTet, asm, condense, solve)
from skfem.models.elasticity import linear_elasticity, lame_parameters

OUT = Path('./results/ASSEMBLY_INTEROP_BENDING')
E, NU = 1.0e6, 0.3          # a declared elastic constant pair; the ratio is designed to be insensitive
PEAK_STRAIN = 0.1           # e0 in the published record


def analytic_ratio() -> float:
    """Ours: a pointwise-linear tension-only response, strain linear through the thickness."""
    # stretch: every fibre at e0 -> mean tensile response e0
    # bending: strain e0 * (2z/h) for z in [-h/2, h/2]; tension only on half, mean e0/2 there
    return 0.25


def solve_strain(mesh: MeshTet, basis: Basis, kind: str, h: float, length: float):
    """Compute the axial strain field under stretch or bending with an independent solver."""
    K = asm(linear_elasticity(*lame_parameters(E, NU)), basis)
    x = basis.doflocs
    dofs = basis.get_dofs(lambda p: np.isclose(p[0], 0.0))
    free_end = basis.get_dofs(lambda p: np.isclose(p[0], length))

    u = np.zeros(basis.N)
    if kind == 'stretch':
        # uniform axial displacement at the far end: u_x = e0 * L, giving uniform strain e0
        u[free_end.nodal['u^1']] = PEAK_STRAIN * length
    else:
        # pure bending imposed kinematically: axial displacement proportional to height, so the
        # surface fibre reaches +e0 and the opposite surface -e0, with the neutral axis at mid-wall
        z = x[2, free_end.nodal['u^1']] - h / 2.0
        u[free_end.nodal['u^1']] = PEAK_STRAIN * length * (2.0 * z / h)

    fixed = np.concatenate([dofs.all(), free_end.nodal['u^1']])
    u = solve(*condense(K, np.zeros(basis.N), x=u, D=fixed))

    # axial strain per element from the displacement gradient, evaluated at element centroids
    eps = basis.interpolate(u).grad[0][0]        # du_x/dx
    return np.asarray(eps).ravel()


def response(eps: np.ndarray) -> float:
    """The declared response: pointwise linear in tensile strain, inert in compression."""
    return float(np.mean(np.clip(eps, 0.0, None)))


def main() -> int:
    rows = []
    # Sweep the slenderness, because that is exactly what the analytic number assumes away.
    for length, h in ((20.0, 1.0), (10.0, 1.0), (5.0, 1.0), (2.0, 1.0)):
        mesh = (MeshTet
                .init_tensor(np.linspace(0, length, 25), np.linspace(0, 1.0, 4), np.linspace(0, h, 9)))
        basis = Basis(mesh, ElementVector(ElementTetP1()))
        r_stretch = response(solve_strain(mesh, basis, 'stretch', h, length))
        r_bend = response(solve_strain(mesh, basis, 'bending', h, length))
        ratio = r_bend / r_stretch if r_stretch else float('nan')
        rows.append({
            'span_over_thickness': length / h,
            'dofs': int(basis.N),
            'mean_tensile_response_stretch': r_stretch,
            'mean_tensile_response_bending': r_bend,
            'external_solver_ratio': ratio,
            'deviation_from_published_floor': ratio - 0.25,
        })
        print(f"  L/h {length/h:5.1f}  dofs {basis.N:6d}  ratio {ratio:.6f}  "
              f"deviation {ratio - 0.25:+.6f}")

    interface = [
        'geometry as a volume mesh, not a thickness scalar — our cell carries no mesh',
        'an elastic constant pair (E, nu); our cell declares a modulus but not a Poisson ratio',
        'boundary conditions as displacement fields on named faces; our cell has none',
        'a kinematic definition of "pure bending" — ours is implicit in a formula',
        'the response functional stated pointwise over a field rather than on a summary',
        'units for every input; the solver does not carry them, so a mismatch is silent',
    ]
    unsuppliable = [
        'a Poisson ratio for the tissue in question, which our cell does not contain',
        'the end conditions of the real specimen, which the published record does not state',
    ]

    summary = {
        'quantity': 'ratio (bending/stretch at matched peak tensile strain)',
        'unit': '1',
        'ours_analytic': analytic_ratio(),
        'published_floor': 0.25,
        'published_scope': 'peak strain 0.1, pure bending, neutral axis at mid-wall, pointwise linear '
                           'tension-only response',
        'external_solver': 'scikit-fem 12.0.2, linear elasticity, tetrahedral P1',
        'rows': rows,
        'interface_items_required': interface,
        'interface_items_we_could_not_supply': unsuppliable,
        'questions_this_makes_askable': [
            'at what span-to-thickness does the analytic quarter stop holding, and by how much',
            'what the ratio becomes under a clamped rather than kinematic end, which the published '
            'record does not specify',
            'whether a tension-only response on a computed field differs from the same response on an '
            'assumed linear field, which is the sufficiency question in mechanical form',
        ],
        'claim_type': 'information_link',
        'control': 'our own analytic value for the same quantity; the published floor is the external '
                   'reference and neither implementation is a facit for the other',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'INTEROP_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"  ours (analytic) {analytic_ratio()}   published floor 0.25")
    print(f"  interface items required {len(interface)}, unsuppliable {len(unsuppliable)}")
    print(f"  written: {OUT / 'INTEROP_V1.json'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
