"""I — the adjustable-resolution toolkit's EXECUTABLE END-TO-END DEMO (certificates visible; BUDGET-AWARE, type-agnostic).

The V14-1 shippable was a DOC; this is the RUNNABLE demo: one entrypoint that takes a decision, runs the full workflow
(declare → grid_trust the SPATIAL variable when a convention is given → classify the governor → CLASS-1 certify-with-budget
via the Richardson continuum, OR CLASS-2/3 abstain with the correct lever), and EMITS the CERTIFICATE. Upgraded with the
B-RATIFIED budget-aware refinement: a CLASS-1 (resolvable) decision whose required resolution exceeds the AVAILABLE budget
ABSTAINS-OVER-BUDGET (practically CLASS-2). decide() is now RESOLUTION-TYPE-AGNOSTIC — spatial (voxel via grid_trust) OR
temporal/count (segments via a required≤available budget) — the adjustable-resolution-per-VARIABLE thesis.

Reuses the built units: grid_trust + richardson_bracket (richardson_continuum_bracket.py) + the budget-aware trigger
(decidability_budget_aware_trigger.py). The class rule is the atlas's uniform rule (divergent→3, refinable→1, floored→2).

PRE-REGISTERED GATES:
  G1 the pipeline RUNS end-to-end on 4 decisions (3 spatial-voxel: modal/fatigue/stent + 1 temporal-budget: Z24 mode-dir)
     producing a CERTIFICATE each.
  G2 the intrinsic classes are CORRECT: modal→CLASS-1 CERTIFY-WITH-BUDGET, fatigue→CLASS-2 ABSTAIN(different-measurement),
     stent→CLASS-3 ABSTAIN(regularize) — matching the atlas + the delivered instances.
  G3 each certificate is VISIBLE + honest: shows {decision, grid-trust|budget, governor, class, action, confidence}; the
     spatial cases DECLARE a convention (an undeclared one would raise).
  G4 the BUDGET-AWARE refinement (B-ratified): the Z24 mode-direction is intrinsic CLASS-1 (duration-reducible) but its
     required resolution (n*≈5601 60s-segments) exceeds the available budget (270) → the demo returns ABSTAIN-OVER-BUDGET
     (resolvable in principle, practically buy a different measurement) — not a false CERTIFY.

Run: .venv-newton/bin/python scripts/physics_exp/adjres_end_to_end_demo.py
"""
import json
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from richardson_continuum_bracket import grid_trust, richardson_bracket
K = 3.0

def classify(signal):
    """the atlas uniform rule (divergent→3, refinable→1, floored→2)."""
    if signal.get('divergent'):
        return 3
    if signal.get('refinable'):
        return 1
    if signal.get('floored'):
        return 2
    return 0
ACTION = {1: 'CERTIFY-WITH-BUDGET (buy the stated resolution)', 2: 'ABSTAIN — different measurement needed (a shape/model σ floors it; more of the same resolution is wasted)', 3: 'ABSTAIN + REGULARIZE (fixed material length / percentile; refining diverges — the void-floor trap)'}
OVERBUDGET_ACTION = 'ABSTAIN — OVER-BUDGET (resolvable in principle, but the required resolution exceeds the available budget; buy a DIFFERENT measurement — same lever as CLASS-2, distinct reason)'

def decide(name, governor, signal, feature_size=None, res=None, convention=None, cont_grids=None, cont_r=2.0, required=None, available=None, budget_unit='', margin=None, margin_tol=None, sigma_meas=None, k=K):
    """resolution-type-agnostic end-to-end entrypoint (adjustable-resolution-per-VARIABLE) spanning ALL FOUR axes: intrinsic
    class (1/2/3) · budget (required ≤ available) · margin-TOLERANCE (worst-case k·(σ+δm)<|m| vs RSS k·√(σ²+δm²)<|m|) · the
    σ_min-governor. declare → [grid_trust] → classify → certify / abstain / abstain-OVER-BUDGET → [margin-tolerance WC/RSS
    refinement] → CERTIFICATE. The tolerance layer (when margin+margin_tol+sigma_meas given, for a CLASS-1 within-budget
    decision) is the margin-tolerance axis: a decision RSS-decidable can fail worst-case → TOLERANCE-SENSITIVE."""
    cls = classify(signal)
    over_budget = bool(cls == 1 and required is not None and (available is not None) and (required > available))
    cert = {'decision': name, 'governor': governor, 'class': cls, 'class_name': {1: 'PRECISION', 2: 'CROSS-QUANTITY', 3: 'DIVERGENT'}[cls], 'action': OVERBUDGET_ACTION if over_budget else ACTION[cls]}
    if convention is not None:
        gt = grid_trust(feature_size, res, convention)
        cert.update({'convention': gt['convention'], 'cells_per_radius': round(gt['cells_per_radius'], 2), 'grid_trusted': gt['trusted']})
    if required is not None and available is not None:
        ratio = required / available
        cert['budget'] = f'required {required} {budget_unit} vs available {available} → ' + (f'OVER by {ratio:.1f}x (practically abstain)' if over_budget else f'within budget (headroom {1 / ratio:.1f}x)')
    if over_budget:
        cert['confidence'] = 'UNDECIDABLE at this BUDGET — required resolution exceeds available; see action'
    elif cls == 1 and cont_grids is not None:
        R = richardson_bracket(cont_grids, cont_r)
        cert['certified_QoI'] = f"continuum {R['continuum']:.3f} ± GCI {(R['gci'] or 0) * 100:.1f}% (order p={R['p_observed']:.2f})"
        cert['confidence'] = 'certifiable at the stated resolution budget' if cert.get('grid_trusted', True) else 'grid below trust (sub-grid) — refine first'
    else:
        cert['confidence'] = 'UNDECIDABLE at this measurement — see action'
    if not over_budget and cls == 1 and (margin is not None) and (margin_tol is not None) and (sigma_meas is not None):
        wc = k * (sigma_meas + margin_tol) < margin
        rss = k * (sigma_meas ** 2 + margin_tol ** 2) ** 0.5 < margin
        if wc:
            cert['tolerance'] = f'CERTIFY-ROBUST: worst-case decidable (k·(σ+δm)={k * (sigma_meas + margin_tol):.3g} < margin {margin:.3g})'
        elif rss:
            cert['tolerance'] = f'TOLERANCE-SENSITIVE: RSS-decidable (k·√(σ²+δm²)={k * (sigma_meas ** 2 + margin_tol ** 2) ** 0.5:.3g} < {margin:.3g}) but NOT worst-case (k·(σ+δm)={k * (sigma_meas + margin_tol):.3g} ≥ {margin:.3g}) — the margin tolerance δm risks a WC miss'
            cert['action'] = 'CERTIFY-STATISTICAL (RSS) — worst-case NOT guaranteed; report both (see tolerance)'
            cert['confidence'] = 'decidable STATISTICALLY (RSS) but not worst-case — margin-tolerance-sensitive'
        else:
            cert['tolerance'] = 'ABSTAIN-TOLERANCE: even RSS undecidable — the margin tolerance δm floors it (reduce δm / different criterion)'
            cert['action'] = 'ABSTAIN — margin tolerance too large'
            cert['confidence'] = 'UNDECIDABLE — margin/threshold tolerance too large'
    return cert

def main():
    print('=' * 104)
    print('I — adjustable-resolution toolkit END-TO-END DEMO (all 4 axes: class · budget · margin-tolerance · σ-governor)')
    print('=' * 104)
    demos = [decide('21700 modal-margin (wall)', feature_size=0.197, res=0.032, convention='width', governor='wall thickness (partial-volume, voxel-limited)', signal={'refinable': True}, cont_grids=[0.9564, 0.9755, 1.0518]), decide('Ti-64 fatigue margin', feature_size=0.872, res=0.03, convention='width', governor='sphericity/shape σ (decidability-blocking; voxel-invariant)', signal={'floored': True}), decide('Zilver stent FSF', feature_size=0.12, res=0.03, convention='width', governor='peak-strain singular readout (λ=2.70, diverges under refinement)', signal={'divergent': True}), decide('Z24 mode-direction detectability', governor='5-mode observability σ_min at noise floor (duration-reducible; TEMPORAL variable)', signal={'refinable': True}, required=5601, available=270, budget_unit='60s-segments'), decide('Z24 damage (scenario-4, temperature regime)', governor='1st-mode Δf; margin-tolerance δm≈σ (temperature-driven threshold spread)', signal={'refinable': True}, margin=0.1534, sigma_meas=0.0271, margin_tol=0.0271)]
    for c in demos:
        print(f"\n  ── CERTIFICATE: {c['decision']}")
        print(f"     governor      : {c['governor']}")
        if 'cells_per_radius' in c:
            print(f"     grid-trust    : {c['cells_per_radius']} cells/radius ({c['convention']}) → trusted={c['grid_trusted']}")
        if 'budget' in c:
            print(f"     budget        : {c['budget']}")
        print(f"     class         : CLASS-{c['class']} {c['class_name']}")
        print(f"     action        : {c['action']}")
        if 'certified_QoI' in c:
            print(f"     certified QoI : {c['certified_QoI']}")
        if 'tolerance' in c:
            print(f"     margin-tol    : {c['tolerance']}")
        print(f"     confidence    : {c['confidence']}")
    g1 = bool(len(demos) == 5 and all(('class' in c for c in demos)))
    spatial = demos[:3]
    spatial_classes = [c['class'] for c in spatial]
    g2 = bool(spatial_classes == [1, 2, 3])
    g3_visible = all(('governor' in c and 'action' in c and ('class' in c) and ('confidence' in c) for c in demos))
    g3_convention_spatial = all(('convention' in c for c in spatial))
    raised = False
    try:
        grid_trust(0.1, 0.02, 'unspecified')
    except ValueError:
        raised = True
    g3 = bool(g3_visible and g3_convention_spatial and raised)
    z = demos[3]
    g4 = bool(z['class'] == 1 and 'OVER-BUDGET' in z['action'] and ('budget' in z) and z['confidence'].startswith('UNDECIDABLE at this BUDGET'))
    t = demos[4]
    g5 = bool(t['class'] == 1 and 'tolerance' in t and ('TOLERANCE-SENSITIVE' in t['tolerance']) and ('CERTIFY-STATISTICAL' in t['action']))
    ok = g1 and g2 and g3 and g4 and g5
    print('\n' + '-' * 104)
    print(f'[G1] pipeline runs end-to-end on 5 decisions (3 spatial-voxel + 1 temporal-budget + 1 margin-tolerance), certificate each -> {g1}')
    print(f'[G2] intrinsic classes CORRECT: spatial classes {spatial_classes} == [1,2,3] (certify / abstain-diff-meas / abstain-regularize) -> {g2}')
    print(f'[G3] each certificate VISIBLE+honest (governor/class/action/confidence; spatial cases declare convention); undeclared RAISES={raised} -> {g3}')
    print(f'[G4] budget-aware (B-ratified): the Z24 mode-direction is intrinsic CLASS-1 but required 5601 > available 270 → ABSTAIN-OVER-BUDGET -> {g4}')
    print(f'[G5] margin-tolerance axis: Z24 scenario-4 is CLASS-1 within-budget but RSS-decidable ∧ NOT worst-case → CERTIFY-STATISTICAL (tolerance-sensitive) -> {g5}')
    os.makedirs('reports', exist_ok=True)
    json.dump({'claim': "the adjustable-resolution toolkit's executable END-TO-END DEMO (certificates visible) — now spans ALL FOUR decision axes at RUNTIME (not just documented): intrinsic class (1/2/3) · budget (required ≤ available) · margin-TOLERANCE (worst-case k·(σ+δm)<|m| vs RSS k·√(σ²+δm²)<|m|) · the σ_min-governor. One entrypoint decide() runs declare → [grid_trust] → classify → certify/abstain → [budget-aware over-budget] → [margin-tolerance WC/RSS refinement] → CERTIFICATE. Demonstrated on 5 real decisions: modal→CERTIFY, fatigue→ABSTAIN-different-measurement, stent→ABSTAIN-regularize (class); Z24 mode-direction→ABSTAIN-OVER-BUDGET (budget, n*=5601≫270); Z24 scenario-4→CERTIFY-STATISTICAL/TOLERANCE-SENSITIVE (margin-tolerance: RSS-decidable 0.115<0.153 but worst-case fails 0.163≥0.153 at δm≈σ). The runnable shippable now matches the documented 4-axis procedure — a customer running it gets the WC-vs-RSS check, not just the class/budget verdict. Reuses grid_trust + richardson_bracket; convention enforced (undeclared raises).", 'gates': {'G1_n_demos': len(demos), 'G1': g1, 'G2_spatial_classes': spatial_classes, 'G2': g2, 'G3_convention_enforced': raised, 'G3': g3, 'G4_budget_aware_overbudget': g4, 'G4': g4, 'G5_margin_tolerance_sensitive': g5, 'G5': g5, 'verdict': 'PASS' if ok else 'FAIL'}, 'certificates': demos, 'provenance': 'reuses richardson_continuum_bracket.grid_trust + richardson_bracket; atlas uniform rule + the B-ratified budget-aware trigger + the margin-tolerance axis (margin_tolerance_decidability); 5 demo decisions (modal/fatigue/stent spatial + Z24 mode-direction temporal + Z24 scenario-4 tolerance-sensitive); CPU — the runnable shippable'}, open('reports/adjres_end_to_end_demo.json', 'w'), indent=1)
    print(f"\nVERDICT: {('PASS' if ok else 'FAIL')} — G1(5 decisions)={g1} G2(classes)={g2} G3(visible+convention)={g3} G4(budget-aware)={g4} G5(margin-tolerance)={g5}")
    print('EVIDENCE -> reports/adjres_end_to_end_demo.json')
    return 0 if ok else 1
if __name__ == '__main__':
    raise SystemExit(main())
