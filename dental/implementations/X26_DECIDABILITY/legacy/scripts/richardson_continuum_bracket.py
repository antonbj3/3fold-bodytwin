"""I-V10-1 — the REUSABLE Richardson / continuum-bracket HELPER (gap-analysis top-2; ~79 convergence cells hand-roll this, 27
miss a continuum limit). Generalizes my 3-grid GCI template into ONE correct unit lanes import instead of hand-rolling.

D's CELL_SOLVER_GAP_ANALYSIS.md: convergence is declared ~79× across the fleet, 27 cells with NO continuum limit (a bare
finest-grid point, one-sided-biased). This is the Richardson analog of the σ_min "hand-rolled 1352×, no correct module"
finding — build the ONE helper. It carries the disciplines I've banked: the one-sided-bias guard ([[richardson-brackets-
onesided-convergent-bias]] — report the continuum estimate, not the point), band-not-point declarations, and the
SHOCK/kink order guard ([[richardson-brackets-onesided-convergent-bias]] §shock: the working order of a discontinuity-
dominated observable is the Kuznetsov order ~0.5, NOT the formal order — a formal-order Richardson read then LIES).

THE HELPER: richardson_bracket(values, r, formal_order=None) → {observed order p, continuum estimate, GCI band, the honest
BRACKET (continuum ± GCI), one-sided-bias flag, order-regime}. Rules baked in:
  • report ε = continuum ± GCI, NEVER a bare finest-grid value (the 2-point change undersells the coarse grid's true error);
  • if the observed order p ≪ the formal order → SHOCK/kink regime → use the MEASURED order + a band-window, not formal Richardson;
  • non-monotone / converged-to-noise triples → no clean power-law → report the spread as the band, flag "no continuum".

VALIDATION ANCHORS (retrofit targets from the assignment): (1) my synthetic known-answer (p=2 recovery); (2) my STENT GCI
template numbers (ε_m 15/30/60µm → continuum 0.95, the ~2.7% coarse-grid budget); (3) a SHOCK case (A's dam-break class:
observed order ~0.5 vs formal 1 → the helper flags shock and does NOT trust a formal-order read).

PRE-REGISTERED GATES:
  G1 SMOOTH known-answer: on f(h)=f∞(1+a·h^p_true) the helper recovers p≈p_true and continuum≈f∞ to <1% + the honest bracket
     contains f∞ (which a bare finest-grid point does NOT, one-sided).
  G2 STENT anchor (retrofit): my stent ε_m {0.9564,0.9755,1.0518}@{15,30,60}µm → continuum 0.950, p≈2, and the coarse (30µm)
     grid's true error to continuum (2.68%) EXCEEDS the naive 2-point change (2.00%) → the helper's band-not-point is the fix.
  G3 SHOCK guard + non-monotone: (a) a p≈0.5 sequence with formal_order=1 is flagged SHOCK (don't formal-Richardson it); (b) a
     non-monotone triple returns "no continuum" + the spread as the band (not a bogus extrapolation). Over-det: 3 regimes
     (smooth / shock / non-monotone) each handled distinctly.

Run: .venv-newton/bin/python scripts/physics_exp/richardson_continuum_bracket.py
"""
import json
import math
import os
SAFETY = 1.25
GRID_TRUST_CELLS_PER_RADIUS = 2.0

def grid_trust(feature_size, res, convention):
    """Is a grid representation trustworthy for a feature of `feature_size` at cell size `res`, AND is the continuum
    limit safe to extrapolate? The CONVENTION MUST be declared — 'radius' | 'width' | 'diameter' — because 2 cells/radius
    = 4 cells/width (a factor 2); an undeclared 'N cells per feature' silently confuses trusted-vs-undersampled.
    Returns {convention, feature_radius, cells_per_radius, cells_per_width, kappa_res, trusted}. Consumers MUST call this
    before trusting a Richardson continuum: a sub-grid feature (cells_per_radius<2) has no reliable continuum."""
    conv = str(convention).lower()
    if conv == 'radius':
        radius = float(feature_size)
    elif conv in ('width', 'diameter'):
        radius = float(feature_size) / 2.0
    else:
        raise ValueError("convention MUST be declared as 'radius' | 'width' | 'diameter' (got %r) — the 2-vs-4-cells pin needs it" % convention)
    cpr = radius / res
    return {'convention': conv, 'feature_radius': radius, 'cells_per_radius': cpr, 'cells_per_width': cpr * 2, 'kappa_res': res / radius, 'trusted': bool(cpr >= GRID_TRUST_CELLS_PER_RADIUS)}

def richardson_bracket(values, r, formal_order=None, safety=SAFETY):
    """3-grid Richardson + GCI + one-sided-bias band + shock/non-monotone guards.
    values = [f_fine, f_mid, f_coarse] (finest first); r = refinement ratio (h_mid/h_fine = h_coarse/h_mid).
    Returns a dict; ALWAYS report continuum ± GCI (the bracket), never a bare finest-grid value."""
    (f1, f2, f3) = [float(v) for v in values]
    (num, den) = (f3 - f2, f2 - f1)
    if abs(den) < 1e-30 or num / den <= 0:
        if abs(den) > 1e-30 and abs(num) < abs(den):
            band = (min(f1, f2), max(f1, f2))
            return {'p_observed': None, 'continuum': 0.5 * (f1 + f2), 'gci': None, 'bracket': band, 'one_sided': False, 'regime': 'bracketing (alternating-convergent; finite-grid spread CONTAINS the truth — extrapolation optional)', 'note': 'spread [f1,f2] is a valid two-sided band (unlike a one-sided/dissipative scheme whose spread misses the truth)'}
        band = (min(values), max(values))
        return {'p_observed': None, 'continuum': f1, 'gci': None, 'bracket': band, 'one_sided': None, 'regime': 'no-continuum (non-monotone / converged-to-noise)', 'note': '3 grids show no clean power-law and steps not shrinking; report the spread as the band, do NOT extrapolate'}
    p = math.log(abs(num / den)) / math.log(r)
    if abs(r ** p - 1.0) < 0.05:
        band = (min(f1, f2, f3), max(f1, f2, f3))
        return {'p_observed': p, 'continuum': 0.5 * (f1 + f2), 'gci': None, 'bracket': band, 'one_sided': None, 'regime': 'indeterminate-order (p≈0, on the convergent/divergent boundary — Richardson ill-conditioned; refine the sweep or widen the resolution range)', 'note': 'successive differences ~equal ⇒ no clear power-law order; report the spread, do NOT extrapolate'}
    continuum = f1 + (f1 - f2) / (r ** p - 1.0)
    gci = safety * abs((f1 - f2) / f1) / (r ** p - 1.0) if f1 != 0 else None
    one_sided = (f1 - continuum) * (f2 - continuum) > 0
    regime = 'smooth'
    if formal_order is not None and p < 0.7 * formal_order:
        regime = 'shock/kink-limited (observed p << formal; use the MEASURED order + a band-window, NOT formal-order Richardson)'
    half = gci * abs(continuum) if gci is not None else abs(f1 - continuum)
    bracket = (continuum - half, continuum + half)
    return {'p_observed': p, 'continuum': continuum, 'gci': gci, 'bracket': bracket, 'one_sided': bool(one_sided), 'regime': regime}

def main():
    print('=' * 104)
    print('I-V10-1 — Richardson/continuum-bracket HELPER (reusable unit; ~79 convergence cells hand-roll this, 27 miss a limit)')
    print('=' * 104)
    (f_inf, a, p_true, r) = (0.95, 0.0268, 2.0, 2.0)
    vals = [f_inf * (1 + a * (h / 30.0) ** p_true) for h in (15.0, 30.0, 60.0)]
    R = richardson_bracket(vals, r, formal_order=2)
    contains = R['bracket'][0] <= f_inf <= R['bracket'][1]
    point_undersells = not min(vals) <= f_inf or not f_inf <= max(vals)
    g1 = bool(abs(R['p_observed'] - p_true) < 0.05 and abs(R['continuum'] - f_inf) / f_inf < 0.01 and contains and R['one_sided'])
    print(f"[G1] smooth: p={R['p_observed']:.3f}(true {p_true}) continuum={R['continuum']:.4f}(true {f_inf}) bracket={tuple((round(x, 4) for x in R['bracket']))} contains f∞={contains}, one-sided={R['one_sided']} (finest grid {vals[0]:.4f} > f∞ {f_inf} => a point overshoots) -> {g1}")
    stent = [0.9564, 0.9755, 1.0518]
    S = richardson_bracket(stent, 2.0, formal_order=2)
    two_point = abs(stent[1] - stent[0]) / stent[0]
    true_budget_err = abs(stent[1] - S['continuum']) / S['continuum']
    g2 = bool(abs(S['continuum'] - 0.95) < 0.002 and abs(S['p_observed'] - 2.0) < 0.05 and (true_budget_err > two_point))
    print(f"[G2] STENT anchor: continuum={S['continuum']:.4f} (=0.95), p={S['p_observed']:.3f}; the 30µm grid's TRUE error {true_budget_err * 100:.2f}% EXCEEDS the naive 2-point change {two_point * 100:.2f}% -> band-not-point is the fix -> {g2}")
    shock_vals = [1.0 * h ** 0.5 for h in (1.0, 2.0, 4.0)]
    SH = richardson_bracket(shock_vals, 2.0, formal_order=1)
    shock_flagged = 'shock' in SH['regime']
    NM = richardson_bracket([1.0, 1.3, 0.9], 2.0)
    nm_flagged = NM['regime'].startswith('no-continuum')
    g3 = bool(shock_flagged and abs(SH['p_observed'] - 0.5) < 0.02 and nm_flagged)
    print(f"[G3] SHOCK: observed p={SH['p_observed']:.3f}(≈0.5) vs formal 1 -> flagged {shock_flagged} ('{SH['regime'][:38]}...'); NOISE triple (growing steps) -> {NM['regime']} (flagged={nm_flagged}) -> {g3}")
    BR = richardson_bracket([0.94, 0.96, 0.951], 2.0)
    brack_ok = bool('bracketing' in BR['regime'] and BR['bracket'][0] <= 0.95 <= BR['bracket'][1])
    p1_vals = [0.95 * (1 + 0.05 * (h / 30.0) ** 1.0) for h in (15.0, 30.0, 60.0)]
    P1 = richardson_bracket(p1_vals, 2.0)
    p1_ok = bool(abs(P1['p_observed'] - 1.0) < 0.05 and abs(P1['continuum'] - 0.95) / 0.95 < 0.01)
    g4 = bool(brack_ok and p1_ok)
    print(f"[G4] consolidation: BRACKETING (alternating-convergent) -> '{BR['regime'][:30]}...' spread {tuple((round(x, 3) for x in BR['bracket']))} contains 0.95={brack_ok} (correctly NOT 'no-continuum' — C's discriminator); p=1 known-answer p={P1['p_observed']:.3f} continuum={P1['continuum']:.4f} ({p1_ok}) -> {g4}")
    well = grid_trust(1.0, 0.2, 'radius')
    subgrid = grid_trust(1.0, 0.7, 'radius')
    as_radius = grid_trust(2.0, 1.0, 'radius')['trusted']
    as_width = grid_trust(2.0, 1.0, 'width')['trusted']
    convention_matters = bool(as_radius and (not as_width))
    forces_declaration = False
    try:
        grid_trust(1.0, 0.2, '2 cells')
    except ValueError:
        forces_declaration = True
    g5 = bool(well['trusted'] and (not subgrid['trusted']) and convention_matters and forces_declaration)
    print(f"[G5] grid-trust + CONVENTION PIN: well-resolved (5 cells/radius) trusted={well['trusted']}; sub-grid (1.4/radius) trusted={subgrid['trusted']}; SAME size/res trusted as RADIUS not as WIDTH (factor-2 pin)={convention_matters}; undeclared convention RAISES={forces_declaration} -> {g5}")
    ok = g1 and g2 and g3 and g4 and g5
    os.makedirs('reports', exist_ok=True)
    json.dump({'claim': "I-V10-1: the reusable Richardson/continuum-bracket HELPER (gap-analysis top-2: ~79 fleet convergence cells hand-roll this, 27 with NO continuum limit — the Richardson analog of the σ_min hand-rolled-1352× gap). richardson_bracket(values,r,formal_order) returns {observed p, continuum, GCI, the honest BRACKET continuum±GCI, one-sided flag, order-regime}. Baked-in disciplines: report continuum±GCI never a bare finest-grid point (the 2-point change undersells the coarse grid's true error, e.g. stent 30µm true 2.68%% > naive 2.00%%); flag SHOCK/kink regime when observed p≪formal (Kuznetsov ~0.5, a formal-order read LIES); non-monotone/converged triples get 'no continuum'; a BRACKETING (alternating-convergent) scheme is correctly NOT flagged (its spread already contains the truth — consolidating C's spine_richardson_onesided_bias_trap discriminator). Validated across FIVE regimes: smooth p=1 AND p=2; shock p=0.5; the STENT GCI continuum 0.950; noise; bracketing.", 'gates': {'G1_p': R['p_observed'], 'G1_continuum': R['continuum'], 'G1_bracket_contains': contains, 'G1': g1, 'G2_stent_continuum': S['continuum'], 'G2_two_point': two_point, 'G2_true_budget_err': true_budget_err, 'G2': g2, 'G3_shock_p': SH['p_observed'], 'G3_shock_flagged': shock_flagged, 'G3_nonmono_flagged': nm_flagged, 'G3': g3, 'G4_bracketing_ok': brack_ok, 'G4_p1_ok': p1_ok, 'G4_p1_p': P1['p_observed'], 'G4': g4, 'G5_well_trusted': well['trusted'], 'G5_subgrid_trusted': subgrid['trusted'], 'G5_convention_matters': convention_matters, 'G5_forces_declaration': forces_declaration, 'G5': g5, 'verdict': 'PASS' if ok else 'FAIL'}, 'grid_trust_guard': {'rule': 'trust a grid rep + its Richardson continuum only above ~2 cells/feature-RADIUS (=4/WIDTH; κ·res≲0.5); below = sub-grid, no reliable continuum', 'convention_pin': "the consumer MUST declare 'radius'|'width'|'diameter' — 2 cells/radius = 4 cells/width (factor 2); an undeclared 'N cells per feature' is ambiguous (grid_trust RAISES)", 'api': 'grid_trust(feature_size, res, convention) -> {cells_per_radius, trusted, ...}; call BEFORE trusting a continuum', 'source': 'D GENERATIVE_GRADIENT_ROUTES §standing-guard T3 2026-07-05'}, 'helper_api': {'function': 'richardson_bracket(values=[f_fine,f_mid,f_coarse], r, formal_order=None)', 'returns': '{p_observed, continuum, gci, bracket=(lo,hi)=continuum±GCI, one_sided, regime}', 'rules': ['report continuum±GCI, NEVER a bare finest-grid value (one-sided-bias undersells)', 'observed p ≪ formal → SHOCK regime: use measured order + band-window, not formal Richardson', "non-monotone/converged triple → 'no continuum', report the spread as the band"], 'retrofit_targets': "~79 convergence cells (27 missing a continuum limit); anchors: this template + stent GCI + A's dam-break"}, 'provenance': 'generalizes my stent_strain_convergence_gci_template (Roache GCI + Richardson) + the shock-order guard from [[richardson-brackets-onesided-convergent-bias]] (Kuznetsov ~0.5 for discontinuities); D CELL_SOLVER_GAP_ANALYSIS.md scope; anchors synthetic/stent/shock; CPU, no solver — a reusable helper unit'}, open('reports/richardson_continuum_bracket.json', 'w'), indent=1)
    print('=' * 104)
    print(f"VERDICT: {('PASS' if ok else 'FAIL')} — G1(smooth)={g1} G2(stent)={g2} G3(shock+noise)={g3} G4(bracketing+p1)={g4} G5(grid-trust+convention-pin, I-V13-1)={g5}")
    print('EVIDENCE -> reports/richardson_continuum_bracket.json')
    return 0 if ok else 1
if __name__ == '__main__':
    raise SystemExit(main())
