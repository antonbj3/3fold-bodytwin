"""I — ATLAS v2 (I_eff-reachability): the decidability CLASS-1/2 boundary REFINED to Fisher-information reachability, per D's
SEED W2 hand-off ("CLASS-1/2 = I_eff-reachability, not n-vs-m — handed to I as atlas v2"; frontier Φ(−Δ√I_eff); D's m*=9.69
reproduces I's n_req≈10).

Atlas v1 classified phenomenologically by dσ_M/d[resolution] (reducible/floored/divergent). Atlas v2 grounds the CLASS-1/2
split in the DECISION FRONTIER P_error = Φ(−Δ·√I_eff), where Δ = the decision margin and I_eff = the EFFECTIVE Fisher
information reachable by the certified variable. With a per-sample measurement variance σ_meas² and an IRREDUCIBLE
cross-quantity/model floor σ_floor² (a shape/model/seasonal term the certified resolution cannot touch):
    I_eff(n) = 1 / (σ_meas²/n + σ_floor²)          →  as n→∞,  I_eff → 1/σ_floor²  (floored).
    decidable at k  ⟺  Δ·√I_eff ≥ k.
So the class is I_eff-REACHABILITY: CLASS-1 iff the ceiling Δ·√I_eff_max = Δ/σ_floor ≥ k (reachable by buying resolution n;
σ_floor→0 ⇒ unbounded ⇒ always reachable), CLASS-2 iff Δ/σ_floor < k (the floor caps I_eff BELOW the decidability threshold —
a DIFFERENT quantity must raise it). ★This UNIFIES three of my axes: the CLASS axis (reachable vs floored), the margin-
TOLERANCE axis (δm IS the σ_floor that caps I_eff), and the σ_min-governor (I_eff = 1/σ² = identifiability = σ_min²).

PRE-REGISTERED GATES:
  G1 the FRONTIER reproduces atlas v1 + D's number: for Z24 magnitude (CLASS-1, σ_floor≈0) the reachability n_req from
     Δ√I_eff(n)=k matches the committed z24_damage n_req≈10 AND D's SEED-W2 m*=9.69 (CONSISTENCY not over-det: the three forms —
     my phenomenological cell, D's Φ-frontier, this I_eff form — are MATHEMATICALLY EQUIVALENT / same computation, so agreement is
     by construction, NOT an independent over-determination; per D's flagship-lesson + the fleet over-det-provenance sweep X6).
  G2 the I_eff-reachability CLASS matches atlas v1 on the real instances (Z24 magnitude→CLASS-1 reachable; Ti-64 fatigue→CLASS-2
     floored, σ_floor = the material δm) — the Fisher framing REPRODUCES the phenomenological classes (refinement, not contradiction).
  G3 the split is by REACHABILITY + a known-bad flip: CLASS-1 ⟺ Δ/σ_floor ≥ k; flipping σ_floor (0 ↔ large) flips the class
     (fatigue with σ_floor→0 becomes CLASS-1; Z24 with a large seasonal σ_floor becomes CLASS-2 — D's L3 second-floor finding).

Run: .venv-newton/bin/python scripts/physics_exp/atlas_v2_ieff_reachability.py
"""
import json
import os
import numpy as np
K = 3.0

def i_eff(n, sigma_meas, sigma_floor):
    """effective Fisher information reachable at resolution n: 1/(σ_meas²/n + σ_floor²) — floored at 1/σ_floor²."""
    return 1.0 / (sigma_meas ** 2 / n + sigma_floor ** 2)

def n_req(delta, sigma_meas, sigma_floor, k=K):
    """resolution n to reach decidability Δ√I_eff = k; None (∞ / CLASS-2) if the floor caps I_eff below the threshold."""
    denom = (delta / k) ** 2 - sigma_floor ** 2
    return sigma_meas ** 2 / denom if denom > 0 else None

def reachable(delta, sigma_floor, k=K):
    """I_eff-reachability: Δ·√I_eff_max = Δ/σ_floor ≥ k (σ_floor→0 ⇒ always reachable)."""
    return bool(sigma_floor <= 1e-12 or delta / sigma_floor >= k)

def main():
    print('=' * 104)
    print('I — ATLAS v2 (I_eff-reachability): CLASS-1/2 = Fisher-information reachability (D SEED-W2 hand-off)')
    print('=' * 104)
    s_meas = 0.0271
    sigma_meas_z = float(np.sqrt(2) * s_meas)
    delta_z = 0.0369
    nreq_z = n_req(delta_z, sigma_meas_z, 0.0)
    committed_nreq = 9.69
    g1 = bool(nreq_z is not None and abs(nreq_z - committed_nreq) < 1.0)
    print(f"\n[G1] Z24 magnitude n_req from Δ√I_eff=k: {nreq_z:.2f} ≈ committed z24_damage {committed_nreq} (== D's SEED-W2 m*=9.69) — CONSISTENCY (the 3 forms are mathematically EQUIVALENT: my σ_Δf(n), D's Φ-frontier, this I_eff — agreement confirms the reframing is the same reachability, NOT an independent over-det; D's flagship-lesson label) -> {g1}")
    (sigma_meas_f, sigma_floor_f, delta_f) = (5.5, 19.2, 40.0)
    z24_class = 1 if reachable(delta_z, 0.0) else 2
    fat_class = 1 if reachable(delta_f, sigma_floor_f) else 2
    g2 = bool(z24_class == 1 and fat_class == 2)
    print(f"[G2] I_eff-reachability class REPRODUCES atlas v1 (consistency, not over-det): Z24 magnitude (σ_floor≈0 ⇒ CLASS-1 by construction) → CLASS-{z24_class}; Ti-64 fatigue (Δ/σ_floor={delta_f / sigma_floor_f:.2f}<k ⇒ criterion eval) → CLASS-{fat_class} — refinement matches v1's phenomenological classes -> {g2}")
    fat_if_no_floor = 1 if reachable(delta_f, 0.0) else 2
    z24_if_seasonal = 1 if reachable(delta_z, 0.03) else 2
    flips = bool(fat_if_no_floor == 1 and z24_if_seasonal == 2)
    g3 = bool(flips)
    print(f"[G3] REACHABILITY split (σ_floor drives the class): fatigue σ_floor→0 ⇒ CLASS-{fat_if_no_floor} (reachable); Z24 + seasonal σ_floor=0.03Hz (Δ/σ_floor={delta_z / 0.03:.2f}<k) ⇒ CLASS-{z24_if_seasonal} (D's L3 second-floor) — the class FOLLOWS I_eff-reachability -> {g3}")
    ok = g1 and g2 and g3
    os.makedirs('reports', exist_ok=True)
    json.dump({'claim': "ATLAS v2 (I_eff-reachability) — the decidability CLASS-1/2 boundary REFINED to Fisher-information reachability, per D's SEED-W2 hand-off. The decision frontier is P_error = Φ(−Δ·√I_eff) with I_eff(n) = 1/(σ_meas²/n + σ_floor²): CLASS-1 iff the ceiling Δ·√I_eff_max = Δ/σ_floor ≥ k (reachable by buying resolution n; σ_floor→0 ⇒ unbounded), CLASS-2 iff Δ/σ_floor < k (an irreducible cross-quantity/model floor σ_floor caps I_eff BELOW the decidability threshold). CONSISTENCY (not over-det, per D's flagship-lesson label): the Z24-magnitude n_req from Δ√I_eff=k gives ~9.7, matching my committed z24_damage n_req AND D's Φ-frontier m*=9.69 — the three forms are MATHEMATICALLY EQUIVALENT (my σ_Δf(n), D's Φ-frontier, this I_eff), so agreement confirms the reframing is the SAME reachability, not an independent validation. The I_eff class REPRODUCES atlas v1 (Z24 magnitude→CLASS-1 reachable; Ti-64 fatigue→CLASS-2 floored, σ_floor = the material δm) — a refinement, not a contradiction. Known-bad flip: σ_floor drives the class (fatigue σ_floor→0 ⇒ CLASS-1; Z24 + a 0.03Hz seasonal σ_floor ⇒ CLASS-2 = D's L3 second-floor). ★UNIFIES three axes via Fisher information: the CLASS axis (reachable vs floored), the margin-TOLERANCE axis (the δm IS the σ_floor capping I_eff), and the σ_min-governor (I_eff = 1/σ² = identifiability = σ_min²). The phenomenological 'swing' of atlas v1 is d I_eff/dn.", 'gates': {'G1_nreq_consistency_not_overdet': {'nreq_ieff': round(nreq_z, 3), 'committed_and_D_mstar': committed_nreq}, 'G1': g1, 'G2_class_reproduces_v1': {'z24_magnitude': z24_class, 'ti64_fatigue': fat_class, 'fatigue_delta_over_floor': round(delta_f / sigma_floor_f, 3)}, 'G2': g2, 'G3_reachability_flip': {'fatigue_if_no_floor': fat_if_no_floor, 'z24_if_seasonal_floor': z24_if_seasonal}, 'G3': g3, 'verdict': 'PASS' if ok else 'FAIL'}, 'honest_scope': "σ_meas/σ_floor/Δ per instance are consumed from my committed cells (z24_damage s_meas=0.0271, fatigue σ_meas=5.5/material δm=19.2 MPa) + the subtlest real |Δf|; the design margin Δ_f and the seasonal σ_floor=0.03Hz are illustrative of the regimes (D's L3 measured the 30% seasonal leak). The I_eff form + the reachability criterion are the content; the frontier Φ(−Δ√I_eff) is D's SEED-W2 (this reproduces its n_req path, not a re-derivation of D's cell). ★D's L4 correctly labelled m*=9.69 as CONSISTENCY-not-overdet vs my n_req, and that label holds here too: my σ_Δf(n) form, D's Φ-frontier, and this I_eff form are the SAME reachability computation in three notations (same inputs) — their agreement confirms the reframing is equivalent, it is NOT an independent over-determination (a same-computation reproduction is consistency, per the overdet-tightness guard). The genuine NEW content is the I_eff-reachability REFRAMING + the unification of the CLASS / margin-tolerance / σ_min-governor axes, not a validation number.", 'provenance': 'D SEED-W2 hand-off (atlas v2 = I_eff-reachability, Φ(−Δ√I_eff), m*=9.69); real σ_meas/σ_floor from committed z24_damage + fatigue_margin_tolerance cells; CPU — a Fisher-information refinement of the atlas'}, open('reports/atlas_v2_ieff_reachability.json', 'w'), indent=1)
    print('=' * 104)
    print(f"VERDICT: {('PASS' if ok else 'FAIL')} — G1(n_req CONSISTENCY vs D m*=9.69)={g1} G2(class reproduces v1)={g2} G3(reachability flip)={g3}")
    print('EVIDENCE -> reports/atlas_v2_ieff_reachability.json')
    return 0 if ok else 1
if __name__ == '__main__':
    raise SystemExit(main())
