"""
Fresh reproduction of node COUNT-MOTOR-UNIT-NUMBER-FORCE-FACTORIZATION.

Pre-registered target (from the node's verify prose, computed BEFORE running):
  F_muscle_predicted = N x mean_single_MU_twitch_torque
  vs MVC = 35.0 N.m (tibialis anterior, young adults)
  Feinstein  N=445 -> 11.35 N.m (32.4% of MVC, -67.6% undershoot)
  McNeil     N=150 -> 3.83  N.m (10.9% of MVC, -89.1% undershoot)
  Van Cutsem N=300 -> 7.65  N.m (21.9% of MVC, -78.1% undershoot)
  mechanism check: generic twitch:tetanus ratio ~0.25 (1:4) should ALONE predict
  an ~75% undershoot (independent of any counting error).
Tolerance (pre-registered): within 1% relative on each N.m/% figure (closed-form
arithmetic from stated inputs); the "~75% undershoot" mechanism claim is checked
as 1 - twitch:tetanus_ratio landing in [70%, 80%] (stated as an approximate,
order-of-magnitude mechanism claim, not to 3 sig figs).

Stated inputs:
  single_MU_twitch_torque = 25.5 +/- 21.5 mN.m  (Van Cutsem et al 1997, PMID9415831)
  MVC (TA, young adults)  = 35.0 N.m            (triangulated literature anchor, ~34-36 N.m)
  N (Feinstein 1955, anatomical axon count, TA)      = 445
  N (McNeil et al 2005, PMID15685623, MUNE)          = 150
  N (Van Cutsem et al 1997, PMID9415831, self-est.)  = 300
  generic mammalian twitch:tetanus ratio ~= 0.25 (1:4, textbook)
"""

SINGLE_MU_TORQUE_MNM = 25.5e-3  # N.m (25.5 mN.m)
MVC_NM = 35.0

counts = {
    "Feinstein": (445, 11.35, 32.4, -67.6),
    "McNeil": (150, 3.83, 10.9, -89.1),
    "VanCutsem": (300, 7.65, 21.9, -78.1),
}

print(f"single-MU twitch torque = {SINGLE_MU_TORQUE_MNM*1e3:.1f} mN.m, MVC = {MVC_NM} N.m\n")

worst_rel = 0.0
for name, (N, prose_nm, prose_pct, prose_undershoot) in counts.items():
    F = N * SINGLE_MU_TORQUE_MNM
    pct = 100 * F / MVC_NM
    undershoot = pct - 100.0
    rel_nm = abs(F - prose_nm) / prose_nm
    rel_pct = abs(pct - prose_pct) / prose_pct
    rel_und = abs(undershoot - prose_undershoot) / abs(prose_undershoot)
    worst_rel = max(worst_rel, rel_nm, rel_pct, rel_und)
    print(f"{name}: N={N} -> F={F:.4f} N.m (prose {prose_nm}), "
          f"{pct:.2f}% of MVC (prose {prose_pct}%), undershoot {undershoot:.2f}% (prose {prose_undershoot}%)")
    print(f"   rel_diff: F={rel_nm*100:.3f}% pct={rel_pct*100:.3f}% undershoot={rel_und*100:.3f}%")

# mechanism check: twitch:tetanus ~0.25 -> predicted undershoot from THIS mechanism alone
twitch_tetanus_ratio = 0.25
mechanism_undershoot_pct = (1 - twitch_tetanus_ratio) * 100  # = 75%
print(f"\nmechanism check: twitch:tetanus={twitch_tetanus_ratio} -> "
      f"predicts {mechanism_undershoot_pct:.1f}% undershoot from fusion alone "
      f"(prose claims ~75%, 'ALONE-sufficient')")
mech_ok = 70.0 <= mechanism_undershoot_pct <= 80.0

TOL = 0.01
verdict_numbers = "REPRODUCES" if worst_rel <= TOL else "DIVERGES"
verdict_mechanism = "REPRODUCES" if mech_ok else "DIVERGES"
print(f"\nWorst relative diff on stated N.m/%/undershoot figures = {worst_rel*100:.4f}% (tol 1%) -> {verdict_numbers}")
print(f"Mechanism (twitch:tetanus alone predicts ~75% undershoot) -> {verdict_mechanism}")
