"""
Fresh reproduction of node COUNT-ALVEOLAR-NUMBER-SURFACE-FACTORIZATION.

Pre-registered target (from the node's verify prose, computed BEFORE running):
  empirical per-alveolus area = 2.9792e+05 um^2
  equal-volume-sphere prediction: r = 100.1 um, S_sphere = 1.2589e+05 um^2
  ratio (real/sphere) = 2.37x
Tolerance (pre-registered): within 1% relative (this is closed-form arithmetic from
stated inputs, not a noisy empirical fit -- any material difference means the prose
number itself does not close).

Stated inputs (from the node claim, cited primary sources):
  N_alveoli   = 480e6           (Ochs et al. 2004, PMID14512270, n=6 human lungs)
  V_alveolus  = 4.2e6 um^3      (Ochs et al. 2004, same paper, per-alveolus volume)
  S_total     = 143 m^2         (Gehr, Bachofen & Weibel 1978, PMID644146, n=8 human lungs)
"""
import math

N_ALVEOLI = 480e6
V_ALVEOLUS_UM3 = 4.2e6
S_TOTAL_M2 = 143.0

S_TOTAL_UM2 = S_TOTAL_M2 * 1e12  # 1 m^2 = 1e12 um^2

# empirical per-alveolus area
s_empirical = S_TOTAL_UM2 / N_ALVEOLI

# equal-volume sphere: V = 4/3 pi r^3 -> r = (3V/4pi)^(1/3); S = 4 pi r^2
r_sphere = (3 * V_ALVEOLUS_UM3 / (4 * math.pi)) ** (1 / 3)
s_sphere = 4 * math.pi * r_sphere ** 2

ratio = s_empirical / s_sphere

print(f"S_total(um^2)        = {S_TOTAL_UM2:.6e}")
print(f"s_empirical (um^2)   = {s_empirical:.6e}  (prose: 2.9792e+05)")
print(f"r_sphere (um)        = {r_sphere:.4f}     (prose: 100.1)")
print(f"s_sphere (um^2)      = {s_sphere:.6e}  (prose: 1.2589e+05)")
print(f"ratio (real/sphere)  = {ratio:.4f}       (prose: 2.37)")

targets = {
    "s_empirical": (2.9792e5, s_empirical),
    "r_sphere": (100.1, r_sphere),
    "s_sphere": (1.2589e5, s_sphere),
    "ratio": (2.37, ratio),
}
print()
worst_rel = 0.0
for name, (prose, mine) in targets.items():
    rel = abs(mine - prose) / prose
    worst_rel = max(worst_rel, rel)
    print(f"{name}: prose={prose} mine={mine:.6g} rel_diff={rel*100:.4f}%")

TOL = 0.01
verdict = "REPRODUCES" if worst_rel <= TOL else "DIVERGES"
print(f"\nWorst relative diff = {worst_rel*100:.4f}%  (tol {TOL*100:.1f}%) -> {verdict}")
