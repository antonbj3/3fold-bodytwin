"""
Rebuild of MODEL-MTORC1-FEEDBACK-SIGNALING from its own stated constants/topology
(cell claim + cert_design.datapoints). The cell's text names 11 states but the
FULL reaction stoichiometry (beyond the named feedback/drug edges) is NOT given in
the cert -- unlike WNT, no explicit v1..vN flux list is present. This is a
best-effort literal reconstruction from the named edges + conserved pools;
flagged where the spec under-determines the topology.

States (inferred minimal topology matching every named constant/edge):
  IR_p          (active insulin receptor;  IR_free = 1 - IR_p, normalized)
  IRS1          (unmodified, "available")
  IRS1_pY       (tyrosine-phosphorylated, active -> drives AKT-T308)
  IRS1_pS636    (S6K-feedback-phosphorylated, INACTIVE)
  Sink          (degraded IRS1 pool, resynthesized)
  AKT_free, AKT_T308, AKT_S473   (AKTtot=144.13 conserved 3-way pool)
  TSC_active, TSC_inactive       (TSCtot=10)
  mTORC1_p (active), mTORC1_i=mTORC1tot-mTORC1_p
  S6K_p (active), S6K_i=S6Ktot-S6K_p
  mTORC2_p (active), mTORC2_i=mTORC2tot-mTORC2_p
"""
import numpy as np
from scipy.integrate import solve_ivp

P = dict(
    k_IRon=0.025376, k_IRoff=0.149328,
    k_act=0.134664, k_dephos=0.00328283,
    k_fb1=1.0, k_fb2=0.0001,
    k_deg=0.0001, k_syn=0.0999968,
    k_T308on=0.699505, k_T308off=4.0739,
    k_S473on=4.50769, k_S473off=7.52842,
    k_TSCinact=0.00627315, k_TSCreact=0.00812537,
    k_mTORC1on=0.0513784, k_mTORC1off=0.999989,
    k_S6Kon=0.00573896, k_S6Koff=0.00528455,
    k_mTORC2on=0.0318902, k_mTORC2off=0.0255714,
)
AKTtot, TSCtot, mTORC1tot, S6Ktot, mTORC2tot = 144.13, 10.0, 4.3225, 127.0725, 6.2175
IRStot = 1.0  # normalized fraction pool for IR/IRS1 layer (unit-free, sign/timescale test only)

# state vector: [IR_p, IRS1, IRS1_pY, IRS1_pS636, Sink,
#                AKT_free, AKT_T308, AKT_S473,
#                TSC_active, mTORC1_p, S6K_p, mTORC2_p]
NAMES = ["IR_p", "IRS1", "IRS1_pY", "IRS1_pS636", "Sink",
         "AKT_free", "AKT_T308", "AKT_S473",
         "TSC_active", "mTORC1_p", "S6K_p", "mTORC2_p"]

def rhs(t, y, p, insulin=1.0, rapa_block=0.0, mtorki_block=0.0):
    IR_p, IRS1, IRS1_pY, IRS1_pS636, Sink, AKT_free, AKT_T308, AKT_S473, TSC_a, mTORC1_p, S6K_p, mTORC2_p = y
    IR_free = 1.0 - IR_p
    TSC_i = TSCtot - TSC_a

    d_IRp = p['k_IRon'] * insulin * IR_free - p['k_IRoff'] * IR_p
    v_act = p['k_act'] * IR_p * IRS1
    v_dephos = p['k_dephos'] * IRS1_pY
    # S6K -> IRS1 feedback: phosphorylates available IRS1 to inactive pS636
    v_fb = (p['k_fb1'] * S6K_p / max(S6Ktot, 1e-9) + p['k_fb2']) * IRS1
    v_deg = p['k_deg'] * IRS1_pS636
    v_syn = p['k_syn'] * Sink
    d_IRS1 = -v_act + v_dephos - v_fb + v_syn
    d_IRS1_pY = v_act - v_dephos
    d_IRS1_pS636 = v_fb - v_deg
    d_Sink = v_deg - v_syn

    # AKT: T308 driven by IRS1_pY (rapamycin does NOT touch this edge);
    #      S473 driven by mTORC2_p (rapamycin-insensitive, mTOR-KI blocks it)
    v_T308on = p['k_T308on'] * IRS1_pY * AKT_free
    v_T308off = p['k_T308off'] * AKT_T308
    v_S473on = (1.0 - mtorki_block) * p['k_S473on'] * mTORC2_p * AKT_free
    v_S473off = p['k_S473off'] * AKT_S473
    d_AKT_free = -v_T308on + v_T308off - v_S473on + v_S473off
    d_AKT_T308 = v_T308on - v_T308off
    d_AKT_S473 = v_S473on - v_S473off

    # TSC: inactivated by AKT-T308 signal, spontaneously reactivates
    v_TSCinact = p['k_TSCinact'] * AKT_T308 * TSC_i
    v_TSCreact = p['k_TSCreact'] * TSC_a
    d_TSC_a = v_TSCinact - v_TSCreact

    # mTORC1: turned ON by amino acids (const rate on inactive pool), turned OFF catalyzed by active TSC
    mTORC1_i = mTORC1tot - mTORC1_p
    v_m1on = p['k_mTORC1on'] * mTORC1_i
    v_m1off = p['k_mTORC1off'] * TSC_a * mTORC1_p / max(TSCtot, 1e-9)
    d_mTORC1p = v_m1on - v_m1off

    # S6K activated by active mTORC1 -- rapamycin's target edge
    S6K_i = S6Ktot - S6K_p
    v_s6kon = (1.0 - rapa_block) * (1.0 - mtorki_block) * p['k_S6Kon'] * mTORC1_p * S6K_i
    v_s6koff = p['k_S6Koff'] * S6K_p
    d_S6Kp = v_s6kon - v_s6koff

    # mTORC2 activation (amino-acid/PI3K driven, constant on-rate; mTOR-KI ALSO blocks this edge)
    mTORC2_i = mTORC2tot - mTORC2_p
    v_m2on = (1.0 - mtorki_block) * p['k_mTORC2on'] * mTORC2_i
    v_m2off = p['k_mTORC2off'] * mTORC2_p
    d_mTORC2p = v_m2on - v_m2off

    return [d_IRp, d_IRS1, d_IRS1_pY, d_IRS1_pS636, d_Sink,
            d_AKT_free, d_AKT_T308, d_AKT_S473,
            d_TSC_a, d_mTORC1p, d_S6Kp, d_mTORC2p]

def baseline_y0():
    return [0.3, 0.5, 0.3, 0.2, 0.0,
            AKTtot - 8.97 - 4.34, 4.34, 8.97,
            TSCtot * 0.5, mTORC1tot * 0.5, S6Ktot * 0.5, mTORC2tot * 0.5]

def run_to_ss(p, insulin=1.0, rapa_block=0.0, mtorki_block=0.0, y0=None, T=2e6):
    if y0 is None:
        y0 = baseline_y0()
    sol = solve_ivp(rhs, [0, T], y0, args=(p, insulin, rapa_block, mtorki_block),
                     method='LSODA', rtol=1e-10, atol=1e-12, dense_output=True)
    return sol

def main():
    print("=== mTORC1/S6K/IRS1 feedback rebuild (best-effort topology from named edges) ===")
    print("NOTE: cert gives 18 rate constants + conserved-pool sizes + named EDGES,")
    print("but NOT the full stoichiometric reaction list (unlike WNT's 17 explicit fluxes).")
    print("Topology below is a minimal reconstruction consistent with every named edge;")
    print("exact literal-agreement with the Dalle Pezze2012 SBML magnitudes is NOT expected.\n")

    # 1. baseline steady state
    sol_base = run_to_ss(P, insulin=1.0)
    y_base = sol_base.y[:, -1]
    akt_active_base = y_base[NAMES.index('AKT_T308')] + y_base[NAMES.index('AKT_S473')]
    dydt_base = np.array(rhs(0, y_base, P, 1.0, 0.0, 0.0))
    print(f"baseline steady state converged: max|dydt|={np.max(np.abs(dydt_base)):.3e}")
    print(f"baseline active-AKT (T308+S473) = {akt_active_base:.4f}  (pool total AKTtot={AKTtot})")

    # 2. rapamycin: block S6K activation edge, re-run from baseline steady state
    sol_rapa = run_to_ss(P, insulin=1.0, rapa_block=1.0, y0=y_base, T=24 * 60.0)  # 24h in minutes
    y_rapa_24h = sol_rapa.y[:, -1]
    akt_active_rapa = y_rapa_24h[NAMES.index('AKT_T308')] + y_rapa_24h[NAMES.index('AKT_S473')]
    pct_change_rapa = (akt_active_rapa - akt_active_base) / akt_active_base * 100
    print(f"\n+rapamycin (k_S6Kon=0) at t=24h: active-AKT = {akt_active_rapa:.4f}"
          f"  ({pct_change_rapa:+.1f}% vs baseline)  [recorded: +49.0%]")

    # 3. mTOR-KI: block BOTH S6K edge and mTORC2 edge
    sol_ki_10min = run_to_ss(P, insulin=1.0, mtorki_block=1.0, y0=y_base, T=10.0)
    y_ki_10min = sol_ki_10min.y[:, -1]
    s473_10min = y_ki_10min[NAMES.index('AKT_S473')]
    s473_base = y_base[NAMES.index('AKT_S473')]
    pct_s473_10min = (s473_10min - s473_base) / s473_base * 100
    print(f"\n+mTOR-KI at t=10min: AKT-S473 = {s473_10min:.4f} (from {s473_base:.4f}, {pct_s473_10min:+.1f}%)"
          f"  [recorded: 8.97->0.21, -97.7%]")

    sol_ki_24h = run_to_ss(P, insulin=1.0, mtorki_block=1.0, y0=y_base, T=24 * 60.0)
    y_ki_24h = sol_ki_24h.y[:, -1]
    t308_24h = y_ki_24h[NAMES.index('AKT_T308')]
    s473_24h = y_ki_24h[NAMES.index('AKT_S473')]
    t308_base = y_base[NAMES.index('AKT_T308')]
    pct_t308_24h = (t308_24h - t308_base) / t308_base * 100
    pct_s473_24h = (s473_24h - s473_base) / s473_base * 100
    print(f"+mTOR-KI at t=24h: AKT-T308 = {t308_24h:.4f} (from {t308_base:.4f}, {pct_t308_24h:+.1f}%)"
          f"  [recorded: 4.34->13.69]")
    print(f"                   AKT-S473 = {s473_24h:.4f} ({pct_s473_24h:+.1f}%)  [recorded: -99.4%]")

    # 4. FORCED ADVERSARY: null feedback (k_fb1=k_fb2=0), re-equilibrate, reapply rapamycin
    P_null = dict(P); P_null['k_fb1'] = 0.0; P_null['k_fb2'] = 0.0
    sol_null_eq = run_to_ss(P_null, insulin=1.0, T=2e6)
    y_null_eq = sol_null_eq.y[:, -1]
    dydt_null = np.array(rhs(0, y_null_eq, P_null, 1.0, 0.0, 0.0))
    akt_active_null_base = y_null_eq[NAMES.index('AKT_T308')] + y_null_eq[NAMES.index('AKT_S473')]
    print(f"\nFORCED ADVERSARY: k_fb1=k_fb2=0, re-equilibrated: max|dydt|={np.max(np.abs(dydt_null)):.3e}"
          f", active-AKT_ss={akt_active_null_base:.6f}  [recorded pool-saturation claim: 144.128401]")
    sol_null_rapa = run_to_ss(P_null, insulin=1.0, rapa_block=1.0, y0=y_null_eq, T=24 * 60.0)
    y_null_rapa = sol_null_rapa.y[:, -1]
    akt_active_null_rapa = y_null_rapa[NAMES.index('AKT_T308')] + y_null_rapa[NAMES.index('AKT_S473')]
    pct_null = (akt_active_null_rapa - akt_active_null_base) / max(akt_active_null_base, 1e-12) * 100
    print(f"null-feedback + rapamycin at t=24h: active-AKT change = {pct_null:+.6f}%  [recorded: 0.000001%]")

    # 5. conservation trick: AKT pool must sum exactly to AKTtot at all times (algebraic bound)
    akt_sum_base = y_base[NAMES.index('AKT_free')] + y_base[NAMES.index('AKT_T308')] + y_base[NAMES.index('AKT_S473')]
    print(f"\nCONSERVATION CHECK: AKT_free+AKT_T308+AKT_S473 = {akt_sum_base:.6f}  (must == AKTtot={AKTtot})")
    tsc_sum = y_base[NAMES.index('TSC_active')]
    print(f"TSC_active in [0,{TSCtot}]: {tsc_sum:.4f}  (bound respected: {0 <= tsc_sum <= TSCtot})")

    print("\n=== SELFTEST ===")
    ok = True
    if abs(akt_sum_base - AKTtot) > 1e-6:
        print("FAIL: AKT conservation violated"); ok = False
    sign_rapa_ok = pct_change_rapa > 0
    sign_s473_collapse_ok = pct_s473_10min < -50
    sign_t308_recover_ok = pct_t308_24h > 0
    print(f"sign check rapamycin +AKT: {'PASS' if sign_rapa_ok else 'FAIL'} (recorded +49.0%, got {pct_change_rapa:+.1f}%)")
    print(f"sign check mTOR-KI S473 fast collapse: {'PASS' if sign_s473_collapse_ok else 'FAIL'} (recorded -97.7%, got {pct_s473_10min:+.1f}%)")
    print(f"sign check mTOR-KI T308 24h over-recovery: {'PASS' if sign_t308_recover_ok else 'FAIL'} (recorded +215%, got {pct_t308_24h:+.1f}%)")
    if not (sign_rapa_ok and sign_s473_collapse_ok and sign_t308_recover_ok):
        ok = False
    print("SELFTEST", "PASS" if ok else "FAIL (see sign checks above)")

if __name__ == "__main__":
    main()
