"""
Rebuild of MODEL-WNT-BETACATENIN from its OWN stated equations/constants
(cell claim text, cert_design.datapoints). Lee2003/BIOMD0000000658-lineage,
15-species mass-action destruction-complex model. NOT tuned to match the
recorded headline numbers -- built from the stated fluxes/constants only,
then compared.

Species order: [Dsh_i, Dsh_a, APC, Axin, GSK3, APC_axin, C1, C1p, Bcat0,
                BcatC1p, BcatC1p_ph, Bcat_ph, Bcat_APC, TCF, Bcat_TCF]
"""
import numpy as np
from scipy.integrate import solve_ivp

K = dict(k1=0.182, k2=0.0182, k3=0.05, k4=0.267, k5=0.133,
         k6=0.0909, k_6=0.909, k7=500.0, k_7=25000.0, k8=500.0, k_8=60000.0,
         k9=206.0, k10=206.0, k11=0.417, k12=0.423, k13=0.000257,
         k14=8.22e-05, k15=0.167, k16=500.0, k_16=15000.0, k17=500.0, k_17=600000.0)

IDX = dict(Dsh_i=0, Dsh_a=1, APC=2, Axin=3, GSK3=4, APC_axin=5, C1=6, C1p=7,
           Bcat0=8, BcatC1p=9, BcatC1p_ph=10, Bcat_ph=11, Bcat_APC=12, TCF=13, Bcat_TCF=14)

def rhs(t, y, W, k):
    Dsh_i, Dsh_a, APC, Axin, GSK3, APC_axin, C1, C1p, Bcat0, BcatC1p, BcatC1p_ph, Bcat_ph, Bcat_APC, TCF, Bcat_TCF = y
    v1 = k['k1'] * Dsh_i * W
    v2 = k['k2'] * Dsh_a
    v3 = k['k3'] * Dsh_a * C1
    v4 = k['k4'] * C1
    v5 = k['k5'] * C1p
    v6 = k['k6'] * GSK3 * APC_axin - k['k_6'] * C1
    v7 = k['k7'] * APC * Axin - k['k_7'] * APC_axin
    v8 = k['k8'] * C1p * Bcat0 - k['k_8'] * BcatC1p
    v9 = k['k9'] * BcatC1p
    v10 = k['k10'] * BcatC1p_ph
    v11 = k['k11'] * Bcat_ph
    v12 = k['k12']
    v13 = k['k13'] * Bcat0
    v14 = k['k14']
    v15 = k['k15'] * Axin
    v16 = k['k16'] * Bcat0 * TCF - k['k_16'] * Bcat_TCF
    v17 = k['k17'] * APC * Bcat0 - k['k_17'] * Bcat_APC

    dDsh_i = -v1 + v2
    dDsh_a = v1 - v2
    dAPC = -v7 - v17
    dAxin = -v7 + v14 - v15
    dGSK3 = -v6 + v3
    dAPC_axin = v7 - v6 + v3
    dC1 = v6 - v3 - v4 + v5
    dC1p = v4 - v5 - v8 + v10
    dBcat0 = -v8 - v16 - v17 + v12 - v13
    dBcatC1p = v8 - v9
    dBcatC1p_ph = v9 - v10
    dBcat_ph = v10 - v11
    dBcat_APC = v17
    dTCF = -v16
    dBcat_TCF = v16
    return [dDsh_i, dDsh_a, dAPC, dAxin, dGSK3, dAPC_axin, dC1, dC1p, dBcat0,
            dBcatC1p, dBcatC1p_ph, dBcat_ph, dBcat_APC, dTCF, dBcat_TCF]

def y0():
    y = np.zeros(15)
    y[IDX['Dsh_i']] = 100.0
    y[IDX['APC']] = 100.0
    y[IDX['GSK3']] = 50.0
    y[IDX['Axin']] = 0.02
    y[IDX['TCF']] = 15.0
    return y

def steady_state(W, k=K, T=2e6):
    sol = solve_ivp(rhs, [0, T], y0(), args=(W, k), method='LSODA',
                     rtol=1e-10, atol=1e-14, dense_output=False)
    yf = sol.y[:, -1]
    dydt = np.array(rhs(sol.t[-1], yf, W, k))
    return yf, np.max(np.abs(dydt))

def main():
    print("=== WNT/beta-catenin rebuild (literal 17-flux/15-species reconstruction) ===")
    y0_ss, res0 = steady_state(0.0)
    y1_ss, res1 = steady_state(1.0)
    bcat0_W0 = y0_ss[IDX['Bcat0']]
    bcat0_W1 = y1_ss[IDX['Bcat0']]
    fold = bcat0_W1 / bcat0_W0
    print(f"max|dydt| at W=0: {res0:.3e}  (should be <1e-9 per cell claim)")
    print(f"max|dydt| at W=1: {res1:.3e}")
    print(f"Bcat0_ss(W=0) = {bcat0_W0:.4f} nM  [recorded: 25.144]")
    print(f"Bcat0_ss(W=1) = {bcat0_W1:.4f} nM  [recorded: 153.423]")
    print(f"fold = {fold:.4f}x  [recorded: 6.102x, PAPER Table2: 6.12x]")
    rel_err_w0 = abs(bcat0_W0 - 25.144) / 25.144 * 100
    rel_err_w1 = abs(bcat0_W1 - 153.423) / 153.423 * 100
    rel_err_fold = abs(fold - 6.102) / 6.102 * 100
    print(f"rel err vs recorded: W0 {rel_err_w0:.3f}%  W1 {rel_err_w1:.3f}%  fold {rel_err_fold:.3f}%")

    # conservation-law check: free-Axin_ss = k14/k15 EXACTLY (claimed algebraic identity)
    axin_pred = K['k14'] / K['k15']
    print(f"\nfree-Axin_ss predicted (k14/k15) = {axin_pred:.6e} nM")
    print(f"free-Axin_ss simulated (W=0)     = {y0_ss[IDX['Axin']]:.6e} nM")
    print(f"free-Axin_ss simulated (W=1)     = {y1_ss[IDX['Axin']]:.6e} nM")
    axin_err = abs(y1_ss[IDX['Axin']] - axin_pred) / axin_pred * 100
    print(f"rel err vs k14/k15 identity: {axin_err:.4f}%")

    # forced Axin-clamp adversary: clamp Axin at k14/k15, verify bit-identical result
    def rhs_clamped(t, y, W, k):
        d = rhs(t, y, W, k)
        d[IDX['Axin']] = 0.0
        return d

    y_clamp = y0()
    y_clamp[IDX['Axin']] = axin_pred
    sol_c = solve_ivp(rhs_clamped, [0, 2e6], y_clamp, args=(1.0, K), method='LSODA',
                       rtol=1e-10, atol=1e-14)
    bcat0_clamped = sol_c.y[IDX['Bcat0'], -1]
    diff_clamp = abs(bcat0_clamped - bcat0_W1)
    print(f"\nAxin-clamp adversary: Bcat0_ss(clamped) = {bcat0_clamped:.6f} vs unclamped {bcat0_W1:.6f}"
          f"  |diff|={diff_clamp:.2e}  (claim: bit-for-bit identical)")

    # conservation trick: does APC moiety (APC+APC_axin+C1+C1p+BcatC1p+BcatC1p_ph+Bcat_APC) stay fixed at 100?
    def apc_moiety(y):
        return y[IDX['APC']] + y[IDX['APC_axin']] + y[IDX['C1']] + y[IDX['C1p']] + \
               y[IDX['BcatC1p']] + y[IDX['BcatC1p_ph']] + y[IDX['Bcat_APC']]
    print(f"\nAPC moiety total (should be ~100nM, conserved): W=0 {apc_moiety(y0_ss):.4f}  W=1 {apc_moiety(y1_ss):.4f}")

    def gsk3_moiety(y):
        return y[IDX['GSK3']] + y[IDX['APC_axin']] + y[IDX['C1']] + y[IDX['C1p']] + \
               y[IDX['BcatC1p']] + y[IDX['BcatC1p_ph']]
    print(f"GSK3 moiety total (should be ~50nM, conserved): W=0 {gsk3_moiety(y0_ss):.4f}  W=1 {gsk3_moiety(y1_ss):.4f}")

    def tcf_moiety(y):
        return y[IDX['TCF']] + y[IDX['Bcat_TCF']]
    print(f"TCF moiety total (should be ~15nM, conserved): W=0 {tcf_moiety(y0_ss):.4f}  W=1 {tcf_moiety(y1_ss):.4f}")

    # APC dysfunction predictions
    K_apc_abund = dict(K);
    y_apc_abund0 = y0(); y_apc_abund0[IDX['APC']] = 100.0 * 0.05  # -95% APC
    sol_a = solve_ivp(rhs, [0, 2e6], y_apc_abund0, args=(0.0, K_apc_abund), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_apc_abund = sol_a.y[IDX['Bcat0'], -1]
    fold_apc_abund = bcat0_apc_abund / bcat0_W0
    print(f"\n-95% APC abundance (W=0): Bcat0_ss = {bcat0_apc_abund:.2f} nM, fold vs WT-W0 = {fold_apc_abund:.2f}x  [recorded: 364.5nM/14.5x]")

    K_apc_bind = dict(K); K_apc_bind['k7'] = K['k7'] * 0.001
    sol_b = solve_ivp(rhs, [0, 2e6], y0(), args=(0.0, K_apc_bind), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_apc_bind = sol_b.y[IDX['Bcat0'], -1]
    fold_apc_bind = bcat0_apc_bind / bcat0_W0
    print(f"k7 x0.001 (Axin-binding collapse, W=0): Bcat0_ss = {bcat0_apc_bind:.2f} nM, fold = {fold_apc_bind:.2f}x  [recorded: 1600.8nM/63.7x]")

    # tankyrase rescue at k7=1%WT (NOTE: distinct background from the 0.1%WT
    # k7-collapse test above -- claim text specifies "k7=1%WT" for this leg)
    K_apc_bind_1pct = dict(K); K_apc_bind_1pct['k7'] = K['k7'] * 0.01
    sol_base1pct = solve_ivp(rhs, [0, 2e6], y0(), args=(0.0, K_apc_bind_1pct), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_base1pct = sol_base1pct.y[IDX['Bcat0'], -1]
    print(f"k7=1%WT baseline (no tankyrase): Bcat0_ss = {bcat0_base1pct:.4f} nM  [recorded: 1243.8nM]")
    K_tank = dict(K_apc_bind_1pct); K_tank['k15'] = K['k15'] * 0.001
    sol_t = solve_ivp(rhs, [0, 2e6], y0(), args=(0.0, K_tank), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_tank = sol_t.y[IDX['Bcat0'], -1]
    print(f"tankyrase rescue (k7=1%WT, k15 x0.001): Bcat0_ss = {bcat0_tank:.4f} nM  [recorded: 2.5nM, from 1243.8nM baseline]")

    # k7=0 exactly, tankyrase rescue should give 0.0% change
    K_k7zero = dict(K_apc_bind); K_k7zero['k7'] = 0.0
    sol_z1 = solve_ivp(rhs, [0, 2e6], y0(), args=(0.0, K_k7zero), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_k7zero_base = sol_z1.y[IDX['Bcat0'], -1]
    K_k7zero_tank = dict(K_k7zero); K_k7zero_tank['k15'] = K['k15'] * 0.001
    sol_z2 = solve_ivp(rhs, [0, 2e6], y0(), args=(0.0, K_k7zero_tank), method='LSODA', rtol=1e-10, atol=1e-14)
    bcat0_k7zero_tank = sol_z2.y[IDX['Bcat0'], -1]
    print(f"k7=0 exactly: baseline Bcat0_ss={bcat0_k7zero_base:.4f}, with tankyrase k15x0.001 Bcat0_ss={bcat0_k7zero_tank:.4f}"
          f"  (claim: bit-for-bit identical, 0.0% rescue)  diff={abs(bcat0_k7zero_base-bcat0_k7zero_tank):.2e}")

    print("\n=== SELFTEST ===")
    ok = True
    if not (res0 < 1e-6 and res1 < 1e-6):
        print("FAIL: steady states not converged"); ok = False
    if rel_err_fold > 5.0:
        print(f"NOTE: fold rel err {rel_err_fold:.2f}% (vs recorded 6.102x)")
    if diff_clamp > 1e-3:
        print(f"FAIL: Axin-clamp adversary NOT bit-identical (diff={diff_clamp:.4e}), contradicts claim"); ok = False
    if abs(bcat0_k7zero_base - bcat0_k7zero_tank) > 1e-3:
        print("FAIL: k7=0 tankyrase-rescue not 0%, contradicts claim"); ok = False
    print("SELFTEST", "PASS" if ok else "FAIL")

if __name__ == "__main__":
    main()
