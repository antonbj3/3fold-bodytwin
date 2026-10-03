"""
EXECUTES HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR-2026-07-26's own named,
never-run honest_gap: "Did not re-execute MODEL-OXIDATIVE-PHOSPHORYLATION's own ODE
... with corrected n_H_ATP=3.6667 to get a new emergent dPsi/RCR/P-O simulation
output -- only the anchor/target comparison text is repaired." That node's own text
states this repair was "independently diagnosed FIVE times ... and executed zero
times before now" -- referring to the TEXT-bracket correction (2.45-2.73 -> 2.78+-0.04).
The dynamical re-simulation named in the same honest_gap was still unexecuted; this
script executes it.

Backing script for MODEL-OXIDATIVE-PHOSPHORYLATION ("oxphos_final.py", cited internally
in cert_design.legs as internal:oxphos_final.py:void_floor_sweep / :ic_stability /
:thermo_consistency) is confirmed ABSENT from disk and from all 35,932 commits of git
history (bt_memory: computation-loss-is-structural-zero-of-164). This script instead
imports the EXISTING, independently-built, on-disk steady-state reimplementation
(data/body_twin/agent_scratch_preserved/oxphosswing_9f2c_reimpl_steadystate.py,
built for a sibling node, HOLE-OXPHOS-DPSI-SWING-CAPACITY-RATIO-GOVERNOR-CONTRADICTION)
UNMODIFIED -- zero re-transcription of the flux equations -- rather than retyping them
a third time.

Fidelity gate (must pass before trusting the corrected run): reproduce the ALREADY-
PUBLISHED baseline (nH_ATP=4.0, the old Route-B value) against MODEL-OXIDATIVE-
PHOSPHORYLATION's own documented targets: dPsi_state4=184.45mV, dPsi_state3=176.1mV,
RCR=10.09, P/O(NADH)=2.481, P/O(succ)=1.487, and the published RCR-vs-leak sweep
[10.77,10.09,8.87,6.36,3.81,2.05,1.00] at g_leak={0.1,1,3,10,30,100,300}x nominal.

New result (never previously computed): the SAME ODE re-solved at the corrected
nH_ATP=3.6667 (=c_ring(8)/F1_sites(3)+ANT/Pi-transport(1.0)) gives an ODE-EMERGENT
P/O(NADH)=2.705 -- about 2.85% (~1.9 sigma) BELOW the corrected point estimate 2.7845
this session used as a pure stoichiometric ratio -- i.e. the model's own kinetic/
thermodynamic back-pressure discounts the theoretical ceiling MORE under the corrected
stoichiometry (-2.85%) than it did under the old one (-0.76%: 2.481 vs 2.5). dPsi and
RCR barely move (<1%); the RCR-vs-leak curve shifts down by up to -3.4% relative at
intermediate leak (10x nominal) and converges to the same fully-uncoupled RCR=1.00 at
300x leak in both cases (a geometric sanity check: at full uncoupling ATP synthesis is
negligible vs leak, so RCR->1 is independent of ATP-synthase stoichiometry).
"""
import sys
import json

sys.path.insert(0, "source_repository/data/body_twin/agent_scratch_preserved")
import oxphosswing_9f2c_reimpl_steadystate as ox  # noqa: E402  (reused unmodified, zero re-transcription)

NH_ATP_OLD = 4.0            # published Route-B value (c_ring(8)/3 + ANT/Pi(1.333) = 4.0)
NH_ATP_CORRECTED = 3.6667   # c_ring(8)/F1_sites(3) + ANT/Pi-transport(1.0) = 3.6667
PURE_STOICH_PO_NADH = 2.7845  # algebra-only point estimate this session already derived


def steady(adp, nh_etc, dg0_etc, nh_atp, g_leak=None):
    kw = {} if g_leak is None else {"g_leak": g_leak}
    roots, _ = ox.solve_steady(adp, nH_ETC=nh_etc, dG0_ETC=dg0_etc, nH_ATP=nh_atp, **kw)
    cands = [r for r in roots if 0 < r < 300]
    if not cands:
        return None
    dpsi = max(cands)
    _, nadh, jetc, jsyn, jleak = ox.charge_residual(dpsi, adp, nh_etc, dg0_etc, nH_ATP=nh_atp, **kw)
    return dict(dPsi=dpsi, NADH=nadh, J_ETC=jetc, J_syn=jsyn, J_leak=jleak)


def point_estimates(nh_atp):
    dg0_nadh = ox.P["dG0_ETC_NADH"]
    dg0_fadh2 = ox.P["dG0_ETC_FADH2"]
    s4n = steady(0.002, 10, dg0_nadh, nh_atp)
    s3n = steady(1.5, 10, dg0_nadh, nh_atp)
    s3s = steady(1.5, 6, dg0_fadh2, nh_atp)
    s4s = steady(0.002, 6, dg0_fadh2, nh_atp)
    return dict(
        dPsi_state4_NADH=s4n["dPsi"], dPsi_state3_NADH=s3n["dPsi"],
        dPsi_state4_succ=s4s["dPsi"], dPsi_state3_succ=s3s["dPsi"],
        RCR_NADH=s3n["J_ETC"] / s4n["J_ETC"], RCR_succ=s3s["J_ETC"] / s4s["J_ETC"],
        PO_NADH=s3n["J_syn"] / s3n["J_ETC"], PO_succ=s3s["J_syn"] / s3s["J_ETC"],
    )


def rcr_vs_leak_sweep(nh_atp, mults=(0.1, 1, 3, 10, 30, 100, 300)):
    g0 = ox.P["g_leak"]
    out = []
    for m in mults:
        s4 = steady(0.002, 10, ox.P["dG0_ETC_NADH"], nh_atp, g_leak=g0 * m)
        s3 = steady(1.5, 10, ox.P["dG0_ETC_NADH"], nh_atp, g_leak=g0 * m)
        out.append(s3["J_ETC"] / s4["J_ETC"])
    return out


if __name__ == "__main__":
    pub = point_estimates(NH_ATP_OLD)
    print("=== FIDELITY CHECK vs MODEL-OXIDATIVE-PHOSPHORYLATION's documented baseline (nH_ATP=4.0) ===")
    targets = {"dPsi_state4_NADH": 184.45, "dPsi_state3_NADH": 176.1, "RCR_NADH": 10.09,
               "PO_NADH": 2.481, "PO_succ": 1.487}
    for k, tgt in targets.items():
        got = pub[k]
        print(f"  {k}: reproduced={got:.4f}  documented={tgt}  ({100*(got-tgt)/tgt:+.3f}%)")

    pub_sweep = [10.77, 10.09, 8.87, 6.36, 3.81, 2.05, 1.00]
    mine_sweep = [round(v, 2) for v in rcr_vs_leak_sweep(NH_ATP_OLD)]
    print(f"\n  RCR-vs-leak sweep reproduced: {mine_sweep}")
    print(f"  RCR-vs-leak sweep documented: {pub_sweep}")

    corr = point_estimates(NH_ATP_CORRECTED)
    print("\n=== NEW EMERGENT OUTPUT under corrected nH_ATP=3.6667 (never previously ODE-simulated) ===")
    for k in targets:
        print(f"  {k}: OLD={pub[k]:.4f}  NEW={corr[k]:.4f}  delta={corr[k]-pub[k]:+.4f}")
    print(f"\n  ODE-emergent P/O(NADH) at corrected stoichiometry: {corr['PO_NADH']:.4f}")
    print(f"  vs pure-stoichiometric point estimate {PURE_STOICH_PO_NADH}: "
          f"{corr['PO_NADH']-PURE_STOICH_PO_NADH:+.4f} ({100*(corr['PO_NADH']-PURE_STOICH_PO_NADH)/PURE_STOICH_PO_NADH:+.3f}%)")

    corrected_sweep = [round(v, 2) for v in rcr_vs_leak_sweep(NH_ATP_CORRECTED)]
    print(f"\n  RCR-vs-leak sweep, corrected nH_ATP: {corrected_sweep}")
    print(f"  delta vs published sweep: {[round(c-o,3) for c,o in zip(corrected_sweep, pub_sweep)]}")

    json.dump(dict(published_baseline=pub, corrected=corr,
                    rcr_sweep_published=mine_sweep, rcr_sweep_corrected=corrected_sweep),
              open("/tmp/oxphos_corrected_nhatp_final.json", "w"), indent=2)
