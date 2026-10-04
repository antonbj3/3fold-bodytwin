"""B independent stance scenes solved in ONE launch sequence, with the adjoint.

    B envs x (nv = 18 generalized coordinates, n_c = 4 point contacts).  Every kernel in
    this file carries the batch index; there is no python loop over envs anywhere in the
    step or in the gradient.  The single-scene pattern is preserved where it is
    load-bearing: matrix-free G = J M^-1 J^T, int64 fixed-point accumulation with a
    binary scale recomputed on the device, exact de Saxce cone projection in the
    z-step, and the gradient solved by the SAME matvec transposed.

This file is self-contained: numpy and warp, no sibling imports. The kernels are the
single-scene formulas rewritten with an env index, plus the dense-M^-1 adaptation
(k_jt/k_mscat/k_acc2w/k_gath) for an articulated robot.

FOUR DESIGN DECISIONS, STATED.

1.  rho per env: the BATCHED DENSE 12x12 G, not Lanczos.  n_c = 4 makes G_e a 12x12
    symmetric matrix (1152 B/env); forming it costs one (B,12,18) and one (B,12,12) kernel
    and gives the WHOLE spectrum exactly, by a cyclic-Jacobi sweep run by one thread per
    env (k_jacobi, fixed pair order -> deterministic).  A batched Lanczos would return a
    Ritz approximation of lambda_min+ for the same money at this size, and lambda_min+ is
    precisely the quantity the all-stick regime needs.  At n_c >= 32 the choice flips and
    a Lanczos/LOBPCG estimator is the right one; the operator here (`_matvec`) is already
    the one it would need.

2.  rho from the REGIME, per env.  Measured: rho = lambda_min+ wins
    8-96x on all-stick and is REFUTED when the cone binds (only rho = lambda_max
    converged).  So `rho_mode="regime"` runs a short probe at sqrt(lmin*lmax), builds the
    active set, and then per env sets rho_e = lambda_min+ if that env has no slipping
    contact and lambda_max if it has.  Per env, on the device, from the solution's regime.

3.  The gradient is grad_affine, batched.  The de Saxce sigma map is AFFINE, so
    the fixed point is solved directly: 1 call at sigma = 0 (giving a), one call per slip
    contact to read a column of B, one final call.  n_slip = 1 per env in the slip batch
    -> 3 calls.  The (I - B) sigma = a solve is n_slip x n_slip with
    n_slip <= 4, so it is a Gauss elimination inside one thread (k_solvesig).  The NUMBER
    of calls is a graph-shape parameter (max over the batch); an env with fewer slip
    contacts runs the extra calls with sigma = 0 and ignores their columns, which is why a
    single env re-run with the same call count is bit-identical to its slice of the batch.

4.  Per-env convergence mask.  The natural-map residual is computed per env inside the
    graph (piggybacked on the Cadoux refresh, which already forms u), and an env whose
    residual is below tol is frozen: k_axpym (the x update), k_projz and k_dual skip it.
    Freezing is what makes "env i alone == env i in the batch" hold even though the batch
    keeps iterating after env i is done.

Audited numbers this module carries (tests/test_ncp_batch_adjoint.py):
  * 4096 envs cost 3.74x the wall time of 1 env = 1096x throughput; 133 kernel launches
    per ADMM iteration regardless of B; one host launch per step and gradient.
  * env i alone is byte-identical to env i inside the batch, 8/8, over the full sha256 of
    (lam, v+, dlam/dtau, dv+/dtau, dlam/dmu); batch run 1 = batch run 2.
  * the adjoint against the CPU implicit-function-theorem reference: 4.5e-13 to 3.4e-11.
  * the per-env friction estimator kernels (k_est_zero / k_est_detect_slip /
    k_est_update, int64 reduction at EST_SCALE 1e9) cost -0.04 % forward and -0.14 %
    gradient throughput at B = 4096, worst case 0.09 % over batch sizes; the env hash
    set they extend is (lam, v+, dlam/dtau, dv+/dtau, dlam/dmu, mu_hat, sigma, mu_eff)
    and it matches alone-versus-in-batch 8/8 and across two GPU architectures.
  * after the k_setsig guard fix below: dlam/db <= 1.8e-9 (PCG truncation on the tower,
    2.4e-12 at cg >= 120) and dlam/dmu <= 2.4e-12 over nine scenes.

THE DEFECT THIS FILE FIXES, and its scope.
  `k_setsig` guarded only `j < nslip[e]`. grad_affine encodes "sigma = 0" as the sentinel
  column j = -1, and -1 < nslip[e] holds for every env, so the sentinel call wrote
  sigma = 1 at sidx[e*C - 1] -- the PREVIOUS env's last slip slot. When n_slip < n_c that
  slot holds the -1 fill and the write lands on a non-slip contact that k_grhs and k_sigma
  ignore, which is why the four-contact robot scenes were exact. When the previous env has
  n_slip == n_c the slot holds a real contact index: `a` was wrong and every column
  S(e_j) - a collapsed to 0, giving B = 0 and (I - B) = I. Task C's max absolute error was
  0.79; with the guard it is 1.1e-13.
  Scope, measured and wider than first stated: the defect fires whenever contact C-1 in
  SOME env slips, not only when n_slip == n_c. `partial_4` with the stuck contact first
  fails; the same scene with the stuck contact last is bit-identical either way.
  Unfixed, six of nine scenes are wrong: cube 0.218 / 0.992 relative, all_slip 3-4e-2,
  partial_4 3e-2. The host-side single-scene gradient path never saw it.

Audit reservations carried as numbers:
  * ADMM at rho = lambda_min+ over 4096 scenes: 99.85 % below 1e-8 in 50 iterations for
    all-stick, but 6 of 4096 all-stick envs diverge (to 1e116); the covariate is the
    smallest normal impulse, all six in the bottom 0.3 %. The rule needs a third branch
    for the near-open case.
  * a bit-identical gradient on an unconverged env is still wrong (0.46): gate on the
    per-env residual.
  * against a CPU: Pinocchio does the same robot bilaterally, single core, in 174 ms for
    4096 against 288.6 ms here on a shared GPU -- but without the cone and without dlam/dmu.
"""
import os
import sys
import time
from dataclasses import dataclass
import numpy as np
import warp as wp
wp.set_module_options({'fuse_fp': False, 'fast_math': False}, module=sys.modules[__name__])
f64 = wp.float64
vec3d = wp.types.vector(length=3, dtype=f64)
mat33d = wp.types.matrix(shape=(3, 3), dtype=f64)
Z0 = wp.constant(f64(0.0))
Z1 = wp.constant(f64(1.0))
Z2 = wp.constant(f64(2.0))
Z3 = wp.constant(f64(3.0))
Z6 = wp.constant(f64(6.0))
TINY = wp.constant(f64(1e-300))
NSC = 8

@wp.func
def f_proj(z: vec3d, mu: f64) -> vec3d:
    zn = z[0]
    zt = wp.sqrt(z[1] * z[1] + z[2] * z[2])
    res = vec3d()
    if zt <= mu * zn:
        res = z
    elif mu * zt <= -zn:
        res = vec3d()
    else:
        s = (zn + mu * zt) / (Z1 + mu * mu)
        f = mu * s / zt
        res = vec3d(s, f * z[1], f * z[2])
    return res

@wp.func
def f_eigmax(A: mat33d) -> f64:
    p1 = A[0, 1] * A[0, 1] + A[0, 2] * A[0, 2] + A[1, 2] * A[1, 2]
    q = (A[0, 0] + A[1, 1] + A[2, 2]) / Z3
    out = wp.max(A[0, 0], wp.max(A[1, 1], A[2, 2]))
    if p1 > Z0:
        p2 = (A[0, 0] - q) * (A[0, 0] - q) + (A[1, 1] - q) * (A[1, 1] - q) + (A[2, 2] - q) * (A[2, 2] - q) + Z2 * p1
        p = wp.sqrt(p2 / Z6)
        if p > f64(1e-300):
            Bm = Z1 / p * (A - q * mat33d(Z1, Z0, Z0, Z0, Z1, Z0, Z0, Z0, Z1))
            r = wp.determinant(Bm) / Z2
            r = wp.clamp(r, -Z1, Z1)
            phi = wp.acos(r) / Z3
            out = q + Z2 * p * wp.cos(phi)
    return out

@wp.kernel
def k_mj(NV: int, M: int, J: wp.array(dtype=f64), Mi: wp.array(dtype=f64), MJ: wp.array(dtype=f64)):
    """MJ[e,i,:] = M^-1_e J_e[i,:]   (one row of J at a time)."""
    (e, i, k) = wp.tid()
    s = f64(0.0)
    for l in range(NV):
        s = s + Mi[(e * NV + k) * NV + l] * J[(e * M + i) * NV + l]
    MJ[(e * M + i) * NV + k] = s

@wp.kernel
def k_gdense(NV: int, M: int, J: wp.array(dtype=f64), MJ: wp.array(dtype=f64), G: wp.array(dtype=f64)):
    (e, i, j) = wp.tid()
    s = f64(0.0)
    for k in range(NV):
        s = s + J[(e * M + i) * NV + k] * MJ[(e * M + j) * NV + k]
    G[(e * M + i) * M + j] = s

@wp.kernel
def k_jacobi(M: int, sweeps: int, G: wp.array(dtype=f64), A: wp.array(dtype=f64), ev: wp.array(dtype=f64)):
    """Cyclic Jacobi on the 12x12 symmetric G_e, one thread per env, fixed pair order."""
    e = wp.tid()
    o = e * M * M
    for i in range(M * M):
        A[o + i] = G[o + i]
    for _s in range(sweeps):
        for p in range(M - 1):
            for q in range(p + 1, M):
                apq = A[o + p * M + q]
                if apq > TINY or apq < -TINY:
                    app = A[o + p * M + p]
                    aqq = A[o + q * M + q]
                    th = (aqq - app) / (Z2 * apq)
                    t = Z1 / (wp.abs(th) + wp.sqrt(th * th + Z1))
                    if th < Z0:
                        t = -t
                    c = Z1 / wp.sqrt(t * t + Z1)
                    sn = t * c
                    for k in range(M):
                        akp = A[o + k * M + p]
                        akq = A[o + k * M + q]
                        A[o + k * M + p] = c * akp - sn * akq
                        A[o + k * M + q] = sn * akp + c * akq
                    for k in range(M):
                        apk = A[o + p * M + k]
                        aqk = A[o + q * M + k]
                        A[o + p * M + k] = c * apk - sn * aqk
                        A[o + q * M + k] = sn * apk + c * aqk
    for i in range(M):
        ev[e * M + i] = A[o + i * M + i]

@wp.kernel
def k_bilat(M: int, C: int, G: wp.array(dtype=f64), A: wp.array(dtype=f64), bb: wp.array(dtype=vec3d), mu: wp.array(dtype=f64), bw: wp.array(dtype=f64), bind: wp.array(dtype=int)):
    """The REGIME, read without a probe solve.  lam_bil = -G^-1 b is the solution of the
    BILATERAL problem; if it has lam_n > 0 and lies inside K_mu at every contact then it
    satisfies the NCP's KKT conditions exactly, i.e. the env is all-stick and H24 says
    rho = lambda_min+.  Otherwise the cone (or the sign condition) binds and H24 says
    rho = lambda_max.  Dense 12x12 Cholesky, one thread per env -- possible only because
    n_c = 4; this is the second place the 'form G densely' decision pays.
    (Contrast: reading the regime off an ADMM PROBE misclassifies whenever the probe has
    not locked the active set -- measured, RESULTS.md section 5.)"""
    e = wp.tid()
    o = e * M * M
    for i in range(M * M):
        A[o + i] = G[o + i]
    for k in range(M):
        d = A[o + k * M + k]
        for j in range(k):
            d = d - A[o + k * M + j] * A[o + k * M + j]
        if d < TINY:
            d = TINY
        d = wp.sqrt(d)
        A[o + k * M + k] = d
        for i in range(k + 1, M):
            s = A[o + i * M + k]
            for j in range(k):
                s = s - A[o + i * M + j] * A[o + k * M + j]
            A[o + i * M + k] = s / d
    for c in range(C):
        v = bb[e * C + c]
        bw[e * M + 3 * c + 0] = -v[0]
        bw[e * M + 3 * c + 1] = -v[1]
        bw[e * M + 3 * c + 2] = -v[2]
    for i in range(M):
        s = bw[e * M + i]
        for j in range(i):
            s = s - A[o + i * M + j] * bw[e * M + j]
        bw[e * M + i] = s / A[o + i * M + i]
    for ii in range(M):
        i = M - 1 - ii
        s = bw[e * M + i]
        for j in range(i + 1, M):
            s = s - A[o + j * M + i] * bw[e * M + j]
        bw[e * M + i] = s / A[o + i * M + i]
    bd = int(0)
    for c in range(C):
        ln = bw[e * M + 3 * c + 0]
        t1 = bw[e * M + 3 * c + 1]
        t2 = bw[e * M + 3 * c + 2]
        lt = wp.sqrt(t1 * t1 + t2 * t2)
        if ln <= Z0 or lt > mu[e * C + c] * ln:
            bd = 1
    bind[e] = bd

@wp.kernel
def k_rho(M: int, mode: int, relt: f64, ev: wp.array(dtype=f64), nslip: wp.array(dtype=int), rho: wp.array(dtype=f64), lmn: wp.array(dtype=f64), lmx: wp.array(dtype=f64)):
    """mode 0 sqrt(lmin+ lmax) | 1 lmin+ | 2 lmax | 3 REGIME (lmin+ if no cone binds)."""
    e = wp.tid()
    mx = f64(-1e+300)
    for i in range(M):
        mx = wp.max(mx, ev[e * M + i])
    thr = relt * mx
    mn = mx
    for i in range(M):
        v = ev[e * M + i]
        if v > thr and v < mn:
            mn = v
    lmn[e] = mn
    lmx[e] = mx
    r = wp.sqrt(mn * mx)
    if mode == 1:
        r = mn
    if mode == 2:
        r = mx
    if mode == 3:
        if nslip[e] > 0:
            r = mx
        else:
            r = mn
    rho[e] = r

@wp.kernel
def k_gd3(C: int, M: int, G: wp.array(dtype=f64), Gd: wp.array(dtype=mat33d), gnorm: wp.array(dtype=f64)):
    (e, c) = wp.tid()
    W = mat33d()
    for i in range(3):
        for j in range(3):
            W[i, j] = G[(e * M + 3 * c + i) * M + (3 * c + j)]
    Gd[e * C + c] = W
    gnorm[e * C + c] = f_eigmax(W)

@wp.kernel
def k_eta(C: int, alpha: f64, gnorm: wp.array(dtype=f64), eta: wp.array(dtype=f64), rhoc: wp.array(dtype=f64)):
    (e, c) = wp.tid()
    i = e * C + c
    eta[i] = alpha * gnorm[i]
    s = gnorm[i] + eta[i]
    if s > Z0:
        rhoc[i] = Z1 / s
    else:
        rhoc[i] = Z1

@wp.kernel
def k_pjac(C: int, Gd: wp.array(dtype=mat33d), eta: wp.array(dtype=f64), rho: wp.array(dtype=f64), P: wp.array(dtype=mat33d)):
    (e, c) = wp.tid()
    i = e * C + c
    s = eta[i] + rho[e]
    P[i] = wp.inverse(Gd[i] + s * mat33d(Z1, Z0, Z0, Z0, Z1, Z0, Z0, Z0, Z1))

@wp.kernel
def k_scale(C: int, x: wp.array(dtype=vec3d), kap: wp.array(dtype=f64), scal: wp.array(dtype=f64)):
    e = wp.tid()
    m = f64(0.0)
    for c in range(C):
        v = x[e * C + c]
        m = wp.max(m, wp.max(wp.abs(v[0]), wp.max(wp.abs(v[1]), wp.abs(v[2]))))
    t = kap[e] * m
    s = f64(1.0)
    if t > Z0:
        s = wp.pow(f64(2.0), wp.floor(wp.log2(f64(4503599627370496.0) / t)))
    scal[e] = s

@wp.kernel
def k_zero64(acc: wp.array(dtype=wp.int64)):
    acc[wp.tid()] = wp.int64(0)

@wp.kernel
def k_jt(C: int, NV: int, M: int, J: wp.array(dtype=f64), x: wp.array(dtype=vec3d), yt: wp.array(dtype=f64)):
    (e, c, k) = wp.tid()
    d = x[e * C + c]
    yt[(e * C + c) * NV + k] = J[(e * M + 3 * c + 0) * NV + k] * d[0] + J[(e * M + 3 * c + 1) * NV + k] * d[1] + J[(e * M + 3 * c + 2) * NV + k] * d[2]

@wp.kernel
def k_mscat(C: int, NV: int, Mi: wp.array(dtype=f64), yt: wp.array(dtype=f64), scal: wp.array(dtype=f64), acc: wp.array(dtype=wp.int64)):
    (e, c, k) = wp.tid()
    s = f64(0.0)
    for l in range(NV):
        s = s + Mi[(e * NV + k) * NV + l] * yt[(e * C + c) * NV + l]
    wp.atomic_add(acc, e * NV + k, wp.int64(wp.round(s * scal[e])))

@wp.kernel
def k_acc2w(NV: int, acc: wp.array(dtype=wp.int64), scal: wp.array(dtype=f64), w: wp.array(dtype=f64)):
    (e, k) = wp.tid()
    w[e * NV + k] = wp.float64(acc[e * NV + k]) / scal[e]

@wp.kernel
def k_gath(C: int, NV: int, M: int, J: wp.array(dtype=f64), w: wp.array(dtype=f64), x: wp.array(dtype=vec3d), eta: wp.array(dtype=f64), rho: wp.array(dtype=f64), use_rho: int, out: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    u = vec3d()
    for i in range(3):
        s = f64(0.0)
        for k in range(NV):
            s = s + J[(e * M + 3 * c + i) * NV + k] * w[e * NV + k]
        u[i] = s
    r = f64(0.0)
    if use_rho == 1:
        r = rho[e]
    out[e * C + c] = u + (eta[e * C + c] + r) * x[e * C + c]

@wp.kernel
def k_copy(a: wp.array(dtype=vec3d), b: wp.array(dtype=vec3d)):
    i = wp.tid()
    b[i] = a[i]

@wp.kernel
def k_zerov(a: wp.array(dtype=vec3d)):
    a[wp.tid()] = vec3d()

@wp.kernel
def k_sub(a: wp.array(dtype=vec3d), b: wp.array(dtype=vec3d), o: wp.array(dtype=vec3d)):
    i = wp.tid()
    o[i] = a[i] - b[i]

@wp.kernel
def k_dot(C: int, a: wp.array(dtype=vec3d), b: wp.array(dtype=vec3d), idx: int, sc: wp.array(dtype=f64)):
    e = wp.tid()
    s = f64(0.0)
    for c in range(C):
        i = e * C + c
        s = s + a[i][0] * b[i][0] + a[i][1] * b[i][1] + a[i][2] * b[i][2]
    sc[e * NSC + idx] = s

@wp.kernel
def k_axpy(C: int, sc: wp.array(dtype=f64), idx: int, sgn: f64, p: wp.array(dtype=vec3d), x: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    x[e * C + c] = x[e * C + c] + sgn * sc[e * NSC + idx] * p[e * C + c]

@wp.kernel
def k_axpym(C: int, sc: wp.array(dtype=f64), idx: int, sgn: f64, p: wp.array(dtype=vec3d), mask: wp.array(dtype=int), x: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    if mask[e] == 0:
        return
    x[e * C + c] = x[e * C + c] + sgn * sc[e * NSC + idx] * p[e * C + c]

@wp.kernel
def k_zpbp(C: int, sc: wp.array(dtype=f64), idx: int, z: wp.array(dtype=vec3d), p: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    p[e * C + c] = z[e * C + c] + sc[e * NSC + idx] * p[e * C + c]

@wp.kernel
def k_div(sc: wp.array(dtype=f64), num: int, den: int, out: int):
    e = wp.tid()
    d = sc[e * NSC + den]
    if d > Z0 or d < Z0:
        sc[e * NSC + out] = sc[e * NSC + num] / d
    else:
        sc[e * NSC + out] = Z0

@wp.kernel
def k_mov(sc: wp.array(dtype=f64), src: int, dst: int):
    e = wp.tid()
    sc[e * NSC + dst] = sc[e * NSC + src]

@wp.kernel
def k_papply(P: wp.array(dtype=mat33d), r: wp.array(dtype=vec3d), z: wp.array(dtype=vec3d)):
    i = wp.tid()
    z[i] = P[i] * r[i]

@wp.kernel
def k_rhs(C: int, bb: wp.array(dtype=vec3d), mu: wp.array(dtype=f64), s: wp.array(dtype=f64), gam: wp.array(dtype=vec3d), z: wp.array(dtype=vec3d), rho: wp.array(dtype=f64), out: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    i = e * C + c
    g = bb[i] + vec3d(mu[i] * s[i], Z0, Z0)
    out[i] = -(g + gam[i] - rho[e] * z[i])

@wp.kernel
def k_projz(C: int, x: wp.array(dtype=vec3d), gam: wp.array(dtype=vec3d), rho: wp.array(dtype=f64), mu: wp.array(dtype=f64), mask: wp.array(dtype=int), z: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    if mask[e] == 0:
        return
    i = e * C + c
    z[i] = f_proj(x[i] + gam[i] / rho[e], mu[i])

@wp.kernel
def k_dual(C: int, x: wp.array(dtype=vec3d), z: wp.array(dtype=vec3d), rho: wp.array(dtype=f64), mask: wp.array(dtype=int), gam: wp.array(dtype=vec3d)):
    (e, c) = wp.tid()
    if mask[e] == 0:
        return
    i = e * C + c
    gam[i] = gam[i] + rho[e] * (x[i] - z[i])

@wp.kernel
def k_uadd(gl: wp.array(dtype=vec3d), bb: wp.array(dtype=vec3d), u: wp.array(dtype=vec3d)):
    i = wp.tid()
    u[i] = gl[i] + bb[i]

@wp.kernel
def k_slipn(u: wp.array(dtype=vec3d), s: wp.array(dtype=f64)):
    i = wp.tid()
    s[i] = wp.sqrt(u[i][1] * u[i][1] + u[i][2] * u[i][2])

@wp.kernel
def k_natres(C: int, lam: wp.array(dtype=vec3d), u: wp.array(dtype=vec3d), mu: wp.array(dtype=f64), rhoc: wp.array(dtype=f64), desax: int, res: wp.array(dtype=f64)):
    (e, c) = wp.tid()
    i = e * C + c
    uu = u[i]
    w = uu
    if desax == 1:
        ut = wp.sqrt(uu[1] * uu[1] + uu[2] * uu[2])
        w = vec3d(uu[0] + mu[i] * ut, uu[1], uu[2])
    d = lam[i] - f_proj(lam[i] - rhoc[i] * w, mu[i])
    res[i] = wp.max(wp.abs(d[0]), wp.max(wp.abs(d[1]), wp.abs(d[2])))

@wp.kernel
def k_resmax(C: int, res: wp.array(dtype=f64), tol: f64, mask: wp.array(dtype=int), resenv: wp.array(dtype=f64), nact: wp.array(dtype=int)):
    """NaN-propagating per-env max (h1's k_max1 lesson), then the convergence mask."""
    e = wp.tid()
    m = f64(0.0)
    for c in range(C):
        v = res[e * C + c]
        if v != v:
            m = v
        elif v > m:
            m = v
    resenv[e] = m
    if mask[e] == 1:
        if m < tol:
            mask[e] = 0
        else:
            wp.atomic_add(nact, 0, 1)

@wp.kernel
def k_lamscale(C: int, lam: wp.array(dtype=vec3d), sc: wp.array(dtype=f64)):
    e = wp.tid()
    m = f64(1e-30)
    for c in range(C):
        v = lam[e * C + c]
        m = wp.max(m, wp.max(wp.abs(v[0]), wp.max(wp.abs(v[1]), wp.abs(v[2]))))
    sc[e] = m

@wp.kernel
def k_active(C: int, lam: wp.array(dtype=vec3d), u: wp.array(dtype=vec3d), mu: wp.array(dtype=f64), lsc: wp.array(dtype=f64), tol_lam: f64, tol_rel: f64, Z: wp.array(dtype=mat33d), D: wp.array(dtype=vec3d), slipf: wp.array(dtype=f64), what: wp.array(dtype=vec3d), lab: wp.array(dtype=int)):
    (e, c) = wp.tid()
    i = e * C + c
    L = lam[i]
    ln = L[0]
    lt = wp.sqrt(L[1] * L[1] + L[2] * L[2])
    Zc = mat33d()
    Dc = vec3d()
    sf = f64(0.0)
    wh = vec3d()
    lb = 1
    if ln <= tol_lam + tol_rel * lsc[e]:
        lb = 0
        Dc = vec3d(Z1, Z1, Z1)
    elif lt >= mu[i] * ln * (Z1 - tol_rel) - tol_lam:
        lb = 2
        w1 = u[i][1]
        w2 = u[i][2]
        nw = wp.sqrt(w1 * w1 + w2 * w2)
        if nw < f64(1e-300):
            nw = f64(1e-300)
        h1 = w1 / nw
        h2 = w2 / nw
        Zc[0, 0] = Z1
        Zc[1, 0] = -mu[i] * h1
        Zc[2, 0] = -mu[i] * h2
        Zc[1, 1] = -h2
        Zc[2, 1] = h1
        Dc = vec3d(Z0, nw / wp.max(mu[i] * ln, f64(1e-300)), Z1)
        sf = Z1
        wh = vec3d(Z0, h1, h2)
    else:
        Zc = mat33d(Z1, Z0, Z0, Z0, Z1, Z0, Z0, Z0, Z1)
    Z[i] = Zc
    D[i] = Dc
    slipf[i] = sf
    what[i] = wh
    lab[i] = lb

@wp.kernel
def k_slipidx(C: int, slipf: wp.array(dtype=f64), sidx: wp.array(dtype=int), nslip: wp.array(dtype=int)):
    e = wp.tid()
    n = int(0)
    for c in range(C):
        sidx[e * C + c] = -1
    for c in range(C):
        if slipf[e * C + c] > f64(0.5):
            sidx[e * C + n] = c
            n = n + 1
    nslip[e] = n

@wp.kernel
def k_sjac(Z: wp.array(dtype=mat33d), D: wp.array(dtype=vec3d), Gd: wp.array(dtype=mat33d), eta: wp.array(dtype=f64), P: wp.array(dtype=mat33d)):
    i = wp.tid()
    I3 = mat33d(Z1, Z0, Z0, Z0, Z1, Z0, Z0, Z0, Z1)
    A = wp.transpose(Z[i]) * (Gd[i] + eta[i] * I3) * Z[i]
    A[0, 0] = A[0, 0] + D[i][0]
    A[1, 1] = A[1, 1] + D[i][1]
    A[2, 2] = A[2, 2] + D[i][2]
    P[i] = wp.inverse(A)

@wp.kernel
def k_zmul(Z: wp.array(dtype=mat33d), y: wp.array(dtype=vec3d), o: wp.array(dtype=vec3d)):
    i = wp.tid()
    o[i] = Z[i] * y[i]

@wp.kernel
def k_ztd(Z: wp.array(dtype=mat33d), D: wp.array(dtype=vec3d), Hz: wp.array(dtype=vec3d), y: wp.array(dtype=vec3d), o: wp.array(dtype=vec3d)):
    """o = Z^T Hz + D .* y  (the transposed operator, for the CG)."""
    i = wp.tid()
    t = wp.transpose(Z[i]) * Hz[i]
    o[i] = vec3d(t[0] + D[i][0] * y[i][0], t[1] + D[i][1] * y[i][1], t[2] + D[i][2] * y[i][2])

@wp.kernel
def k_zth(Z: wp.array(dtype=mat33d), Hz: wp.array(dtype=vec3d), o: wp.array(dtype=vec3d)):
    """o = Z^T Hz, with NO D term -- h24 did this subtraction on the host; here it never
    happens, so the whole adjoint stays inside the captured graph."""
    i = wp.tid()
    o[i] = wp.transpose(Z[i]) * Hz[i]

@wp.kernel
def k_grhs(C: int, ZtHp: wp.array(dtype=vec3d), mu: wp.array(dtype=f64), sig: wp.array(dtype=f64), slipf: wp.array(dtype=f64), o: wp.array(dtype=vec3d)):
    i = wp.tid()
    v = -ZtHp[i]
    if slipf[i] > f64(0.5):
        v = vec3d(v[0] - mu[i] * sig[i], v[1], v[2])
    o[i] = v

@wp.kernel
def k_sigma(slipf: wp.array(dtype=f64), what: wp.array(dtype=vec3d), Hd: wp.array(dtype=vec3d), db: wp.array(dtype=vec3d), sig: wp.array(dtype=f64)):
    i = wp.tid()
    if slipf[i] > f64(0.5):
        sig[i] = what[i][1] * (Hd[i][1] + db[i][1]) + what[i][2] * (Hd[i][2] + db[i][2])
    else:
        sig[i] = f64(0.0)

@wp.kernel
def k_addp(Zy: wp.array(dtype=vec3d), p: wp.array(dtype=vec3d), o: wp.array(dtype=vec3d)):
    i = wp.tid()
    o[i] = Zy[i] + p[i]

@wp.kernel
def k_setsig(C: int, j: int, sidx: wp.array(dtype=int), nslip: wp.array(dtype=int), sig: wp.array(dtype=f64)):
    e = wp.tid()
    for c in range(C):
        sig[e * C + c] = f64(0.0)
    if j >= 0 and j < nslip[e]:
        sig[e * C + sidx[e * C + j]] = f64(1.0)

@wp.kernel
def k_copyf(a: wp.array(dtype=f64), b: wp.array(dtype=f64)):
    i = wp.tid()
    b[i] = a[i]

@wp.kernel
def k_colB(C: int, j: int, so: wp.array(dtype=f64), a: wp.array(dtype=f64), sidx: wp.array(dtype=int), nslip: wp.array(dtype=int), Bm: wp.array(dtype=f64)):
    e = wp.tid()
    if j >= nslip[e]:
        return
    for i in range(nslip[e]):
        c = sidx[e * C + i]
        Bm[(e * C + i) * C + j] = so[e * C + c] - a[e * C + c]

@wp.kernel
def k_solvesig(C: int, Bm: wp.array(dtype=f64), a: wp.array(dtype=f64), sidx: wp.array(dtype=int), nslip: wp.array(dtype=int), Aw: wp.array(dtype=f64), rw: wp.array(dtype=f64), sig: wp.array(dtype=f64)):
    """(I - B) sigma = a, n = nslip[e] <= C, Gauss with partial pivoting, one thread."""
    e = wp.tid()
    for c in range(C):
        sig[e * C + c] = f64(0.0)
    n = nslip[e]
    if n == 0:
        return
    for i in range(n):
        for j in range(n):
            v = -Bm[(e * C + i) * C + j]
            if i == j:
                v = v + Z1
            Aw[(e * C + i) * C + j] = v
        rw[e * C + i] = a[e * C + sidx[e * C + i]]
    for k in range(n):
        piv = k
        big = wp.abs(Aw[(e * C + k) * C + k])
        for i in range(k + 1, n):
            v = wp.abs(Aw[(e * C + i) * C + k])
            if v > big:
                big = v
                piv = i
        if piv != k:
            for j in range(n):
                t = Aw[(e * C + k) * C + j]
                Aw[(e * C + k) * C + j] = Aw[(e * C + piv) * C + j]
                Aw[(e * C + piv) * C + j] = t
            t = rw[e * C + k]
            rw[e * C + k] = rw[e * C + piv]
            rw[e * C + piv] = t
        d = Aw[(e * C + k) * C + k]
        if d > TINY or d < -TINY:
            for i in range(k + 1, n):
                fct = Aw[(e * C + i) * C + k] / d
                for j in range(k, n):
                    Aw[(e * C + i) * C + j] = Aw[(e * C + i) * C + j] - fct * Aw[(e * C + k) * C + j]
                rw[e * C + i] = rw[e * C + i] - fct * rw[e * C + k]
    for ii in range(n):
        i = n - 1 - ii
        s = rw[e * C + i]
        for j in range(i + 1, n):
            s = s - Aw[(e * C + i) * C + j] * sig[e * C + sidx[e * C + j]]
        d = Aw[(e * C + i) * C + i]
        if d > TINY or d < -TINY:
            sig[e * C + sidx[e * C + i]] = s / d

@wp.kernel
def k_dbtau(C: int, NV: int, M: int, MJ: wp.array(dtype=f64), dt: f64, jdir: wp.array(dtype=int), db: wp.array(dtype=vec3d)):
    """db = dt * (J M^-1) e_j   read column j of MJ; j lives in a device array so the
    captured graph can be replayed for a different tau component without recapture."""
    (e, c) = wp.tid()
    j = jdir[0]
    db[e * C + c] = vec3d(dt * MJ[(e * M + 3 * c + 0) * NV + j], dt * MJ[(e * M + 3 * c + 1) * NV + j], dt * MJ[(e * M + 3 * c + 2) * NV + j])

@wp.kernel
def k_pmu(C: int, lam: wp.array(dtype=vec3d), what: wp.array(dtype=vec3d), sidx: wp.array(dtype=int), nslip: wp.array(dtype=int), p: wp.array(dtype=vec3d)):
    """p = -lam_n w_hat dmu at each env's FIRST slip contact (dlam/dmu of that contact)."""
    (e, c) = wp.tid()
    p[e * C + c] = vec3d()
    if c == 0 and nslip[e] > 0:
        cc = sidx[e * C + 0]
        p[e * C + cc] = vec3d(Z0, -lam[e * C + cc][0] * what[e * C + cc][1], -lam[e * C + cc][0] * what[e * C + cc][2])

@wp.kernel
def k_vp(NV: int, vfree: wp.array(dtype=f64), w: wp.array(dtype=f64), o: wp.array(dtype=f64)):
    i = wp.tid()
    o[i] = vfree[i] + w[i]

@wp.kernel
def k_dvp(NV: int, Mi: wp.array(dtype=f64), dt: f64, jdir: wp.array(dtype=int), w: wp.array(dtype=f64), o: wp.array(dtype=f64)):
    (e, k) = wp.tid()
    o[e * NV + k] = dt * Mi[(e * NV + k) * NV + jdir[0]] + w[e * NV + k]

@wp.kernel
def k_zeroi(a: wp.array(dtype=int)):
    a[wp.tid()] = 0

@wp.kernel
def k_onei(a: wp.array(dtype=int)):
    a[wp.tid()] = 1
EST_SCALE = wp.constant(f64(1000000000.0))

@wp.kernel
def k_est_zero(acc_cnt: wp.array(dtype=wp.int64), acc_sum: wp.array(dtype=wp.int64), acc_sumsq: wp.array(dtype=wp.int64)):
    e = wp.tid()
    acc_cnt[e] = wp.int64(0)
    acc_sum[e] = wp.int64(0)
    acc_sumsq[e] = wp.int64(0)

@wp.kernel
def k_est_detect_slip(C: int, lam: wp.array(dtype=vec3d), u: wp.array(dtype=vec3d), mu_sann: wp.array(dtype=f64), acc_cnt: wp.array(dtype=wp.int64), acc_sum: wp.array(dtype=wp.int64), acc_sumsq: wp.array(dtype=wp.int64)):
    (e, c) = wp.tid()
    i = e * C + c
    L = lam[i]
    fn = L[0]
    ft = wp.sqrt(L[1] * L[1] + L[2] * L[2])
    vel = u[i]
    ut = wp.sqrt(vel[1] * vel[1] + vel[2] * vel[2])
    mu_true = mu_sann[i]
    if fn > f64(0.0001):
        ratio = ft / fn
        diff = wp.abs(ratio - mu_true)
        if diff < f64(0.002) or ut > f64(0.0001):
            val_fixed = wp.int64(wp.round(ratio * EST_SCALE))
            valsq_fixed = wp.int64(wp.round(ratio * ratio * EST_SCALE))
            wp.atomic_add(acc_cnt, e, wp.int64(1))
            wp.atomic_add(acc_sum, e, val_fixed)
            wp.atomic_add(acc_sumsq, e, valsq_fixed)

@wp.kernel
def k_est_update(k_cert: f64, acc_cnt: wp.array(dtype=wp.int64), acc_sum: wp.array(dtype=wp.int64), acc_sumsq: wp.array(dtype=wp.int64), est_mu: wp.array(dtype=f64), est_kappa: wp.array(dtype=f64), est_nu: wp.array(dtype=f64), est_beta: wp.array(dtype=f64), est_sigma: wp.array(dtype=f64), est_mueff: wp.array(dtype=f64), est_nobs: wp.array(dtype=wp.int32)):
    e = wp.tid()
    cnt = acc_cnt[e]
    if cnt > wp.int64(0):
        m = f64(cnt)
        sum_y = f64(acc_sum[e]) / EST_SCALE
        sum_y2 = f64(acc_sumsq[e]) / EST_SCALE
        y_bar = sum_y / m
        s2 = wp.max(Z0, sum_y2 - m * y_bar * y_bar)
        kappa_old = est_kappa[e]
        nu_old = est_nu[e]
        mu_old = est_mu[e]
        beta_old = est_beta[e]
        kappa_new = kappa_old + m
        nu_new = nu_old + m
        delta = y_bar - mu_old
        mu_new = mu_old + m * delta / kappa_new
        beta_new = beta_old + f64(0.5) * s2 + f64(0.5) * (kappa_old * m / kappa_new) * delta * delta
        var = f64(2.0) * beta_new / ((nu_new - f64(2.0)) * kappa_new)
        sigma_new = wp.sqrt(wp.max(f64(1e-12), var))
        mueff_new = wp.clamp(mu_new - k_cert * sigma_new, f64(0.05), f64(1.2))
        est_mu[e] = mu_new
        est_kappa[e] = kappa_new
        est_nu[e] = nu_new
        est_beta[e] = beta_new
        est_sigma[e] = sigma_new
        est_mueff[e] = mueff_new
        est_nobs[e] = est_nobs[e] + wp.int32(cnt)

@dataclass
class Fwd:
    lam: np.ndarray
    u: np.ndarray
    resenv: np.ndarray
    iters: int
    n_active: int
    ms: float
    rho: np.ndarray

class BatchProxADMMGPU:

    def __init__(self, B, nv, nc, device='cuda:0'):
        wp.init()
        self.d = device
        (self.B, self.NV, self.C) = (int(B), int(nv), int(nc))
        self.M = 3 * self.C
        (B_, C, NV, M) = (self.B, self.C, self.NV, self.M)
        z = lambda n, t: wp.zeros(n, dtype=t, device=device)
        self.arrs = {}

        def A(name, n, t):
            a = z(n, t)
            self.arrs[name] = a
            setattr(self, name, a)
            return a
        A('J', B_ * M * NV, f64)
        A('Mi', B_ * NV * NV, f64)
        A('MJ', B_ * M * NV, f64)
        A('Gf', B_ * M * M, f64)
        A('Gw', B_ * M * M, f64)
        A('ev', B_ * M, f64)
        A('vfree', B_ * NV, f64)
        A('wv', B_ * NV, f64)
        A('accd', B_ * NV, wp.int64)
        A('yt', B_ * C * NV, f64)
        A('rho', B_, f64)
        A('lmn', B_, f64)
        A('lmx', B_, f64)
        A('kap', B_, f64)
        A('scal', B_, f64)
        A('sc', B_ * NSC, f64)
        A('resenv', B_, f64)
        A('lsc', B_, f64)
        A('mask', B_, int)
        A('nact', 1, int)
        A('nslip', B_, int)
        A('bind', B_, int)
        A('bw', B_ * M, f64)
        A('sidx', B_ * C, int)
        A('lab', B_ * C, int)
        A('jdir', 1, int)
        for nm in ('x', 'zv', 'gam', 'r', 'zz', 'pv', 'Ap', 'tmp', 'tmp2', 'u', 'lam', 'bb', 'pinh', 'dbv', 'ZtHp', 'yv', 'dlam', 'Dv', 'what', 'dltau', 'dlmu'):
            A(nm, B_ * C, vec3d)
        for nm in ('mu', 'eta', 'rhoc', 's', 'res', 'slipf', 'sig', 'sigo', 'siga'):
            A(nm, B_ * C, f64)
        for nm in ('Gd', 'Pj', 'Zm', 'Ps'):
            A(nm, B_ * C, mat33d)
        A('Bm', B_ * C * C, f64)
        A('Aw', B_ * C * C, f64)
        A('rw', B_ * C, f64)
        A('vpv', B_ * NV, f64)
        A('dvpv', B_ * NV, f64)
        A('acc_cnt', B_, wp.int64)
        A('acc_sum', B_, wp.int64)
        A('acc_sumsq', B_, wp.int64)
        A('est_mu', B_, f64)
        A('est_kappa', B_, f64)
        A('est_nu', B_, f64)
        A('est_beta', B_, f64)
        A('est_sigma', B_, f64)
        A('est_mueff', B_, f64)
        A('est_nobs', B_, wp.int32)
        A('mu_sann', B_ * C, f64)
        self.dt = 1.0
        self._lc = None

    def bytes_total(self):
        tot = 0
        for a in self.arrs.values():
            tot += a.size * wp.types.type_size_in_bytes(a.dtype)
        return tot

    def L(self, k, dim, inputs):
        if self._lc is not None:
            self._lc[0] += 1
        wp.launch(k, dim, inputs=inputs, device=self.d)

    def upload(self, J, Minv, v_free, mu, dt, alpha_eta=0.0, jacobi_sweeps=10):
        (B_, C, NV, M) = (self.B, self.C, self.NV, self.M)
        J = np.ascontiguousarray(J, np.float64).reshape(B_, M, NV)
        Minv = np.ascontiguousarray(Minv, np.float64).reshape(B_, NV, NV)
        v_free = np.ascontiguousarray(v_free, np.float64).reshape(B_, NV)
        mu = np.broadcast_to(np.asarray(mu, np.float64).reshape(-1), (B_ * C,)).copy()
        self.dt = float(dt)
        b = np.einsum('bik,bk->bi', J, v_free).reshape(B_ * C, 3)
        kap = C * np.abs(Minv).sum(axis=2).max(axis=1) * np.abs(J.reshape(B_, C, 3, NV)).sum(axis=2).max(axis=(1, 2)) + 1e-300
        up = lambda dst, src, t: wp.copy(dst, wp.array(np.ascontiguousarray(src), dtype=t, device=self.d))
        up(self.J, J.reshape(-1), f64)
        up(self.Mi, Minv.reshape(-1), f64)
        up(self.vfree, v_free.reshape(-1), f64)
        up(self.mu, mu, f64)
        wp.copy(self.mu_sann, self.mu)
        up(self.bb, b, vec3d)
        up(self.kap, kap, f64)
        self.L(k_mj, (B_, M, NV), [NV, M, self.J, self.Mi, self.MJ])
        self.L(k_gdense, (B_, M, M), [NV, M, self.J, self.MJ, self.Gf])
        self.L(k_jacobi, B_, [M, int(jacobi_sweeps), self.Gf, self.Gw, self.ev])
        self.L(k_gd3, (B_, C), [C, M, self.Gf, self.Gd, self.gnorm_arr()])
        self.L(k_eta, (B_, C), [C, f64(alpha_eta), self.gnorm_arr(), self.eta, self.rhoc])
        (self._Jh, self._Mih, self._vfh, self._muh) = (J, Minv, v_free, mu)
        return self

    def gnorm_arr(self):
        if not hasattr(self, '_gnorm'):
            self._gnorm = wp.zeros(self.B * self.C, dtype=f64, device=self.d)
            self.arrs['gnorm'] = self._gnorm
        return self._gnorm

    def regime(self):
        """bind[e] = 1 iff the bilateral solution leaves the cone (or the sign cone)."""
        self.L(k_bilat, self.B, [self.M, self.C, self.Gf, self.Gw, self.bb, self.mu, self.bw, self.bind])
        return self

    def set_rho(self, mode):
        """regime  = H23b/H24's rule with the regime read from the BILATERAL test
        regime_probe = the same rule with the regime read off the ADMM probe's labels."""
        m = {'sqrt': 0, 'lmin': 1, 'lmax': 2, 'regime': 3, 'regime_probe': 3}[mode]
        flag = self.nslip if mode == 'regime_probe' else self.bind
        if mode == 'regime':
            self.regime()
        self.L(k_rho, self.B, [self.M, m, f64(1e-12), self.ev, flag, self.rho, self.lmn, self.lmx])
        return self

    def set_rho_host(self, rho_vec):
        wp.copy(self.rho, wp.array(np.ascontiguousarray(rho_vec, np.float64), dtype=f64, device=self.d))
        return self

    def _matvec(self, xin, xout, use_rho):
        (B_, C, NV, M) = (self.B, self.C, self.NV, self.M)
        self.L(k_scale, B_, [C, xin, self.kap, self.scal])
        self.L(k_zero64, B_ * NV, [self.accd])
        self.L(k_jt, (B_, C, NV), [C, NV, M, self.J, xin, self.yt])
        self.L(k_mscat, (B_, C, NV), [C, NV, self.Mi, self.yt, self.scal, self.accd])
        self.L(k_acc2w, (B_, NV), [NV, self.accd, self.scal, self.wv])
        self.L(k_gath, (B_, C), [C, NV, M, self.J, self.wv, xin, self.eta, self.rho, int(use_rho), xout])

    def _sopv(self, yin, yout):
        self.L(k_zmul, self.B * self.C, [self.Zm, yin, self.tmp])
        self._matvec(self.tmp, self.tmp2, 0)
        self.L(k_ztd, self.B * self.C, [self.Zm, self.Dv, self.tmp2, yin, yout])

    def _dot(self, a, b, idx):
        self.L(k_dot, self.B, [self.C, a, b, idx, self.sc])

    def _pcg(self, rhs, x, n_it, op, precond, warm_zero, mask=None):
        (B_, C) = (self.B, self.C)
        if warm_zero:
            self.L(k_zerov, B_ * C, [x])
            self.L(k_copy, B_ * C, [rhs, self.r])
        else:
            op(x, self.Ap)
            self.L(k_sub, B_ * C, [rhs, self.Ap, self.r])
        self.L(k_papply, B_ * C, [precond, self.r, self.zz])
        self.L(k_copy, B_ * C, [self.zz, self.pv])
        self._dot(self.r, self.zz, 0)
        for _ in range(int(n_it)):
            op(self.pv, self.Ap)
            self._dot(self.pv, self.Ap, 1)
            self.L(k_div, B_, [self.sc, 0, 1, 2])
            if mask is None:
                self.L(k_axpy, (B_, C), [C, self.sc, 2, f64(1.0), self.pv, x])
            else:
                self.L(k_axpym, (B_, C), [C, self.sc, 2, f64(1.0), self.pv, mask, x])
            self.L(k_axpy, (B_, C), [C, self.sc, 2, f64(-1.0), self.Ap, self.r])
            self.L(k_papply, B_ * C, [precond, self.r, self.zz])
            self._dot(self.r, self.zz, 3)
            self.L(k_div, B_, [self.sc, 3, 0, 4])
            self.L(k_zpbp, (B_, C), [C, self.sc, 4, self.zz, self.pv])
            self.L(k_mov, B_, [self.sc, 3, 0])

    def _admm_iter(self, cg_iters):
        (B_, C) = (self.B, self.C)
        self.L(k_rhs, (B_, C), [C, self.bb, self.mu, self.s, self.gam, self.zv, self.rho, self.tmp])
        self._pcg(self.tmp, self.x, cg_iters, lambda a, b: self._matvec(a, b, 1), self.Pj, False, mask=self.mask)
        self.L(k_projz, (B_, C), [C, self.x, self.gam, self.rho, self.mu, self.mask, self.zv])
        self.L(k_dual, (B_, C), [C, self.x, self.zv, self.rho, self.mask, self.gam])

    def _refresh_and_check(self, tol, desax=True, check=True):
        (B_, C) = (self.B, self.C)
        self._matvec(self.zv, self.tmp2, 0)
        self.L(k_uadd, B_ * C, [self.tmp2, self.bb, self.u])
        self.L(k_slipn, B_ * C, [self.u, self.s])
        if check:
            self.L(k_natres, (B_, C), [C, self.zv, self.u, self.mu, self.rhoc, 1 if desax else 0, self.res])
            self.L(k_resmax, B_, [C, self.res, f64(tol), self.mask, self.resenv, self.nact])

    def _fwd_block(self, n_iter, cg_iters, refresh_every, tol, desax=True):
        for i in range(int(n_iter)):
            self._admm_iter(cg_iters)
            if refresh_every and (i + 1) % refresh_every == 0:
                self._refresh_and_check(tol, desax, check=True)

    def reset_state(self, reset_est=False):
        (B_, C) = (self.B, self.C)
        self.L(k_zerov, B_ * C, [self.gam])
        self.L(k_zerov, B_ * C, [self.zv])
        self.L(k_zerov, B_ * C, [self.x])
        self.s.zero_()
        self.L(k_onei, B_, [self.mask])
        self.nact.zero_()
        if reset_est and hasattr(self, 'k_cert'):
            self.init_estimator(k_cert=self.k_cert)

    def forward(self, n_iter=100, cg_iters=8, refresh_every=4, tol=1e-08, chunk=25, desax=True, capture=True, readback=True):
        (B_, C) = (self.B, self.C)
        self.L(k_pjac, (B_, C), [C, self.Gd, self.eta, self.rho, self.Pj])
        self.reset_state()
        graph = None
        cw = min(int(chunk), int(n_iter))
        if capture and self.d.startswith('cuda'):
            wp.load_module(device=self.d)
            with wp.ScopedCapture(device=self.d) as cap:
                self._fwd_block(cw, cg_iters, refresh_every, tol, desax)
            graph = cap.graph
        wp.synchronize_device(self.d)
        t0 = time.perf_counter()
        done = 0
        while done < n_iter:
            if graph is not None:
                wp.capture_launch(graph)
            else:
                self._fwd_block(cw, cg_iters, refresh_every, tol, desax)
            done += cw
        wp.synchronize_device(self.d)
        ms = (time.perf_counter() - t0) * 1000.0
        self.L(k_copy, B_ * C, [self.zv, self.lam])
        if not readback:
            return Fwd(None, None, None, done, -1, ms, None)
        return Fwd(lam=self.lam.numpy().reshape(B_, C, 3).copy(), u=self.u.numpy().reshape(B_, C, 3).copy(), resenv=self.resenv.numpy().copy(), iters=done, n_active=int(self.mask.numpy().sum()), ms=ms, rho=self.rho.numpy().copy())

    def build_active(self, tol_lam=1e-10, tol_rel=1e-06):
        (B_, C) = (self.B, self.C)
        self.L(k_lamscale, B_, [C, self.zv, self.lsc])
        self.L(k_active, (B_, C), [C, self.zv, self.u, self.mu, self.lsc, f64(tol_lam), f64(tol_rel), self.Zm, self.Dv, self.slipf, self.what, self.lab])
        self.L(k_slipidx, B_, [C, self.slipf, self.sidx, self.nslip])
        self.L(k_sjac, B_ * C, [self.Zm, self.Dv, self.Gd, self.eta, self.Ps])
        return self

    def _grad_setup(self):
        """ZtHp = Z^T (H p_inh + db), independent of sigma -> computed once per direction."""
        (B_, C) = (self.B, self.C)
        self._matvec(self.pinh, self.tmp, 0)
        self.L(k_uadd, B_ * C, [self.tmp, self.dbv, self.tmp2])
        self.L(k_zth, B_ * C, [self.Zm, self.tmp2, self.ZtHp])

    def _grad_once(self, cg_iters):
        """sig (in) -> dlam, sigo (out).  One application of the AFFINE sigma map."""
        (B_, C) = (self.B, self.C)
        self.L(k_grhs, B_ * C, [C, self.ZtHp, self.mu, self.sig, self.slipf, self.tmp])
        self._pcg(self.tmp, self.yv, cg_iters, self._sopv, self.Ps, True)
        self.L(k_zmul, B_ * C, [self.Zm, self.yv, self.tmp])
        self.L(k_addp, B_ * C, [self.tmp, self.pinh, self.dlam])
        self._matvec(self.dlam, self.tmp2, 0)
        self.L(k_sigma, B_ * C, [self.slipf, self.what, self.tmp2, self.dbv, self.sigo])

    def grad_affine(self, n_calls, cg_iters):
        """n_calls = n_slip_max + 2 (graph-shape parameter; see module docstring)."""
        (B_, C) = (self.B, self.C)
        K = int(n_calls) - 2
        self._grad_setup()
        self.L(k_setsig, B_, [C, -1, self.sidx, self.nslip, self.sig])
        self._grad_once(cg_iters)
        self.L(k_copyf, B_ * C, [self.sigo, self.siga])
        for j in range(K):
            self.L(k_setsig, B_, [C, j, self.sidx, self.nslip, self.sig])
            self._grad_once(cg_iters)
            self.L(k_colB, B_, [C, j, self.sigo, self.siga, self.sidx, self.nslip, self.Bm])
        self.L(k_solvesig, B_, [C, self.Bm, self.siga, self.sidx, self.nslip, self.Aw, self.rw, self.sig])
        self._grad_once(cg_iters)

    def _fwd_part(self, n_iter, cg_iters, refresh_every, tol):
        (B_, C, NV) = (self.B, self.C, self.NV)
        self._fwd_block(n_iter, cg_iters, refresh_every, tol)
        self.L(k_copy, B_ * C, [self.zv, self.lam])
        self._matvec(self.lam, self.tmp, 0)
        self.L(k_vp, B_ * NV, [NV, self.vfree, self.wv, self.vpv])

    def _grad_part(self, n_calls, grad_cg, do_mu):
        (B_, C, NV) = (self.B, self.C, self.NV)
        self.build_active()
        self.L(k_zerov, B_ * C, [self.pinh])
        self.L(k_dbtau, (B_, C), [C, NV, self.M, self.MJ, f64(self.dt), self.jdir, self.dbv])
        self.grad_affine(n_calls, grad_cg)
        self._matvec(self.dlam, self.tmp, 0)
        self.L(k_dvp, (B_, NV), [NV, self.Mi, f64(self.dt), self.jdir, self.wv, self.dvpv])
        self.L(k_copy, B_ * C, [self.dlam, self.dltau])
        self.L(k_zerov, B_ * C, [self.dlmu])
        if do_mu:
            self.L(k_pmu, (B_, C), [C, self.lam, self.what, self.sidx, self.nslip, self.pinh])
            self.L(k_zerov, B_ * C, [self.dbv])
            self.grad_affine(n_calls, grad_cg)
            self.L(k_copy, B_ * C, [self.dlam, self.dlmu])

    def _step_and_grad(self, n_iter, cg_iters, refresh_every, tol, n_calls, grad_cg, do_mu):
        self._fwd_part(n_iter, cg_iters, refresh_every, tol)
        self._grad_part(n_calls, grad_cg, do_mu)

    def init_estimator(self, mu0=0.8, kappa0=0.35, nu0=4.0, sigma0=0.18, k_cert=1.0):
        B_ = self.B
        beta0 = 0.5 * (nu0 - 2.0) * kappa0 * sigma0 ** 2
        mueff0 = float(np.clip(mu0 - k_cert * sigma0, 0.05, 1.2))
        up = lambda a, v, dt: wp.copy(a, wp.array(np.full(B_, v, dtype=np.float64), dtype=dt, device=self.d))
        up(self.est_mu, mu0, f64)
        up(self.est_kappa, kappa0, f64)
        up(self.est_nu, nu0, f64)
        up(self.est_beta, beta0, f64)
        up(self.est_sigma, sigma0, f64)
        up(self.est_mueff, mueff0, f64)
        wp.copy(self.est_nobs, wp.array(np.zeros(B_, dtype=np.int32), dtype=wp.int32, device=self.d))
        self.k_cert = float(k_cert)

    def upload_mu_sann(self, mu_sann_arr):
        mu_sann_arr = np.ascontiguousarray(mu_sann_arr, np.float64).reshape(-1)
        if len(mu_sann_arr) == self.B:
            mu_sann_arr = np.repeat(mu_sann_arr, self.C)
        wp.copy(self.mu_sann, wp.array(mu_sann_arr, dtype=f64, device=self.d))

    def _est_part(self, k_cert=1.0):
        self.L(k_est_zero, self.B, [self.acc_cnt, self.acc_sum, self.acc_sumsq])
        self.L(k_est_detect_slip, (self.B, self.C), [self.C, self.zv, self.u, self.mu_sann, self.acc_cnt, self.acc_sum, self.acc_sumsq])
        self.L(k_est_update, self.B, [f64(k_cert), self.acc_cnt, self.acc_sum, self.acc_sumsq, self.est_mu, self.est_kappa, self.est_nu, self.est_beta, self.est_sigma, self.est_mueff, self.est_nobs])

    def _fwd_part_with_est(self, n_iter, cg_iters, refresh_every, tol, k_cert=1.0):
        self._fwd_part(n_iter, cg_iters, refresh_every, tol)
        self._est_part(k_cert)

    def _step_and_grad_with_est(self, n_iter, cg_iters, refresh_every, tol, n_calls, grad_cg, do_mu, k_cert=1.0):
        self._fwd_part_with_est(n_iter, cg_iters, refresh_every, tol, k_cert)
        self._grad_part(n_calls, grad_cg, do_mu)

    def capture_step_grad_with_est(self, n_iter, cg_iters=8, refresh_every=4, tol=1e-08, n_calls=3, grad_cg=60, do_mu=True, k_cert=1.0):
        self.L(k_pjac, (self.B, self.C), [self.C, self.Gd, self.eta, self.rho, self.Pj])
        self.reset_state()
        return self._capture(self._step_and_grad_with_est, n_iter, cg_iters, refresh_every, tol, n_calls, grad_cg, do_mu, k_cert)

    def capture_fwd_with_est(self, n_iter, cg_iters=8, refresh_every=4, tol=1e-08, k_cert=1.0):
        self.L(k_pjac, (self.B, self.C), [self.C, self.Gd, self.eta, self.rho, self.Pj])
        self.reset_state()
        return self._capture(self._fwd_part_with_est, n_iter, cg_iters, refresh_every, tol, k_cert)

    def count_launches(self, fn, *a, **kw):
        self._lc = [0]
        fn(*a, **kw)
        n = self._lc[0]
        self._lc = None
        return n

    def _capture(self, fn, *a):
        wp.load_module(device=self.d)
        wp.synchronize_device(self.d)
        t0 = time.perf_counter()
        with wp.ScopedCapture(device=self.d) as cap:
            fn(*a)
        wp.synchronize_device(self.d)
        return (cap.graph, (time.perf_counter() - t0) * 1000.0)

    def capture_step_grad(self, n_iter, cg_iters=8, refresh_every=4, tol=1e-08, n_calls=3, grad_cg=60, do_mu=True):
        self.L(k_pjac, (self.B, self.C), [self.C, self.Gd, self.eta, self.rho, self.Pj])
        self.reset_state()
        return self._capture(self._step_and_grad, n_iter, cg_iters, refresh_every, tol, n_calls, grad_cg, do_mu)

    def capture_fwd(self, n_iter, cg_iters=8, refresh_every=4, tol=1e-08):
        self.L(k_pjac, (self.B, self.C), [self.C, self.Gd, self.eta, self.rho, self.Pj])
        self.reset_state()
        return self._capture(self._fwd_part, n_iter, cg_iters, refresh_every, tol)

    def capture_grad(self, n_calls=3, grad_cg=60, do_mu=True):
        return self._capture(self._grad_part, n_calls, grad_cg, do_mu)

    def replay(self, graph, reps=1, reset=True, reset_est=False):
        if reset:
            self.L(k_pjac, (self.B, self.C), [self.C, self.Gd, self.eta, self.rho, self.Pj])
            self.reset_state(reset_est=reset_est)
        wp.synchronize_device(self.d)
        t0 = time.perf_counter()
        for _ in range(reps):
            wp.capture_launch(graph)
        wp.synchronize_device(self.d)
        return (time.perf_counter() - t0) * 1000.0 / reps

    def set_jdir(self, j):
        wp.copy(self.jdir, wp.array(np.array([int(j)], np.int32), dtype=int, device=self.d))
