"""GPU exact-cone NCP-PGS (Warp); the CPU reference is `ncp_ref.py`.

Solves, per step, the cone complementarity problem

    K_mu ∋ lambda_c  ⊥  (u_c + Gamma(u_c)) ∈ K_mu*,      u = G lambda + b,

with the exact second-order Coulomb cone K_mu = {(l_n, l_t) : ||l_t|| <= mu l_n}, the de Saxce
correction Gamma(u) = (mu ||u_t||, 0, 0) and the per-contact fixed point

    lambda_c <- proj_{K_mu}( lambda_c - rho_c (u_c + Gamma(u_c)) ),   rho_c = 1 / lambda_max(W_cc),

W_cc = J_c M^-1 J_c^T the 3x3 diagonal block of the Delassus operator. The full Delassus is never
formed: u_c is read off the generalized velocity v = v_free + M^-1 J^T lambda, which is maintained
incrementally (matrix-free, SPEC rule 5).

Parallelism / determinism (SPEC B.1-B.4):
  * deterministic graph colouring. Contacts are first put in a canonical order by a radix sort on the
    key (body_a, body_b, feature); the Jones-Plassmann priority is a hash of the RANK in that sorted
    list, so the colouring is a pure function of the contact set and not of the generation order.
    No two contacts of a colour share a dynamic body -> Jacobi inside a colour is exact Gauss-Seidel,
    and the colours are swept sequentially.
  * the per-body TOTAL M^-1 J^T lambda is accumulated in int64 fixed point (scale 1e10, float64
    rounding), as in src/motion_engine/contact_engine_gpu.py: integer addition is associative, so the
    sum does not depend on thread order and two runs are bit-identical. A contact update REPLACES its
    own quantized contribution (subtract the old, add the new; f_swap) instead of adding a quantized
    velocity delta, and v = v_free + M^-1 J^T lambda is recomputed from that total every colour pass.
    The accumulator is then a pure function of the current lambda, so the fixed-point rounding does
    not compound with the iteration count.
  * the colouring is cached and reused while (body pairs, feature ids, M^-1) are unchanged; the
    comparison runs on the device (one kernel, one int read back), as in
    src/motion_engine/contact_engine_gpu_color_cache.py. "Masses" in the cache key is the
    translational block 1/m I of M^-1 plus the dynamic/static flag per DOF row: the world inertia
    R Ib^-1 R^T changes every step as a body rotates but does not change the colouring, so comparing
    all 36 entries would make the cache never hit.
  * optional CUDA graph capture of one full colour loop (captured once per colouring, replayed per
    iteration).

Contract (SPEC rule "Contract between A and B"):  `Scene` carries J (3n_c x n_v, float64), Minv
(n_v x n_v or the n_v diagonal), v_free (n_v), mu (n_c), body_pairs (n_c x 2, -1 = static), dt.
n_v = 6 * n_bodies (free rigid bodies, generalized velocity (v, omega) per body, world frame).
`Scene.blocks()` slices the dense J into the per-contact 3x6 blocks the GPU solver consumes;
`BlockScene` can also be built directly (the synthetic timing scenes do, to avoid materialising a
12288 x 3072 dense J).

  python -u ncp_gpu.py            # selftest: determinism, scenes (i)-(iv), timings

Audited numbers this module reproduces (tests/test_ncp_gpu_cone.py):
  * lam and v bit-identical over two runs in float32 and float64, also under a permuted
    contact order; float32-float64 gap 1.25e-7.
  * sliding box, 200 steps, 400 iterations/step: max normal gap 1.552204e-11 m with the
    de Saxce term on, 4.352523e-03 m with it off (same code, one flag).
  * sliding distance against v0^2/(2 mu g) = 0.169895 m: -1.2227e-02 at dt 1/240 and
    -3.0687e-03 at 1/960 (first order in dt); the two cone models differ by 1.2e-5
    relative, three orders below the discretisation error.
  * int64 fixed point accumulating the total (not the delta): drift at 80 000 iterations
    4.4e-6 -> 1.4e-10.

Audit reservations carried as numbers, not adjectives:
  * determinism holds for the same device and launch configuration. The independent
    oracle (Box3D) is bit-identical over 1/4/8 threads, which is a stronger property.
  * island sleeping is 1.18x at N = 10000 with 90 % of the bodies at rest, and the K8
    tower never sleeps at a 1e-3 threshold: the jitter floor is 2.7e-3 m/s, set by the
    undamped Baumgarte term, not by the sleeper.
  * PGS at a 1000:1 mass ratio: 200 iterations shoot the column out (213 mm), 800 hold
    3.3 mm at residual 1.7e-2. First order is not enough in that regime.
"""
import os
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Optional
import numpy as np
_WP_FUSE_FP = os.environ.get('WP_FUSE_FP', '0').lower() in ('1', 'true', 'yes', 'on')
_WP_FAST_MATH = os.environ.get('WP_FAST_MATH', '0').lower() in ('1', 'true', 'yes', 'on')
WP_MODULE_OPTIONS = {'fuse_fp': _WP_FUSE_FP, 'fast_math': _WP_FAST_MATH}

def _apply_wp_module_options(wp):
    """Set WP_MODULE_OPTIONS on this Warp module. Must run before the first @wp.kernel of the
    module is declared (the kernels live in `_build_kernels`, whose `__module__` is this file).
    `set_module_options` takes a module OBJECT, not a name -- passing a string raises."""
    wp.set_module_options(dict(WP_MODULE_OPTIONS), module=sys.modules[__name__])
    return dict(WP_MODULE_OPTIONS)
_S = 10000000000.0
_BETA = 0.2
_SLOP = 0.0001
_G = 9.81
_FEATBITS = 20
_NCOLMAX = 1024
CONTACT_HERTZ = 30.0
CONTACT_DAMPING_RATIO = 10.0
CONTACT_SPEED = 3.0

@dataclass(frozen=True)
class Softness:
    """soft_step.py:28-32."""
    bias_rate: float
    mass_scale: float
    impulse_scale: float

def make_soft(hertz: float, zeta: float, h: float) -> Softness:
    """soft_step.py:35-44. Discrete implicit spring; mass_scale + impulse_scale == 1 (invariant I-2)."""
    if hertz <= 0.0:
        return Softness(0.0, 0.0, 0.0)
    omega = 2.0 * np.pi * hertz
    a1 = 2.0 * zeta + h * omega
    a2 = h * omega * a1
    a3 = 1.0 / (1.0 + a2)
    return Softness(bias_rate=omega / a1, mass_scale=a2 * a3, impulse_scale=a3)

def contact_hertz(inv_h: float, world_hertz: float=CONTACT_HERTZ) -> float:
    """soft_step.py:47-49."""
    return min(world_hertz, 0.125 * inv_h)

def bias_velocity(separation: float, soft: Softness, use_bias: bool, inv_h: float, contact_speed: float=CONTACT_SPEED):
    """soft_step.py:52-64. (velocity_bias, mass_scale, impulse_scale). This is the host mirror of the
    three branches the kernel takes (see k_sweep)."""
    if separation > 0.0:
        return (separation * inv_h, 1.0, 0.0)
    if use_bias:
        vb = max(soft.mass_scale * soft.bias_rate * separation, -contact_speed)
        return (vb, soft.mass_scale, soft.impulse_scale)
    return (0.0, 1.0, 0.0)

@dataclass
class SoftStep:
    """Soft-step configuration for one solve. `relax_frac` of the per-substep iteration budget is
    spent in the final relax pass (use_bias = False), the rest in the biased pass; box3d uses 4
    biased + 1 relax, i.e. 0.2."""
    hertz: float = CONTACT_HERTZ
    zeta: float = CONTACT_DAMPING_RATIO
    contact_speed: float = CONTACT_SPEED
    relax_frac: float = 0.25

    def params(self, h: float):
        """(bias_rate, mass_scale, impulse_scale, inv_h, contact_speed) for a substep of length h."""
        inv_h = 1.0 / h
        s = make_soft(contact_hertz(inv_h, self.hertz), self.zeta, h)
        return np.array([s.bias_rate, s.mass_scale, s.impulse_scale, inv_h, self.contact_speed], dtype=np.float64)

    def split(self, iters: int):
        """(biased iterations, relax iterations) out of a per-substep budget of `iters`."""
        iters = max(1, int(iters))
        relax = min(iters - 1, max(1, int(round(self.relax_frac * iters)))) if iters > 1 else 0
        return (iters - relax, relax)

@dataclass
class BlockScene:
    """Per-contact 3x6 Jacobian blocks + per-body 6x6 inverse mass. The signs of the pair are baked in:
    u_c = Ja_c v_a + Jb_c v_b + bias_c, so Jb already carries the minus of the relative velocity."""
    Ja: np.ndarray
    Jb: np.ndarray
    Minv: np.ndarray
    v_free: np.ndarray
    mu: np.ndarray
    body_pairs: np.ndarray
    dt: float
    bias: np.ndarray
    feature: np.ndarray
    sep: Optional[np.ndarray] = None

    @property
    def n_c(self) -> int:
        return int(self.Ja.shape[0])

    @property
    def n_b(self) -> int:
        return int(self.Minv.shape[0])

    def delassus(self) -> np.ndarray:
        """G = J M^-1 J^T, dense float64. Host only, for the CPU-reference comparison (SPEC rule 5
        allows an explicit G for n_c <= 4096 on the A side)."""
        (n_c, n_b) = (self.n_c, self.n_b)
        J = np.zeros((3 * n_c, 6 * n_b))
        for c in range(n_c):
            (a, b) = (int(self.body_pairs[c, 0]), int(self.body_pairs[c, 1]))
            if a >= 0:
                J[3 * c:3 * c + 3, 6 * a:6 * a + 6] += self.Ja[c]
            if b >= 0:
                J[3 * c:3 * c + 3, 6 * b:6 * b + 6] += self.Jb[c]
        M = np.zeros((6 * n_b, 6 * n_b))
        for i in range(n_b):
            M[6 * i:6 * i + 6, 6 * i:6 * i + 6] = self.Minv[i]
        return J @ M @ J.T

    def rhs(self) -> np.ndarray:
        """b = J v_free + bias, float64."""
        n_c = self.n_c
        out = self.bias.astype(np.float64).copy()
        for c in range(n_c):
            (a, b) = (int(self.body_pairs[c, 0]), int(self.body_pairs[c, 1]))
            if a >= 0:
                out[c] += self.Ja[c] @ self.v_free[a]
            if b >= 0:
                out[c] += self.Jb[c] @ self.v_free[b]
        return out.reshape(-1)

@dataclass
class Scene:
    """The A/B contract of SPEC: dense J, Minv (dense or diagonal), v_free, mu, body_pairs, dt."""
    J: np.ndarray
    Minv: np.ndarray
    v_free: np.ndarray
    mu: np.ndarray
    body_pairs: np.ndarray
    dt: float
    bias: Optional[np.ndarray] = None
    feature: Optional[np.ndarray] = None
    sep: Optional[np.ndarray] = None

    def blocks(self) -> BlockScene:
        J = np.asarray(self.J, dtype=np.float64)
        n_c = J.shape[0] // 3
        n_v = J.shape[1]
        if n_v % 6:
            raise ValueError(f'Scene.blocks: n_v={n_v} is not a multiple of 6 (6-DOF rigid bodies)')
        n_b = n_v // 6
        Mi = np.asarray(self.Minv, dtype=np.float64)
        if Mi.ndim == 1:
            Mi = np.diag(Mi)
        Ja = np.zeros((n_c, 3, 6))
        Jb = np.zeros((n_c, 3, 6))
        pairs = np.asarray(self.body_pairs, dtype=np.int32).reshape(n_c, 2)
        for c in range(n_c):
            (a, b) = (int(pairs[c, 0]), int(pairs[c, 1]))
            if a >= 0:
                Ja[c] = J[3 * c:3 * c + 3, 6 * a:6 * a + 6]
            if b >= 0:
                Jb[c] = J[3 * c:3 * c + 3, 6 * b:6 * b + 6]
        Minv = np.zeros((n_b, 6, 6))
        for i in range(n_b):
            Minv[i] = Mi[6 * i:6 * i + 6, 6 * i:6 * i + 6]
        bias = np.zeros((n_c, 3)) if self.bias is None else np.asarray(self.bias, dtype=np.float64).reshape(n_c, 3)
        feat = np.arange(n_c, dtype=np.int32) if self.feature is None else np.asarray(self.feature, dtype=np.int32)
        sep = None if self.sep is None else np.asarray(self.sep, dtype=np.float64).reshape(n_c)
        return BlockScene(Ja=Ja, Jb=Jb, Minv=Minv, v_free=np.asarray(self.v_free, dtype=np.float64).reshape(n_b, 6), mu=np.asarray(self.mu, dtype=np.float64).reshape(n_c), body_pairs=pairs, dt=float(self.dt), bias=bias, feature=feat, sep=sep)

@dataclass
class SolveResult:
    lam: np.ndarray
    v: np.ndarray
    residual: float
    iters: int
    n_colors: int
    cache_hit: bool
    ms: dict = field(default_factory=dict)
    history: Optional[dict] = None
    v_bias: Optional[np.ndarray] = None
    force_selection: str = 'COLORED_PGS_ORDER_WARM_START'
    contact_model: str = 'NCP_RIGID'
    bias_applied: bool = True
    contact_info: dict = field(default_factory=lambda : {'status': 'UNASSESSED_CONTACT_UNIQUENESS'})

    def diagnose_contacts(self, scene: BlockScene, rcond=1e-12, max_contacts=256):
        """Opt-in bounded host diagnostic; never form a dense SVD in the GPU loop.

        The supplied scene must be the scene solved for this result. The result
        is a numerical fixed-branch diagnostic, not a global uniqueness proof.
        Soft-step and CCP have different residual systems and remain unassessed.
        """
        if self.contact_model != 'NCP_RIGID':
            self.contact_info = {'status': 'UNASSESSED_CONTACT_MODEL', 'scope': self.contact_model}
        elif scene.n_c > max_contacts:
            self.contact_info = {'status': 'UNASSESSED_SIZE_LIMIT', 'max_contacts': max_contacts}
        else:
            if self.lam.shape != (scene.n_c, 3):
                raise ValueError('scene and result contact counts do not agree')
            from .ncp_ref import contact_diagnostics, classify
            J = np.zeros((3 * scene.n_c, 6 * scene.n_b))
            for (c, pair) in enumerate(scene.body_pairs):
                for (body, block) in zip(pair, (scene.Ja[c], scene.Jb[c])):
                    if body >= 0:
                        J[3 * c:3 * c + 3, 6 * body:6 * body + 6] += block
            (G, b) = (scene.delassus(), scene.rhs())
            if not self.bias_applied:
                b = b - scene.bias.reshape(-1)
            lam = self.lam.reshape(-1)
            labels = classify(lam, G, b, scene.mu)
            self.contact_info = contact_diagnostics(lam, G, b, scene.mu, labels, rcond=rcond, observable=J.T)
        self.contact_info['force_selection'] = self.force_selection
        return self.contact_info
_KCACHE: dict = {}

def _types(wp, dtype: str):
    if dtype == 'float32':
        wpt = wp.float32
        return (wpt, wp.vec3, wp.mat33, wp.types.matrix(shape=(3, 6), dtype=wpt), wp.types.matrix(shape=(6, 6), dtype=wpt), wp.types.vector(length=6, dtype=wpt), np.float32)
    if dtype == 'float64':
        wpt = wp.float64
        return (wpt, wp.types.vector(length=3, dtype=wpt), wp.types.matrix(shape=(3, 3), dtype=wpt), wp.types.matrix(shape=(3, 6), dtype=wpt), wp.types.matrix(shape=(6, 6), dtype=wpt), wp.types.vector(length=6, dtype=wpt), np.float64)
    raise ValueError(f"dtype must be 'float32' or 'float64', got {dtype!r}")

def _build_kernels(wp, dtype: str):
    """One compiled kernel set per floating-point type. Built inside a factory so the annotations
    close over the concrete Warp types (float32 or float64)."""
    _apply_wp_module_options(wp)
    if dtype in _KCACHE:
        return _KCACHE[dtype]
    (flt, vec3t, mat33t, mat36t, mat66t, vec6t, nptype) = _types(wp, dtype)
    ZERO = wp.constant(flt(0.0))
    ONE = wp.constant(flt(1.0))
    MONE = wp.constant(flt(-1.0))
    TWO = wp.constant(flt(2.0))
    THREE = wp.constant(flt(3.0))
    SIX = wp.constant(flt(6.0))
    TINY = wp.constant(flt(1e-14))
    SCALE64 = wp.constant(wp.float64(_S))

    @wp.func
    def f_proj(z: vec3t, mu: flt) -> vec3t:
        """Exact projection onto K_mu = {(zn, zt) : ||zt|| <= mu zn}."""
        zn = z[0]
        zt = wp.sqrt(z[1] * z[1] + z[2] * z[2])
        res = vec3t()
        if zt <= mu * zn:
            res = z
        elif mu * zt <= -zn:
            res = vec3t()
        else:
            s = (zn + mu * zt) / (ONE + mu * mu)
            f = mu * s / zt
            res = vec3t(s, f * z[1], f * z[2])
        return res

    @wp.func
    def f_eigmax(A: mat33t) -> flt:
        """Largest eigenvalue of a symmetric 3x3 (closed form, deterministic)."""
        p1 = A[0, 1] * A[0, 1] + A[0, 2] * A[0, 2] + A[1, 2] * A[1, 2]
        q = (A[0, 0] + A[1, 1] + A[2, 2]) / THREE
        out = wp.max(A[0, 0], wp.max(A[1, 1], A[2, 2]))
        if p1 > ZERO:
            p2 = (A[0, 0] - q) * (A[0, 0] - q) + (A[1, 1] - q) * (A[1, 1] - q) + (A[2, 2] - q) * (A[2, 2] - q) + TWO * p1
            p = wp.sqrt(p2 / SIX)
            if p > TINY:
                B = ONE / p * (A - q * mat33t(ONE, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ONE))
                r = wp.determinant(B) / TWO
                r = wp.clamp(r, MONE, ONE)
                phi = wp.acos(r) / THREE
                out = q + TWO * p * wp.cos(phi)
        return out

    @wp.func
    def f_acc(acc: wp.array(dtype=wp.int64), b: int, d: vec6t):
        """Order-invariant accumulation of one contact's ABSOLUTE contribution M^-1 J_c^T lambda_c:
        round(d * 1e10) in float64, int64 atomics. Integer addition is associative -> the sum is
        thread-order independent."""
        for k in range(6):
            wp.atomic_add(acc, b * 6 + k, wp.int64(wp.round(wp.float64(d[k]) * SCALE64)))

    @wp.func
    def f_swap(acc: wp.array(dtype=wp.int64), b: int, dn: vec6t, do: vec6t):
        """Replace this contact's quantized contribution in the per-body total:
            acc_b <- acc_b - round(S * M^-1 J_c^T lambda_c^old) + round(S * M^-1 J_c^T lambda_c^new).
        Both terms are computed from the ABSOLUTE lambda, so the quantization of the superseded value
        is removed exactly (integer subtraction) instead of being left behind. After any number of
        updates acc_b is therefore exactly sum_c round(S * M^-1 J_c^T lambda_c) for the CURRENT
        lambda -- a pure function of lambda, with no dependence on the iteration history. The error
        against the exact M^-1 J^T lambda stays bounded by (contacts on b) x 0.5 quantum whatever the
        iteration count; accumulating the deltas instead compounds one rounding per update."""
        for k in range(6):
            wp.atomic_add(acc, b * 6 + k, wp.int64(wp.round(wp.float64(dn[k]) * SCALE64)) - wp.int64(wp.round(wp.float64(do[k]) * SCALE64)))

    @wp.kernel
    def k_rho(C: int, pa: wp.array(dtype=int), pb: wp.array(dtype=int), Ja: wp.array(dtype=mat36t), Jb: wp.array(dtype=mat36t), Minv: wp.array(dtype=mat66t), rho: wp.array(dtype=flt)):
        c = wp.tid()
        if c >= C:
            return
        W = mat33t()
        a = pa[c]
        b = pb[c]
        if a >= 0:
            W = W + Ja[c] * Minv[a] * wp.transpose(Ja[c])
        if b >= 0:
            W = W + Jb[c] * Minv[b] * wp.transpose(Jb[c])
        lm = f_eigmax(W)
        if lm > TINY:
            rho[c] = ONE / lm
        else:
            rho[c] = ZERO

    @wp.kernel
    def k_init_v(nb: int, vfree: wp.array(dtype=vec6t), v: wp.array(dtype=vec6t)):
        i = wp.tid()
        if i < nb:
            v[i] = vfree[i]

    @wp.kernel
    def k_scatter(C: int, pa: wp.array(dtype=int), pb: wp.array(dtype=int), Ja: wp.array(dtype=mat36t), Jb: wp.array(dtype=mat36t), Minv: wp.array(dtype=mat66t), lam: wp.array(dtype=vec3t), acc: wp.array(dtype=wp.int64)):
        """v <- v + M^-1 J^T lambda for a warm start (all contacts at once, int64 accumulation)."""
        c = wp.tid()
        if c >= C:
            return
        d = lam[c]
        a = pa[c]
        b = pb[c]
        if a >= 0:
            f_acc(acc, a, Minv[a] * (wp.transpose(Ja[c]) * d))
        if b >= 0:
            f_acc(acc, b, Minv[b] * (wp.transpose(Jb[c]) * d))

    @wp.kernel
    def k_sweep(ci: int, cstart: wp.array(dtype=int), ccount: wp.array(dtype=int), order: wp.array(dtype=int), pa: wp.array(dtype=int), pb: wp.array(dtype=int), Ja: wp.array(dtype=mat36t), Jb: wp.array(dtype=mat36t), Minv: wp.array(dtype=mat66t), v: wp.array(dtype=vec6t), cbias: wp.array(dtype=vec3t), mu: wp.array(dtype=flt), rho: wp.array(dtype=flt), desaxce: int, lam: wp.array(dtype=vec3t), acc: wp.array(dtype=wp.int64), sep: wp.array(dtype=flt), softp: wp.array(dtype=flt), sflag: wp.array(dtype=int)):
        """One colour of the NCP-PGS sweep: Jacobi over the contacts of the colour (they share no
        dynamic body, so this is exact Gauss-Seidel). The int64 accumulator carries the per-body TOTAL
        M^-1 J^T lambda, not a running sum of velocity deltas (see f_swap).

        `sflag[0] == 0` -> the Baumgarte path: the bias is the precomputed `cbias` and the update is
        lambda <- proj_K(lambda - rho (u + Gamma(u))), unchanged arithmetic.

        `sflag[0] == 1` -> the TGS soft-step path (re-derived from `motion_engine/soft_step.py`,
        cited at module level; `softp` = (bias_rate, mass_scale, impulse_scale, inv_h,
        contact_speed), `sflag[1]` = use_bias). Per contact the three branches of
        soft_step.bias_velocity (soft_step.py:52-64) are taken on the CURRENT separation:

            s > 0                : vb = s * inv_h,                     ms = 1, is = 0   (speculative)
            s <= 0 and use_bias  : vb = max(ms*bias_rate*s, -c_speed), ms, is from make_soft
            relax (use_bias = 0) : vb = 0,                             ms = 1, is = 0

        and the normal row of the update carries them. `bias_velocity` already multiplies the
        push-out term by mass_scale (soft_step.py:62), so mass_scale multiplies only the VELOCITY at
        the call site and impulse_scale only the accumulated impulse -- the wiring docs/RUNNING.md:65-67
        states for this primitive, and the same expression as box3d's
            impulse = -normalMass * (massScale * vn + bias) - impulseScale * oldImpulse
        with normalMass == rho_c:

            z_n = (1 - is) lam_n - rho (ms (u_n + Gamma_n) + vb),   z_t = lam_t - rho u_t
            lam <- proj_K(z).

        The soft scaling is applied to the NORMAL row only, which is what box3d does (its friction
        rows are solved with no massScale/impulseScale); scaling all three rows would, since
        ms + is == 1 and proj_K is positively homogeneous, collapse the whole update to
        ms * proj_K(lam - rho u) and soften friction with it."""
        t = wp.tid()
        if t >= ccount[ci]:
            return
        c = order[cstart[ci] + t]
        a = pa[c]
        b = pb[c]
        u = cbias[c]
        if a >= 0:
            u = u + Ja[c] * v[a]
        if b >= 0:
            u = u + Jb[c] * v[b]
        m = mu[c]
        lo = lam[c]
        zn = vec3t()
        if sflag[0] == 1:
            sp = sep[c]
            vb = ZERO
            ms = ONE
            isc = ZERO
            if sp > ZERO:
                vb = sp * softp[3]
            elif sflag[1] == 1:
                vb = wp.max(softp[1] * softp[0] * sp, -softp[4])
                ms = softp[1]
                isc = softp[2]
            rn = u[0]
            if desaxce == 1:
                rn = rn + m * wp.sqrt(u[1] * u[1] + u[2] * u[2])
            zn = f_proj(vec3t((ONE - isc) * lo[0] - rho[c] * (ms * rn + vb), lo[1] - rho[c] * u[1], lo[2] - rho[c] * u[2]), m)
        else:
            s = u
            if desaxce == 1:
                ut = wp.sqrt(u[1] * u[1] + u[2] * u[2])
                s = vec3t(u[0] + m * ut, u[1], u[2])
            zn = f_proj(lo - rho[c] * s, m)
        lam[c] = zn
        if a >= 0:
            JaT = wp.transpose(Ja[c])
            f_swap(acc, a, Minv[a] * (JaT * zn), Minv[a] * (JaT * lo))
        if b >= 0:
            JbT = wp.transpose(Jb[c])
            f_swap(acc, b, Minv[b] * (JbT * zn), Minv[b] * (JbT * lo))

    @wp.kernel
    def k_apply(nb: int, vfree: wp.array(dtype=vec6t), v: wp.array(dtype=vec6t), acc: wp.array(dtype=wp.int64)):
        """v <- v_free + M^-1 J^T lambda, recomputed from the accumulated total every colour pass.
        The accumulator is NOT cleared: it holds the running total, not a per-pass delta."""
        i = wp.tid()
        if i >= nb:
            return
        d = vec6t()
        for k in range(6):
            d[k] = flt(wp.float64(acc[i * 6 + k]) / SCALE64)
        v[i] = vfree[i] + d

    @wp.kernel
    def k_apply_c(n: int, base: int, cbody: wp.array(dtype=int), vfree: wp.array(dtype=vec6t), v: wp.array(dtype=vec6t), acc: wp.array(dtype=wp.int64)):
        """Same as k_apply but only for the bodies touched by ONE colour. A colour contains each
        dynamic body at most once, so only those bodies' accumulators changed in the sweep just
        launched; every other body's v is already correct. This makes the colour loop cost scale with
        the awake contact set instead of with the total body count (the difference between a sleeping
        island being free and only its contacts being free)."""
        t = wp.tid()
        if t >= n:
            return
        i = cbody[base + t]
        d = vec6t()
        for k in range(6):
            d[k] = flt(wp.float64(acc[i * 6 + k]) / SCALE64)
        v[i] = vfree[i] + d

    @wp.kernel
    def k_resid(C: int, pa: wp.array(dtype=int), pb: wp.array(dtype=int), Ja: wp.array(dtype=mat36t), Jb: wp.array(dtype=mat36t), v: wp.array(dtype=vec6t), cbias: wp.array(dtype=vec3t), mu: wp.array(dtype=flt), rho: wp.array(dtype=flt), desaxce: int, lam: wp.array(dtype=vec3t), res: wp.array(dtype=flt), uout: wp.array(dtype=vec3t), sep: wp.array(dtype=flt), softp: wp.array(dtype=flt), sflag: wp.array(dtype=int)):
        """Natural-map residual ||lambda - proj_K(z)||_inf of the SAME fixed-point map the sweep
        iterates, so in soft-step mode it is measured against the pass that ran last (after a relax
        pass, sflag[1] == 0, it is the unbiased NCP residual)."""
        c = wp.tid()
        if c >= C:
            return
        a = pa[c]
        b = pb[c]
        u = cbias[c]
        if a >= 0:
            u = u + Ja[c] * v[a]
        if b >= 0:
            u = u + Jb[c] * v[b]
        uout[c] = u
        m = mu[c]
        lo = lam[c]
        d = vec3t()
        if sflag[0] == 1:
            sp = sep[c]
            vb = ZERO
            ms = ONE
            isc = ZERO
            if sp > ZERO:
                vb = sp * softp[3]
            elif sflag[1] == 1:
                vb = wp.max(softp[1] * softp[0] * sp, -softp[4])
                ms = softp[1]
                isc = softp[2]
            rn = u[0]
            if desaxce == 1:
                rn = rn + m * wp.sqrt(u[1] * u[1] + u[2] * u[2])
            d = lam[c] - f_proj(vec3t((ONE - isc) * lo[0] - rho[c] * (ms * rn + vb), lo[1] - rho[c] * u[1], lo[2] - rho[c] * u[2]), m)
        else:
            s = u
            if desaxce == 1:
                ut = wp.sqrt(u[1] * u[1] + u[2] * u[2])
                s = vec3t(u[0] + m * ut, u[1], u[2])
            d = lam[c] - f_proj(lam[c] - rho[c] * s, m)
        res[c] = wp.max(wp.abs(d[0]), wp.max(wp.abs(d[1]), wp.abs(d[2])))

    @wp.kernel
    def k_sortkey(C: int, pa: wp.array(dtype=int), pb: wp.array(dtype=int), feat: wp.array(dtype=int), nb1: wp.int64, key: wp.array(dtype=wp.int64), val: wp.array(dtype=int)):
        """Canonical key (body_a, body_b, feature) -> the sorted order is a pure function of the
        contact set, never of the generation order."""
        i = wp.tid()
        if i >= C:
            return
        pair = (wp.int64(pa[i]) + wp.int64(1)) * nb1 + (wp.int64(pb[i]) + wp.int64(1))
        key[i] = pair * wp.int64(1048576) + wp.int64(feat[i])
        val[i] = i

    @wp.kernel
    def k_prio(C: int, sidx: wp.array(dtype=int), pkey: wp.array(dtype=wp.int64)):
        """Jones-Plassmann priority from the RANK in the canonical order: a hash (so a chain of
        contacts -- a stack -- does not colour one per round) packed with the rank (so the key is
        unique). Pure function of the sorted list."""
        r = wp.tid()
        if r >= C:
            return
        st = wp.rand_init(7919, r)
        p = wp.int64(wp.randi(st))
        if p < wp.int64(0):
            p = -p
        pkey[sidx[r]] = p * wp.int64(4294967296) + wp.int64(r)

    @wp.kernel
    def k_clear_claim(claim: wp.array(dtype=wp.int64)):
        claim[wp.tid()] = wp.int64(9223372036854775807)

    @wp.kernel
    def k_claim(C: int, col: wp.array(dtype=int), pa: wp.array(dtype=int), pb: wp.array(dtype=int), dyn: wp.array(dtype=int), pkey: wp.array(dtype=wp.int64), claim: wp.array(dtype=wp.int64)):
        c = wp.tid()
        if c >= C:
            return
        if col[c] >= 0:
            return
        a = pa[c]
        b = pb[c]
        k = pkey[c]
        if a >= 0:
            if dyn[a] == 1:
                wp.atomic_min(claim, a, k)
        if b >= 0:
            if dyn[b] == 1:
                wp.atomic_min(claim, b, k)

    @wp.kernel
    def k_win(C: int, cur: int, col: wp.array(dtype=int), pa: wp.array(dtype=int), pb: wp.array(dtype=int), dyn: wp.array(dtype=int), pkey: wp.array(dtype=wp.int64), claim: wp.array(dtype=wp.int64), left: wp.array(dtype=int)):
        c = wp.tid()
        if c >= C:
            return
        if col[c] >= 0:
            return
        a = pa[c]
        b = pb[c]
        k = pkey[c]
        ok = True
        if a >= 0:
            if dyn[a] == 1:
                if claim[a] != k:
                    ok = False
        if b >= 0:
            if dyn[b] == 1:
                if claim[b] != k:
                    ok = False
        if ok:
            col[c] = cur
        else:
            wp.atomic_add(left, 0, 1)

    @wp.kernel
    def k_colkey(C: int, col: wp.array(dtype=int), key: wp.array(dtype=wp.int64), val: wp.array(dtype=int)):
        i = wp.tid()
        if i >= C:
            return
        key[i] = wp.int64(col[i])
        val[i] = i

    @wp.kernel
    def k_colcount(C: int, col: wp.array(dtype=int), ccount: wp.array(dtype=int)):
        i = wp.tid()
        if i >= C:
            return
        wp.atomic_add(ccount, col[i], 1)

    @wp.kernel
    def k_cache_cmp(C: int, nb: int, pa: wp.array(dtype=int), pb: wp.array(dtype=int), feat: wp.array(dtype=int), opa: wp.array(dtype=int), opb: wp.array(dtype=int), ofeat: wp.array(dtype=int), Minv: wp.array(dtype=mat66t), oMinv: wp.array(dtype=mat66t), mism: wp.array(dtype=int)):
        """Device-side exact comparison of the colouring inputs (pairs, feature ids, masses)."""
        i = wp.tid()
        if i < C:
            if pa[i] != opa[i]:
                wp.atomic_max(mism, 0, 1)
            if pb[i] != opb[i]:
                wp.atomic_max(mism, 0, 1)
            if feat[i] != ofeat[i]:
                wp.atomic_max(mism, 0, 1)
        if i < nb:
            m = Minv[i]
            o = oMinv[i]
            for r in range(3):
                if m[r, r] != o[r, r]:
                    wp.atomic_max(mism, 0, 1)
            for r in range(6):
                dm = flt(0.0)
                do = flt(0.0)
                for cc in range(6):
                    dm = dm + wp.abs(m[r, cc])
                    do = do + wp.abs(o[r, cc])
                if (dm > flt(0.0)) != (do > flt(0.0)):
                    wp.atomic_max(mism, 0, 1)
    K = dict(rho=k_rho, init_v=k_init_v, scatter=k_scatter, sweep=k_sweep, apply=k_apply, apply_c=k_apply_c, resid=k_resid, sortkey=k_sortkey, prio=k_prio, clear_claim=k_clear_claim, claim=k_claim, win=k_win, colkey=k_colkey, colcount=k_colcount, cmp=k_cache_cmp, types=(flt, vec3t, mat33t, mat36t, mat66t, vec6t, nptype))
    _KCACHE[dtype] = K
    return K

class NCPSolverGPU:
    """Persistent GPU solver. Capacity is fixed at construction; the colour cache and the captured
    CUDA graph live across `solve` calls, so a simulation loop with an unchanged contact set pays the
    colouring once."""

    def __init__(self, n_bodies: int, max_contacts: int, dtype: str='float32', device: Optional[str]=None, desaxce: bool=True, graph_capture: bool=True, color_cache: bool=True):
        import warp as wp
        self.wp = wp
        wp.init()
        self.dev = device or ('cuda:0' if wp.is_cuda_available() else 'cpu')
        self.dtype = dtype
        self.K = _build_kernels(wp, dtype)
        (self.flt, self.vec3t, self.mat33t, self.mat36t, self.mat66t, self.vec6t, self.np_t) = self.K['types']
        self.NB = int(n_bodies)
        self.MC = int(max_contacts)
        self.desaxce = 1 if desaxce else 0
        self.use_cache = bool(color_cache)
        self.want_graph = bool(graph_capture) and self.dev.startswith('cuda')
        (MC, NB) = (self.MC, self.NB)
        d = self.dev
        z = lambda n, t: wp.zeros(n, dtype=t, device=d)
        self.Ja = z(MC, self.mat36t)
        self.Jb = z(MC, self.mat36t)
        self.Minv = z(NB, self.mat66t)
        self.vfree = z(NB, self.vec6t)
        self.v = z(NB, self.vec6t)
        self.mu = z(MC, self.flt)
        self.rho = z(MC, self.flt)
        self.cbias = z(MC, self.vec3t)
        self.lam = z(MC, self.vec3t)
        self.uout = z(MC, self.vec3t)
        self.res = z(MC, self.flt)
        self.sep = z(MC, self.flt)
        self.softp = z(5, self.flt)
        self.sflag = z(2, int)
        self.pa = z(MC, int)
        self.pb = z(MC, int)
        self.feat = z(MC, int)
        self.dyn = z(NB, int)
        self.acc = z(NB * 6, wp.int64)
        self.col = z(MC, int)
        self.order = z(MC, int)
        self.skey = z(2 * MC, wp.int64)
        self.sval = z(2 * MC, int)
        self.pkey = z(MC, wp.int64)
        self.claim = z(NB, wp.int64)
        self.left = z(1, int)
        self.cstart = z(_NCOLMAX, int)
        self.ccount = z(_NCOLMAX, int)
        self.cbody = z(2 * MC, int)
        self._bstart = None
        self._bcount = None
        self.c_pa = z(MC, int)
        self.c_pb = z(MC, int)
        self.c_feat = z(MC, int)
        self.c_Minv = z(NB, self.mat66t)
        self.mism = z(1, int)
        self._cached_C = -1
        self.n_colors = 0
        self.cache_hits = 0
        self.cache_misses = 0
        self._graph = None
        self._graph_key = None
        self._counts_host = None
        self.prof = {}

    @contextmanager
    def _t(self, name):
        self.wp.synchronize_device(self.dev)
        t0 = time.perf_counter()
        yield
        self.wp.synchronize_device(self.dev)
        self.prof[name] = (time.perf_counter() - t0) * 1000.0

    def upload(self, bs: BlockScene):
        wp = self.wp
        (C, NB) = (bs.n_c, bs.n_b)
        if C > self.MC:
            raise ValueError(f'NCPSolverGPU: {C} contacts exceeds capacity {self.MC}')
        if NB != self.NB:
            raise ValueError(f'NCPSolverGPU: scene has {NB} bodies, solver built for {self.NB}')
        nt = self.np_t
        cp = lambda dst, src, t, n: wp.copy(dst, wp.array(src, dtype=t, device=self.dev), count=n)
        cp(self.Ja, np.ascontiguousarray(bs.Ja, nt), self.mat36t, C)
        cp(self.Jb, np.ascontiguousarray(bs.Jb, nt), self.mat36t, C)
        cp(self.mu, np.ascontiguousarray(bs.mu, nt), self.flt, C)
        cp(self.cbias, np.ascontiguousarray(bs.bias, nt), self.vec3t, C)
        sep = np.zeros(C) if bs.sep is None else np.asarray(bs.sep, np.float64).reshape(C)
        cp(self.sep, np.ascontiguousarray(sep, nt), self.flt, C)
        cp(self.pa, np.ascontiguousarray(bs.body_pairs[:, 0], np.int32), int, C)
        cp(self.pb, np.ascontiguousarray(bs.body_pairs[:, 1], np.int32), int, C)
        cp(self.feat, np.ascontiguousarray(bs.feature, np.int32), int, C)
        cp(self.Minv, np.ascontiguousarray(bs.Minv, nt), self.mat66t, NB)
        self._pairs_host = np.ascontiguousarray(bs.body_pairs, np.int32)
        cp(self.vfree, np.ascontiguousarray(bs.v_free, nt), self.vec6t, NB)
        dyn = (np.abs(bs.Minv).reshape(NB, 36).max(axis=1) > 0).astype(np.int32)
        cp(self.dyn, dyn, int, NB)
        self.C = C

    def _color(self, C: int) -> int:
        """Deterministic colouring: canonical radix sort on (body_a, body_b, feature), rank-hashed
        Jones-Plassmann priority, atomic-min claim per body, counting sort by colour."""
        (wp, K, d) = (self.wp, self.K, self.dev)
        wp.launch(K['sortkey'], C, inputs=[C, self.pa, self.pb, self.feat, wp.int64(self.NB + 1), self.skey, self.sval], device=d)
        wp.utils.radix_sort_pairs(self.skey, self.sval, C)
        wp.launch(K['prio'], C, inputs=[C, self.sval, self.pkey], device=d)
        wp.copy(self.col, wp.array(np.full(C, -1, np.int32), dtype=int, device=d), count=C)
        cur = 0
        while True:
            wp.launch(K['clear_claim'], self.NB, inputs=[self.claim], device=d)
            self.left.zero_()
            wp.launch(K['claim'], C, inputs=[C, self.col, self.pa, self.pb, self.dyn, self.pkey, self.claim], device=d)
            wp.launch(K['win'], C, inputs=[C, cur, self.col, self.pa, self.pb, self.dyn, self.pkey, self.claim, self.left], device=d)
            cur += 1
            if int(self.left.numpy()[0]) == 0:
                break
            if cur > C:
                raise RuntimeError('ncp_gpu._color: colouring did not terminate')
        if cur > _NCOLMAX:
            raise RuntimeError(f'ncp_gpu._color: {cur} colours exceeds capacity {_NCOLMAX}')
        wp.launch(K['colkey'], C, inputs=[C, self.col, self.skey, self.sval], device=d)
        wp.utils.radix_sort_pairs(self.skey, self.sval, C)
        wp.copy(self.order, self.sval, count=C)
        self.ccount.zero_()
        wp.launch(K['colcount'], C, inputs=[C, self.col, self.ccount], device=d)
        wp.utils.array_scan(self.ccount[:cur], self.cstart[:cur], False)
        self._counts_host = self.ccount.numpy()[:cur].copy()
        order = self.order.numpy()[:C]
        starts = self.cstart.numpy()[:cur]
        pr = self._pairs_host
        (lists, bstart, bcount) = ([], np.zeros(cur, np.int32), np.zeros(cur, np.int32))
        pos = 0
        for ci in range(cur):
            idx = order[starts[ci]:starts[ci] + int(self._counts_host[ci])]
            bd = np.concatenate([pr[idx, 0], pr[idx, 1]])
            bd = bd[bd >= 0]
            lists.append(bd)
            bstart[ci] = pos
            bcount[ci] = len(bd)
            pos += len(bd)
        flat = np.concatenate(lists).astype(np.int32) if lists else np.zeros(0, np.int32)
        wp.copy(self.cbody, wp.array(flat, dtype=int, device=d), count=max(len(flat), 1))
        (self._bstart, self._bcount) = (bstart, bcount)
        return cur

    def _color_cached(self, C: int) -> bool:
        """True = the cached colouring is reused (pairs, feature ids and masses all unchanged)."""
        (wp, d) = (self.wp, self.dev)
        hit = False
        if self.use_cache and C == self._cached_C:
            self.mism.zero_()
            wp.launch(self.K['cmp'], max(C, self.NB), inputs=[C, self.NB, self.pa, self.pb, self.feat, self.c_pa, self.c_pb, self.c_feat, self.Minv, self.c_Minv, self.mism], device=d)
            hit = int(self.mism.numpy()[0]) == 0
        if hit:
            self.cache_hits += 1
            return True
        self.n_colors = self._color(C)
        wp.copy(self.c_pa, self.pa, count=C)
        wp.copy(self.c_pb, self.pb, count=C)
        wp.copy(self.c_feat, self.feat, count=C)
        wp.copy(self.c_Minv, self.Minv, count=self.NB)
        self._cached_C = C
        self.cache_misses += 1
        self._graph = None
        return False

    def _iter_launches(self):
        (wp, K, d) = (self.wp, self.K, self.dev)
        for ci in range(self.n_colors):
            n = int(self._counts_host[ci])
            if n == 0:
                continue
            wp.launch(K['sweep'], n, inputs=[ci, self.cstart, self.ccount, self.order, self.pa, self.pb, self.Ja, self.Jb, self.Minv, self.v, self.cbias, self.mu, self.rho, self.desaxce, self.lam, self.acc, self.sep, self.softp, self.sflag], device=d)
            wp.launch(K['apply_c'], int(self._bcount[ci]), inputs=[int(self._bcount[ci]), int(self._bstart[ci]), self.cbody, self.vfree, self.v, self.acc], device=d)

    def _resid_now(self, C: int) -> float:
        (wp, K, d) = (self.wp, self.K, self.dev)
        wp.launch(K['resid'], C, inputs=[C, self.pa, self.pb, self.Ja, self.Jb, self.v, self.cbias, self.mu, self.rho, self.desaxce, self.lam, self.res, self.uout, self.sep, self.softp, self.sflag], device=d)
        return float(np.max(self.res.numpy()[:C]))

    def _set_flags(self, soft_on: int, use_bias: int):
        self.wp.copy(self.sflag, self.wp.array(np.array([soft_on, use_bias], np.int32), dtype=int, device=self.dev), count=2)

    def solve(self, bs: BlockScene, iters: int=200, warm: Optional[np.ndarray]=None, measure_residual: bool=True, history_every: int=0, soft: Optional[np.ndarray]=None, relax_iters: int=0) -> SolveResult:
        """`soft` = the five scalars of `SoftStep.params(h)` (bias_rate, mass_scale, impulse_scale,
        inv_h, contact_speed); None = the Baumgarte path of `bs.bias`. With `relax_iters > 0` the
        solve runs `iters` biased iterations, records the velocity at that point in
        `SolveResult.v_bias` (the one soft-step integrates positions with) and then `relax_iters`
        further iterations with use_bias = False, continuing from the same lambda and accumulator."""
        (wp, K, d) = (self.wp, self.K, self.dev)
        with self._t('upload'):
            self.upload(bs)
        C = self.C
        if C == 0:
            return SolveResult(np.zeros((0, 3)), bs.v_free.copy(), 0.0, 0, 0, False, dict(self.prof), v_bias=bs.v_free.copy())
        with self._t('setup'):
            if soft is None:
                self._set_flags(0, 1)
            else:
                wp.copy(self.softp, wp.array(np.ascontiguousarray(soft, self.np_t), dtype=self.flt, device=d), count=5)
                self._set_flags(1, 1)
            wp.launch(K['rho'], C, inputs=[C, self.pa, self.pb, self.Ja, self.Jb, self.Minv, self.rho], device=d)
            self.acc.zero_()
            if warm is None:
                self.lam.zero_()
            else:
                wp.copy(self.lam, wp.array(np.ascontiguousarray(warm, self.np_t), dtype=self.vec3t, device=d), count=C)
                wp.launch(K['scatter'], C, inputs=[C, self.pa, self.pb, self.Ja, self.Jb, self.Minv, self.lam, self.acc], device=d)
            wp.launch(K['apply'], self.NB, inputs=[self.NB, self.vfree, self.v, self.acc], device=d)
        with self._t('color'):
            hit = self._color_cached(C)
        with self._t('iter'):
            if self.want_graph:
                if self._graph is None:
                    try:
                        wp.load_module(device=d)
                        with wp.ScopedCapture(device=d) as cap:
                            self._iter_launches()
                        self._graph = cap.graph
                    except Exception:
                        self.want_graph = False
                        self._graph = None

            def _run(n):
                if self._graph is not None:
                    for _ in range(n):
                        wp.capture_launch(self._graph)
                else:
                    for _ in range(n):
                        self._iter_launches()
            if history_every > 0:
                hist = {'iter': [], 'res': [], 'lam': []}
                hist['iter'].append(0)
                hist['res'].append(self._resid_now(C))
                hist['lam'].append(self.lam.numpy()[:C].astype(np.float64).copy())
                done = 0
                while done < iters:
                    n = min(history_every, iters - done)
                    _run(n)
                    done += n
                    hist['iter'].append(done)
                    hist['res'].append(self._resid_now(C))
                    hist['lam'].append(self.lam.numpy()[:C].astype(np.float64).copy())
            else:
                hist = None
                _run(iters)
            v_bias = None
            if relax_iters > 0:
                v_bias = self.v.numpy().astype(np.float64).copy()
                self._set_flags(1 if soft is not None else 0, 0)
                _run(relax_iters)
        r = 0.0
        if measure_residual:
            with self._t('resid'):
                wp.launch(K['resid'], C, inputs=[C, self.pa, self.pb, self.Ja, self.Jb, self.v, self.cbias, self.mu, self.rho, self.desaxce, self.lam, self.res, self.uout, self.sep, self.softp, self.sflag], device=d)
                r = float(np.max(self.res.numpy()[:C])) if C else 0.0
        return SolveResult(lam=self.lam.numpy()[:C].astype(np.float64).copy(), v=self.v.numpy().astype(np.float64).copy(), residual=r, iters=iters + int(relax_iters), n_colors=self.n_colors, cache_hit=hit, ms=dict(self.prof), history=hist, v_bias=v_bias, contact_model='SOFT_STEP' if soft is not None else 'NCP_RIGID' if self.desaxce else 'CCP', bias_applied=relax_iters <= 0)

    def coloring(self):
        """(colour per contact, colour-grouped order, per-colour offset, per-colour count)."""
        C = self.C
        return (self.col.numpy()[:C].copy(), self.order.numpy()[:C].copy(), self.cstart.numpy()[:self.n_colors].copy(), self.ccount.numpy()[:self.n_colors].copy())

    def raw_bytes(self):
        """Raw device bytes of lambda and v, for the bit-identity gate."""
        C = self.C
        return (self.lam.numpy()[:C].tobytes(), self.v.numpy().tobytes())

    def contact_velocity(self):
        return self.uout.numpy()[:self.C].astype(np.float64).copy()

def solve_scene_gpu(bs: BlockScene, iters: int=200, dtype: str='float32', **kw) -> SolveResult:
    s = NCPSolverGPU(bs.n_b, max(bs.n_c, 1), dtype=dtype, **kw)
    return s.solve(bs, iters=iters)

def _quat_to_R(q):
    (x, y, z, w) = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])

@dataclass
class Box:
    half: np.ndarray
    mass: float
    x: np.ndarray
    q: np.ndarray = field(default_factory=lambda : np.array([0.0, 0.0, 0.0, 1.0]))
    v: np.ndarray = field(default_factory=lambda : np.zeros(3))
    w: np.ndarray = field(default_factory=lambda : np.zeros(3))

    @property
    def Ib(self):
        h = self.half * 2.0
        return self.mass / 12.0 * np.diag([h[1] ** 2 + h[2] ** 2, h[0] ** 2 + h[2] ** 2, h[0] ** 2 + h[1] ** 2])
_CORNER = np.array([[sx, sy, sz] for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)])

def _corners(b: Box):
    R = _quat_to_R(b.q)
    return b.x + _CORNER * b.half @ R.T

def _face_corners(half, x, R, axis, sgn):
    """The 4 corners of the (axis, sgn) face, in cyclic order."""
    (i, j, k) = (axis, (axis + 1) % 3, (axis + 2) % 3)
    n = R[:, i] * sgn
    c = x + n * half[i]
    (u, v) = (R[:, j] * half[j], R[:, k] * half[k])
    return (np.array([c + u + v, c - u + v, c - u - v, c + u - v]), n, c)

def _clip(poly, c, ax, ha):
    """Sutherland-Hodgman clip of `poly` against |(p-c).ax| <= ha."""
    for s in (1.0, -1.0):
        out = []
        n = len(poly)
        for i in range(n):
            (p, q) = (poly[i], poly[(i + 1) % n])
            dp = s * np.dot(p - c, ax) - ha
            dq = s * np.dot(q - c, ax) - ha
            if dp <= 0:
                out.append(p)
            if dp * dq < 0:
                out.append(p + (q - p) * (dp / (dp - dq)))
        poly = out
        if not poly:
            return []
    return poly

def box_plane_contacts(bi, b: Box, margin=0.002):
    """Vertex-plane contacts against the static ground z = 0. normal points from ground to box."""
    out = []
    for (k, p) in enumerate(_corners(b)):
        if p[2] < margin:
            out.append((bi, -1, p.copy(), np.array([0.0, 0.0, 1.0]), float(p[2]), k))
    return out

def box_box_contacts(ia, A: Box, ib, B: Box, margin=0.002):
    """SAT + reference/incident face clipping; up to 4 points per face pair. The normal points from B
    to A. Face axes are preferred over edge-cross axes by an absolute tolerance: for two nearly
    parallel boxes an edge axis is numerically a hair deeper than the true face axis, and taking it
    collapses a 4-point face manifold to a single point (measured: a 3-box column loses one interface
    to a 1-point manifold within 20 steps and falls over)."""
    (Ra, Rb) = (_quat_to_R(A.q), _quat_to_R(B.q))
    d = B.x - A.x

    def sep_of(ax):
        ra = float(np.sum(np.abs(Ra.T @ ax) * A.half))
        rb = float(np.sum(np.abs(Rb.T @ ax) * B.half))
        proj = float(np.dot(d, ax))
        return (abs(proj) - (ra + rb), ax if proj >= 0 else -ax)
    bestF = None
    for (kind, R) in (('A', Ra), ('B', Rb)):
        for i in range(3):
            (sp, ax) = sep_of(R[:, i])
            if sp > margin:
                return []
            better = 0.0001 if bestF is not None and kind == 'B' and (bestF[1] == 'A') else 1e-12
            if bestF is None or sp > bestF[0] + better:
                bestF = (sp, kind, i, ax)
    bestE = None
    for i in range(3):
        for j in range(3):
            cr = np.cross(Ra[:, i], Rb[:, j])
            nn = np.linalg.norm(cr)
            if nn <= 1e-06:
                continue
            (sp, ax) = sep_of(cr / nn)
            if sp > margin:
                return []
            if bestE is None or sp > bestE[0] + 1e-12:
                bestE = (sp, 'E', 3 * i + j, ax)
    if bestE is not None and bestE[0] > bestF[0] + 0.001:
        (sep, kind, idx, nAB) = bestE
        (ca, cb) = (_corners(A), _corners(B))
        pa = ca[int(np.argmax(ca @ nAB))]
        pb = cb[int(np.argmin(cb @ nAB))]
        return [(ia, ib, 0.5 * (pa + pb), -nAB, float(sep), 100 + idx)]
    (sep, kind, idx, nAB) = bestF
    if kind == 'A':
        (rhalf, rx, rR, ihalf, ix, iR) = (A.half, A.x, Ra, B.half, B.x, Rb)
        (raxis, rsgn) = (idx, 1.0 if np.dot(Ra[:, idx], nAB) >= 0 else -1.0)
        nBA = -nAB
    else:
        (rhalf, rx, rR, ihalf, ix, iR) = (B.half, B.x, Rb, A.half, A.x, Ra)
        (raxis, rsgn) = (idx, 1.0 if np.dot(Rb[:, idx], -nAB) >= 0 else -1.0)
        nBA = nAB
    (_, rn, rc) = _face_corners(rhalf, rx, rR, raxis, rsgn)
    dots = [float(np.dot(iR[:, a] * sgn, rn)) for a in range(3) for sgn in (1.0, -1.0)]
    mmin = int(np.argmin(dots))
    (iaxis, isgn) = (mmin // 2, 1.0 if mmin % 2 == 0 else -1.0)
    (poly, _, _) = _face_corners(ihalf, ix, iR, iaxis, isgn)
    (j, k) = ((raxis + 1) % 3, (raxis + 2) % 3)
    poly = _clip(list(poly), rc, rR[:, j], rhalf[j])
    if not poly:
        return []
    poly = _clip(poly, rc, rR[:, k], rhalf[k])
    if not poly:
        return []
    pts = [(float(np.dot(p - rc, rn)), p) for p in poly]
    pts = [pq for pq in pts if pq[0] < margin]
    if not pts:
        return []
    if len(pts) > 4:
        pts = sorted(pts, key=lambda t: t[0])[:4]
    base = (raxis * 2 + (0 if rsgn > 0 else 1)) * 8
    return [(ia, ib, p.copy(), nBA.copy(), g, base + n) for (n, (g, p)) in enumerate(pts)]

def _basis(n):
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = a - np.dot(a, n) * n
    t1 /= np.linalg.norm(t1)
    return (t1, np.cross(n, t1))

class BoxWorld:
    """Minimal rigid-body world: free boxes + a static ground plane at z = 0. Self-contained contact
    generation (vertex-plane + SAT box-box), velocity-level step with the GPU NCP solver."""

    def __init__(self, dt=1.0 / 240.0, mu=0.5, gravity=True, dtype='float32', solver_kw=None, iters=200, margin=0.002, desaxce=True, sleep: Optional['SleepConfig']=None, substeps=1, warm_start=None, soft: Optional['SoftStep']=None):
        self.bodies: list[Box] = []
        self.soft = soft
        self.warm_start = soft is not None if warm_start is None else bool(warm_start)
        self._warm = {}
        self.dt = float(dt)
        self.substeps = int(substeps)
        self.sleep_cfg = sleep
        self.sleeper: Optional[Sleeper] = None
        self.n_asleep = 0
        self.mu = float(mu)
        self.gravity = gravity
        self.iters = int(iters)
        self.margin = float(margin)
        self.dtype = dtype
        self.desaxce = desaxce
        self.solver_kw = solver_kw or {}
        self._solver = None
        self.last = None
        self.n_contacts = 0
        self.n_colors = 0

    def add(self, half, mass, x, v=None):
        self.bodies.append(Box(np.asarray(half, float), float(mass), np.asarray(x, float), v=np.zeros(3) if v is None else np.asarray(v, float)))
        return len(self.bodies) - 1

    def contacts(self):
        out = []
        for (i, b) in enumerate(self.bodies):
            out += box_plane_contacts(i, b, self.margin)
        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)):
                out += box_box_contacts(i, self.bodies[i], j, self.bodies[j], self.margin)
        return out

    def _Minv(self):
        n = len(self.bodies)
        Mi = np.zeros((n, 6, 6))
        for (i, b) in enumerate(self.bodies):
            R = _quat_to_R(b.q)
            Mi[i, :3, :3] = np.eye(3) / b.mass
            Mi[i, 3:, 3:] = R @ np.linalg.inv(b.Ib) @ R.T
        return Mi

    def build_scene(self, f_ext=None, dt=None, contacts=None, asleep=None) -> BlockScene:
        n = len(self.bodies)
        dt = self.dt if dt is None else float(dt)
        Mi = self._Minv()
        vf = np.zeros((n, 6))
        for (i, b) in enumerate(self.bodies):
            f = np.zeros(3)
            if self.gravity:
                f = f + np.array([0.0, 0.0, -_G]) * b.mass
            if f_ext is not None:
                f = f + np.asarray(f_ext(i), float)
            vf[i, :3] = b.v + dt * Mi[i, :3, :3] @ f
            vf[i, 3:] = b.w
            if asleep is not None and asleep[i]:
                vf[i] = 0.0
        cs = self.contacts() if contacts is None else contacts
        nc = len(cs)
        Ja = np.zeros((nc, 3, 6))
        Jb = np.zeros((nc, 3, 6))
        pairs = np.zeros((nc, 2), np.int32)
        bias = np.zeros((nc, 3))
        feat = np.zeros(nc, np.int32)
        mu = np.full(nc, self.mu)
        sep = np.zeros(nc)
        for (c, (ia, ib, p, nrm, g, fid)) in enumerate(cs):
            (t1, t2) = _basis(nrm)
            D = np.stack([nrm, t1, t2])
            ra = p - self.bodies[ia].x
            Ja[c, :, :3] = D
            Ja[c, :, 3:] = np.cross(ra[None, :], D)
            if ib >= 0:
                rb = p - self.bodies[ib].x
                Jb[c, :, :3] = -D
                Jb[c, :, 3:] = -np.cross(rb[None, :], D)
            pairs[c] = (ia, ib)
            sep[c] = g
            if self.soft is None:
                bias[c, 0] = g / dt if g > 0 else -_BETA * max(-g - _SLOP, 0.0) / dt
            feat[c] = ia * 4096 + (ib + 1) * 256 + fid
        return BlockScene(Ja=Ja, Jb=Jb, Minv=Mi, v_free=vf, mu=mu, body_pairs=pairs, dt=dt, bias=bias, feature=feat, sep=sep)

    def step(self, f_ext=None, history_every=0):
        """One macro step. With `substeps` > 1 the macro step is split into S sub-steps of dt/S, each
        one a full contact solve of `iters` iterations (Catto small-steps: the iteration budget moves
        from the solver into the time discretisation)."""
        r = None
        for _ in range(max(1, self.substeps)):
            r = self._substep(f_ext, self.dt / max(1, self.substeps), history_every)
        return r

    def _substep(self, f_ext, dt, history_every=0):
        n = len(self.bodies)
        cs = self.contacts()
        asleep = None
        if self.sleep_cfg is not None and self.sleep_cfg.enabled:
            if self.sleeper is None:
                self.sleeper = Sleeper(n, self.sleep_cfg)
            pairs_all = np.array([[c[0], c[1]] for c in cs], np.int64).reshape(-1, 2)
            vel = np.array([b.v for b in self.bodies]).reshape(n, 3)
            omg = np.array([b.w for b in self.bodies]).reshape(n, 3)
            forced = np.zeros(n, bool)
            if f_ext is not None:
                forced = np.array([np.any(np.asarray(f_ext(i), float) != 0.0) for i in range(n)])
            asleep = self.sleeper.begin(pairs_all, vel, omg, forced)
            for i in range(n):
                if asleep[i]:
                    self.bodies[i].v = np.zeros(3)
                    self.bodies[i].w = np.zeros(3)
            cs = [c for c in cs if not (c[0] >= 0 and asleep[c[0]] or (c[1] >= 0 and asleep[c[1]]))]
            self.n_asleep = int(asleep.sum())
        bs = self.build_scene(f_ext, dt=dt, contacts=cs, asleep=asleep)
        if self._solver is None:
            self._solver = NCPSolverGPU(n, max(4096, len(cs)), dtype=self.dtype, desaxce=self.desaxce, **self.solver_kw)
        if bs.n_c:
            warm = None
            if self.warm_start and self._warm:
                warm = np.array([self._warm.get(int(f), (0.0, 0.0, 0.0)) for f in bs.feature])
            if self.soft is None:
                r = self._solver.solve(bs, iters=self.iters, warm=warm, history_every=history_every)
            else:
                (nb_it, nr_it) = self.soft.split(self.iters)
                r = self._solver.solve(bs, iters=nb_it, warm=warm, history_every=history_every, soft=self.soft.params(dt), relax_iters=nr_it)
                self.last_split = (nb_it, nr_it)
            if self.warm_start:
                self._warm = {int(f): tuple(l) for (f, l) in zip(bs.feature, r.lam)}
            v = r.v
            self.n_colors = r.n_colors
        else:
            v = bs.v_free
            self.n_colors = 0
            r = SolveResult(np.zeros((0, 3)), v, 0.0, 0, 0, False, v_bias=v)
        self.n_contacts = bs.n_c
        self.last = r
        vp = v if r.v_bias is None else r.v_bias
        for (i, b) in enumerate(self.bodies):
            if asleep is not None and asleep[i]:
                continue
            b.v = v[i, :3].copy()
            b.w = v[i, 3:].copy()
            b.x = b.x + vp[i, :3] * dt
            wq = np.array([vp[i, 3], vp[i, 4], vp[i, 5], 0.0])
            qq = b.q + 0.5 * _qmul(wq, b.q) * dt
            b.q = qq / np.linalg.norm(qq)
        if self.sleeper is not None:
            self.sleeper.end(np.array([b.v for b in self.bodies]).reshape(n, 3), np.array([b.w for b in self.bodies]).reshape(n, 3))
        return r

    def penetration(self):
        """Deepest penetration in the current configuration (m)."""
        cs = self.contacts()
        return max([0.0] + [-g for (_, _, _, _, g, _) in cs])

    def kinetic_energy(self):
        return float(sum((0.5 * b.mass * b.v @ b.v + 0.5 * b.w @ (b.Ib @ b.w) for b in self.bodies)))

def _qmul(a, b):
    (x1, y1, z1, w1) = a
    (x2, y2, z2, w2) = b
    return np.array([w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2, w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2, w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2])

def contact_labels(lam: np.ndarray, mu: np.ndarray, tol_lam=1e-10, tol_rel=1e-06) -> np.ndarray:
    """Active-set label per contact as int: 0 = open, 1 = stick, 2 = slip. Same rule as
    `ncp_ref.classify` (scale = max |lambda| over the problem)."""
    L = np.asarray(lam, np.float64).reshape(-1, 3)
    mu = np.asarray(mu, np.float64).reshape(-1)
    scale = max(float(np.max(np.abs(L))), 1e-30)
    ln = L[:, 0]
    lt = np.linalg.norm(L[:, 1:3], axis=1)
    out = np.ones(len(L), np.int8)
    out[lt >= mu * ln * (1.0 - tol_rel) - tol_lam] = 2
    out[ln <= tol_lam + tol_rel * scale] = 0
    return out

@dataclass
class SleepConfig:
    enabled: bool = True
    v_sleep: float = 0.001
    w_sleep: float = 0.001
    k_steps: int = 8

def connected_islands(n_b: int, pairs: np.ndarray) -> np.ndarray:
    """Union-find (min-label propagation + pointer jumping) over the DYNAMIC-body contact pairs.
    A static endpoint (-1) never unions, so two boxes resting on the same ground plane are two
    islands. Returns the island label per body (the smallest body index of its component), a pure
    function of the pair set -- the label does not depend on the order the pairs arrive in."""
    lab = np.arange(n_b, dtype=np.int64)
    pairs = np.asarray(pairs, np.int64).reshape(-1, 2)
    if len(pairs) == 0:
        return lab
    m = (pairs[:, 0] >= 0) & (pairs[:, 1] >= 0)
    (a, b) = (pairs[m, 0], pairs[m, 1])
    if len(a) == 0:
        return lab
    while True:
        prev = lab
        mn = np.minimum(lab[a], lab[b])
        lab = lab.copy()
        np.minimum.at(lab, a, mn)
        np.minimum.at(lab, b, mn)
        for _ in range(64):
            nl = lab[lab]
            if np.array_equal(nl, lab):
                break
            lab = nl
        if np.array_equal(lab, prev):
            return lab

class Sleeper:
    """Island sleeping.

    An island (see `connected_islands`) is put to sleep once EVERY one of its bodies has stayed below
    both (v_sleep, w_sleep) for `k_steps` consecutive steps. A sleeping island's contacts are dropped
    before colouring, so they are neither coloured nor swept, and its bodies skip integration
    entirely (they do not even see gravity).

    Wake conditions:
      * a contact pair that was not in the previous step's pair set resets the still-counter of both
        its endpoints. A new pair to an awake or moving body also puts that body in the island, so
        the island's minimum counter is that body's and the whole island wakes with it;
      * a body carrying a nonzero external force is reset;
      * a body moving faster than the thresholds is reset.

    Every input is a deterministic function of the step state (the velocities come back bit-identical
    from the GPU, the pair set and the island labels are order-independent), so the sleep mask is
    deterministic too."""

    def __init__(self, n_b: int, cfg: SleepConfig):
        self.cfg = cfg
        self.n_b = int(n_b)
        self.still = np.zeros(n_b, np.int64)
        self.asleep = np.zeros(n_b, bool)
        self.labels = np.arange(n_b, dtype=np.int64)
        self.prev_keys = None
        self._keys_raw = None

    def begin(self, pairs: np.ndarray, v: np.ndarray, w: np.ndarray, forced: Optional[np.ndarray]=None) -> np.ndarray:
        cfg = self.cfg
        pairs = np.asarray(pairs, np.int64).reshape(-1, 2)
        keys = pairs[:, 0] * (self.n_b + 1) + pairs[:, 1] + 1 if len(pairs) else np.zeros(0, np.int64)
        same = self._keys_raw is not None and np.array_equal(keys, self._keys_raw)
        if not same:
            if self.prev_keys is not None and len(keys):
                fresh = ~np.isin(keys, self.prev_keys)
                if fresh.any():
                    bd = np.concatenate([pairs[fresh, 0], pairs[fresh, 1]])
                    self.still[bd[bd >= 0]] = 0
            self.prev_keys = np.unique(keys)
            self._keys_raw = keys.copy()
        if forced is not None:
            self.still[np.asarray(forced, bool)] = 0
        moving = (np.linalg.norm(v, axis=1) > cfg.v_sleep) | (np.linalg.norm(w, axis=1) > cfg.w_sleep)
        self.still[moving] = 0
        if not same:
            self.labels = connected_islands(self.n_b, pairs)
        isl = np.full(self.n_b, np.iinfo(np.int64).max, np.int64)
        np.minimum.at(isl, self.labels, self.still)
        self.asleep = isl[self.labels] >= cfg.k_steps
        return self.asleep

    def end(self, v: np.ndarray, w: np.ndarray):
        cfg = self.cfg
        quiet = (np.linalg.norm(v, axis=1) <= cfg.v_sleep) & (np.linalg.norm(w, axis=1) <= cfg.w_sleep)
        awake = ~self.asleep
        self.still[awake & quiet] += 1
        self.still[awake & ~quiet] = 0

    def n_islands_asleep(self) -> int:
        return len(np.unique(self.labels[self.asleep])) if self.asleep.any() else 0

def _quat_to_R_batch(q):
    (x, y, z, w) = (q[:, 0], q[:, 1], q[:, 2], q[:, 3])
    R = np.empty((len(q), 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z)
    R[:, 0, 1] = 2 * (x * y - z * w)
    R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w)
    R[:, 1, 1] = 1 - 2 * (x * x + z * z)
    R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w)
    R[:, 2, 1] = 2 * (y * z + x * w)
    R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R

class LatticeWorld:
    """Two-layer N-body box lattice, the geometry of the motion engine's N10000 fixture
    (`docs/RT_BODY_PAIR_EXPERIMENT.md`: 18-point box cloud, horizontal spacing 0.35, lower centre
    0.0999, vertical spacing 0.1999, one ground pair per lower body and one vertical pair per upper
    body; N = 10000 -> 5041 lower + 4959 upper). The pair set is structural, so the host work is O(n)
    vectorised numpy -- BoxWorld's O(n^2) SAT is unusable at n = 10000. Each pair is expanded to a
    4-point corner manifold, so n_c = 4 N.

    Simplification recorded: the contact normals are the world +z axis and the interface plane of a
    stacked pair is taken through the lower box centre plus R_b (0,0,half). The fixture's bodies rest
    flat, so this is exact at rest and O(theta * half) under the small rotations a pushed body picks
    up; it is a timing/sleeping fixture, not a trajectory reference."""

    def __init__(self, n=10000, half=0.1, spacing=0.35, z0=0.0999, dz=0.1999, mass=1.0, dt=1.0 / 240.0, mu=0.5, iters=40, dtype='float32', sleep=None, solver_kw=None, points_per_pair=4):
        g = int(np.ceil(np.sqrt(n / 2.0)))
        self.n_low = min(g * g, n)
        self.n_up = n - self.n_low
        if self.n_up > self.n_low:
            raise ValueError('lattice: more upper than lower bodies')
        self.n = int(n)
        self.half = float(half)
        self.dt = float(dt)
        self.mu = float(mu)
        self.iters = int(iters)
        self.dtype = dtype
        self.mass = np.full(n, float(mass))
        (ix, iy) = np.divmod(np.arange(self.n_low), g)
        x = np.zeros((n, 3))
        x[:self.n_low, 0] = (ix - (g - 1) / 2.0) * spacing
        x[:self.n_low, 1] = (iy - (g - 1) / 2.0) * spacing
        x[:self.n_low, 2] = z0
        x[self.n_low:, :2] = x[:self.n_up, :2]
        x[self.n_low:, 2] = z0 + dz
        self.x = x
        self.v = np.zeros((n, 3))
        self.w = np.zeros((n, 3))
        self.q = np.tile(np.array([0.0, 0.0, 0.0, 1.0]), (n, 1))
        h = float(half)
        Ib = mass / 12.0 * ((2 * h) ** 2 + (2 * h) ** 2)
        self.Ib = np.full(3, Ib)
        if points_per_pair == 4:
            off = np.array([[-h, -h, -h], [h, -h, -h], [-h, h, -h], [h, h, -h]])
        else:
            off = np.array([[0.0, 0.0, -h]])
        K = len(off)
        ca = np.concatenate([np.repeat(np.arange(self.n_low), K), np.repeat(np.arange(self.n_low, n), K)])
        cb = np.concatenate([np.full(self.n_low * K, -1), np.repeat(np.arange(self.n_up), K)])
        (self.ca, self.cb) = (ca.astype(np.int64), cb.astype(np.int64))
        self.off = np.tile(off, (self.n_low + self.n_up, 1))
        self.feat = (np.arange(len(ca)) % K).astype(np.int32) + 8 * (self.cb + 1).astype(np.int32)
        self.pairs = np.stack([self.ca, self.cb], axis=1)
        self.sleep_cfg = sleep
        self.sleeper = None if sleep is None or not sleep.enabled else Sleeper(n, sleep)
        self.solver = NCPSolverGPU(n, len(ca), dtype=dtype, **solver_kw or {})
        self.n_contacts = 0
        self.n_colors = 0
        self.n_asleep = 0
        self.last = None

    def _gaps_points(self, idx=None):
        ca = self.ca if idx is None else self.ca[idx]
        cb = self.cb if idx is None else self.cb[idx]
        off = self.off if idx is None else self.off[idx]
        R = _quat_to_R_batch(self.q)
        pa = self.x[ca] + np.einsum('nij,nj->ni', R[ca], off)
        gap = np.empty(len(ca))
        gnd = cb < 0
        gap[gnd] = pa[gnd, 2]
        up = ~gnd
        plane = self.x[cb[up], 2] + R[cb[up], 2, 2] * self.half
        gap[up] = pa[up, 2] - plane
        return (gap, pa)

    def build_scene(self, f_ext=None, keep=None, asleep=None) -> BlockScene:
        (n, dt, h) = (self.n, self.dt, self.half)
        idx = np.arange(len(self.ca)) if keep is None else np.flatnonzero(keep)
        (gap, pa) = self._gaps_points(idx)
        nc = len(idx)
        Minv = np.zeros((n, 6, 6))
        Minv[:, 0, 0] = Minv[:, 1, 1] = Minv[:, 2, 2] = 1.0 / self.mass
        R = _quat_to_R_batch(self.q)
        Iinv = R @ (np.diag(1.0 / self.Ib)[None] @ np.transpose(R, (0, 2, 1)))
        Minv[:, 3:, 3:] = Iinv
        vf = np.zeros((n, 6))
        f = np.zeros((n, 3))
        f[:, 2] = -_G * self.mass
        if f_ext is not None:
            f = f + f_ext
        vf[:, :3] = self.v + dt * f / self.mass[:, None]
        vf[:, 3:] = self.w
        if asleep is not None:
            vf[asleep] = 0.0
        D = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        (a, b) = (self.ca[idx], self.cb[idx])
        ra = pa - self.x[a]
        Ja = np.zeros((nc, 3, 6))
        Jb = np.zeros((nc, 3, 6))
        Ja[:, :, :3] = D
        Ja[:, :, 3:] = np.cross(ra[:, None, :], D[None, :, :])
        hasb = b >= 0
        rb = np.zeros((nc, 3))
        rb[hasb] = pa[hasb] - self.x[b[hasb]]
        Jb[hasb, :, :3] = -D
        Jb[hasb, :, 3:] = -np.cross(rb[hasb][:, None, :], D[None, :, :])
        g = gap
        bias = np.zeros((nc, 3))
        bias[:, 0] = np.where(g > 0, g / dt, -_BETA * np.maximum(-g - _SLOP, 0.0) / dt)
        return BlockScene(Ja=Ja, Jb=Jb, Minv=Minv, v_free=vf, mu=np.full(nc, self.mu), body_pairs=np.stack([a, b], axis=1).astype(np.int32), dt=dt, bias=bias, feature=self.feat[idx])

    def step(self, f_ext=None):
        n = self.n
        asleep = None
        keep = None
        if self.sleeper is not None:
            forced = None if f_ext is None else np.any(np.asarray(f_ext) != 0.0, axis=1)
            asleep = self.sleeper.begin(self.pairs, self.v, self.w, forced)
            self.v[asleep] = 0.0
            self.w[asleep] = 0.0
            keep = ~(asleep[self.ca] | np.where(self.cb >= 0, asleep[np.maximum(self.cb, 0)], False))
            self.n_asleep = int(asleep.sum())
        bs = self.build_scene(f_ext, keep=keep, asleep=asleep)
        self.n_contacts = bs.n_c
        if bs.n_c:
            r = self.solver.solve(bs, iters=self.iters, measure_residual=False)
            v6 = r.v
            self.n_colors = r.n_colors
        else:
            v6 = bs.v_free
            self.n_colors = 0
            r = SolveResult(np.zeros((0, 3)), v6, 0.0, 0, 0, False)
        self.last = r
        mv = np.ones(n, bool) if asleep is None else ~asleep
        self.v[mv] = v6[mv, :3]
        self.w[mv] = v6[mv, 3:]
        self.x[mv] = self.x[mv] + self.v[mv] * self.dt
        wq = np.zeros((n, 4))
        wq[:, :3] = self.w
        qd = np.stack([wq[:, 3] * self.q[:, 0] + wq[:, 0] * self.q[:, 3] + wq[:, 1] * self.q[:, 2] - wq[:, 2] * self.q[:, 1], wq[:, 3] * self.q[:, 1] - wq[:, 0] * self.q[:, 2] + wq[:, 1] * self.q[:, 3] + wq[:, 2] * self.q[:, 0], wq[:, 3] * self.q[:, 2] + wq[:, 0] * self.q[:, 1] - wq[:, 1] * self.q[:, 0] + wq[:, 2] * self.q[:, 3], wq[:, 3] * self.q[:, 3] - wq[:, 0] * self.q[:, 0] - wq[:, 1] * self.q[:, 1] - wq[:, 2] * self.q[:, 2]], axis=1)
        qq = self.q.copy()
        qq[mv] = self.q[mv] + 0.5 * qd[mv] * self.dt
        self.q = qq / np.linalg.norm(qq, axis=1, keepdims=True)
        if self.sleeper is not None:
            self.sleeper.end(self.v, self.w)
        return r

def lattice_pushed_forces(w: LatticeWorld, frac=0.1, fx=40.0):
    """A deterministic `frac` of the bodies carry a constant horizontal force. Whole ISLANDS are
    pushed (a stacked pair and its lower body together), so the pushed fraction is also the awake
    fraction: pushing single bodies would drag their stack partner through friction and leave the
    awake set larger than `frac`."""
    f = np.zeros((w.n, 3))
    stride = int(round(1.0 / frac))
    j = np.arange(0, w.n_up, stride)
    f[j, 0] = fx
    f[w.n_low + j, 0] = fx
    return f

def scene_pushed_cube(dtype='float32', steps=500, iters=200, solver_kw=None):
    """(i) ODYNSim pushed cube: 0.20 m, 15 kg, mu = 0.5, dt = 0.004 s, horizontal 147.15 N
    alternating x/y every 0.5 s, ramped over 0.3 s."""
    (dt, F, per, ramp) = (0.004, 147.15, 0.5, 0.3)
    w = BoxWorld(dt=dt, mu=0.5, dtype=dtype, iters=iters, solver_kw=solver_kw)
    w.add([0.1, 0.1, 0.1], 15.0, [0.0, 0.0, 0.1])
    (traj, pen, res) = ([], [], [])
    for s in range(steps):
        t = s * dt
        k = int(t // per)
        amp = F * min(1.0, (t - k * per) / ramp)
        f = np.array([amp, 0.0, 0.0]) if k % 2 == 0 else np.array([0.0, amp, 0.0])
        r = w.step(f_ext=lambda i, f=f: f)
        traj.append(w.bodies[0].x.copy())
        pen.append(w.penetration())
        res.append(r.residual)
    return (w, np.array(traj), np.array(pen), np.array(res))

def scene_sliding_box(dtype='float32', steps=200, iters=200, desaxce=True, v0=1.0, mu=0.3, dt=1.0 / 240.0, substeps=1, soft=None):
    """(ii) sliding box, v0 = 1 m/s, mu = 0.3: normal gap during the slide."""
    w = BoxWorld(dt=dt, mu=mu, dtype=dtype, iters=iters, desaxce=desaxce, substeps=substeps, soft=soft)
    w.add([0.1, 0.1, 0.1], 1.0, [0.0, 0.0, 0.1], v=[v0, 0.0, 0.0])
    (gap, vx) = ([], [])
    for _ in range(steps):
        w.step()
        gap.append(w.bodies[0].x[2] - 0.1)
        vx.append(w.bodies[0].v[0])
    return (w, np.array(gap), np.array(vx))

def sliding_distance(dtype='float32', iters=400, desaxce=True, v0=1.0, mu=0.3, dt=1.0 / 240.0, max_steps=4000, substeps=1, soft=None):
    """(ii) rider: total sliding distance until the box stops, against the closed form
    d = v0^2 / (2 mu g). Integration stops the first step vx <= 0. With `substeps` > 1 the macro
    step of `dt` is split into S sub-steps of dt/S, each with `iters` iterations."""
    w = BoxWorld(dt=dt, mu=mu, dtype=dtype, iters=iters, desaxce=desaxce, substeps=substeps, soft=soft)
    w.add([0.1, 0.1, 0.1], 1.0, [0.0, 0.0, 0.1], v=[v0, 0.0, 0.0])
    x0 = w.bodies[0].x[0]
    n = 0
    for n in range(1, max_steps + 1):
        w.step()
        if w.bodies[0].v[0] <= 0.0:
            break
    d = float(w.bodies[0].x[0] - x0)
    ref = v0 * v0 / (2.0 * mu * _G)
    return (d, ref, n)

def scene_massratio(ratio, dtype='float32', steps=400, iters=200, substeps=1, soft=None):
    """(iii) 3-box column, H = 0.20, top body heavy by `ratio` (geometry of
    an independent mass-ratio study script)."""
    H = 0.2
    w = BoxWorld(dt=1.0 / 240.0, mu=0.5, dtype=dtype, iters=iters, substeps=substeps, soft=soft)
    for k in range(3):
        m = 1.0 * (ratio if k == 2 else 1.0)
        w.add([H / 2, H / 2, H / 2], m, [0.0, 0.0, (k + 0.5) * H])
    pen = []
    for _ in range(steps):
        w.step()
        pen.append(w.penetration())
    return (w, np.array(pen))

def scene_tower(K=8, dtype='float32', steps=400, iters=200, sleep=None, substeps=1, soft=None):
    """(iv) K-box tower. With `sleep` the tower's single island goes to sleep once it has settled.
    Records the per-step jitter (max body speed / spin over the bodies) so the floor of RESULTS B.e
    can be read off the same run: `w.vmax_history`, `w.wmax_history`."""
    H = 0.2
    w = BoxWorld(dt=1.0 / 240.0, mu=0.5, dtype=dtype, iters=iters, sleep=sleep, substeps=substeps, soft=soft)
    for k in range(K):
        w.add([H / 2, H / 2, H / 2], 1.0, [0.0, 0.0, (k + 0.5) * H])
    (pen, asleep, vmax, wmax) = ([], [], [], [])
    for _ in range(steps):
        w.step()
        pen.append(w.penetration())
        asleep.append(w.n_asleep)
        vmax.append(max((float(np.linalg.norm(b.v)) for b in w.bodies)))
        wmax.append(max((float(np.linalg.norm(b.w)) for b in w.bodies)))
    w.asleep_history = np.array(asleep)
    w.vmax_history = np.array(vmax)
    w.wmax_history = np.array(wmax)
    return (w, np.array(pen))

def synthetic_scene(n_c, n_b=512, seed=0, mu=0.5):
    """Random contacts over `n_b` free bodies, built directly as blocks (a dense J at n_c = 4096 would
    be 12288 x 3072)."""
    rng = np.random.default_rng(seed)
    a = rng.integers(0, n_b, n_c)
    b = (a + 1 + rng.integers(0, n_b - 1, n_c)) % n_b
    nrm = rng.normal(size=(n_c, 3))
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    Ja = np.zeros((n_c, 3, 6))
    Jb = np.zeros((n_c, 3, 6))
    for c in range(n_c):
        (t1, t2) = _basis(nrm[c])
        D = np.stack([nrm[c], t1, t2])
        ra = rng.uniform(-0.1, 0.1, 3)
        rb = rng.uniform(-0.1, 0.1, 3)
        Ja[c, :, :3] = D
        Ja[c, :, 3:] = np.cross(ra[None, :], D)
        Jb[c, :, :3] = -D
        Jb[c, :, 3:] = -np.cross(rb[None, :], D)
    Ib = np.linalg.inv(1.0 / 12.0 * np.diag([0.08, 0.08, 0.08]))
    Mi = np.tile(np.block([[np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), Ib]]), (n_b, 1, 1))
    vf = rng.normal(scale=0.2, size=(n_b, 6))
    pairs = np.stack([a, b], axis=1).astype(np.int32)
    return BlockScene(Ja=Ja, Jb=Jb, Minv=Mi, v_free=vf, mu=np.full(n_c, mu), body_pairs=pairs, dt=1.0 / 240.0, bias=np.zeros((n_c, 3)), feature=np.arange(n_c, dtype=np.int32))

def random_wellposed_scene(n_c=8, n_b=8, seed=1, mu=0.4):
    """A random scene with 3 n_c <= 6 n_b and full-rank J, so the Delassus is SPD and the NCP solution
    is unique -- the case where a GPU-vs-CPU comparison in lambda itself is meaningful."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_b)
    a = perm[np.arange(n_c) % n_b]
    b = perm[(np.arange(n_c) + 3) % n_b]
    Ja = rng.normal(size=(n_c, 3, 6))
    Jb = rng.normal(size=(n_c, 3, 6))
    Mi = np.tile(np.eye(6), (n_b, 1, 1)) * rng.uniform(0.5, 2.0, (n_b, 1, 1))
    vf = rng.normal(scale=0.5, size=(n_b, 6))
    bias = np.zeros((n_c, 3))
    bias[:, 0] = -rng.uniform(0.1, 1.0, n_c)
    return BlockScene(Ja=Ja, Jb=Jb, Minv=Mi, v_free=vf, mu=np.full(n_c, mu), body_pairs=np.stack([a, b], axis=1).astype(np.int32), dt=1.0 / 240.0, bias=bias, feature=np.arange(n_c, dtype=np.int32))

def _selftest():
    print('ncp_gpu selftest')
    bs = random_wellposed_scene()
    s = NCPSolverGPU(bs.n_b, bs.n_c)
    r1 = s.solve(bs, iters=2000)
    s2 = NCPSolverGPU(bs.n_b, bs.n_c)
    r2 = s2.solve(bs, iters=2000)
    same = r1.lam.tobytes() == r2.lam.tobytes() and r1.v.tobytes() == r2.v.tobytes()
    print(f"  determinism (raw bytes): {('identical' if same else 'DIFFER')}  residual={r1.residual:.3e} colours={r1.n_colors}")
    (w, traj, pen, res) = scene_pushed_cube(steps=250)
    print(f'  (i) pushed cube: end xy=({traj[-1, 0]:.4f},{traj[-1, 1]:.4f}) max_pen={pen.max() * 1000.0:.4f} mm max_res={res.max():.3e} contacts={w.n_contacts} colours={w.n_colors}')
    for ratio in (1, 100, 1000):
        (w3, p3) = scene_massratio(ratio, steps=240)
        print(f'  (iii) mass ratio {ratio}: final_pen={p3[-1] * 1000.0:.4f} mm  KE={w3.kinetic_energy():.3e}')
    for n in (4, 32, 256, 4096):
        sc = synthetic_scene(n)
        sol = NCPSolverGPU(sc.n_b, n)
        sol.solve(sc, iters=50, measure_residual=False)
        t0 = time.perf_counter()
        for _ in range(5):
            sol.solve(sc, iters=50, measure_residual=False)
        dt = (time.perf_counter() - t0) / 5 * 1000.0
        print(f'  timing n_c={n:5d}: {dt:8.3f} ms/solve(50 it)  colours={sol.n_colors} hits={sol.cache_hits}')
    return True
if __name__ == '__main__':
    import sys
    sys.exit(0 if _selftest() else 1)
