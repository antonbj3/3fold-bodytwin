"""Scene assembly for object system identification (numpy only) -- shared by the
generator and the estimator so both sides see bit-identical (J, Minv, v_free, mu,
b_offset). That sharing is the point: generator and estimator are one model, and the
cost at the true parameters (2.42e-06 over 512 x 100 x 32 residuals at zero noise) is
the floor it implies, not a second model's disagreement.

Coupled system:  nv = 32 (fixed-base Talos, every dof actuated) + 6 (box) = 38.
nc = 5 contacts = 1 hand-box point contact + 4 box-table corner contacts, 15 rows.
Conventions and formulas are `ncp_ref`'s (build_scene / Box / _tangents /
_quat_integrate); they are re-expressed batched over envs, not changed.
"""
import numpy as np
G_VEC = np.array([0.0, 0.0, -9.81])
NVR = 32
NVB = 6
NVREAL = NVR + NVB
NVIRT = 15
NV = NVREAL + NVIRT
NC = 5
MROWS = 3 * NC
B_FLOOR = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
CORNERS = np.array([[1, 1, -1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1]], dtype=np.float64)

def quat_to_R(q):
    (w, x, y, z) = (q[..., 0], q[..., 1], q[..., 2], q[..., 3])
    R = np.empty(q.shape[:-1] + (3, 3))
    R[..., 0, 0] = 1 - 2 * (y * y + z * z)
    R[..., 0, 1] = 2 * (x * y - w * z)
    R[..., 0, 2] = 2 * (x * z + w * y)
    R[..., 1, 0] = 2 * (x * y + w * z)
    R[..., 1, 1] = 1 - 2 * (x * x + z * z)
    R[..., 1, 2] = 2 * (y * z - w * x)
    R[..., 2, 0] = 2 * (x * z - w * y)
    R[..., 2, 1] = 2 * (y * z + w * x)
    R[..., 2, 2] = 1 - 2 * (x * x + y * y)
    return R

def quat_integrate(q, w, dt):
    (qw, qx, qy, qz) = (q[..., 0], q[..., 1], q[..., 2], q[..., 3])
    (vx, vy, vz) = (w[..., 0], w[..., 1], w[..., 2])
    dq = 0.5 * np.stack([-qx * vx - qy * vy - qz * vz, qw * vx + qy * vz - qz * vy, qw * vy - qx * vz + qz * vx, qw * vz + qx * vy - qy * vx], axis=-1)
    out = q + dt * dq
    return out / np.linalg.norm(out, axis=-1, keepdims=True)

def skew(r):
    B = r.shape[:-1]
    S = np.zeros(B + (3, 3))
    S[..., 0, 1] = -r[..., 2]
    S[..., 0, 2] = r[..., 1]
    S[..., 1, 0] = r[..., 2]
    S[..., 1, 2] = -r[..., 0]
    S[..., 2, 0] = -r[..., 1]
    S[..., 2, 1] = r[..., 0]
    return S

def tangents_hand(n):
    """ncp_ref._tangents, batched.  n is the box rear-face outward normal (~ -x)."""
    a = np.where(np.abs(n[..., 0:1]) < 0.9, np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    t1 = np.cross(n, a)
    t1 /= np.linalg.norm(t1, axis=-1, keepdims=True)
    t2 = np.cross(n, t1)
    return (t1, t2)

def inertia_body(m, I_zz, L):
    """diag body inertia; I_xx, I_yy from mass and the known box shape L=(Lx,Ly,Lz),
    I_zz free.  Matches ncp_ref.Box.I_body = m/12 diag(Ly^2+Lz^2, Lx^2+Lz^2, Lx^2+Ly^2)."""
    Ixx = m / 12.0 * (L[..., 1] ** 2 + L[..., 2] ** 2)
    Iyy = m / 12.0 * (L[..., 0] ** 2 + L[..., 2] ** 2)
    return np.stack([Ixx, Iyy, I_zz], axis=-1)

def assemble(Minv_r, nle_r, J_h, hand_pos, v_r, tau, pos, quat, vel, omega, m, I_zz, mu_floor, L, mu_hand, plane_z, dt, eps=0.0, pos_corr=True):
    """Return J (B,15,NV), Minv (B,NV,NV), v_free (B,NV), mu (B,5), b_off (B,15).

    Contact regularization.  Four coplanar corner contacts on a rigid box make the
    Delassus operator rank deficient (15 rows, 9 independent), so lam is not unique and
    two solvers return different corner load splits.  A linear contact compliance eps
    (u = G lam + b + eps lam) makes it unique.  It is realized WITHOUT touching the h1
    Scene contract or any solver: NVIRT extra generalized coordinates with Minv = eps*I
    and J = I on those columns, v_free = 0.  Then J Minv J^T = G_real + eps I exactly,
    b is unchanged and v^+ on the real dofs is unchanged.  Generator and estimator use
    the same eps, so this is a model choice, not a mismatch.
    pos_corr: b_offset = gap/dt for every contact (position-level Signorini), which is
    ncp_ref.build_scene's speculative term extended to negative gaps so the compliant
    contact settles at a fixed penetration instead of sinking."""
    B = pos.shape[0]
    R = quat_to_R(quat)
    Ib = inertia_body(m, I_zz, L)
    Iw_inv = np.einsum('bij,bj,bkj->bik', R, 1.0 / Ib, R)
    Iw = np.einsum('bij,bj,bkj->bik', R, Ib, R)
    Minv = np.zeros((B, NV, NV))
    Minv[:, :NVR, :NVR] = Minv_r
    idx = np.arange(NVR, NVR + 3)
    Minv[:, idx, idx] = (1.0 / m)[:, None]
    Minv[:, NVR + 3:NVREAL, NVR + 3:NVREAL] = Iw_inv
    vi = np.arange(NVREAL, NV)
    Minv[:, vi, vi] = eps
    v_free = np.zeros((B, NV))
    v_free[:, :NVR] = v_r + dt * np.einsum('bij,bj->bi', Minv_r, tau - nle_r)
    v_free[:, NVR:NVR + 3] = vel + dt * G_VEC
    gyro = -np.cross(omega, np.einsum('bij,bj->bi', Iw, omega))
    v_free[:, NVR + 3:NVREAL] = omega + dt * np.einsum('bij,bj->bi', Iw_inv, gyro)
    J = np.zeros((B, MROWS, NV))
    J[:, np.arange(MROWS), NVREAL + np.arange(MROWS)] = 1.0
    b_off = np.zeros((B, MROWS))
    half = 0.5 * L
    n0 = -R[:, :, 0]
    (t1, t2) = tangents_hand(n0)
    Bm0 = np.stack([n0, t1, t2], axis=1)
    r0 = hand_pos - pos
    J[:, 0:3, :NVR] = np.einsum('bij,bjk->bik', Bm0, J_h)
    J[:, 0:3, NVR:NVR + 3] = -Bm0
    J[:, 0:3, NVR + 3:NVREAL] = np.einsum('bij,bjk->bik', Bm0, skew(r0))
    gap0 = np.einsum('bi,bi->b', n0, r0) - half[:, 0]
    if pos_corr:
        b_off[:, 0] = gap0 / dt
    gaps = np.zeros((B, 4))
    for k in range(4):
        rk = np.einsum('bij,bj->bi', R, CORNERS[k] * half)
        row = 3 + 3 * k
        J[:, row:row + 3, NVR:NVR + 3] = B_FLOOR
        J[:, row:row + 3, NVR + 3:NVREAL] = -B_FLOOR @ skew(rk)
        g = pos[:, 2] + rk[:, 2] - plane_z
        gaps[:, k] = g
        if pos_corr:
            b_off[:, row] = g / dt
    mu = np.empty((B, NC))
    mu[:, 0] = mu_hand
    mu[:, 1:] = mu_floor[:, None]
    return (J, Minv, v_free, mu, b_off, gap0, gaps)

def delassus_b(J, Minv, v_free, b_off):
    MJt = np.einsum('bij,bkj->bik', Minv, J)
    G = np.einsum('bij,bjk->bik', J, MJt)
    b = np.einsum('bij,bj->bi', J, v_free) + b_off
    return (G, b)

def apply_impulse(Minv, J, v_free, lam):
    imp = np.einsum('bji,bj->bi', J, lam)
    return v_free + np.einsum('bij,bj->bi', Minv, imp)
