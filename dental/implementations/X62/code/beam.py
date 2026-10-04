"""Conventional 3D Euler--Bernoulli frame. mm, N, MPa. No biological law."""
import numpy as np
from scipy.linalg import solve
from scipy.integrate import quad
from fractions import Fraction

def moments(w=1.5, t=1.0):
    a = w / 2
    b = t
    area = np.pi * a * b / 2
    cz = 4 * b / (3 * np.pi)
    Iy = a * b ** 3 * (np.pi / 8 - 8 / (9 * np.pi))
    Iz = np.pi * a ** 3 * b / 8
    return (area, Iy, Iz, cz)

def integrated_moments(w=1.5, t=1.0):
    a = w / 2
    b = t
    A = quad(lambda z: 2 * a * np.sqrt(max(0, 1 - (z / b) ** 2)), 0, b, epsabs=1e-13)[0]
    cz = quad(lambda z: z * 2 * a * np.sqrt(max(0, 1 - (z / b) ** 2)), 0, b, epsabs=1e-13)[0] / A
    Iy = quad(lambda z: (z - cz) ** 2 * 2 * a * np.sqrt(max(0, 1 - (z / b) ** 2)), 0, b, epsabs=1e-13)[0]
    Iz = quad(lambda z: 2 / 3 * (a * np.sqrt(max(0, 1 - (z / b) ** 2))) ** 3, 0, b, epsabs=1e-13)[0]
    return (A, Iy, Iz, cz)

def stiffness(points, vertical, E=4800.0, nu=0.4, EI=None):
    (A, Iy, Iz, _) = moments()
    (EIy, EIz) = (E * Iy, E * Iz) if EI is None else EI
    GJ = E / (2 * (1 + nu)) * (Iy + Iz)
    K = np.zeros((6 * len(points), 6 * len(points)))
    lengths = []
    for (i, (p, q)) in enumerate(zip(points[:-1], points[1:])):
        L = np.linalg.norm(q - p)
        x = (q - p) / L
        y = vertical - x * np.dot(vertical, x)
        y /= np.linalg.norm(y)
        z = np.cross(x, y)
        R = np.array([x, y, z])
        T = np.zeros((12, 12))
        for j in range(4):
            T[j * 3:j * 3 + 3, j * 3:j * 3 + 3] = R
        k = np.zeros((12, 12))
        for (a, b, s) in [(0, 6, E * A / L), (3, 9, GJ / L)]:
            k[np.ix_([a, b], [a, b])] += s * np.array([[1, -1], [-1, 1]])
        B = lambda ei: ei / L ** 3 * np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L * L, -6 * L, 2 * L * L], [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L * L, -6 * L, 4 * L * L]])
        idx = [1, 5, 7, 11]
        k[np.ix_(idx, idx)] += B(EIz)
        idx = [2, 4, 8, 10]
        signs = np.diag([1, -1, 1, -1])
        k[np.ix_(idx, idx)] += signs @ B(EIy) @ signs
        ids = np.arange(6 * i, 6 * i + 12)
        K[np.ix_(ids, ids)] += T.T @ k @ T
        lengths.append(float(L))
    return (K, np.array(lengths))

def condense(K):
    n = K.shape[0] // 6
    u = np.array([6 * i + j for i in range(n) for j in range(3)])
    r = np.array([6 * i + j for i in range(n) for j in range(3, 6)])
    S = K[np.ix_(u, u)] - K[np.ix_(u, r)] @ solve(K[np.ix_(r, r)], K[np.ix_(r, u)], assume_a='pos')
    return ((S + S.T) / 2, u, r)

def direct(K, displacement):
    (_, u, r) = condense(K)
    q = np.zeros(len(K))
    q[u] = displacement.reshape(-1)
    q[r] = solve(K[np.ix_(r, r)], -K[np.ix_(r, u)] @ q[u], assume_a='pos')
    return ((K @ q)[u].reshape(-1, 3), q)

def affine_box(S, center, radius):
    """Exact-real matrix enclosure, rounded outward; floating S error is NOT certified."""
    c = S @ np.ravel(center)
    rad = np.abs(S) @ np.ravel(radius)
    return (np.nextafter(c - rad, -np.inf), np.nextafter(c + rad, np.inf))

def exact_binary_box(S, center, radius):
    """Rigorous enclosure of the affine map whose coefficients are frozen binary floats.

    This covers evaluation roundoff. It does not enclose error in building S, nor
    the physical constitutive remainder. Those are separately UNKNOWN.
    """
    u = [Fraction(float(x)) for x in np.ravel(center)]
    r = [Fraction(float(x)) for x in np.ravel(radius)]
    lows = []
    highs = []

    def out(x, up):
        f = float(x)
        if up and Fraction(f) < x or (not up and Fraction(f) > x):
            f = float(np.nextafter(f, np.inf if up else -np.inf))
        return f
    for row in S:
        a = [Fraction(float(x)) for x in row]
        c = sum((x * y for (x, y) in zip(a, u)))
        d = sum((abs(x) * y for (x, y) in zip(a, r)))
        lows.append(out(c - d, False))
        highs.append(out(c + d, True))
    return (np.array(lows), np.array(highs))
