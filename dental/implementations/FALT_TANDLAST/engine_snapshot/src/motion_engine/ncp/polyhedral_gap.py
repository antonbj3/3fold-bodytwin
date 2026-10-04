"""Complete finite face cover for a rational 3D polygon friction law.

Includes edge slip: vertex-only enumeration is incomplete. Cone nesting
does not imply nesting of coupled dynamic solution sets. No circle claim.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from .gap_box import GapProblem, dot, rational_json
from .observable_cert import _q
DIAMOND = ((1, 0), (0, 1), (-1, 0), (0, -1))
SQUARE = ((1, 1), (-1, 1), (-1, -1), (1, -1))

@dataclass(frozen=True)
class PolyhedralGapProblem(GapProblem):
    polygon: tuple

    @classmethod
    def make(cls, J, inverse_mass, free_velocity, gaps, step, *, mu, compliance=None, initial_velocity=None, polygon=DIAMOND):
        q = tuple((tuple((_q(t) for t in v)) for v in polygon))
        if len(q) < 3 or any((len(v) != 2 for v in q)):
            raise ValueError('invalid polygon')
        faces = []
        for (a, b) in zip(q, q[1:] + q[:1]):
            normal = (b[1] - a[1], a[0] - b[0])
            rho = dot(normal, a)
            if rho <= 0 or any((dot(normal, v) > rho for v in q)):
                raise ValueError('polygon must be convex CCW with origin interior')
            other_vertices = [v for v in q if v not in (a, b)]
            if not other_vertices:
                raise ValueError('at least three distinct polygon vertices required')
            if any((dot(normal, v) == rho for v in other_vertices)):
                raise ValueError('redundant or collinear polygon vertex')
            faces.append((normal, rho))
        n = len(gaps)
        if len(J) != 3 * n:
            raise ValueError('3D stride')
        normal = GapProblem.make(J[::3], inverse_mass, free_velocity, gaps, step, compliance=compliance, initial_velocity=initial_velocity)
        full = tuple((tuple((_q(t) for t in row)) for row in J))
        mus = tuple((_q(t) for t in mu))
        if any((len(row) != len(normal.inverse_mass) for row in full)) or len(mus) != n or any((t < 0 for t in mus)):
            raise ValueError('dimensions or friction')
        return cls(full, normal.inverse_mass, normal.free_velocity, normal.initial_velocity, normal.gaps, normal.step, mus, normal.compliance, 'polyhedral3d', q)

    @property
    def stride(self):
        return 3

    @property
    def modes(self):
        return ('open', 'stick') + tuple((f'vertex:{i}' for i in range(len(self.polygon)))) + tuple((f'edge:{i}' for i in range(len(self.polygon))))

    def system(self, prefix):
        n = len(self.gaps)
        m = len(self.J)
        nd = m + n
        if len(prefix) > n or any((p not in self.modes for p in prefix)):
            raise ValueError('invalid prefix')
        A = [[sum((self.J[i][k] * self.inverse_mass[k] * self.J[j][k] for k in range(len(self.inverse_mass)))) for j in range(m)] for i in range(m)]
        for (i, d) in enumerate(self.compliance):
            A[3 * i][3 * i] += d
        b = [dot(row, self.free_velocity) for row in self.J]
        H = []
        r = []
        E = []
        e = []

        def unit(i, s=Q(1)):
            row = [Q(0)] * nd
            row[i] = s
            return row

        def le(row, rhs=Q(0)):
            H.append(row)
            r.append(rhs)

        def eq(row, rhs=Q(0)):
            E.append(row)
            e.append(rhs)
        for (i, (lo, hi)) in enumerate(self.gaps):
            ni = 3 * i
            ti = ni + 1
            mu = self.mu[i]
            le(unit(m + i), hi)
            le(unit(m + i, -1), -lo)
            le(unit(ni, -1))
            wn = A[ni][:] + [Q(0)] * n
            wn[m + i] = 1 / self.step
            le([-t for t in wn], b[ni])
            ut = [A[ti][:] + [Q(0)] * n, A[ti + 1][:] + [Q(0)] * n]
            for (a, z) in zip(self.polygon, self.polygon[1:] + self.polygon[:1]):
                normal = (z[1] - a[1], a[0] - z[0])
                rho = dot(normal, a)
                row = unit(ti, normal[0])
                row[ti + 1] = normal[1]
                row[ni] = -mu * rho
                le(row)
            if i >= len(prefix):
                continue
            mode = prefix[i]
            if mode == 'open':
                for j in [ni, ti, ti + 1]:
                    eq(unit(j))
                continue
            eq(wn, -b[ni])
            if mode == 'stick':
                for j in range(2):
                    eq(ut[j], -b[ti + j])
            elif mode.startswith('vertex:'):
                a = self.polygon[int(mode.split(':')[1])]
                for j in range(2):
                    row = unit(ti + j)
                    row[ni] = -mu * a[j]
                    eq(row)
                for z in self.polygon:
                    c = (a[0] - z[0], a[1] - z[1])
                    le([c[0] * x + c[1] * y for (x, y) in zip(*ut)], -dot(c, b[ti:ti + 2]))
            else:
                k = int(mode.split(':')[1])
                a = self.polygon[k]
                z = self.polygon[(k + 1) % len(self.polygon)]
                normal = (z[1] - a[1], a[0] - z[0])
                rho = dot(normal, a)
                row = unit(ti, normal[0])
                row[ti + 1] = normal[1]
                row[ni] = -mu * rho
                eq(row)
                c = (normal[1], -normal[0])
                eq([c[0] * x + c[1] * y for (x, y) in zip(*ut)], -dot(c, b[ti:ti + 2]))
                le([normal[0] * x + normal[1] * y for (x, y) in zip(*ut)], -dot(normal, b[ti:ti + 2]))
        return (H, r, E, e)

def cone_sandwich():
    """Exact diamond subset disk subset square, with radial distortion bound."""
    return dict(inner=DIAMOND, outer=SQUARE, inner_vertex_norm_squared=[dot(v, v) for v in DIAMOND], outer_contains_disk='abs(t_x)<=r and abs(t_y)<=r whenever t_x²+t_y²<=r²', inner_contains_disk_radius_squared=Q(1, 2), outer_radius_squared=Q(2), dynamic_solution_sandwich=False)

def solve_polyhedral_gap(J, inverse_mass, free_velocity, gaps, step, *, mu, compliance=None, initial_velocity=None, polygon=DIAMOND, max_nodes=2048, exhaustive=False):
    """Fail closed before converting large sparse scenes or creating modes."""
    from .gap_box import enclose_gap_box
    try:
        row_count = J.shape[0] if hasattr(J, 'shape') else len(J)
        if row_count + len(gaps) > 128:
            return dict(status='UNCERTAIN', complete=False, reason='EXACT_ARITHMETIC_DIMENSION_BUDGET', bounds={}, records=[])
        if hasattr(J, 'toarray'):
            J = J.toarray().tolist()
        pb = PolyhedralGapProblem.make(J, inverse_mass, free_velocity, gaps, step, mu=mu, compliance=compliance, initial_velocity=initial_velocity, polygon=polygon)
        return enclose_gap_box(pb, max_nodes=max_nodes, exhaustive=exhaustive)
    except (ValueError, TypeError, ArithmeticError) as exc:
        return dict(status='UNCERTAIN', complete=False, reason=str(exc), bounds={}, records=[])
