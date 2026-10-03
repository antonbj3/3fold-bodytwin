from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np
from scipy.special import factorial, lpmv


LAYER_NAMES = ("white_matter", "gray_matter", "csf", "skull", "scalp")
MODEL_EQUATIONS = {
    "poisson": "nabla . (sigma(x) grad(phi(x))) = -nabla . (p delta(x-r_s))",
    "current_dipole_discretization": "J = I [delta(x-r_s+d/2) - delta(x-r_s-d/2)]",
    "spherical_harmonic_radial_form": "phi_l(r) = A_l r^l + B_l r^(-(l+1))",
    "interface_conditions": "phi_left=phi_right and sigma_left dphi_left/dr=sigma_right dphi_right/dr",
    "infinity_condition": "sigma dphi/dr + (l+1) sigma phi/r = 0 at the outer air interface",
    "electrode_reference": "u_i = phi(electrode_i) - mean_j(phi(electrode_j))",
    "rdm": "RDM = ||u_test/||u_test||_2 - u_ref/||u_ref||_2||_2",
    "inverse_fit": "argmin_depth ||u_measured/||u_measured||_2 - alpha L_depth||_2",
}


@dataclass(frozen=True)
class HeadGeometry:
    r_wm: float = 0.050
    r_gm: float = 0.080
    csf_thickness: float = 0.003
    skull_thickness: float = 0.007
    scalp_thickness: float = 0.005

    @property
    def r_csf(self) -> float:
        return self.r_gm + self.csf_thickness

    @property
    def r_skull(self) -> float:
        return self.r_csf + self.skull_thickness

    @property
    def r_scalp(self) -> float:
        return self.r_skull + self.scalp_thickness

    @property
    def interfaces(self) -> tuple[float, ...]:
        return (0.0, self.r_wm, self.r_gm, self.r_csf, self.r_skull, self.r_scalp)

    def validate(self) -> None:
        values = np.array(
            [
                self.r_wm,
                self.r_gm,
                self.csf_thickness,
                self.skull_thickness,
                self.scalp_thickness,
            ],
            dtype=float,
        )
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("head geometry must contain finite positive lengths")
        if not self.r_wm < self.r_gm < self.r_csf < self.r_skull < self.r_scalp:
            raise ValueError("head interfaces must be strictly increasing")

    def with_skull_thickness(self, value: float) -> "HeadGeometry":
        value = float(value)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("skull thickness must be positive")
        outer_scalp = self.r_scalp
        new_outer_skull = self.r_csf + value
        new_scalp_thickness = outer_scalp - new_outer_skull
        if new_scalp_thickness <= 0.0:
            raise ValueError("skull perturbation would move the scalp boundary")
        result = replace(
            self,
            skull_thickness=value,
            scalp_thickness=new_scalp_thickness,
        )
        result.validate()
        return result

    def with_csf_thickness(self, value: float) -> "HeadGeometry":
        value = float(value)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("CSF thickness must be positive")
        outer_skull = self.r_skull
        new_r_csf = self.r_gm + value
        new_skull_thickness = outer_skull - new_r_csf
        if new_skull_thickness <= 0.0:
            raise ValueError("CSF perturbation would move the skull boundary")
        result = replace(
            self,
            csf_thickness=value,
            skull_thickness=new_skull_thickness,
        )
        result.validate()
        return result


@dataclass(frozen=True)
class Conductivity:
    white_matter: float = 0.140
    gray_matter: float = 0.330
    csf: float = 1.790
    skull: float = 0.010
    scalp: float = 0.430

    def as_array(self) -> np.ndarray:
        return np.array(
            [self.white_matter, self.gray_matter, self.csf, self.skull, self.scalp],
            dtype=float,
        )

    def validate(self) -> None:
        values = self.as_array()
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("conductivities must be finite positive values in S/m")

    def with_skull(self, value: float) -> "Conductivity":
        result = replace(self, skull=float(value))
        result.validate()
        return result

    def with_csf(self, value: float) -> "Conductivity":
        result = replace(self, csf=float(value))
        result.validate()
        return result

    def with_scalp(self, value: float) -> "Conductivity":
        result = replace(self, scalp=float(value))
        result.validate()
        return result


@dataclass(frozen=True)
class Electrode:
    name: str
    theta: float
    phi: float

    def position(self, radius: float) -> np.ndarray:
        return np.array(
            [
                radius * np.sin(self.theta) * np.cos(self.phi),
                radius * np.sin(self.theta) * np.sin(self.phi),
                radius * np.cos(self.theta),
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class CurrentDipole:
    center: tuple[float, float, float]
    orientation: tuple[float, float, float]
    current: float = 1.0e-9
    separation: float = 0.006

    @property
    def direction(self) -> np.ndarray:
        vector = np.asarray(self.orientation, dtype=float)
        norm = np.linalg.norm(vector)
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError("source orientation must be a finite non-zero vector")
        return vector / norm

    @property
    def moment(self) -> np.ndarray:
        return self.current * self.separation * self.direction

    def endpoints(self) -> tuple[np.ndarray, np.ndarray]:
        center = np.asarray(self.center, dtype=float)
        delta = 0.5 * self.separation * self.direction
        return center + delta, center - delta

    def validate(self) -> None:
        center = np.asarray(self.center, dtype=float)
        if center.shape != (3,) or not np.all(np.isfinite(center)):
            raise ValueError("source center must be a finite three-vector")
        if not np.isfinite(self.current) or self.current <= 0.0:
            raise ValueError("source current must be positive in A")
        if not np.isfinite(self.separation) or self.separation <= 0.0:
            raise ValueError("source separation must be positive in m")
        self.direction
        for point in self.endpoints():
            if not np.all(np.isfinite(point)):
                raise ValueError("source endpoints must be finite")


@dataclass(frozen=True)
class ForwardResult:
    raw_potential_v: np.ndarray
    referenced_potential_v: np.ndarray
    reference_potential_v: float
    electrode_names: tuple[str, ...]
    electrode_positions_m: np.ndarray

    @property
    def topography_norm_v(self) -> float:
        return float(np.linalg.norm(self.referenced_potential_v))

    def finite(self) -> bool:
        return bool(
            np.all(np.isfinite(self.raw_potential_v))
            and np.all(np.isfinite(self.referenced_potential_v))
            and np.isfinite(self.reference_potential_v)
        )


def default_electrodes() -> tuple[Electrode, ...]:
    specs = (
        ("Fp1", 72.0, 0.0),
        ("Fp2", 72.0, 180.0),
        ("F7", 54.0, -60.0),
        ("F8", 54.0, 60.0),
        ("F3", 36.0, -30.0),
        ("F4", 36.0, 30.0),
        ("Fz", 18.0, 0.0),
        ("T3", 90.0, -90.0),
        ("T4", 90.0, 90.0),
        ("C3", 72.0, -30.0),
        ("C4", 72.0, 30.0),
        ("Cz", 36.0, 0.0),
        ("T7", 108.0, -90.0),
        ("T8", 108.0, 90.0),
        ("P3", 108.0, -30.0),
        ("P4", 108.0, 30.0),
        ("Pz", 144.0, 0.0),
        ("O1", 162.0, 0.0),
        ("O2", 162.0, 180.0),
        ("M1", 90.0, -135.0),
        ("M2", 90.0, 135.0),
    )
    return tuple(
        Electrode(name, np.deg2rad(theta), np.deg2rad(phi))
        for name, theta, phi in specs
    )


def default_source(
    depth_m: float = 0.020,
    geometry: HeadGeometry | None = None,
) -> CurrentDipole:
    geometry = geometry or HeadGeometry()
    return CurrentDipole(
        center=(0.0, 0.0, geometry.r_csf - depth_m),
        orientation=(1.0, 0.0, 0.0),
    )


def _real_harmonic_values(ell: int, theta: float, phi: float) -> np.ndarray:
    x = float(np.cos(theta))
    normalization = np.sqrt((2.0 * ell + 1.0) / (4.0 * np.pi))
    if ell == 0:
        return np.array([normalization], dtype=float)
    values = [normalization * float(lpmv(0, ell, x))]
    for m in range(1, ell + 1):
        n = normalization * np.sqrt(
            float(factorial(ell - m) / factorial(ell + m))
        )
        p = float(lpmv(m, ell, x))
        if m % 2:
            p = -p
        scale = np.sqrt(2.0) * n * p
        values.extend((scale * np.cos(m * phi), scale * np.sin(m * phi)))
    return np.asarray(values, dtype=float)


def _angular_dot(
    ell: int,
    source_spherical: tuple[float, float],
    target_spherical: tuple[float, float],
) -> float:
    source_values = _real_harmonic_values(ell, source_spherical[0], source_spherical[1])
    target_values = _real_harmonic_values(ell, target_spherical[0], target_spherical[1])
    return float(np.dot(source_values, target_values))


def _spherical(position: Sequence[float]) -> tuple[float, float, float]:
    vector = np.asarray(position, dtype=float)
    radius = float(np.linalg.norm(vector))
    if radius <= 1.0e-15:
        raise ValueError("spherical source must be away from the center")
    theta = float(np.arccos(np.clip(vector[2] / radius, -1.0, 1.0)))
    phi = float(np.arctan2(vector[1], vector[0]) % (2.0 * np.pi))
    return radius, theta, phi


@dataclass
class _Segment:
    layer: int
    sigma: float
    lower: float
    upper: float
    a_index: int
    b_index: int | None


class SphericalHeadModel:
    def __init__(
        self,
        geometry: HeadGeometry = HeadGeometry(),
        conductivity: Conductivity = Conductivity(),
        lmax: int = 28,
        air_sigma: float = 1.0e-6,
    ) -> None:
        geometry.validate()
        conductivity.validate()
        if not np.isfinite(air_sigma) or air_sigma <= 0.0:
            raise ValueError("air conductivity must be positive in S/m")
        if not isinstance(lmax, (int, np.integer)) or lmax < 1:
            raise ValueError("lmax must be a positive integer")
        self.geometry = geometry
        self.conductivity = conductivity
        self.lmax = int(lmax)
        self.air_sigma = float(air_sigma)
        self.interfaces = geometry.interfaces
        self.layer_sigma = conductivity.as_array()
        self._radial_cache: dict[tuple[int, int, float], tuple[list[_Segment], np.ndarray]] = {}

    def _source_layer(self, radius: float) -> int:
        if radius >= self.geometry.r_scalp:
            raise ValueError("source must be inside the scalp boundary")
        boundaries = self.interfaces[1:-1]
        for boundary in boundaries:
            if np.isclose(radius, boundary, rtol=0.0, atol=1.0e-12):
                raise ValueError("source cannot lie on a spherical interface")
        return int(np.searchsorted(boundaries, radius, side="right"))

    def _radial_solution(
        self,
        ell: int,
        source_layer: int,
        source_radius: float,
    ) -> tuple[list[_Segment], np.ndarray]:
        key = (ell, source_layer, round(source_radius, 14))
        if key in self._radial_cache:
            return self._radial_cache[key]
        segments: list[_Segment] = []
        unknown_count = 0
        for layer in range(len(self.layer_sigma)):
            if layer < source_layer:
                has_inverse_term = self.interfaces[layer] > 0.0
                b_index = unknown_count + 1 if has_inverse_term else None
                segments.append(
                    _Segment(
                        layer=layer,
                        sigma=float(self.layer_sigma[layer]),
                        lower=self.interfaces[layer],
                        upper=self.interfaces[layer + 1],
                        a_index=unknown_count,
                        b_index=b_index,
                    )
                )
                unknown_count += 2 if has_inverse_term else 1
            elif layer == source_layer:
                has_inverse_term = self.interfaces[layer] > 0.0
                b_index = unknown_count + 1 if has_inverse_term else None
                segments.append(
                    _Segment(
                        layer=layer,
                        sigma=float(self.layer_sigma[layer]),
                        lower=self.interfaces[layer],
                        upper=source_radius,
                        a_index=unknown_count,
                        b_index=b_index,
                    )
                )
                unknown_count += 2 if has_inverse_term else 1
                segments.append(
                    _Segment(
                        layer=layer,
                        sigma=float(self.layer_sigma[layer]),
                        lower=source_radius,
                        upper=self.interfaces[layer + 1],
                        a_index=unknown_count,
                        b_index=unknown_count + 1,
                    )
                )
                unknown_count += 2
            else:
                segments.append(
                    _Segment(
                        layer=layer,
                        sigma=float(self.layer_sigma[layer]),
                        lower=self.interfaces[layer],
                        upper=self.interfaces[layer + 1],
                        a_index=unknown_count,
                        b_index=unknown_count + 1,
                    )
                )
                unknown_count += 2
        matrix = np.zeros((unknown_count, unknown_count), dtype=float)
        rhs = np.zeros(unknown_count, dtype=float)
        source_index = source_layer
        inner = segments[source_index]
        outer = segments[source_index + 1]
        rows: list[np.ndarray] = []
        values: list[float] = []

        def state(segment: _Segment, radius: float) -> tuple[np.ndarray, np.ndarray]:
            value_row = np.zeros(unknown_count, dtype=float)
            derivative_row = np.zeros(unknown_count, dtype=float)
            if radius == 0.0:
                if ell == 0:
                    value_row[segment.a_index] = 1.0
                return value_row, segment.sigma * derivative_row
            power = radius**ell
            inverse = radius ** (-(ell + 1))
            value_row[segment.a_index] += power
            if ell > 0:
                derivative_row[segment.a_index] += ell * radius ** (ell - 1)
            if segment.b_index is not None:
                value_row[segment.b_index] += inverse
                derivative_row[segment.b_index] += -(ell + 1.0) * inverse / radius
            return value_row, segment.sigma * derivative_row

        source_u_in, source_q_in = state(inner, source_radius)
        source_u_out, source_q_out = state(outer, source_radius)
        rows.append(source_u_out - source_u_in)
        values.append(0.0)
        rows.append(source_q_out - source_q_in)
        values.append(-1.0 / (source_radius**2))
        for left, right in zip(segments[:-1], segments[1:]):
            if left is inner and right is outer:
                continue
            boundary = right.lower
            left_u, left_q = state(left, boundary)
            right_u, right_q = state(right, boundary)
            rows.append(right_u - left_u)
            values.append(0.0)
            rows.append(right_q - left_q)
            values.append(0.0)
        last = segments[-1]
        outer_u, outer_q = state(last, self.geometry.r_scalp)
        robin = outer_q + (ell + 1.0) * self.air_sigma * outer_u / self.geometry.r_scalp
        rows.append(robin)
        values.append(0.0)
        if len(rows) != unknown_count:
            raise RuntimeError("radial boundary system has the wrong number of equations")
        matrix = np.asarray(rows, dtype=float)
        rhs = np.asarray(values, dtype=float)
        coefficients = np.linalg.solve(matrix, rhs)
        result = (segments, coefficients)
        self._radial_cache[key] = result
        return result

    @staticmethod
    def _radial_value(
        segments: list[_Segment],
        coefficients: np.ndarray,
        ell: int,
        radius: float,
    ) -> float:
        selected = None
        for segment in segments:
            if radius >= segment.lower - 1.0e-13 and radius <= segment.upper + 1.0e-13:
                selected = segment
                break
        if selected is None:
            raise ValueError("radial evaluation point is outside the modeled head")
        if radius == 0.0:
            return float(coefficients[selected.a_index] if ell == 0 else 0.0)
        value = coefficients[selected.a_index] * radius**ell
        if selected.b_index is not None:
            value += coefficients[selected.b_index] * radius ** (-(ell + 1))
        return float(value)

    def green_monomial(
        self,
        position: Sequence[float],
        current: float,
        electrode_positions_m: np.ndarray,
    ) -> np.ndarray:
        if not np.isfinite(current) or current <= 0.0:
            raise ValueError("monopole current must be positive in A")
        source_spherical = _spherical(position)
        source_radius = source_spherical[0]
        source_layer = self._source_layer(source_radius)
        targets = [_spherical(point) for point in electrode_positions_m]
        result = np.zeros(len(targets), dtype=float)
        for ell in range(self.lmax + 1):
            segments, coefficients = self._radial_solution(ell, source_layer, source_radius)
            angular_values = np.asarray(
                [
                    _angular_dot(ell, (source_spherical[1], source_spherical[2]), (theta, phi))
                    for _, theta, phi in targets
                ],
                dtype=float,
            )
            radial_values = np.asarray(
                [
                    self._radial_value(segments, coefficients, ell, target[0])
                    for target in targets
                ],
                dtype=float,
            )
            result += current * radial_values * angular_values
        return result

    def forward(
        self,
        source: CurrentDipole,
        electrodes: Sequence[Electrode] = default_electrodes(),
        reference: str = "average",
    ) -> ForwardResult:
        source.validate()
        positions = np.asarray(
            [electrode.position(self.geometry.r_scalp) for electrode in electrodes],
            dtype=float,
        )
        first, second = source.endpoints()
        raw = self.green_monomial(first, source.current, positions)
        raw -= self.green_monomial(second, source.current, positions)
        names = tuple(electrode.name for electrode in electrodes)
        normalized_reference = reference.lower().replace("-", "_")
        if normalized_reference in {"average", "average_reference"}:
            reference_value = float(np.mean(raw))
        elif normalized_reference in {"mastoid", "linked_mastoid", "linked_mastoids"}:
            try:
                left = names.index("M1")
                right = names.index("M2")
            except ValueError as error:
                raise ValueError("mastoid reference requires M1 and M2 electrodes") from error
            reference_value = float(0.5 * (raw[left] + raw[right]))
        else:
            raise ValueError(f"unknown reference: {reference}")
        referenced = raw - reference_value
        result = ForwardResult(
            raw_potential_v=raw,
            referenced_potential_v=referenced,
            reference_potential_v=reference_value,
            electrode_names=names,
            electrode_positions_m=positions,
        )
        if not result.finite():
            raise FloatingPointError("non-finite EEG forward result")
        return result

    def leadfields(
        self,
        source_template: CurrentDipole,
        depths_m: Sequence[float],
        electrodes: Sequence[Electrode] = default_electrodes(),
        reference: str = "average",
    ) -> tuple[np.ndarray, list[ForwardResult]]:
        leadfields: list[np.ndarray] = []
        results: list[ForwardResult] = []
        for depth in depths_m:
            source = replace(
                source_template,
                center=(0.0, 0.0, self.geometry.r_csf - float(depth)),
            )
            result = self.forward(source, electrodes=electrodes, reference=reference)
            results.append(result)
            norm = np.linalg.norm(result.referenced_potential_v)
            if norm <= 0.0 or not np.isfinite(norm):
                raise FloatingPointError("zero or non-finite leadfield")
            leadfields.append(result.referenced_potential_v / norm)
        return np.asarray(leadfields), results


def analytic_dipole_potential(
    source: CurrentDipole,
    electrode_positions_m: np.ndarray,
    sigma: float,
) -> np.ndarray:
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("homogeneous conductivity must be positive in S/m")
    first, second = source.endpoints()
    first_term = source.current / (4.0 * np.pi * sigma)
    second_term = source.current / (4.0 * np.pi * sigma)
    return first_term / np.linalg.norm(electrode_positions_m - first, axis=1) - second_term / np.linalg.norm(
        electrode_positions_m - second, axis=1
    )


def relative_metrics(reference: ForwardResult, test: ForwardResult) -> dict[str, float]:
    reference_vector = reference.referenced_potential_v
    test_vector = test.referenced_potential_v
    reference_norm = np.linalg.norm(reference_vector)
    test_norm = np.linalg.norm(test_vector)
    if reference_norm <= 0.0 or test_norm <= 0.0:
        raise ValueError("relative metrics require non-zero referenced vectors")
    rdm = float(
        np.linalg.norm(test_vector / test_norm - reference_vector / reference_norm)
    )
    return {
        "rdm": rdm,
        "mag": float(test_norm / reference_norm),
        "max_abs_delta_v": float(np.max(np.abs(test_vector - reference_vector))),
    }


def fit_depth(
    measured: ForwardResult,
    leadfields: np.ndarray,
    depths_m: Sequence[float],
) -> dict[str, float | int | list[float]]:
    measured_vector = measured.referenced_potential_v
    measured_norm = np.linalg.norm(measured_vector)
    if measured_norm <= 0.0:
        raise ValueError("cannot fit a zero measured vector")
    scores: list[float] = []
    amplitudes: list[float] = []
    for leadfield in leadfields:
        amplitude = float(np.dot(leadfield, measured_vector) / measured_norm)
        residual = measured_vector / measured_norm - amplitude * leadfield
        scores.append(float(np.linalg.norm(residual)))
        amplitudes.append(amplitude)
    index = int(np.argmin(scores))
    return {
        "index": index,
        "depth_m": float(depths_m[index]),
        "score": float(scores[index]),
        "amplitudes": [float(value) for value in amplitudes],
    }


def unit_check() -> dict[str, float | str | bool]:
    geometry = HeadGeometry()
    conductivity = Conductivity()
    source = default_source(geometry=geometry)
    electrodes = default_electrodes()
    positions = np.asarray(
        [electrode.position(geometry.r_scalp) for electrode in electrodes],
        dtype=float,
    )
    homogeneous = Conductivity(
        white_matter=0.2,
        gray_matter=0.2,
        csf=0.2,
        skull=0.2,
        scalp=0.2,
    )
    analytic = analytic_dipole_potential(source, positions, homogeneous.white_matter)
    return {
        "length_unit": "m",
        "conductivity_unit": "S/m",
        "source_current_unit": "A",
        "source_separation_unit": "m",
        "source_moment": float(np.linalg.norm(source.moment)),
        "source_moment_unit": "A m",
        "potential_unit": "V",
        "air_conductivity": 1.0e-6,
        "air_conductivity_unit": "S/m",
        "homogeneous_peak_potential_v": float(np.max(np.abs(analytic))),
        "all_dimensions_finite": bool(
            np.all(np.isfinite(positions))
            and np.all(np.isfinite(analytic))
            and np.isfinite(np.linalg.norm(source.moment))
        ),
    }


def parameter_table() -> list[dict[str, object]]:
    geometry = HeadGeometry()
    conductivity = Conductivity()
    source = default_source(geometry=geometry)
    rows: list[dict[str, object]] = [
        {
            "parameter": "r_wm",
            "value_m": geometry.r_wm,
            "unit": "m",
            "source_or_assumption": "spherical first-run assumption",
        },
        {
            "parameter": "r_gm",
            "value_m": geometry.r_gm,
            "unit": "m",
            "source_or_assumption": "spherical first-run assumption",
        },
        {
            "parameter": "csf_thickness",
            "value_m": geometry.csf_thickness,
            "unit": "m",
            "source_or_assumption": "spherical first-run assumption",
        },
        {
            "parameter": "skull_thickness",
            "value_m": geometry.skull_thickness,
            "unit": "m",
            "source_or_assumption": "spherical first-run assumption",
        },
        {
            "parameter": "scalp_thickness",
            "value_m": geometry.scalp_thickness,
            "unit": "m",
            "source_or_assumption": "spherical first-run assumption",
        },
    ]
    rows.append(
        {
            "parameter": "sigma_air",
            "value": 1.0e-6,
            "unit": "S/m",
            "source_or_assumption": "quasi-static exterior-air assumption; positive numerical value",
        }
    )
    for name, value in zip(LAYER_NAMES, conductivity.as_array()):
        rows.append(
            {
                "parameter": f"sigma_{name}",
                "value": float(value),
                "unit": "S/m",
                "source_or_assumption": "Vorwerk et al. 2024 Table 1 standard, DOI 10.3389/fnhum.2024.1335212",
            }
        )
    rows.extend(
        [
            {
                "parameter": "source_current",
                "value": source.current,
                "unit": "A",
                "source_or_assumption": "frozen synthetic finite-difference dipole",
            },
            {
                "parameter": "source_separation",
                "value": source.separation,
                "unit": "m",
                "source_or_assumption": "frozen synthetic finite-difference dipole",
            },
            {
                "parameter": "source_depth_below_inner_skull",
                "value": 0.020,
                "unit": "m",
                "source_or_assumption": "frozen preregistration",
            },
            {
                "parameter": "reference",
                "value": "average",
                "unit": "categorical",
                "source_or_assumption": "frozen preregistration; linked mastoid alternative implemented",
            },
        ]
    )
    return rows


if __name__ == "__main__":
    check = unit_check()
    baseline = SphericalHeadModel().forward(default_source())
    print(check)
    print(
        {
            "reference_v": baseline.reference_potential_v,
            "peak_abs_referenced_v": float(np.max(np.abs(baseline.referenced_potential_v))),
            "topography_norm_v": baseline.topography_norm_v,
        }
    )
