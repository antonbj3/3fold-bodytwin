"""Empirical femur priors from Imperial and VSD; dimensions are mm and degrees.

These are initialization distributions, not measurements of an individual. The
intervals are source-cohort leave-person-out residual calibrations. A returned
``in_range`` of False means extrapolation beyond the training covariates.
"""
from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class ScalarPrior:
    mean: float
    lower: float
    upper: float
    unit: str
    source: str
    in_range: bool


def _finite(**values):
    if not all(isfinite(float(v)) for v in values.values()):
        raise ValueError("prior inputs must be finite")


def prior_from_anthropometry(height_m, weight_kg, sex):
    """VSD femur length, 26 complete persons; 90% LOSO residual interval.

    ``sex`` is the source-table F/M code. The outcome is the mean of both
    femur lengths for a person, so this is not a side-specific measurement.
    """
    _finite(height_m=height_m, weight_kg=weight_kg)
    if sex not in ('F', 'M'):
        raise ValueError("sex must be the source-table code 'F' or 'M'")
    center = (-28.70531864813394 + 287.58649078123136*height_m
              - 0.43002443915241884*weight_kg + 19.437650371069545*(sex == 'M'))
    radius = 26.907281971921236
    valid = 1.52 <= height_m <= 1.87 and 37.4 <= weight_kg <= 90.0
    return ScalarPrior(center, center-radius, center+radius, 'mm', 'VSD DA-034; 90% LOSO', valid)


def contralateral(observed, side, parameter='head_radius_mm'):
    """VSD bilateral location and 95% residual interval for radius or CCD.

    ``side`` is the observed side. Calibration uses 30 paired VSD persons,
    including flagged rows; head-radius fit flags can make uncertainty large.
    """
    _finite(observed=observed)
    if side not in ('L', 'R'):
        raise ValueError("side must be 'L' or 'R'")
    specs = {
        'head_radius_mm': (-0.18970055021132914, 3.2350690346634945,
                           (19.05774484680095, 27.144961931926858), 'mm'),
        'ccd_deg': (-0.16876652337992226, 9.960779194746657,
                    (122.39082596063201, 141.27498466197173), 'deg'),
    }
    if parameter not in specs:
        raise ValueError("supported parameters: head_radius_mm, ccd_deg")
    delta, radius, bounds, unit = specs[parameter]
    center = observed + (delta if side == 'L' else -delta)
    return ScalarPrior(center, center-radius, center+radius, unit,
                       'VSD DA-035/DA-029; 95% paired residual', bounds[0] <= observed <= bounds[1])


def prior_from_tibia(tibial_length_mm, plateau_width_mm):
    """Imperial radius from two tibial measures; 95% LOSO residual interval.

    This simpler two-input model has LOSO RMSE 1.042 mm. DA-030's incremental
    0.741-mm model additionally requires 28 femur features and is not this API.
    """
    _finite(tibial_length_mm=tibial_length_mm, plateau_width_mm=plateau_width_mm)
    center = (-1.9228814746744358 + 0.015740061858935603*tibial_length_mm
              + 0.2847346970357698*plateau_width_mm)
    radius = 2.222820148100059
    valid = 297.0508167703634 <= tibial_length_mm <= 384.6587628975647 and 52.32292714357307 <= plateau_width_mm <= 78.21117348498284
    return ScalarPrior(center, center-radius, center+radius, 'mm', 'Imperial DA-030; 95% LOSO', valid)


def scalar_population_prior():
    """Imperial Gaussian for (radius mm, CCD deg, signed anteversion deg).

    The 95% ellipsoid is (x-mean)' covariance^-1 (x-mean) <= 11.432;
    this covered 51/54 unflagged VSD femurs. Off-diagonal gains were weak.
    """
    return {
        'parameters': ('head_radius_mm', 'ccd_deg', 'anteversion_signed_deg'),
        'mean': (21.606044855861057, 128.77570853852862, -0.847172182642586),
        'covariance': ((4.787193502475404, 2.7029712748462504, -0.8134765757657664),
                       (2.7029712748462504, 17.282834367296918, -5.317804395435018),
                       (-0.8134765757657664, -5.317804395435018, 416.1040169243852)),
        'mahalanobis_squared_95': 11.432164672207174,
        'imperial_observed_ranges': ((17.79043933342422, 26.779072841954545),
                                     (120.2942205355748, 139.43739749633758),
                                     (-38.94477243142709, 38.69155570533006)),
        'source': 'Imperial DA-033; threshold calibrated on Imperial subject folds',
    }
