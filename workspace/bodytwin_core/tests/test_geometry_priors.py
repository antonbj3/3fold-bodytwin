"""Reproduce a held-out source number and check the public prior API."""
import csv
from pathlib import Path

import numpy as np
import pytest

from bodytwin_core import geometry


def test_anthropometry_person_held_out_mae():
    source = Path(__file__).resolve().parents[2] / 'results/BT-DA-034/inputs/VSD_FEMUR.csv'
    with source.open() as f:
        rows = list(csv.DictReader(f))
    people = {}
    for row in rows:
        people.setdefault(row['subject'], {})[row['side']] = row
    eligible = [s for s in people.values() if len(s) == 2 and s['L']['height_m'] and s['L']['weight_kg']]
    x = np.array([[float(s['L']['height_m']), float(s['L']['weight_kg']), s['L']['sex'] == 'M'] for s in eligible])
    y = np.array([np.mean([float(s[k]['femur_length_mm']) for k in ('L', 'R')]) for s in eligible])
    predictions = []
    for i in range(len(y)):
        fit = np.arange(len(y)) != i
        coef = np.linalg.lstsq(np.c_[np.ones(fit.sum()), x[fit]], y[fit], rcond=None)[0]
        predictions.append(np.r_[1, x[i]] @ coef)
    assert np.mean(abs(y-predictions)) == pytest.approx(9.7776870266, abs=1e-8)
    prior = geometry.prior_from_anthropometry(1.7, 70, 'F')
    assert prior.in_range and prior.lower < prior.mean < prior.upper
    assert prior.upper-prior.lower == pytest.approx(53.8145639438, abs=1e-8)
    assert not geometry.prior_from_anthropometry(1.95, 70, 'F').in_range


def test_bilateral_and_tibia_bounds():
    bilateral = geometry.contralateral(22.0, 'L')
    assert bilateral.mean == pytest.approx(21.81029944978867)
    assert bilateral.lower < 22 < bilateral.upper
    assert geometry.contralateral(22.0, 'R').mean > 22
    tibia = geometry.prior_from_tibia(350.0, 65.0)
    assert tibia.in_range and tibia.lower < tibia.mean < tibia.upper
    assert not geometry.prior_from_tibia(500.0, 65.0).in_range
    pop = geometry.scalar_population_prior()
    assert pop['parameters'] == ('head_radius_mm', 'ccd_deg', 'anteversion_signed_deg')
    assert np.linalg.eigvalsh(pop['covariance']).min() > 0
