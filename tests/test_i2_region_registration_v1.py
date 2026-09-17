#!/usr/bin/env python3
"""Tests for the opt-in I2 region-registration cross-check (``i2_region_registration_v1``).

Synthetic fixture only (the repo's own layered box example); no physiological claims. Plain
python3 runnable + pytest compatible, mirroring the other I2 test modules' style.

The real frozen example binds its synthetic geometry labels (``SYNTH-BOX-*``) to the audited
``MSK-MUSCLE`` roster row via an explicit ``material_region_id``; the cross-check must find no
issue. Each mutation probe below makes the check fail for exactly one reason.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLE_PAYLOAD = REPO / 'examples' / 'geometry' / 'region_example.json'


from bodytwin.geometry.i2_region_registration_v1 import (  # noqa: E402
    assert_regions_registered, registration_issues)


def _example():
    return json.loads(EXAMPLE_PAYLOAD.read_text(encoding='utf-8'))


def _structural(issues):
    return [issue for issue in issues if issue.startswith('structural: ')]


def test_example_passes_with_empty_issue_list():
    """The real frozen example binds explicitly to MSK-MUSCLE (density 1.06 kg/L): no issue."""
    payload = _example()
    assert registration_issues(payload) == [], registration_issues(payload)
    assert assert_regions_registered(payload) is payload
    for region in payload['regions']:
        assert region['region_id'].startswith('SYNTH-BOX-')
        assert region['material_region_id'] == 'MSK-MUSCLE'
        assert region['density']['name'] == 'density_muscle'
        assert region['density']['value'] == 1.06
        assert region['density']['units'] == 'kg/L'


def test_plain_roster_region_id_without_explicit_binding_passes():
    """When material_region_id is absent, a region_id that is itself a roster id is accepted."""
    payload = _example()
    payload['regions'] = payload['regions'][:1]
    del payload['regions'][0]['material_region_id']
    payload['regions'][0]['region_id'] = 'MSK-MUSCLE'
    assert registration_issues(payload) == [], registration_issues(payload)


def test_unregistered_material_region_id_is_reported_by_name():
    """An explicit material_region_id outside the roster must be reported."""
    payload = _example()
    payload['regions'][0]['material_region_id'] = 'NOT-A-REGISTERED-REGION'
    issues = registration_issues(payload)
    assert any('NOT-A-REGISTERED-REGION' in issue for issue in issues), issues
    assert _structural(issues) == [], issues


def test_unregistered_region_id_without_binding_is_reported():
    """A synthetic label with no material_region_id and no roster region_id must be reported."""
    payload = _example()
    del payload['regions'][0]['material_region_id']
    issues = registration_issues(payload)
    assert any('SYNTH-BOX-1' in issue for issue in issues), issues


def test_tendon_density_missing_is_reported_and_opt_out():
    """MSK-TENDON is registered but its density is MISSING; opt-out drops only that issue."""
    payload = _example()
    payload['regions'][0]['material_region_id'] = 'MSK-TENDON'
    strict = registration_issues(payload)
    missing_issue = [issue for issue in strict if 'MSK-TENDON' in issue and 'MISSING' in issue]
    assert missing_issue, strict

    relaxed = registration_issues(payload, require_known_density=False)
    assert not [issue for issue in relaxed if 'MSK-TENDON' in issue and 'MISSING' in issue], relaxed


def test_density_value_disagreement_is_reported():
    """A declared density that disagrees with the KNOWN roster row is reported."""
    payload = _example()
    payload['regions'][0]['density']['value'] = 1.20
    issues = registration_issues(payload)
    assert any('density' in issue and 'disagree' in issue.lower() for issue in issues), issues
    assert any('1.2' in issue for issue in issues), issues
    assert _structural(issues) == [], issues


def test_structurally_invalid_payload_does_not_raise():
    """A bare-number mass still yields structural: issues and never raises."""
    payload = _example()
    payload['regions'][0]['mass'] = 0.03816  # bare number, no units
    issues = registration_issues(payload)
    assert _structural(issues), issues


def test_duplicate_region_ids_do_not_raise():
    """Duplicate region ids are caught structurally and never raise here."""
    payload = _example()
    payload['regions'][1]['region_id'] = payload['regions'][0]['region_id']
    issues = registration_issues(payload)
    assert any(issue.startswith('structural: ') and 'duplicate' in issue for issue in issues), issues


def test_non_dict_payload_does_not_raise():
    """A non-object payload returns a structural issue instead of raising."""
    issues = registration_issues([1, 2, 3])
    assert _structural(issues), issues


def test_assert_regions_registered_raises_on_invalid_payload():
    """assert_regions_registered raises ValueError carrying the issues, or returns the payload."""
    payload = _example()
    payload['regions'][0]['material_region_id'] = 'NOT-A-REGISTERED-REGION'
    try:
        assert_regions_registered(payload)
    except ValueError as error:
        assert 'NOT-A-REGISTERED-REGION' in str(error), error
    else:
        raise AssertionError('assert_regions_registered accepted an unregistered region id')


if __name__ == '__main__':
    functions = [name for name in sorted(globals())
                 if name.startswith('test_') and callable(globals()[name])]
    for name in functions:
        globals()[name]()
        print('PASS', name)
    print(f'{len(functions)} tests passed')
