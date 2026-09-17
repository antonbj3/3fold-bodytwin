"""Optional external reference data for tests.

A few checks compare against reference files that are not part of this repository. Each file
or directory is located by one environment variable:

* unset -> the check is skipped (the reference is optional);
* set, but the path does not exist -> the check FAILS, so a misconfigured reference run cannot
  look green.
"""
import os
from pathlib import Path

import pytest

REFERENCE_VARIABLES = {
    "BODYTWIN_REF_AREA_DESCRIPTOR": (
        "file", "earlier area-descriptor prototype (geometry_descriptor_v1.py)"),
    "BODYTWIN_REF_PROPERTY_REGISTRY_DIR": (
        "dir", "audited property registry (property_registry.csv, "
               "registry_regions_resolved.csv)"),
    "BODYTWIN_REF_ADAPTED_PRODUCER": (
        "file", "adapted public metabolic-cost producer script (metabolic_cost_adapted.py)"),
    "BODYTWIN_REF_SLS_SCRIPT": (
        "file", "reference exact SLS propagation script exposing simulate()"),
    "BODYTWIN_REF_UNIPD_PACK": (
        "dir", "parsed UniPD tendon pack (data/experimental_data.xlsx, data/unipd_series.npz, "
               "evidence/unipd_series_meta.json)"),
}


RETIRED_VARIABLE = "BODYTWIN_RESEARCH_ROOT"
RETIRED_MESSAGE = (
    "BODYTWIN_RESEARCH_ROOT is no longer used: set the individual reference variables instead ("
    + ", ".join(sorted(REFERENCE_VARIABLES)) + "; see tests/research_data.py)"
)


def check_retired_variable():
    """Fail loudly when the retired BODYTWIN_RESEARCH_ROOT is set (it would otherwise be ignored)."""
    if os.environ.get(RETIRED_VARIABLE):
        pytest.fail(RETIRED_MESSAGE, pytrace=False)


def optional_path(variable):
    """Return the configured reference path, skip when unset, fail when set but missing."""
    check_retired_variable()
    kind, what = REFERENCE_VARIABLES[variable]
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} not set: {what} not available")
    path = Path(value)
    exists = path.is_dir() if kind == "dir" else path.is_file()
    if not exists:
        pytest.fail(f"{variable}={value!r} is set but is not an existing {kind} ({what})")
    return path
