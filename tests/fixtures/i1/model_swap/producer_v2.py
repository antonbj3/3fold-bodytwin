#!/usr/bin/env python
"""Deterministic I1 model-swap producer stub, version 2 (contract mechanics only)."""
PRODUCER_VERSION = 2
MODEL_ID = "bodytwin-model-swap-stub-v2"
CHANGED = "v2 changes the model identity, so its sha256 differs from v1"


def produce():
    """Return the stub's (non-physiological) result marker."""
    return {"model_id": MODEL_ID, "version": PRODUCER_VERSION, "note": CHANGED}


if __name__ == "__main__":
    print(produce())
