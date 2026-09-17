#!/usr/bin/env python
"""Deterministic I1 model-swap producer stub, version 1 (contract mechanics only)."""
PRODUCER_VERSION = 1
MODEL_ID = "bodytwin-model-swap-stub-v1"


def produce():
    """Return the stub's (non-physiological) result marker."""
    return {"model_id": MODEL_ID, "version": PRODUCER_VERSION}


if __name__ == "__main__":
    print(produce())
