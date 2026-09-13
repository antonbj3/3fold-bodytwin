PYTHON ?= python

.PHONY: verify verify-detectors verify-cells verify-selected verify-inventory
verify:
	$(PYTHON) examples/verification/verify_bodytwin.py

verify-detectors:
	$(PYTHON) examples/verification/verify_bodytwin.py --selected --family detectors

verify-cells:
	$(PYTHON) examples/verification/verify_bodytwin.py --selected --family cells

verify-selected:
	$(PYTHON) examples/verification/verify_bodytwin.py --selected

verify-inventory:
	$(PYTHON) examples/verification/verify_bodytwin.py --inventory
