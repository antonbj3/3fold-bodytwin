"""R3 remaps immutable artifacts into the unchanged R2 validator."""
from pathlib import Path
from preparation import evaluate

class Remap:

    def __init__(self, bundle):
        self.base = bundle
        self.payload = bundle.payload / 'round3'

    def name(self, name):
        return name.replace('FROZEN_PREPARATIONS.json', 'FROZEN_PREPARATIONS_R3.json').replace('PREREG_R2.json', 'PREREG_R3.json').replace('payload/preparations/', 'payload/preparations_r3/')

    def json(self, name):
        return self.base.json(self.name(name))

    def npz(self, name):
        return self.base.npz(self.name(name))

    def bytes(self, name):
        return self.base.bytes(self.name(name))

def run(bundle):
    return evaluate(Remap(bundle))
