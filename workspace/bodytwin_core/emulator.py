"""A303 active-set emulator from BT-XF4, on its audited K1 scope."""
import sys
from .vendor.emulator import xf4_model
sys.modules.setdefault('xf4_model',xf4_model)
from .vendor.emulator import xf4_core
sys.modules.setdefault('xf4_core',xf4_core)
from .vendor.emulator.xf4_emul import XF4Emulator, rel_error
from .config import paths

def from_source(config=None):
    """Construct the source model with its original relative K1 asset path."""
    import os
    from .vendor.emulator.xf4_core import RegionEmulator
    source=paths(config)['xf4_results']
    before=os.getcwd()
    try:
        os.chdir(source)
        return RegionEmulator()
    finally:
        os.chdir(before)

def emulator_from_source(config=None, K=6, knn=1):
    import os
    source=paths(config)['xf4_results']
    before=os.getcwd()
    try:
        os.chdir(source)
        return XF4Emulator(K=K,knn=knn)
    finally:
        os.chdir(before)
