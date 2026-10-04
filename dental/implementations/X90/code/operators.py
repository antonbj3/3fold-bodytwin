"""Call reviewed functions; never mutate their trees or recalibrate their models."""
import importlib.util
import sys
from pathlib import Path

def module(name, path):
    path = Path(path)
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj

def load(package, nerve):
    package = Path(package)
    r4 = package / 'batch27/demos/PROOF_LANE_FULL_CROWN_R4/code'
    sys.modules.pop('common', None)
    module('common', r4 / 'common.py')
    geom = module('geometry', r4 / 'geometry.py')
    module('generate', r4 / 'generate.py')
    generator = module('generate_c', r4 / 'generate_c.py')
    adapt = module('generate_d', r4 / 'generate_d.py')
    dep = package / 'batch9/demos/PATIENT360/dependencies/results'
    source = module('x90_source_geometry', dep / 'LANE_X18_CROWN_ANTAGONIST/code/geometry.py')
    sys.path.insert(0, str(dep / 'LANE_X34_STL_DESIGN_GATE/code'))
    from designgate.gate import check
    force = module('x90_height_force', package / 'batch23/demos/X82/code/crown_corridor.py')
    force_bridge = module('x90_force_bridge', package / 'batch23/demos/X82/code/consumer_bridge.py')
    pulpa = module('x90_pulp_consumer', package / 'demos/X12/code/consumer_ports.py')
    nerve = module('x90_nerve_science', Path(nerve) / 'science.py')
    cap = module('x90_cap_plugin', package / 'batch19/demos/FULL_CROWN/code/homolog_plugin.py').cap
    return (geom, generator, adapt, source, check, force, pulpa, nerve, cap, force_bridge)
