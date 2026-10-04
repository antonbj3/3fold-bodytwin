"""Guard coefficient roundoff while preserving frozen model values."""
import numpy as np

def models_match(old, new):
    if old.keys() != new.keys():
        return False
    for (device, model) in old.items():
        if model.keys() != new[device].keys():
            return False
        for (key, value) in model.items():
            if key in ('x', 'y'):
                if value != new[device][key]:
                    return False
            elif not np.allclose(value, new[device][key], rtol=1e-12, atol=1e-10):
                return False
    return True
