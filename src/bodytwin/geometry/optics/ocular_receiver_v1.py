"""Normalize integer receiver energy; no implicit artifact writes."""
import numbers
import numpy as np


def normalize(histogram, edges, input_units, overflow_units=0):
    counts = np.asarray(histogram)
    edges = np.asarray(edges, dtype=np.float64)
    if (counts.ndim != 2 or counts.dtype.kind not in 'iu' or
            np.any(counts < 0) or edges.ndim != 1 or
            len(edges) < 2 or counts.shape != (len(edges)-1, len(edges)-1) or
            not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0)):
        raise ValueError('Expected nonnegative integer square map and increasing finite edges')
    if (not isinstance(input_units, numbers.Integral) or isinstance(input_units, bool) or input_units <= 0 or
            not isinstance(overflow_units, numbers.Integral) or isinstance(overflow_units, bool) or overflow_units < 0):
        raise ValueError('Expected positive integer input energy and nonnegative overflow')
    # Python integer summation avoids silently overflowing the energy contract.
    total = sum(map(int, counts.flat)) + int(overflow_units)
    if total > input_units:
        raise ValueError('Receiver energy exceeds launched energy')
    area = np.multiply.outer(np.diff(edges), np.diff(edges))
    if not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError('Bin area must be finite and positive')
    density = counts.astype(np.float64) / int(input_units) / area
    if not np.isfinite(density).all():
        raise ValueError('Normalized density is not finite')
    return density, area, int(overflow_units)/int(input_units), total/int(input_units)


if __name__ == '__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[4]/'probes/optics/ocular_receiver_probe.py'), run_name='__main__')
