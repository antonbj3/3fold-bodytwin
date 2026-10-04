import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from packing import *
from scipy.sparse import csr_matrix
records = []
for p in sorted((PAYLOAD / 'public/scenes').glob('*.npz')):
    before = p.stat().st_size
    with np.load(p, allow_pickle=False) as z:
        a = dict(z)
    A = csr_matrix((a['a_data'], a['a_indices'], a['a_indptr']), shape=(int(a['a_rows']), len(a['xy'])))
    (B, b) = reduce_equal_rows(A, a['b'])
    a.update(a_data=B.data, a_indices=B.indices, a_indptr=B.indptr, a_rows=B.shape[0], b=b)
    temp = p.with_suffix('.npz.tmp')
    savez(temp, **a)
    with np.load(temp, allow_pickle=False) as z:
        if any((not np.array_equal(z[k], np.asarray(v), equal_nan=True) for (k, v) in a.items())):
            raise RuntimeError('packing roundtrip failed')
    temp.replace(p)
    records.append(dict(file=p.name, before_bytes=before, after_bytes=p.stat().st_size, original_rows=A.shape[0], retained_rows=B.shape[0]))
dump(ROOT / 'raw/COMPACTION.json', dict(records=records, original_bytes=sum((r['before_bytes'] for r in records)), packed_bytes=sum((r['after_bytes'] for r in records)), scope='Exactly equal stored coefficient rows reduced to min RHS; lossless compression; no numerical tolerance or gate changed'))
