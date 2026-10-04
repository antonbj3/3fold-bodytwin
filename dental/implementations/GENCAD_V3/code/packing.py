"""Exactly equal inequality rows share their smallest RHS; lossless NPZ/LZMA."""
import numpy as np, io, zipfile
from scipy.sparse import csr_matrix

def reduce_equal_rows(A, b):
    A = csr_matrix(A)
    if not len(b) or not np.all(np.diff(A.indptr) == 3):
        return (A, b)
    key = np.c_[A.indices.reshape(-1, 3), A.data.reshape(-1, 3)]
    (unique, first, inverse) = np.unique(key, axis=0, return_index=True, return_inverse=True)
    rhs = np.full(len(first), np.inf)
    np.minimum.at(rhs, inverse, b)
    if not np.array_equal(unique[inverse], key) or np.any(rhs[inverse] > b):
        raise ValueError('lossless inequality packing failed')
    order = np.argsort(first)
    u = unique[order]
    return (csr_matrix((u[:, 3:].ravel(), u[:, :3].astype(np.int32).ravel(), 3 * np.arange(len(u) + 1, dtype=np.int32)), shape=(len(u), A.shape[1])), rhs[order])

def savez(path, **arrays):
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_LZMA) as z:
        for (name, array) in arrays.items():
            stream = io.BytesIO()
            np.save(stream, np.asarray(array), allow_pickle=False)
            info = zipfile.ZipInfo(name + '.npy', date_time=(2026, 10, 2, 0, 0, 0))
            info.compress_type = zipfile.ZIP_LZMA
            z.writestr(info, stream.getvalue())
