"ToothFairy2 I/O: reads MetaImage (.mha, okomprimerad) directly out of zip without unpacking.\nIn: zip path, case-id (e.g. 'ToothFairy2F_001'). Output : numpy volume i (z,y,x) order + spacing ( mm , (z,y,x)).\nNo external ITK required; ElementType MET_UCHAR / MET_INT / MET_DOUBLE is supported.\n"
from dental_release.paths import expand as _release_expand
import zipfile, numpy as np
ZIP = _release_expand('@DENTAL_CORPUS_ROOT@/geometry/ToothFairy2/ToothFairy2_Dataset.zip')
ROOT = 'Dataset112_ToothFairy2'
_DT = {'MET_UCHAR': np.uint8, 'MET_CHAR': np.int8, 'MET_SHORT': np.int16, 'MET_USHORT': np.uint16, 'MET_INT': np.int32, 'MET_UINT': np.uint32, 'MET_FLOAT': np.float32, 'MET_DOUBLE': np.float64}

def read_mha_bytes(buf):
    i = buf.index(b'ElementDataFile')
    j = buf.index(b'\n', i) + 1
    hdr = dict((l.split(' = ', 1) for l in buf[:j].decode().strip().split('\n') if ' = ' in l))
    assert hdr.get('CompressedData', 'False') == 'False', "compressed MHA not supported"
    assert hdr['ElementDataFile'].strip() == 'LOCAL'
    dims = [int(v) for v in hdr['DimSize'].split()]
    sp = [float(v) for v in hdr['ElementSpacing'].split()]
    dt = np.dtype(_DT[hdr['ElementType']]).newbyteorder('>' if hdr.get('BinaryDataByteOrderMSB') == 'True' else '<')
    n = dims[0] * dims[1] * dims[2]
    arr = np.frombuffer(buf, dtype=dt, count=n, offset=j).reshape(dims[2], dims[1], dims[0])
    return (arr, tuple(sp[::-1]), hdr)

def case_ids(zpath=ZIP):
    z = zipfile.ZipFile(zpath)
    return sorted((n.split('/')[-1][:-4] for n in z.namelist() if '/labelsTr/' in n and n.endswith('.mha')))

def load(case, kind='label', zpath=ZIP):
    z = zipfile.ZipFile(zpath)
    name = f'{ROOT}/labelsTr/{case}.mha' if kind == 'label' else f'{ROOT}/imagesTr/{case}_0000.mha'
    return read_mha_bytes(z.read(name))
