"""Anonymous, byte-range ZIP extraction. Never downloads whole archive."""
from dental_release.paths import expand as _release_expand
import io, json, pathlib, struct, zipfile, requests, time, hashlib, zlib, shutil, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_sol_night/X69'))

class RemoteZIP(io.RawIOBase):

    def __init__(self, url, size):
        self.url = url
        self.size = size
        self.pos = 0
        self.transferred = 0
        self.logs = []

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        return self.pos

    def read(self, n=-1):
        n = min(n if n >= 0 else self.size - self.pos, self.size - self.pos)
        if n <= 0:
            return b''
        a = self.pos
        b = a + n - 1
        for attempt in range(3):
            r = requests.get(self.url, headers={'Range': f'bytes={a}-{b}'}, timeout=90)
            if r.status_code == 206 and r.headers.get('Content-Range', '').startswith(f'bytes {a}-{b}/'):
                break
            if attempt == 2:
                raise RuntimeError((r.status_code, dict(r.headers)))
        if len(r.content) != n:
            raise RuntimeError('short range')
        self.transferred += n
        self.logs.append({'start': a, 'end': b, 'bytes': n, 'status': r.status_code})
        self.pos += n
        return r.content

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    zenodo = '--zenodo' in sys.argv
    if zenodo:
        sys.argv.remove('--zenodo')
        meta = json.loads((ROOT / 'raw/zenodo_8027553.json').read_text())
        zf = meta['files'][0]
        f = {'download_url': zf['links']['self'], 'size': zf['size']}
    else:
        meta = json.loads((ROOT / 'raw/figshare_26965903.json').read_text())
        f = meta['files'][0]
    t = time.monotonic()
    remote = RemoteZIP(f['download_url'], f['size'])
    z = zipfile.ZipFile(remote)
    listing = [{'name': i.filename, 'bytes': i.file_size, 'compressed': i.compress_size, 'crc32': f'{i.CRC:08x}'} for i in z.infolist()]
    (ROOT / 'raw' / ('zenodo_zip_directory.json' if zenodo else 'remote_zip_directory.json')).write_text(json.dumps(listing, indent=2))
    print('members', len(listing))
    if '--list' in sys.argv:
        print(json.dumps(listing[:100], indent=2))
    if '--list' in sys.argv:
        return
    selected = sys.argv[1:]
    if not selected:
        raise SystemExit('Pass exact member names, or --list')
    manifest = []
    for name in selected:
        info = z.getinfo(name)
        if info.file_size > 1500000000:
            raise RuntimeError('member too large')
        used = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
        if used + info.file_size + info.compress_size > 3000000000:
            raise RuntimeError('lane intermediate cap')
        if shutil.disk_usage(DATA).free - info.file_size - info.compress_size < 5000000000:
            raise RuntimeError('disk reserve')
        dst = DATA / 'hao_demo' / name if zenodo else DATA / pathlib.PurePosixPath(name).name
        dst.parent.mkdir(parents=True, exist_ok=True)
        with z.open(name) as src, dst.open('wb') as out:
            sha = hashlib.sha256()
            n = 0
            while True:
                b = src.read(8 * 1024 * 1024)
                if not b:
                    break
                out.write(b)
                sha.update(b)
                n += len(b)
        manifest.append({'archive_url': f['download_url'], 'archive_size': f['size'], 'member': name, 'local_path': str(dst), 'bytes': n, 'sha256': sha.hexdigest(), 'zip_crc32_verified': True})
        print(manifest[-1], flush=True)
    (ROOT / 'raw' / ('HAO_ACQUISITION_MANIFEST.json' if zenodo else 'ACQUISITION_MANIFEST.json')).write_text(json.dumps({'files': manifest, 'transferred_bytes': remote.transferred, 'seconds': time.monotonic() - t, 'range_requests': remote.logs}, indent=2))
if __name__ == '__main__':
    main()
