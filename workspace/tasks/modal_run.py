"""BodyTwin cloud run on Modal (replaces OVH when it is unreachable). Usage:
    modal run tasks/modal_run.py --job-id <ID> --src <local directory> --cmd "<command>" [--cpu 4|8|16] [--mem-gb 8|16|32]  (size s/m/l is selected; max 3 h)
The directory is uploaded, the command runs in it (bash -c), the whole directory is retrieved back to <src>/modal_out/<ID>/.
Image: debian_slim + numpy 2.2.6, scipy, h5py. Receipt is placed in tasks/cloud_receipts.jsonl. Write results in the working directory.
"""
import io, json, os, pathlib, tarfile, time
import modal

app = modal.App("bodytwin-jobs")
img = modal.Image.debian_slim().pip_install("numpy==2.2.6", "scipy", "h5py")
REC = pathlib.Path("tasks/cloud_receipts.jsonl")


def _tar(src: str) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        for p in pathlib.Path(src).rglob("*"):
            if "__pycache__" in p.parts or "modal_out" in p.parts:
                continue
            if p.is_file() and p.stat().st_size <= 200 * 1024 * 1024:
                t.add(p, arcname=str(p.relative_to(src)))
    return buf.getvalue()


def _run(payload: bytes, cmd: str, threads: int):
    import subprocess
    work = pathlib.Path("/work"); work.mkdir(exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as t:
        t.extractall(work)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    p = subprocess.run(["bash", "-c", cmd], cwd=work, env=env, capture_output=True, text=True)
    (work / "run.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    (work / "EXIT_CODE").write_text(str(p.returncode))
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        t.add(work, arcname=".")
    return buf.getvalue()


@app.function(image=img, cpu=4, memory=8192, timeout=3 * 3600)
def run_s(payload: bytes, cmd: str, threads: int):
    return _run(payload, cmd, threads)


@app.function(image=img, cpu=8, memory=16384, timeout=3 * 3600)
def run_m(payload: bytes, cmd: str, threads: int):
    return _run(payload, cmd, threads)


@app.function(image=img, cpu=16, memory=32768, timeout=3 * 3600)
def run_l(payload: bytes, cmd: str, threads: int):
    return _run(payload, cmd, threads)


def _fn(cpu, mem_gb, timeout_min):
    return run_s if cpu <= 4 and mem_gb <= 8 else run_m if cpu <= 8 and mem_gb <= 16 else run_l


@app.local_entrypoint()
def main(job_id: str, src: str, cmd: str, cpu: float = 8, mem_gb: float = 16, timeout_min: int = 60):
    t0 = time.time()
    out = run_remote(job_id, src, cmd, cpu, mem_gb, timeout_min)
    print(json.dumps({"job": job_id, "out": out, "seconds": round(time.time() - t0, 1)}))


def run_remote(job_id, src, cmd, cpu, mem_gb, timeout_min):
    f = _fn(cpu, mem_gb, timeout_min)
    res = f.remote(_tar(src), cmd, int(cpu))
    out = pathlib.Path(src) / "modal_out" / job_id
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(res), mode="r:gz") as t:
        t.extractall(out)
    with REC.open("a") as r:
        r.write(json.dumps({"t": time.strftime("%Y-%m-%dT%H:%M:%S"), "job": job_id, "machine": "modal", "dir": src,
                            "cpus": cpu, "mem_gb": mem_gb, "timeout_min": timeout_min, "cmd": cmd,
                            "exit_code": (out / "EXIT_CODE").read_text().strip()}) + "\n")
    return str(out)
