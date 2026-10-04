"""DENT-SOLVER-HIGHRES — CalculiX 2.23 with MKL PARDISO (and SPOOLES/ITERATIVE for comparison) for large static FE.

Why: the stock ccx (~/.local/opt/squashfs-root/usr/bin/ccx, FreeCAD AppImage) only has SPOOLES; it segfaulted (signal 11)
at 1.35 M DOF / 13.6 GB in F1 (C3D10 h/2, P1 UNN 24). This cell provides a PARDISO build and a runner.

BUILD (reproducible, nothing system-wide; everything under @DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/opt, 2026-09-23):
  O=@DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/opt; mkdir -p $O/src && cd $O/src
  curl -O http://www.dhondt.de/ccx_2.23.src.tar.bz2 && tar xjf ccx_2.23.src.tar.bz2
  curl -L -o spooles.2.2.tgz https://www.netlib.org/linalg/spooles/spooles.2.2.tgz
  mkdir SPOOLES.2.2 && tar xzf spooles.2.2.tgz -C SPOOLES.2.2 && cd SPOOLES.2.2
  sed -i 's#^  CC = /usr/lang-4.0/bin/cc#  CC = gcc#; s#^  OPTLEVEL = -O#  OPTLEVEL = -O3#' Make.inc
  make lib && (cd MT/src && make -f makeGlobalLib)          # MT objects are appended to spooles.a
  conda create -y -p $O/mklenv -c conda-forge mkl-devel=2024 mkl-include llvm-openmp   # MKL 2024 (PARDISO)
  cd $O/src/CalculiX/ccx_2.23/src && <write Makefile_F2 as below> && make -j10 -f Makefile_F2
  Makefile_F2 = stock Makefile_MT with
    CFLAGS += -I$(MKL)/include -DSPOOLES -DARPACK -DPARDISO -DMATRIXSTORAGE -DUSE_MT=1
    LIBS = spooles.a -larpack (system libarpack2) -L$(MKL)/lib -Wl,-rpath,$(MKL)/lib -Wl,--no-as-needed
           -Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -lm -ldl
  (ARPACK is required at link time because ccx_2.23.c/robustdesign.c/sensi_coor.c reference ARPACK-guarded symbols.)
  Result: $O/src/CalculiX/ccx_2.23/src/ccx_2.23_pardiso (LP64 integers: ≤ 2^31 matrix non-zeros in the upper triangle;
  PARDISO keeps its factor with 64-bit internal indexing; out-of-core not enabled).

USE: select the solver in the deck's *STATIC line (SOLVER=PARDISO | SPOOLES | ITERATIVE CHOLESKY | ITERATIVE SCALING).
  set_solver(inp_path, 'PARDISO') patches an existing deck; run(job, solver_bin=...) goes through tasks/heavy_run.sh
  and records /usr/bin/time -v maxRSS.  Units are whatever the deck uses (dental: mm, N, MPa).
Measured (F2, results/F2_pdl_nonlinear/solver_runs.jsonl): see REPORT.md there.
"""
from dental_release.paths import expand as _release_expand
import os
import re
import subprocess
import time
from pathlib import Path
OPT = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/opt'))
CCX_PARDISO = OPT / 'src/CalculiX/ccx_2.23/src/ccx_2.23_pardiso'
CCX_STOCK = Path(os.path.expanduser('~/.local/opt/squashfs-root/usr/bin/ccx'))
HEAVY = _release_expand('@DENTAL_INPUT_ROOT@/workspace/tasks/heavy_run.sh')

def set_solver(inp, solver):
    """Rewrite every '*STATIC' keyword line so that it carries SOLVER=<solver> (keeps other parameters)."""
    p = Path(inp)
    out = []
    for line in p.read_text().split('\n'):
        if line.upper().startswith('*STATIC'):
            parts = [x for x in line.split(',') if not x.strip().upper().startswith('SOLVER')]
            line = ','.join(parts + [f'SOLVER={solver}'])
        out.append(line)
    p.write_text('\n'.join(out))

def run(job, solver_bin=CCX_PARDISO, ram_gb=6, threads=4):
    """Run ccx on <job>.inp via heavy_run.sh (ram_gb >= 10 takes the bigmem lock). Returns ok, seconds, maxrss_GB, solver line."""
    job = Path(job)
    env = dict(os.environ, THREADS=str(threads))
    for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'CCX_NPROC_EQUATION_SOLVER', 'CCX_NPROC_STIFFNESS', 'CCX_NPROC_RESULTS'):
        env[k] = str(threads)
    t = time.time()
    r = subprocess.run([HEAVY, str(ram_gb), '0', '/usr/bin/time', '-v', '-o', f'{job.name}.time', str(solver_bin), '-i', job.name], cwd=job.parent, capture_output=True, text=True, env=env)
    (job.parent / f'{job.name}.log').write_text(r.stdout + r.stderr)
    maxrss = None
    tf = job.parent / f'{job.name}.time'
    if tf.exists():
        m = re.search('Maximum resident set size \\(kbytes\\): (\\d+)', tf.read_text())
        maxrss = int(m.group(1)) / 1000000.0 if m else None
    used = re.findall('(Factoring the system of equations using[^\\n]*|Solving the system of equations using[^\\n]*)', r.stdout)
    return {'ok': r.returncode == 0 and 'Job finished' in r.stdout, 'returncode': r.returncode, 'seconds': round(time.time() - t, 1), 'maxrss_GB': maxrss, 'solver_msgs': sorted(set(used))[:3], 'binary': str(solver_bin)}
