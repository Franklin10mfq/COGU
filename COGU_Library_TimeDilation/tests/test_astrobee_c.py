"""
Astrobee end to end in C: generate, compile, run, compare against Python.

This is the case that exercises the parts of the generator the BIRDS-TM
example never touches, because BIRDS-TM has no path constraints:

    ng > 0          C/D/z Jacobians as per-node CVXPY parameters, the vi
                    slack, and the g() term of the nonlinear cost J
    nu > 1          cpg_update_ou_k takes an index (a size-1 parameter gets a
                    different CVXPYgen signature -- see cpg_compat.h)
    ZOH, fixed tf   the non-dilation branch of every #if in the template

Needs cmake and a C compiler on PATH, plus cvxpygen.  Generating the convex
solver is the slow part (~1 min); an existing `solver/` in the output
directory is reused, so a second run is much faster.

Run from the package root:  python tests/test_astrobee_c.py
"""

import os
import re
import shutil
import subprocess
import sys
import warnings

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'examples'))

import numpy as np

from astrobee_model import END_POS, START_POS, astrobee_problem
from scvx_toolbox import solve_trajectory

OUTPUT_DIR = os.path.join(_ROOT, 'astrobee_c_output')

# The two runs are the same algorithm on the same data, but ECOS is called
# from two different builds, and from iteration ~5 on this problem sits in a
# regime where ECOS reports "solution may be inaccurate" (rho negative, trust
# region collapsing -- see test_astrobee_vs_enhanced.py).
#
# The late bound is not a number picked to make this pass.  It is calibrated
# against an independent measurement of how noisy this tail is: running
# COGU_Library_Enhanced's Python solve against scvx_toolbox's Python solve --
# same language, same solver build, only the library differs -- the two walk
# apart by up to 5.8e-4 in J over iterations 20..25 while L moves by 2e-6.
# So anything under 1e-3 in the tail is that noise, not a code difference.
#
# What is actually load-bearing here is the early band, where the agreement is
# 1e-9..1e-7, and the terminal state check below, which is a hard equality
# constraint the C loop has to satisfy exactly.
TOL_EARLY = 1e-6      # iterations 1..N_EARLY
TOL_LATE = 1e-3       # everything after
N_EARLY = 4


def _binary():
    for candidate in ('build/scvx_cli_demo.exe', 'build/scvx_cli_demo',
                      'build/Release/scvx_cli_demo.exe'):
        path = os.path.join(OUTPUT_DIR, candidate)
        if os.path.isfile(path):
            return path
    return None


def parse_c_log(text):
    """(L, J) per iteration from the generated loop's verbose output."""
    return [(float(L), float(J)) for L, J in
            re.findall(r'L=\s*([\d.eE+-]+)\s*\|\s*J=\s*([\d.eE+-]+)', text)]


def main():
    warnings.filterwarnings('ignore')

    if shutil.which('cmake') is None:
        print("SKIP: cmake not on PATH")
        return 0

    # 1. Python solve + C generation
    print("[1/4] solving in Python and generating C ...")
    kwargs = astrobee_problem(T=30, tf=200.0, time_mode='ZOH',
                              legacy_quaternion=True)
    result = solve_trajectory(size_N=20, max_iter=25,
                              pin_aware_trust_region=False,
                              defect_penalty_space='unscaled',
                              generate_c=True, c_output_dir=OUTPUT_DIR,
                              verbose=False, **kwargs)
    py_hist = [(h['L_SCVx'], h['J_SCVx']) for h in result['history']]
    print(f"      python: {result['status']}, {result['iterations']} iterations")

    for name in ('scvx_main.c', 'dynamics.c', 'dynamics.h', 'cpg_compat.h',
                 'cogu_scvx.h', 'embedded_data.h', 'scvx_cli_demo.c',
                 'CMakeLists.txt'):
        path = os.path.join(OUTPUT_DIR, name)
        if not os.path.isfile(path):
            print(f"FAIL: generator did not emit {name}")
            return 1
    print("      all expected files emitted")

    # 2. Configure
    print("[2/4] cmake configure ...")
    cfg = subprocess.run(
        ['cmake', '-B', 'build', '-DCMAKE_BUILD_TYPE=Release',
         '-DCMAKE_C_FLAGS=-O3'],
        cwd=OUTPUT_DIR, capture_output=True, text=True)
    if cfg.returncode != 0:
        print("FAIL: cmake configure\n" + cfg.stdout[-2000:] + cfg.stderr[-2000:])
        return 1

    # 3. Build
    print("[3/4] cmake build ...")
    bld = subprocess.run(['cmake', '--build', 'build', '--config', 'Release'],
                         cwd=OUTPUT_DIR, capture_output=True, text=True)
    if bld.returncode != 0:
        print("FAIL: build\n" + bld.stdout[-3000:] + bld.stderr[-3000:])
        return 1

    binary = _binary()
    if binary is None:
        print("FAIL: scvx_cli_demo was not produced")
        return 1

    # 4. Run and compare
    print("[4/4] running the generated binary ...")
    run = subprocess.run([binary], cwd=OUTPUT_DIR, capture_output=True,
                         text=True, timeout=900)
    if run.returncode != 0:
        print(f"FAIL: binary exited {run.returncode}\n{run.stdout[-2000:]}")
        return 1

    c_hist = parse_c_log(run.stdout)
    if not c_hist:
        print("FAIL: could not parse the C iteration log\n" + run.stdout[-2000:])
        return 1

    m = re.search(r'iterations=(\d+)', run.stdout)
    if m and int(m.group(1)) != len(c_hist):
        print(f"FAIL: C reported iterations={m.group(1)} but printed "
              f"{len(c_hist)} iterations")
        return 1

    print()
    print(" it |     L Python          L C      rel "
          "|     J Python          J C      rel")
    ok = True
    worst_early = worst_late = 0.0
    n = min(len(py_hist), len(c_hist))
    for i in range(n):
        (Lp, Jp), (Lc, Jc) = py_hist[i], c_hist[i]
        dL = abs(Lc - Lp) / max(abs(Lp), 1e-12)
        dJ = abs(Jc - Jp) / max(abs(Jp), 1e-12)
        tol = TOL_EARLY if i < N_EARLY else TOL_LATE
        if i < N_EARLY:
            worst_early = max(worst_early, dL, dJ)
        else:
            worst_late = max(worst_late, dL, dJ)
        bad = max(dL, dJ) > tol
        ok &= not bad
        print(f"{i + 1:3d} | {Lp:11.5f} {Lc:12.5f} {dL:8.1e} "
              f"| {Jp:11.5f} {Jc:12.5f} {dJ:8.1e}{'  <-- over tol' if bad else ''}")

    print(f"\nworst relative difference, iterations 1..{N_EARLY}: "
          f"{worst_early:.1e}  (tolerance {TOL_EARLY:.0e})")
    if n > N_EARLY:
        print(f"worst relative difference, iterations {N_EARLY + 1}..{n}: "
              f"{worst_late:.1e}  (tolerance {TOL_LATE:.0e})")

    # The terminal state is a hard equality constraint; C must land it too.
    fin = re.search(r'final state:\s+(.*)', run.stdout)
    if fin:
        x_final = np.array([float(v) for v in fin.group(1).split()])
        err = np.linalg.norm(x_final[0:3] - END_POS[0:3, 0])
        print(f"\nC final position error: {err:.2e} m")
        if err > 1e-6:
            print("FAIL: C did not reach the terminal position")
            ok = False

    print("\nPASS" if ok else "\nFAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
