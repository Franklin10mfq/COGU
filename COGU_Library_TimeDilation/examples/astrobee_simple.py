"""
Astrobee via solve_trajectory() -- port of the example that shipped with
`COGU_Library_Enhanced` (sympy_cvx_converter), the library used before the
Waskiri fork.

Same problem, same numbers, same T=30: this is the reference case for checking
that scvx_toolbox did not change any behaviour on the way to adding free final
time.  `tests/test_astrobee_vs_enhanced.py` turns the comparison into a test.

Three flags reproduce the old library exactly:

    pin_aware_trust_region=False    Enhanced charges the pinned x_0 and x_T
                                    against the trust region
    defect_penalty_space='unscaled' Enhanced measures J's defect in physical
                                    units
    legacy_quaternion=True          Enhanced's angular_vel feeds scalar-first
                                    quaternions to scipy's scalar-last
                                    from_quat, so the warm-start angular
                                    velocity is not the rate of the slerp

Set them the other way (see astrobee_improved()) to get scvx_toolbox's
corrected formulation.

Run from the package root:  python examples/astrobee_simple.py
"""

import os
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from astrobee_model import (END_POS, OMEGA_MAX, START_POS, VEL_MAX,
                            astrobee_problem)
from scvx_toolbox import solve_trajectory

T = 30
TF = 200.0

# Enhanced's own output for this exact problem (its C binary,
# astrobee_c_output/build/cogu_scvx.exe, and its Python solve agree):
ENHANCED = {
    1: (4919.371848, 14199.891376),
    2: (896.069943, 4835.872459),
    3: (738.630373, 1058.388600),
    4: (726.894148, 758.281004),
}


def run(legacy=True, generate_c=False, c_output_dir='astrobee_c_output',
        max_iter=25, verbose=True):
    """Solve the Astrobee problem; `legacy=True` reproduces Enhanced."""
    kwargs = astrobee_problem(T=T, tf=TF, time_mode='ZOH',
                              legacy_quaternion=legacy)
    return solve_trajectory(
        size_N=20, max_iter=max_iter,
        pin_aware_trust_region=not legacy,
        defect_penalty_space='unscaled' if legacy else 'scaled',
        generate_c=generate_c, c_output_dir=c_output_dir,
        verbose=verbose, **kwargs)


def report(result, label):
    x = result['x_opt_unscaled']
    print(f"\n--- {label} ---")
    print(f"status     : {result['status']}  ({result['iterations']} iterations)")
    print(f"J final    : {result['history'][-1]['J_SCVx']:.4f}")
    print(f"start pos  : {x[0:3, 0].round(4)}  (err "
          f"{np.linalg.norm(x[0:3, 0:1] - START_POS[0:3]):.2e} m)")
    print(f"end pos    : {x[0:3, -1].round(4)}  (err "
          f"{np.linalg.norm(x[0:3, -1:] - END_POS[0:3]):.2e} m)")
    print(f"max |v|    : {np.max(np.linalg.norm(x[3:6], axis=0)):.4f} m/s "
          f"(limit {VEL_MAX})")
    print(f"max |w|    : {np.max(np.linalg.norm(x[10:13], axis=0)):.6f} rad/s "
          f"(limit {OMEGA_MAX:.6f})")


def main():
    t0 = time.time()
    legacy = run(legacy=True, verbose=False)
    report(legacy, "legacy formulation (reproduces COGU_Library_Enhanced)")

    print("\n it |          L         J |   L Enhanced   J Enhanced |     rel L     rel J")
    for it in sorted(ENHANCED):
        L_ref, J_ref = ENHANCED[it]
        h = legacy['history'][it - 1]
        print(f"{it:3d} | {h['L_SCVx']:10.4f} {h['J_SCVx']:9.4f} |"
              f" {L_ref:12.4f} {J_ref:12.4f} |"
              f" {abs(h['L_SCVx'] - L_ref) / L_ref:9.1e}"
              f" {abs(h['J_SCVx'] - J_ref) / J_ref:9.1e}")

    improved = run(legacy=False, verbose=False)
    report(improved, "corrected formulation (scvx_toolbox defaults)")

    print(f"\ntotal {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
