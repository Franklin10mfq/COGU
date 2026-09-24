"""
Astrobee: scvx_toolbox must reproduce COGU_Library_Enhanced.

Enhanced (`sympy_cvx_converter`) is the library that was in use before the
Waskiri fork, and its Astrobee example is the case the team validated by hand.
This test pins the numbers that library produces so a regression in
scvx_toolbox shows up as a failure rather than as a plot that looks a bit off.

The reference below is Enhanced's own Python run at T=30, tf=200.  Its
generated C binary (astrobee_c_output/build/cogu_scvx.exe) agrees with it, so
one reference covers both.

Why only the first ITERS_COMPARED iterations are checked
--------------------------------------------------------
Neither library converges on this problem: both hit max_iter with rho negative
and the trust region collapsing, and ECOS reports "solution may be inaccurate"
from roughly iteration 5 on.  Past iteration ~19 the two runs walk apart into
different points of an almost-flat basin (J moves by 5e-4 relative while L
moves by 2e-6), which is solver noise, not a formulation difference.  Checking
the early iterations tightly is the meaningful test; iteration 1 agrees to
1e-12.

That non-convergence is a property of the Astrobee example itself and predates
scvx_toolbox -- see examples/astrobee_compare_settings.py for what the
new formulation's corrections do (and do not do) about it.

Run from the package root:  python tests/test_astrobee_vs_enhanced.py
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'examples'))

import numpy as np

from astrobee_model import astrobee_problem
from scvx_toolbox import solve_trajectory

# COGU_Library_Enhanced, Astrobee at T=30 / tf=200: (L_SCVx, J_SCVx)
ENHANCED = {
    1: (4919.371870, 14199.891376),
    2: (896.069943, 4835.872459),
    3: (738.630283, 1058.388493),
    4: (726.894138, 758.281008),
    5: (723.947138, 765.663941),
    6: (723.948460, 765.641915),
    7: (723.949163, 765.636085),
    8: (723.947813, 765.648575),
    9: (723.948217, 765.644169),
    10: (723.948793, 765.635010),
    11: (723.948445, 765.636523),
    12: (723.952708, 764.647583),
}

ITERS_COMPARED = 12
TOL = 1e-5


def main():
    kwargs = astrobee_problem(T=30, tf=200.0, time_mode='ZOH',
                              legacy_quaternion=True)
    result = solve_trajectory(size_N=20, max_iter=ITERS_COMPARED,
                              pin_aware_trust_region=False,
                              defect_penalty_space='unscaled',
                              verbose=False, **kwargs)

    history = result['history']
    print(" it |            L     L Enhanced     rel "
          "|            J     J Enhanced     rel")
    worst = 0.0
    for it in sorted(ENHANCED):
        if it > len(history):
            print(f"{it:3d} | run stopped early ({result['status']})")
            worst = float('inf')
            break
        L_ref, J_ref = ENHANCED[it]
        h = history[it - 1]
        dL = abs(h['L_SCVx'] - L_ref) / abs(L_ref)
        dJ = abs(h['J_SCVx'] - J_ref) / abs(J_ref)
        worst = max(worst, dL, dJ)
        print(f"{it:3d} | {h['L_SCVx']:12.6f} {L_ref:14.6f} {dL:7.1e} "
              f"| {h['J_SCVx']:12.6f} {J_ref:14.6f} {dJ:7.1e}")

    print()
    print(f"worst relative difference: {worst:.2e}  (tolerance {TOL:.0e})")
    ok = worst < TOL
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
