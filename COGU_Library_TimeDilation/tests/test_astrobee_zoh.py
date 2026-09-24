"""
Astrobee, fixed-time ZOH: parity against the reference notebook.

Reference: COGU_algorithms/Astrobee_template/
           SCVx_official_template_cvxpy_cvxpygen_ZOH.ipynb, cell 43.

The notebook runs T = 51-1 = 50, so this test does too.  An earlier Astrobee
example in sympy_cvx_converter compared a T=30 run against these T=50 numbers,
which is why the library looked wrong when it was not.

Three legacy switches are on, because that is what parity means here:

    pin_aware_trust_region=False   the notebook charges the pinned x_0, x_T
                                   against eta
    legacy_quaternion=True         the notebook feeds scalar-first quaternions
                                   to scipy's scalar-last from_quat, so its
                                   warm-start angular velocity is not the
                                   angular velocity of the slerp
    defect_penalty_space='unscaled'

Known difference that does NOT affect this comparison: the notebook adds the
vi (path-constraint slack) penalty to L twice -- cell 36 has the guarded block
followed by an unguarded copy of itself.  The obstacle constraints are
satisfied at the solution, so vi = 0 and the duplicate term vanishes.

Run from the package root:  python tests/test_astrobee_zoh.py
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'examples'))

import numpy as np

from astrobee_model import END_POS, OMEGA_MAX, START_POS, VEL_MAX, astrobee_problem
from scvx_toolbox import solve_trajectory

# notebook cell 43: iteration -> (L_SCVx_opt, J_SCVx_opt, oJ_SCVx)
NOTEBOOK = {
    1: (3267.226174, 8899.070398, 67974.116139),
    2: (947.987305, 3468.934579, 8899.070398),
    3: (739.425863, 959.195131, 3468.934579),
    4: (725.836670, 741.646912, 959.195131),
    5: (723.574405, 733.681643, 741.646912),
    6: (723.543791, 729.496018, 733.681643),
    7: (723.539466, 729.480704, 729.496018),
    8: (723.539309, 729.466957, 729.480704),
}

TOL = 5e-5      # relative; the notebook's numbers are quoted to 6 decimals


def main():
    kwargs = astrobee_problem(T=50, tf=200.0, time_mode='ZOH',
                              legacy_quaternion=True)
    result = solve_trajectory(size_N=20, max_iter=max(NOTEBOOK),
                              pin_aware_trust_region=False,
                              defect_penalty_space='unscaled',
                              verbose=False, **kwargs)

    history = result['history']
    print(" it |        L (lib)      L (nb)  rel "
          "|        J (lib)      J (nb)  rel")
    worst = 0.0
    for it in sorted(NOTEBOOK):
        if it > len(history):
            print(f"{it:3d} | run stopped early ({result['status']})")
            worst = float('inf')
            break
        L_nb, J_nb, _ = NOTEBOOK[it]
        h = history[it - 1]
        dL = abs(h['L_SCVx'] - L_nb) / abs(L_nb)
        dJ = abs(h['J_SCVx'] - J_nb) / abs(J_nb)
        worst = max(worst, dL, dJ)
        print(f"{it:3d} | {h['L_SCVx']:14.5f} {L_nb:11.5f} {dL:7.1e} "
              f"| {h['J_SCVx']:14.5f} {J_nb:11.5f} {dJ:7.1e}")

    x = result['x_opt_unscaled']
    print()
    print(f"err start : {np.linalg.norm(x[0:3, 0:1] - START_POS[0:3]):.2e} m")
    print(f"err end   : {np.linalg.norm(x[0:3, -1:] - END_POS[0:3]):.2e} m")
    print(f"max |v|   : {np.max(np.linalg.norm(x[3:6], axis=0)):.4f} "
          f"(limit {VEL_MAX})")
    print(f"max |w|   : {np.max(np.linalg.norm(x[10:13], axis=0)):.6f} "
          f"(limit {OMEGA_MAX:.6f})")
    print()
    print(f"worst relative difference: {worst:.2e}  (tolerance {TOL:.0e})")

    ok = worst < TOL
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
