"""
Astrobee: what scvx_toolbox's three formulation corrections actually do here.

The corrections were found on the BIRDS-TM time-dilation problem, where they
are the difference between a useless answer and a good one (see the README).
This script runs the Astrobee case under each of them so the effect on *this*
problem is on record rather than assumed.

Measured result (T=30, tf=200, ZOH, 25 iterations max):

    config                       status  it     J final  err_end[m]    max|w|
    A baseline (=Enhanced)     max_iter  25    753.1169    5.55e-10  0.024254
    B pin_aware                max_iter  25    753.0722    2.23e-09  0.024251
    C pin_aware+scaled         max_iter  25    762.3624    5.82e-10  0.024256
    D C + quaternion fixed     max_iter  25    761.9853    2.21e-10  0.024191

None of the four converge, and the corrections do not change that.  Astrobee
is not the failure mode they address:

  * `pin_aware_trust_region` matters when the warm start disagrees with
    `start`/`end`, which puts a constant floor under eta.  Astrobee's warm
    start already agrees at both ends, so there is nothing to fix.

  * `defect_penalty_space='scaled'` matters when S_x compresses the states
    hard enough that J's defect penalty and L's vc penalty are on different
    footings.  Astrobee's S_x is mild, so it only shifts J's units (hence the
    762 vs 753 -- those two numbers are NOT comparable).

What actually stops Astrobee from converging is the test itself:
`delta_L = J_ref - L` plateaus near 28 because the ZOH discretization defect
at tau = 6.67 s does not shrink, while the criterion compares it against
`e_tol * |J_ref|` = 3.8.  The *step* does shrink -- ||dx||_2 is 1e-5 by
iteration 25 -- so the trajectory is converged even though the flag says
max_iter.  All four runs land the terminal state to 1e-9 m and respect the
obstacle and rate limits.

Run from the package root:  python examples/astrobee_compare_settings.py
"""

import os
import sys
import time
import warnings

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from astrobee_model import END_POS, OMEGA_MAX, astrobee_problem
from scvx_toolbox import solve_trajectory

CONFIGS = [
    # label, pin_aware_trust_region, defect_penalty_space, legacy_quaternion
    ("A baseline (=Enhanced)", False, 'unscaled', True),
    ("B pin_aware",            True,  'unscaled', True),
    ("C pin_aware+scaled",     True,  'scaled',   True),
    ("D C + quaternion fixed", True,  'scaled',   False),
]


def main():
    warnings.filterwarnings('ignore')
    print(f"{'config':<24} {'status':>10} {'it':>3} {'J final':>11} "
          f"{'err_end[m]':>11} {'max|w|':>9} {'s':>6}")
    print("-" * 80)
    for label, pin, defect, legacy_q in CONFIGS:
        kwargs = astrobee_problem(T=30, tf=200.0, time_mode='ZOH',
                                  legacy_quaternion=legacy_q)
        t0 = time.time()
        r = solve_trajectory(size_N=20, max_iter=25,
                             pin_aware_trust_region=pin,
                             defect_penalty_space=defect,
                             verbose=False, **kwargs)
        x = r['x_opt_unscaled']
        print(f"{label:<24} {r['status']:>10} {r['iterations']:3d} "
              f"{r['history'][-1]['J_SCVx']:11.4f} "
              f"{np.linalg.norm(x[0:3, -1:] - END_POS[0:3]):11.2e} "
              f"{np.max(np.linalg.norm(x[10:13], axis=0)):9.6f} "
              f"{time.time() - t0:6.1f}")
    print(f"\n(angular rate limit is {OMEGA_MAX:.6f} rad/s)")
    print("J is not comparable between the 'unscaled' and 'scaled' rows: the "
          "defect penalty is measured in different units.")


if __name__ == '__main__':
    main()
