"""
BIRDS-VAT: what does time dilation actually cost, and what does it buy?

The controlled comparison.  One library, one model, one warm start, one set of
SCvx tuning -- the only thing that changes is `time_mode`:

    ZOH             tf pinned at tf_init
    ZOH+dilation    tf free in [tf_min, tf_max]

so any difference in wall time, iteration count or trajectory is attributable
to the dilation formulation and nothing else.

What to expect, mechanically: dilation adds one decision variable (p = tf),
one Jacobian per node (F = df/dp), one more ODE inside the RK4 discretizer
(H' = A H + F), one extra term in the trust region and two bound constraints.
So an *iteration* is strictly more expensive.  What it buys is a better
answer -- the final time stops being a guess.  Whether it also needs fewer
iterations is not something the formulation guarantees; that is what this
script measures.

Note that `J` is NOT comparable between the two rows: `tau` is tf/T in the
fixed-time formulation and the normalized 1/T under dilation, so the running
cost is scaled differently.  Compare tf, the open-loop accuracy, wall time
and iteration count.

Run from the package root:  python examples/birds_vat_timing.py
"""

import os
import sys
import time
import warnings

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from birds_vat_only_dilation import (_build_symbolic_model,
                                     build_vat_problem_data)
from scvx_toolbox import solve_trajectory


def run(mode, data, states, controls, dynamics, dyn_params_sym, verbose=False):
    """Solve the BIRDS-VAT problem in one time mode; return timings + result."""
    T = data['T']
    common = dict(
        states=states, controls=controls, dynamics=dynamics,
        start=data['start_pos'], end=data['end_pos'], T=T,
        dynamic_parameters_sym=dyn_params_sym,
        dynamic_parameters_val=data['dyn_params_val'],
        S_x=data['S_x'], c_x=data['c_x'], S_u=data['S_u'], c_u=data['c_u'],
        warm_start_x=data['warm_x'], warm_start_u=data['warm_u'],
        cost_terms=data['cost_terms'],
        state_bounds=data['state_bounds'],
        control_bounds=data['control_bounds'],
        size_N=20, lamb=data['lamb'], max_iter=25,
        rho0=0.0, rho1=0.1, rho2=0.7,
        eta0=0.001, eta1=10.0, eta_init=1.0,
        beta_sh=2.0, beta_gr=2.0, e_tol=0.005,
        epsilon_stop=1e-3, defect_penalty_space='scaled',
        solver='ECOS', verbose=verbose, generate_c=False,
    )
    if mode == 'ZOH':
        # tf pinned at the dilation run's initial guess: the "before" case.
        extra = dict(time_mode='ZOH', tf=data['tf_init'])
    else:
        extra = dict(time_mode='ZOH+dilation',
                     tf_init=data['tf_init'], tf_min=data['tf_min'],
                     tf_max=data['tf_max'], alpha_p=data['alpha_p'],
                     S_p=data['S_p'], c_p=data['c_p'])

    t0 = time.time()
    result = solve_trajectory(**common, **extra)
    wall = time.time() - t0

    hist = result['history']
    return {
        'mode': mode,
        'wall': wall,
        'iterations': result['iterations'],
        'converged': result['converged'],
        'status': result['status'],
        'tf': result['tf'],
        'J': hist[-1]['J_SCVx'],
        'disc_ms': float(np.mean([h['discretization_ms'] for h in hist])),
        'solve_ms': float(np.mean([h['solve_ms'] for h in hist])),
        'per_iter': wall / max(result['iterations'], 1),
        'result': result,
    }


def main():
    warnings.filterwarnings('ignore')

    data = build_vat_problem_data()
    states, controls, dynamics, dyn_params_sym = _build_symbolic_model()

    print("BIRDS-VAT, same model and tuning, only time_mode differs")
    print(f"T = {data['T']}, tf_init = {data['tf_init']} s, "
          f"tf in [{data['tf_min']}, {data['tf_max']}], "
          f"alpha_p = {data['alpha_p']}")
    print()

    rows = [run(m, data, states, controls, dynamics, dyn_params_sym)
            for m in ('ZOH', 'ZOH+dilation')]

    print(f"{'mode':<14} {'status':>10} {'it':>3} {'tf [s]':>9} "
          f"{'wall [s]':>9} {'s/iter':>8} {'disc ms':>8} {'solve ms':>9}")
    print("-" * 76)
    for r in rows:
        print(f"{r['mode']:<14} {r['status']:>10} {r['iterations']:3d} "
              f"{r['tf']:9.3f} {r['wall']:9.1f} {r['per_iter']:8.2f} "
              f"{r['disc_ms']:8.1f} {r['solve_ms']:9.1f}")

    fixed, dil = rows
    print()
    print(f"tf:              {fixed['tf']:.2f} s fixed  ->  {dil['tf']:.2f} s optimized "
          f"({100 * (fixed['tf'] - dil['tf']) / fixed['tf']:.1f}% shorter manoeuvre)")
    print(f"iterations:      {fixed['iterations']}  ->  {dil['iterations']}")
    print(f"wall time:       {fixed['wall']:.1f} s  ->  {dil['wall']:.1f} s "
          f"({100 * (dil['wall'] - fixed['wall']) / fixed['wall']:+.1f}%)")
    print(f"cost per iter:   {fixed['per_iter']:.2f} s  ->  {dil['per_iter']:.2f} s "
          f"({100 * (dil['per_iter'] - fixed['per_iter']) / fixed['per_iter']:+.1f}%)")
    print()
    print("J is not comparable across the two rows (tau is tf/T vs 1/T).")


if __name__ == '__main__':
    main()
