"""
COGU_C_caller_test.py -- solve the BIRDS-VAT time-dilation problem through
the compiled COGU_SCVx_main C library (via ctypes) instead of the slow
Python/cvxpy solve path, and reproduce the exact same plots and open-loop
validation as birds_vat_only_dilation.run_birds_vat_dilation(...).

This never imports sympy, cvxpy or cvxpygen: the problem's numeric data
comes from birds_vat_only_dilation.build_vat_problem_data(...) (pure numpy),
and the trajectory itself comes from the O3-compiled shared library, so a
solve here takes milliseconds instead of the ~10s+ the cvxpy/ECOS Python path
needs. That is what makes this file (and cogu_bindings.cogu_scvx_main) the
right building block for a Monte Carlo sweep, an interactive UI, or the
closed-loop / NMPC demo in VAT_COGU_C_caller_test.py.

Build the library first (from this directory):
    cmake -B build -DCMAKE_C_FLAGS="-O3"
    cmake --build build --config Release

Run:
    python COGU_C_caller_test.py                      # warm start built in Python (matches Python solve exactly)
    python COGU_C_caller_test.py --warm-mode slerp     # warm start built in C instead
    python COGU_C_caller_test.py --warm-mode linear --no-plot
"""

import argparse
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                    # this dir: cogu_bindings.py
sys.path.insert(0, os.path.dirname(_HERE))   # parent dir: birds_vat_only_dilation.py

from cogu_bindings import load_library, cogu_scvx_main, COGU_WarmMode  # noqa: E402
from birds_vat_only_dilation import (  # noqa: E402
    build_vat_problem_data, simulate_open_loop_u, _plot_birds_vat,
)

_WARM_MODE_MAP = {'user': COGU_WarmMode.USER, 'slerp': COGU_WarmMode.SLERP,
                  'linear': COGU_WarmMode.LINEAR}


def run_cogu_c_caller(
    problem_kwargs=None, solver_kwargs=None,
    warm_mode='user', lib_path=None,
    dt_state=0.1, plot=True, verbose=True,
):
    """
    Solve the same BIRDS-VAT time-dilation problem as
    birds_vat_only_dilation.run_birds_vat_dilation(...), but through the
    compiled C shared library instead of Python/cvxpy.

    Parameters
    ----------
    problem_kwargs : dict, optional -- passed to build_vat_problem_data(...).
    solver_kwargs : dict, optional -- overrides any COGU_SCVx_main tuning
        argument (see cogu_bindings.cogu_scvx_main's keyword list): e.g.
        {'max_iter': 40, 'e_tol': 1e-3}.
    warm_mode : {'user', 'slerp', 'linear'}
        'user'   -- Python builds warm_x/warm_u via numpy (slerp +
                    angular_velocity, no cvxpy) exactly as
                    birds_vat_only_dilation.py does, and passes them to C
                    verbatim (COGU_WARM_USER). Reproduces the Python solve's
                    trajectory bit-for-bit.
        'slerp'  -- C itself builds the SLERP + finite-difference-rate warm
                    start (COGU_WARM_SLERP) from `start`/`end` alone.
        'linear' -- C builds a plain linear-interpolation warm start
                    (COGU_WARM_LINEAR).
    lib_path : str, optional -- explicit path to libcogu_scvx.so / .dll;
        default searches this directory's build/build_final subfolders.
    dt_state : float -- open-loop RK4 validation step [s].
    plot : bool -- show the matplotlib figures (same 5 as birds_vat_only_dilation.py).
    verbose : bool -- print progress/results.

    Returns
    -------
    (result, data) -- `result` has the same shape birds_vat_only_dilation.py's
    caller expects (x_opt_unscaled, u_opt_unscaled, tf, history, open_loop,
    B, iterations, converged, status), plus the raw COGU outputs under
    result['cogu'] (J, L, rho, disc_ms, solve_ms, ...). `data` is
    build_vat_problem_data(...)'s return dict.
    """
    problem_kwargs = dict(problem_kwargs or {})
    solver_kwargs = dict(solver_kwargs or {})

    if warm_mode not in _WARM_MODE_MAP:
        raise ValueError(f"warm_mode must be one of {list(_WARM_MODE_MAP)}, got {warm_mode!r}")
    wm = _WARM_MODE_MAP[warm_mode]

    data = build_vat_problem_data(**problem_kwargs)
    lib = load_library(lib_path)

    if verbose:
        print("=" * 68)
        print("COGU_C_caller_test -- BIRDS-VAT via compiled C (O3) solver")
        print("=" * 68)
        print(f"warm_mode={warm_mode!r}   tf free in [{data['tf_min']}, {data['tf_max']}] s, "
              f"initial guess {data['tf_init']} s")

    call_kwargs = dict(
        defect_penalty_scaled=1, epsilon_stop=1e-3, e_tol=0.005,
        beta_gr=2.0, beta_sh=2.0, eta_init=1.0, eta1=10.0, eta0=0.001,
        rho2=0.7, rho1=0.1, rho0=0.0, max_stall=3, max_iter=15, size_n=10,
        verbose=(1 if verbose else 0),
        tf_fixed=data['tf_init'], lam=data['lamb'],
        Sx=data['S_x'], cx=data['c_x'], Su=data['S_u'], cu=data['c_u'],
        alpha_p=data['alpha_p'], tf_min=data['tf_min'], tf_max=data['tf_max'],
        tf_init=data['tf_init'], Sp=data['S_p'], cp=data['c_p'],
    )

    if wm == COGU_WarmMode.USER:
        call_kwargs.update(warm_x=data['warm_x'], warm_u=data['warm_u'],
                           quat_offset=-1, omega_offset=-1, omega_transform=None)
    elif wm == COGU_WarmMode.SLERP:
        # quat_offset=0 / omega_offset=4 / omega_transform=Jsat is exactly
        # BIRDS-VAT's state layout (quaternion at 0:4, angular momentum at
        # 4:7 = Jsat @ angular_rate) -- see cogu_scvx.h's COGU_WARM_SLERP doc.
        call_kwargs.update(warm_x=None, warm_u=None,
                           quat_offset=0, omega_offset=4, omega_transform=data['Jsat'])
    else:  # LINEAR
        call_kwargs.update(warm_x=None, warm_u=None,
                           quat_offset=-1, omega_offset=-1, omega_transform=None)

    call_kwargs.update(solver_kwargs)

    t0 = time.time()
    res = cogu_scvx_main(
        lib,
        start=data['start_pos'][:, 0], end=data['end_pos'][:, 0],
        dp=data['dyn_params_val'], warm_mode=wm,
        **call_kwargs,
    )
    total_ms = (time.time() - t0) * 1e3

    tf_opt = res['tf']
    x_value, u_value = res['x'], res['u']

    if verbose:
        print("\n" + "=" * 68)
        print("RESULTS (C solver)")
        print("=" * 68)
        print(f"Wall time (Python call overhead incl.): {total_ms:.2f} ms")
        print(f"  C-reported: discretization {res['disc_ms']:.2f} ms, "
              f"convex solve {res['solve_ms']:.2f} ms")
        print(f"Iterations:       {res['iterations']}")
        print(f"Converged:        {res['converged']}   (status: {res['status_name']})")
        print(f"Optimal tf:       {tf_opt:.4f} s     (initial guess {data['tf_init']} s, "
              f"bounds [{data['tf_min']}, {data['tf_max']}])")
        print(f"J = {res['J']:.6f}   L = {res['L']:.6f}   rho(last) = {res['rho']:.4f}")
        print(f"Quat start:       {x_value[0:4, 0].round(6)}")
        print(f"Quat final:       {x_value[0:4, -1].round(6)}   target "
              f"{data['end_quat'][:, 0]}")
        print(f"hb final:         {x_value[4:7, -1].round(9)}")

    # ---- Wrap the COGU outputs into the shape birds_vat_only_dilation's
    # ---- plotter (_plot_birds_vat) expects, so the two scripts share it.
    history = [{'iteration': i, 'tf': res['hist_tf'][i], 'J_SCVx': res['hist_J'][i]}
              for i in range(res['hist_len'])]
    result = dict(
        x_opt_unscaled=x_value, u_opt_unscaled=u_value, tf=tf_opt,
        iterations=res['iterations'], converged=res['converged'],
        status=res['status_name'], history=history, B=data['B'],
        cogu=res,
    )

    ol = simulate_open_loop_u(
        x0_quat=x_value[0:4, 0], x0_omega=data['Jsat_inv'] @ x_value[4:7, 0],
        u_value=u_value, tf_opt=tf_opt, Jsat=data['Jsat'], Jsat_inv=data['Jsat_inv'],
        B=data['B'], end_quat=data['end_quat'][:, 0], dt_state=dt_state, verbose=verbose,
    )
    ol['omega_guidance'] = data['Jsat_inv'] @ x_value[4:7, :]
    result['open_loop'] = ol

    if plot:
        _plot_birds_vat(result, data, ol)

    return result, data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--warm-mode', choices=list(_WARM_MODE_MAP), default='user',
                        help="how COGU_SCVx_main builds its warm start (default: user)")
    parser.add_argument('--no-plot', action='store_true', help='skip matplotlib figures')
    parser.add_argument('--quiet', action='store_true', help='suppress progress printing')
    parser.add_argument('--lib', default=None,
                        help='explicit path to libcogu_scvx.so / cogu_scvx.dll')
    args = parser.parse_args()
    run_cogu_c_caller(warm_mode=args.warm_mode, lib_path=args.lib,
                      plot=not args.no_plot, verbose=not args.quiet)
