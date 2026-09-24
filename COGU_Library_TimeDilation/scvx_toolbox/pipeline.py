"""
`solve_trajectory` -- the one call a user needs.

The user describes the problem in SymPy and numpy; the library generates the
evaluation functions, builds the convex subproblem, runs SCvx, and optionally
emits C code.
"""

import importlib.util
import os
import tempfile

import numpy as np

from .modes import control_node_count, parse_time_mode
from .problem import build_problem
from .scaling import default_scaling, dilation_scaling
from .solver import solve_scvx
from .symbolic import SymbolicModel, generate_numpy_module

_LEGACY_ALIASES = {
    'etta0': 'eta0',
    'etta1': 'eta1',
    'etta_init': 'eta_init',
}


def _load_generated_module(source_path):
    spec = importlib.util.spec_from_file_location("_scvx_generated", source_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _default_warm_start_x(start, end, nx, T):
    warm = np.zeros((nx, T + 1))
    for k in range(T + 1):
        warm[:, k:k + 1] = start + (end - start) * k / T
    return warm


def solve_trajectory(states, controls, dynamics, start, end, T,
                     tf=None,
                     time_mode='ZOH',
                     tf_init=None, tf_min=None, tf_max=None, alpha_p=0.0,
                     nonconvex_constraints=None,
                     dynamic_parameters_sym=None, dynamic_parameters_val=None,
                     S_x=None, c_x=None, S_u=None, c_u=None,
                     S_p=None, c_p=None,
                     warm_start_x=None, warm_start_u=None,
                     state_bounds=None, control_bounds=None,
                     cost_terms=None,
                     size_N=20, lamb=1000.0, max_iter=25,
                     rho0=0.0, rho1=0.1, rho2=0.7,
                     eta0=1e-8, eta1=10.0, eta_init=1.0,
                     beta_sh=2.0, beta_gr=2.0,
                     e_tol=0.005, epsilon_stop=1e-5,
                     solver='ECOS', solver_options=None, verbose=True,
                     defect_substeps=1, defect_penalty_space='unscaled',
                     max_stall=3, pin_aware_trust_region=True,
                     dilation_symbol=None,
                     generate_c=False, c_output_dir='c_output',
                     solver_source_dir=None,
                     **legacy):
    """
    Solve a free or fixed final time trajectory optimization problem via SCvx.

    Required
    --------
    states, controls : list of SymPy symbols
    dynamics : list of SymPy expressions
        dx/dt = f(x, u, t) in *absolute* time.  Under a dilation mode the
        library forms p*f itself; do not pre-multiply it yourself unless you
        also pass `dilation_symbol`.
    start, end : (nx,1) numpy arrays
    T : int
        Number of intervals (so T+1 nodes).
    cost_terms : list of dict
        Running cost, see cost_dsl.

    Time mode
    ---------
    time_mode : 'ZOH' | 'FOH' | 'ZOH+dilation' | 'FOH+dilation'
    tf : float
        Fixed final time.  Required for the non-dilation modes; used as the
        initial guess for the dilation modes if `tf_init` is not given.
    tf_init, tf_min, tf_max : float
        Initial guess and bounds on the free final time.  Defaults:
        tf_init = tf, tf_min = 0.1*tf_init, tf_max = 3*tf_init.
    alpha_p : float
        Weight of the convex terminal cost alpha_p*tf.  0 leaves the final
        time free to whatever the rest of the cost prefers; larger values
        push it towards tf_min (minimum-time).

    Optional
    --------
    nonconvex_constraints : list of SymPy expressions g(x,u) <= 0
    dynamic_parameters_sym / _val : symbols and values of model constants
    S_x, c_x, S_u, c_u : scaling (default identity / zero)
    S_p, c_p : dilation scaling (default derived from tf_min/tf_max)
    warm_start_x : (nx, T+1);  warm_start_u : (nu, T) or (nu, T+1) for FOH
    state_bounds, control_bounds : bound DSL, see problem.build_bound_constraints
    defect_penalty_space : 'unscaled' | 'scaled'
        Units the defect penalty of J is measured in.  'unscaled' reproduces
        the reference notebooks; 'scaled' puts J's defect penalty on exactly
        the same footing as the vc penalty inside L, which is worth trying if
        rho comes out negative on steps that clearly improve J.
    max_stall : int
        Give up after this many consecutive rejected steps taken with the
        trust region already at eta0.
    pin_aware_trust_region : bool
        Keep the boundary-pinned states x_0 and x_T out of the trust region.
        Leave this True unless you are reproducing the legacy formulation:
        when the warm start disagrees with `start`/`end`, including them puts
        a hard floor under eta and the subproblem eventually goes infeasible.
    generate_c : emit C code after a successful Python solve

    Returns
    -------
    dict with x_opt / u_opt (scaled), x_opt_unscaled / u_opt_unscaled, tf,
    t_nodes, converged, iterations and history.
    """
    for old, new in _LEGACY_ALIASES.items():
        if old in legacy and legacy[old] is not None:
            if old == 'etta0':
                eta0 = legacy.pop(old)
            elif old == 'etta1':
                eta1 = legacy.pop(old)
            else:
                eta_init = legacy.pop(old)
    if legacy:
        raise TypeError(f"unexpected keyword arguments: {sorted(legacy)}")

    mode_name, hold, dilation = parse_time_mode(time_mode)

    nx = len(states)
    nu = len(controls)
    ng = len(nonconvex_constraints) if nonconvex_constraints else 0
    u_nodes = control_node_count(hold, T)

    start = np.asarray(start, dtype=float).reshape(nx, 1)
    end = np.asarray(end, dtype=float).reshape(nx, 1)

    # ------------------------------------------------------------------
    # Final time bookkeeping
    # ------------------------------------------------------------------
    if dilation:
        if tf_init is None:
            if tf is None:
                raise ValueError(
                    "a dilation mode needs an initial final time: pass tf_init "
                    "(or tf, which is used as the initial guess)"
                )
            tf_init = float(tf)
        tf_init = float(tf_init)
        tf_min = 0.1 * tf_init if tf_min is None else float(tf_min)
        tf_max = 3.0 * tf_init if tf_max is None else float(tf_max)
        if not (tf_min <= tf_init <= tf_max):
            raise ValueError(
                f"tf_init ({tf_init}) must lie within [tf_min, tf_max] "
                f"([{tf_min}, {tf_max}])"
            )
        if S_p is None or c_p is None:
            S_p, c_p = dilation_scaling(tf_min, tf_max)
        tf_for_config = tf_init
    else:
        if tf is None:
            raise ValueError("a fixed-time mode needs tf")
        tf_for_config = float(tf)
        S_p = c_p = None

    # ------------------------------------------------------------------
    # Symbolic model -> numpy evaluation functions
    # ------------------------------------------------------------------
    model = SymbolicModel(states, controls, dynamics,
                          nonconvex_constraints=nonconvex_constraints,
                          dynamic_parameters=dynamic_parameters_sym,
                          dilation=dilation,
                          dilation_symbol=dilation_symbol)

    handle = tempfile.NamedTemporaryFile(suffix='.py', delete=False,
                                         prefix='scvx_generated_')
    handle.close()
    generated_path = handle.name

    try:
        generate_numpy_module(model, output_file=generated_path)
        generated = _load_generated_module(generated_path)

        funcs = {
            'f': generated.f_dynamics,
            'A': generated.A_jacobian,
            'B': generated.B_jacobian,
            'F': generated.F_jacobian,
            'y': generated.y_affine,
            'g': generated.g_constraints if ng > 0 else None,
            'C': generated.C_jacobian,
            'D': generated.D_jacobian,
            'G': generated.G_jacobian,
            'z': generated.z_affine,
        }

        # --------------------------------------------------------------
        # Scaling and warm start
        # --------------------------------------------------------------
        S_x, c_x, S_u, c_u = default_scaling(nx, nu, S_x, c_x, S_u, c_u)

        if warm_start_x is None:
            warm_start_x = _default_warm_start_x(start, end, nx, T)
        warm_start_x = np.asarray(warm_start_x, dtype=float)

        if warm_start_u is None:
            warm_start_u = np.zeros((nu, u_nodes))
        warm_start_u = np.asarray(warm_start_u, dtype=float)

        # --------------------------------------------------------------
        # Convex subproblem
        # --------------------------------------------------------------
        tau_val = 1.0 / T if dilation else tf_for_config / T

        build_kwargs = dict(
            nx=nx, nu=nu, T=T, ng=ng,
            time_mode=mode_name,
            tau_val=tau_val, lamb_val=lamb,
            S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
            S_p=S_p, c_p=c_p,
            alpha_p=alpha_p, tf_min=tf_min, tf_max=tf_max,
            state_bounds=state_bounds,
            control_bounds=control_bounds,
            cost_terms=cost_terms,
            pin_aware_trust_region=pin_aware_trust_region,
        )
        prob_dict = build_problem(**build_kwargs)

        if verbose:
            print(f"[scvx] mode={mode_name}  nx={nx} nu={nu} ng={ng} "
                  f"T={T} control nodes={u_nodes}  tau={tau_val:.6g}"
                  + (f"  tf in [{tf_min}, {tf_max}] (init {tf_init})"
                     if dilation else f"  tf={tf_for_config}"))

        # --------------------------------------------------------------
        # Solve
        # --------------------------------------------------------------
        dyn_par = (np.asarray(dynamic_parameters_val, dtype=float)
                   if dynamic_parameters_val is not None else np.array([]))

        result = solve_scvx(prob_dict, funcs, config={
            'tf': tf_for_config,
            'size_N': size_N,
            'dynamic_parameters': dyn_par,
            'S_x': S_x, 'c_x': c_x, 'S_u': S_u, 'c_u': c_u,
            'S_p': S_p, 'c_p': c_p,
            'start_pos': start, 'end_pos': end,
            'warm_start_x': warm_start_x,
            'warm_start_u': warm_start_u,
            'warm_start_tf': tf_init if dilation else tf_for_config,
            'rho0': rho0, 'rho1': rho1, 'rho2': rho2,
            'eta0': eta0, 'eta1': eta1, 'eta_init': eta_init,
            'beta_sh': beta_sh, 'beta_gr': beta_gr,
            'lamb': lamb, 'alpha_p': alpha_p,
            'e_tol': e_tol, 'epsilon_stop': epsilon_stop,
            'max_iter': max_iter,
            'solver': solver, 'solver_options': solver_options,
            'verbose': verbose,
            'defect_substeps': defect_substeps,
            'defect_penalty_space': defect_penalty_space,
            'max_stall': max_stall,
        })

        result['model'] = model
        result['problem'] = prob_dict

        # --------------------------------------------------------------
        # C generation (only after a successful Python solve)
        # --------------------------------------------------------------
        if generate_c:
            from .cgen import emit_c_project
            emit_c_project(
                model=model, prob_dict=prob_dict, build_kwargs=build_kwargs,
                output_dir=c_output_dir,
                solver_source_dir=solver_source_dir,
                embed={
                    'time_mode': mode_name,
                    'tf': tf_for_config, 'lam': lamb, 'alpha_p': alpha_p,
                    'tf_min': tf_min, 'tf_max': tf_max, 'tf_init': tf_init if dilation else tf_for_config,
                    'start': start, 'end': end, 'dp': dyn_par,
                    'S_x': S_x, 'c_x': c_x, 'S_u': S_u, 'c_u': c_u,
                    'S_p': S_p, 'c_p': c_p,
                    'warm_x': warm_start_x, 'warm_u': warm_start_u,
                    'solver_params': {
                        'SIZE_N': int(size_N),
                        'MAX_ITER': int(max_iter),
                        'MAX_STALL': int(max_stall),
                        'RHO0': rho0, 'RHO1': rho1, 'RHO2': rho2,
                        'ETA0': eta0, 'ETA1': eta1, 'ETA_INIT': eta_init,
                        'BETA_SH': beta_sh, 'BETA_GR': beta_gr,
                        'E_TOL': e_tol, 'EPSILON_STOP': epsilon_stop,
                        'DEFECT_PENALTY_SCALED':
                            defect_penalty_space == 'scaled',
                    },
                },
                verbose=verbose,
            )

        return result

    finally:
        try:
            os.unlink(generated_path)
        except OSError:
            pass
