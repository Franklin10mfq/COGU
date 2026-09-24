"""
BIRDS-TM (VAT only) with ZOH + TIME DILATION.

Same physical problem as birds_vat_only.py -- pure attitude control,

    nx = 7   quaternion + body angular momentum
    nu = 4   VAT thrust magnitudes

but the final time is now a decision variable instead of a fixed 480 s.

WHY THIS FILE IS PARAMETERIZED
-------------------------------
Every physical/numeric quantity of the problem -- inertia, thruster geometry
and orientation, bounds, boundary conditions, tf settings, cost weight,
warm-start choice, SCvx tuning -- is a keyword argument, not a hardcoded
module-level constant. That is what lets this same file drive a Monte Carlo
study (resample the physical parameters, re-solve, repeat) or feed an
interface, instead of being a fixed one-shot script. Calling every function
below with no arguments reproduces today's example exactly.

The module is split in two:

    build_vat_problem_data(...)   pure numpy/scipy -- assembles every NUMERIC
                                  array the problem needs (inertia, scaling,
                                  boundary conditions, warm start, thruster
                                  torque matrix B, ...) and returns it as a
                                  plain dict. No sympy/cvxpy/cvxpygen import,
                                  so this is cheap and safe to call from
                                  anything, including a Monte Carlo loop or
                                  the fast-C caller scripts that don't want
                                  the slow Python solve path at all.

    run_birds_vat_dilation(...)   the SymPy model, the scvx_toolbox SCvx
                                  solve (optionally emitting C via
                                  generate_c=True), the open-loop RK4
                                  validation, and the plots -- i.e. today's
                                  script, now callable with any problem/
                                  solver configuration.

Run:
    python birds_vat_only_dilation.py
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from scvx_toolbox.utils import angular_velocity, slerp

# =============================================================================
# 1. PURE-NUMPY PROBLEM ASSEMBLY -- no sympy/cvxpy/cvxpygen import
# =============================================================================


def build_vat_problem_data(
    T=30,
    tf_init=480.0, tf_min=60.0, tf_max=900.0, alpha_p=0.05,
    # Inertia [kg m^2]
    J1=0.0333, J2=0.0333, J3=0.0067, Jxy=0.0001, Jxz=0.0001, Jyz=0.0001,
    # Thruster positions [mm] (converted to [m] via rth_scale) and VAT cant
    # angle [deg] -- the classic 4-thruster antagonistic-pair VAT layout.
    rth1=(-26.43, -26.61, -147.785),
    rth2=(26.43, -26.61, -147.785),
    rth3=(26.43, 26.61, -147.785),
    rth4=(-26.43, 26.61, -147.785),
    rth_scale=1.0e-3,
    theta_vat_deg=40.0,
    # Physical limits
    omegab_max_deg_s=1.5, fth_max=25e-6,
    # Boundary conditions: quaternion (scalar-first) + angular rate [deg/s]
    start_quat=(1.0, 0.0, 0.0, 0.0), end_quat=(0.0, 1.0, 0.0, 0.0),
    start_omega_deg_s=(0.2, 0.1, -0.2), end_omega_deg_s=(0.0, 0.0, 0.0),
    # Warm start: 'slerp' (default, matches the original script), 'linear',
    # or an explicit (warm_x, warm_u) pair of numpy arrays.
    warm_start='slerp',
    # Cost / SCvx-visible weights
    cost_weight=0.5, lamb=5000.0,
):
    """
    Assemble every numeric array the BIRDS-VAT time-dilation problem needs.

    Returns a dict with (among others): Jsat, Jsat_inv, dyn_params_val,
    start_pos, end_pos, S_x, c_x, S_u, c_u, S_p, c_p, warm_x, warm_u, B
    (thruster torque matrix), rth_list, fth_list, omegab_max, Fth_max,
    cost_terms, state_bounds, control_bounds, and every scalar argument
    passed in (T, tf_init, tf_min, tf_max, alpha_p, lamb, ...).
    """
    omegab_max = omegab_max_deg_s * np.pi / 180.0
    theta_vat = theta_vat_deg * np.pi / 180.0
    tau_absolute = tf_init / T  # step of the initial guess (warm start only)

    # ---- Inertia ------------------------------------------------------
    Jsat = np.array([[J1, Jxy, Jxz],
                     [Jxy, J2, Jyz],
                     [Jxz, Jyz, J3]])
    Jsat_inv = np.linalg.inv(Jsat)

    # ---- Thruster geometry ---------------------------------------------
    rth_d = [np.array(r, dtype=float).reshape(3, 1) * rth_scale
             for r in (rth1, rth2, rth3, rth4)]

    fth_dirs_unscaled = [
        np.array([[0.0], [np.sin(theta_vat)], [np.cos(theta_vat)]]),
        np.array([[0.0], [np.sin(theta_vat)], [np.cos(theta_vat)]]),
        np.array([[0.0], [-np.sin(theta_vat)], [np.cos(theta_vat)]]),
        np.array([[0.0], [-np.sin(theta_vat)], [np.cos(theta_vat)]]),
    ]
    fth_d = [f / np.linalg.norm(f) for f in fth_dirs_unscaled]

    dyn_params_val = np.array([
        J1, J2, J3, Jxy, Jxz, Jyz,
        *rth_d[0][:, 0], *rth_d[1][:, 0], *rth_d[2][:, 0], *rth_d[3][:, 0],
        *fth_d[0][:, 0], *fth_d[1][:, 0], *fth_d[2][:, 0], *fth_d[3][:, 0],
    ])

    B = np.column_stack([
        np.cross(rth_d[i][:, 0], fth_d[i][:, 0]) for i in range(4)
    ])

    # ---- Boundary conditions -------------------------------------------
    start_quat_a = np.array(start_quat, dtype=float).reshape(4, 1)
    end_quat_a = np.array(end_quat, dtype=float).reshape(4, 1)
    start_hb = Jsat @ (np.array(start_omega_deg_s, dtype=float).reshape(3, 1) * np.pi / 180.0)
    end_hb = Jsat @ (np.array(end_omega_deg_s, dtype=float).reshape(3, 1) * np.pi / 180.0)

    start_pos = np.block([[start_quat_a], [start_hb]]).astype(float)
    end_pos = np.block([[end_quat_a], [end_hb]]).astype(float)

    # ---- Scaling ---------------------------------------------------------
    S_x = np.diag([2, 2, 2, 2, 2 * J1 * omegab_max, 2 * J2 * omegab_max, 2 * J3 * omegab_max])
    c_x = np.array([[-1], [-1], [-1], [-1],
                    [-J1 * omegab_max], [-J2 * omegab_max], [-J3 * omegab_max]])
    S_u = np.diag([2 * fth_max] * 4)
    c_u = np.array([[-fth_max]] * 4)

    # Dilation scaling: hp = S_p*p + c_p, chosen so scaled p sits in [0,1]
    # across [tf_min, tf_max]. solve_trajectory derives this automatically
    # if omitted -- computed here too so build_vat_problem_data is self
    # contained (e.g. for a caller that never touches scvx_toolbox at all).
    S_p = np.array([[tf_max - tf_min]])
    c_p = np.array([[tf_min]])

    # ---- Warm start --------------------------------------------------
    if isinstance(warm_start, tuple):
        warm_x, warm_u = warm_start
        warm_x = np.asarray(warm_x, dtype=float)
        warm_u = np.asarray(warm_u, dtype=float)
    elif warm_start == 'slerp':
        quat = slerp(start_pos[0:4, 0], end_pos[0:4, 0], T + 1).T
        omega = angular_velocity(quat.T, tau_absolute).T
        hb_ws = Jsat @ omega
        warm_x = np.zeros((7, T + 1))
        warm_x[0:4, :] = quat
        warm_x[4:7, :] = hb_ws
        warm_u = np.zeros((4, T))
    elif warm_start == 'linear':
        warm_x = np.zeros((7, T + 1))
        for k in range(T + 1):
            a = k / T
            warm_x[:, k] = (1 - a) * start_pos[:, 0] + a * end_pos[:, 0]
        warm_u = np.zeros((4, T))
    else:
        raise ValueError(f"warm_start must be 'slerp', 'linear', or an "
                         f"(warm_x, warm_u) tuple; got {warm_start!r}")

    # ---- Cost / constraints (DSL dicts, consumed by scvx_toolbox) --------
    cost_terms = [
        {'kind': 'norm1', 'var': 'u', 'slice': slice(0, 4),
         'coeff': 'tau', 'weight': cost_weight},
    ]

    inv_Jb_S = Jsat_inv @ S_x[4:7, 4:7]
    inv_Jb_c = Jsat_inv @ c_x[4:7, 0:1]
    state_bounds = [
        {'kind': 'norm', 'slice': slice(4, 7), 'norm_type': 'norm2',
         'limit': omegab_max, 'S': inv_Jb_S, 'c': inv_Jb_c},
    ]
    control_bounds = [
        {'kind': 'box', 'slice': slice(0, 4), 'lower': 0.0, 'upper': fth_max},
        {'kind': 'box_affine', 'c_affine': np.array([[1.0, 1.0, 0.0, 0.0]]),
         'lower': 0.0, 'upper': fth_max},
        {'kind': 'box_affine', 'c_affine': np.array([[0.0, 0.0, 1.0, 1.0]]),
         'lower': 0.0, 'upper': fth_max},
    ]

    return dict(
        T=T, tf_init=tf_init, tf_min=tf_min, tf_max=tf_max, alpha_p=alpha_p,
        J1=J1, J2=J2, J3=J3, Jxy=Jxy, Jxz=Jxz, Jyz=Jyz,
        Jsat=Jsat, Jsat_inv=Jsat_inv,
        rth_list=rth_d, fth_list=fth_d, theta_vat=theta_vat,
        dyn_params_val=dyn_params_val, B=B,
        omegab_max=omegab_max, Fth_max=fth_max,
        start_quat=start_quat_a, end_quat=end_quat_a,
        start_pos=start_pos, end_pos=end_pos,
        S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u, S_p=S_p, c_p=c_p,
        warm_x=warm_x, warm_u=warm_u,
        cost_terms=cost_terms, cost_weight=cost_weight, lamb=lamb,
        state_bounds=state_bounds, control_bounds=control_bounds,
        tau_absolute=tau_absolute,
    )


# =============================================================================
# 2. OPEN-LOOP RK4 SIMULATION -- shared by every caller (Python-solve,
#    COGU_C_caller_test.py, VAT_COGU_C_caller_test.py's per-step plant model)
# =============================================================================

def u_interpol_zoh(t, u_value, dt_control):
    """Zero-order hold of a control trajectory."""
    T_control = u_value.shape[1]
    index = int(np.floor(t / dt_control))
    index = min(max(index, 0), T_control - 1)
    return u_value[:, index]


def attitude_dynamics_open_loop(t, x, u, params):
    """Original nonlinear rigid-body attitude dynamics, x=[q(4), hb(3)]."""
    q = x[0:4]; w = x[4:7]
    I, I_inv, Bm = params["I"], params["I_inv"], params["B"]
    torque = Bm @ u
    w_dot = I_inv @ (torque - np.cross(w, I @ w))
    wx, wy, wz = w
    q0_, q1_, q2_, q3_ = q
    q_dot = 0.5 * np.array([-wx*q1_ - wy*q2_ - wz*q3_,
                            wx*q0_ + wz*q2_ - wy*q3_,
                            wy*q0_ - wz*q1_ + wx*q3_,
                            wz*q0_ + wy*q1_ - wx*q2_])
    return np.concatenate([q_dot, w_dot])


def rk4_step(f, t, x, u, dt, params):
    k1 = f(t, x, u, params)
    k2 = f(t + dt/2, x + dt/2*k1, u, params)
    k3 = f(t + dt/2, x + dt/2*k2, u, params)
    k4 = f(t + dt, x + dt*k3, u, params)
    return x + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)


def simulate_open_loop_u(x0_quat, x0_omega, u_value, tf_opt, Jsat, Jsat_inv, B,
                         end_quat, dt_state=0.1, verbose=False):
    """
    Same as simulate_open_loop but takes the initial condition directly
    (quaternion + angular RATE, not momentum) and the control trajectory --
    this is the form both C-caller scripts have on hand (COGU_SCVx_main's
    x_out is in the same physical units as the Python solve's x_opt_unscaled,
    but a caller may also want to simulate from an arbitrary perturbed state
    that never went through the solver, e.g. in NMPC).
    """
    T_control = u_value.shape[1]
    dt_control = tf_opt / T_control
    sim_time = tf_opt
    N = int(round(sim_time / dt_state)) + 1
    if verbose and abs(dt_control / dt_state - round(dt_control / dt_state)) > 1e-9:
        print(f"\n[note] dt_state={dt_state}s does not divide the optimized "
              f"dt_control={dt_control:.4f}s; hold edges land between samples.")

    params = {"I": Jsat, "I_inv": Jsat_inv, "B": B}
    x_sim = np.zeros(7)
    x_sim[0:4] = x0_quat
    x_sim[4:7] = x0_omega

    q_history = np.zeros((N, 4))
    w_history = np.zeros((N, 3))
    u_history = np.zeros((N, 4))
    t_history = np.zeros(N)

    t = 0.0
    for i in range(N):
        q_history[i] = x_sim[0:4]
        w_history[i] = x_sim[4:7]
        t_history[i] = t
        u_current = u_interpol_zoh(t, u_value, dt_control)
        u_history[i] = u_current
        if i < N - 1:
            dt = min(dt_state, tf_opt - t)
            if dt > 0:
                x_sim = rk4_step(attitude_dynamics_open_loop, t, x_sim, u_current, dt, params)
                q_norm = np.linalg.norm(x_sim[0:4])
                if q_norm > 0:
                    x_sim[0:4] /= q_norm
                t += dt

    q_final_ol = q_history[-1]
    w_final_ol = w_history[-1]
    q_target = np.asarray(end_quat).reshape(-1)
    q_error_norm = np.linalg.norm(q_final_ol - q_target)
    w_error_norm = np.linalg.norm(w_final_ol)

    if verbose:
        print("\n" + "=" * 68)
        print("OPEN-LOOP VALIDATION  (RK4, ZOH, original nonlinear dynamics)")
        print("=" * 68)
        print(f"Integrated over:     {tf_opt:.4f} s (the optimized final time)")
        print(f"Final quaternion OL: {q_final_ol.round(8)}")
        print(f"Target quaternion:   {q_target}")
        print(f"Max ||omega|| OL:    {np.max(np.linalg.norm(w_history, axis=1)):.8f} rad/s")
        print("\n--- OPEN-LOOP ERRORS ---")
        print(f"||q_final - q_target||: {q_error_norm:.8e}")
        print(f"||omega_final||:        {w_error_norm:.8e} rad/s")

    return dict(t=t_history, q=q_history, w=w_history, u=u_history,
               q_error_norm=q_error_norm, w_error_norm=w_error_norm)


# =============================================================================
# 3. SYMPY MODEL + SOLVE + OPEN-LOOP VALIDATION + PLOTS
# =============================================================================

def _build_symbolic_model():
    """The VAT rigid-body attitude model. The dynamics' functional FORM is
    not something Monte Carlo varies (only the numeric dyn_params_val is) --
    so, unlike every numeric quantity above, this is not a parameter."""
    import sympy as sp

    b1th, b2th, b3th, b4th = sp.symbols('b1th b2th b3th b4th')
    q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
    hbx, hby, hbz = sp.symbols('hbx hby hbz')
    Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
    Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')

    rth1x, rth1y, rth1z = sp.symbols('rth1x rth1y rth1z')
    rth2x, rth2y, rth2z = sp.symbols('rth2x rth2y rth2z')
    rth3x, rth3y, rth3z = sp.symbols('rth3x rth3y rth3z')
    rth4x, rth4y, rth4z = sp.symbols('rth4x rth4y rth4z')
    rth1 = sp.Matrix([rth1x, rth1y, rth1z])
    rth2 = sp.Matrix([rth2x, rth2y, rth2z])
    rth3 = sp.Matrix([rth3x, rth3y, rth3z])
    rth4 = sp.Matrix([rth4x, rth4y, rth4z])

    fth1x, fth1y, fth1z = sp.symbols('fth1x fth1y fth1z')
    fth2x, fth2y, fth2z = sp.symbols('fth2x fth2y fth2z')
    fth3x, fth3y, fth3z = sp.symbols('fth3x fth3y fth3z')
    fth4x, fth4y, fth4z = sp.symbols('fth4x fth4y fth4z')
    fth1 = sp.Matrix([fth1x, fth1y, fth1z])
    fth2 = sp.Matrix([fth2x, fth2y, fth2z])
    fth3 = sp.Matrix([fth3x, fth3y, fth3z])
    fth4 = sp.Matrix([fth4x, fth4y, fth4z])

    J = sp.Matrix([[Ixx, Ixy, Ixz],
                   [Ixy, Iyy, Iyz],
                   [Ixz, Iyz, Izz]])
    hb_sympy = sp.Matrix([hbx, hby, hbz])

    f_quaternion = 0.5 * sp.Matrix([[-q1, -q2, -q3],
                                    [q0, -q3, q2],
                                    [q3, q0, -q1],
                                    [-q2, q1, q0]]) * (J.inv() * hb_sympy)

    torque_th = (b1th * rth1.cross(fth1) + b2th * rth2.cross(fth2)
                + b3th * rth3.cross(fth3) + b4th * rth4.cross(fth4))

    f_rigid_body = sp.simplify(-(J.inv() * hb_sympy).cross(hb_sympy) + torque_th)

    # NOTE: absolute-time dynamics dx/dt = f(x,u); solve_trajectory builds
    # p*f itself under a dilation mode.
    dynamics = [f_quaternion[0], f_quaternion[1], f_quaternion[2], f_quaternion[3],
               f_rigid_body[0], f_rigid_body[1], f_rigid_body[2]]

    states = [q0, q1, q2, q3, hbx, hby, hbz]
    controls = [b1th, b2th, b3th, b4th]
    dyn_params_sym = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz,
                      rth1x, rth1y, rth1z, rth2x, rth2y, rth2z,
                      rth3x, rth3y, rth3z, rth4x, rth4y, rth4z,
                      fth1x, fth1y, fth1z, fth2x, fth2y, fth2z,
                      fth3x, fth3y, fth3z, fth4x, fth4y, fth4z]
    return states, controls, dynamics, dyn_params_sym


def run_birds_vat_dilation(
    problem_kwargs=None,
    solver_kwargs=None,
    generate_c=True, c_output_dir='VAT_c_output',
    dt_state=0.1,
    plot=True, verbose=True,
):
    """
    Solve the BIRDS-VAT time-dilation problem via scvx_toolbox (SymPy model
    + Python SCvx), validate it open-loop with independent RK4 integration,
    and (optionally) plot -- i.e. today's script, fully parameterized.

    Parameters
    ----------
    problem_kwargs : dict, optional
        Passed to `build_vat_problem_data(...)` -- every physical/numeric
        knob (inertia, thruster geometry, bounds, boundary conditions, tf
        settings, warm start, cost weight, ...).
    solver_kwargs : dict, optional
        Passed through to `scvx_toolbox.pipeline.solve_trajectory(...)` and
        overrides its defaults -- SCvx tuning (size_N, max_iter, rho*, eta*,
        beta*, e_tol, epsilon_stop, defect_penalty_space, solver, ...).
    generate_c, c_output_dir : as in `solve_trajectory`.
    dt_state : float
        Integration step of the independent open-loop RK4 validation [s].
    plot : bool
        Show the matplotlib figures.
    verbose : bool
        Print progress and the results summary.

    Returns
    -------
    (result, data) : the scvx_toolbox result dict and the build_vat_problem_data dict.
    """
    from scvx_toolbox.pipeline import solve_trajectory

    problem_kwargs = dict(problem_kwargs or {})
    solver_kwargs = dict(solver_kwargs or {})

    data = build_vat_problem_data(**problem_kwargs)
    T = data['T']
    tf_init, tf_min, tf_max, alpha_p = data['tf_init'], data['tf_min'], data['tf_max'], data['alpha_p']
    Jsat, Jsat_inv, B = data['Jsat'], data['Jsat_inv'], data['B']
    omegab_max, Fth_max = data['omegab_max'], data['Fth_max']

    states, controls, dynamics, dyn_params_sym = _build_symbolic_model()

    if verbose:
        print("=" * 68)
        print("BIRDS-TM (VAT only) -- ZOH + TIME DILATION")
        print("=" * 68)
        print(f"tf is free in [{tf_min}, {tf_max}] s, initial guess {tf_init} s, "
              f"alpha_p = {alpha_p}")
        print()

    t0 = time.time()

    solve_defaults = dict(
        size_N=20, lamb=data['lamb'], max_iter=25,
        rho0=0.0, rho1=0.1, rho2=0.7,
        eta0=0.001, eta1=10.0, eta_init=1.0,
        beta_sh=2.0, beta_gr=2.0, e_tol=0.005,
        epsilon_stop=1e-3, defect_penalty_space='scaled',
        solver='ECOS', verbose=verbose,
    )
    solve_defaults.update(solver_kwargs)

    result = solve_trajectory(
        states=states, controls=controls, dynamics=dynamics,
        start=data['start_pos'], end=data['end_pos'], T=T,
        time_mode='ZOH+dilation',
        tf_init=tf_init, tf_min=tf_min, tf_max=tf_max, alpha_p=alpha_p,
        dynamic_parameters_sym=dyn_params_sym,
        dynamic_parameters_val=data['dyn_params_val'],
        S_x=data['S_x'], c_x=data['c_x'], S_u=data['S_u'], c_u=data['c_u'],
        S_p=data['S_p'], c_p=data['c_p'],
        warm_start_x=data['warm_x'], warm_start_u=data['warm_u'],
        cost_terms=data['cost_terms'],
        state_bounds=data['state_bounds'], control_bounds=data['control_bounds'],
        generate_c=generate_c, c_output_dir=c_output_dir,
        **solve_defaults,
    )

    total = time.time() - t0

    x_value = result['x_opt_unscaled']
    u_value = result['u_opt_unscaled']
    tf_opt = result['tf']
    tau_val = tf_opt / T

    if verbose:
        print("\n" + "=" * 68)
        print("RESULTS")
        print("=" * 68)
        print(f"Wall time:        {total:.1f} s")
        print(f"Iterations:       {result['iterations']}")
        print(f"Converged:        {result['converged']}   (status: {result['status']})")
        print(f"Optimal tf:       {tf_opt:.4f} s     (initial guess {tf_init} s, "
              f"bounds [{tf_min}, {tf_max}])")
        print(f"Absolute step:    {tau_val:.4f} s")
        print(f"Quat start:       {x_value[0:4, 0].round(6)}")
        print(f"Quat final:       {x_value[0:4, -1].round(6)}   target "
              f"{data['end_quat'][:, 0]}")
        print(f"hb final:         {x_value[4:7, -1].round(9)}")

    omega_guidance = Jsat_inv @ x_value[4:7, :]

    # ---- Open-loop validation: RK4 of the original nonlinear dynamics ----
    ol = simulate_open_loop_u(
        x0_quat=x_value[0:4, 0], x0_omega=Jsat_inv @ x_value[4:7, 0],
        u_value=u_value, tf_opt=tf_opt, Jsat=Jsat, Jsat_inv=Jsat_inv, B=B,
        end_quat=data['end_quat'][:, 0], dt_state=dt_state, verbose=verbose,
    )
    ol['omega_guidance'] = omega_guidance
    result['open_loop'] = ol
    result['B'] = B

    if plot:
        _plot_birds_vat(result, data, ol)

    return result, data


def _plot_birds_vat(result, data, ol):
    import matplotlib.pyplot as plt

    x_value = result['x_opt_unscaled']
    u_value = result['u_opt_unscaled']
    tf_opt = result['tf']
    T = data['T']
    tau_val = tf_opt / T
    t_nodes = np.arange(T + 1) * tau_val
    t_history, q_history, w_history, u_history = ol['t'], ol['q'], ol['w'], ol['u']
    omega_guidance = ol['omega_guidance']
    omegab_max = data['omegab_max']

    plt.figure(figsize=(7, 4.5))
    plt.plot(t_history, q_history[:, 0], label='q0')
    plt.plot(t_history, q_history[:, 1], label='q1')
    plt.plot(t_history, q_history[:, 2], label='q2')
    plt.plot(t_history, q_history[:, 3], label='q3')
    plt.plot(t_history, np.linalg.norm(q_history, axis=1), '--', label=r'$\|q\|$')
    plt.plot(t_nodes, x_value[0:4, :].T, 'o', ms=3, alpha=0.45)
    plt.xlabel('Time [s]'); plt.ylabel('Quaternion')
    plt.title(f'Quaternion: Open-Loop vs SCvx   (optimized $t_f$ = {tf_opt:.2f} s)')
    plt.legend(); plt.grid(True); plt.tight_layout()

    plt.figure(figsize=(7, 4.5))
    plt.plot(t_history, w_history[:, 0] * 180/np.pi, label=r'$\omega_x$')
    plt.plot(t_history, w_history[:, 1] * 180/np.pi, label=r'$\omega_y$')
    plt.plot(t_history, w_history[:, 2] * 180/np.pi, label=r'$\omega_z$')
    plt.plot(t_history, np.linalg.norm(w_history, axis=1) * 180/np.pi, '--', label=r'$\|\omega\|$')
    plt.plot(t_nodes, omega_guidance.T * 180/np.pi, 'o', ms=3, alpha=0.45)
    plt.axhline(omegab_max * 180/np.pi, linestyle='--', alpha=0.7, label=r'$\omega_{\max}$')
    plt.axhline(-omegab_max * 180/np.pi, linestyle='--', alpha=0.7)
    plt.xlabel('Time [s]'); plt.ylabel('Angular velocity [deg/s]')
    plt.title('Angular Velocity: Open-Loop vs SCvx')
    plt.legend(); plt.grid(True); plt.tight_layout()

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(7, 7))
    for i in range(4):
        axs[i].step(t_history, u_history[:, i] * 1e6, where='post', label=f'VAT{i+1}')
        axs[i].set_ylabel(f'VAT{i+1} [uN]'); axs[i].grid(True); axs[i].legend(loc='upper right')
    axs[-1].set_xlabel('Time [s]')
    fig.suptitle(f'VAT Thrust Sequence -- SCvx ZOH + dilation ($t_f$ = {tf_opt:.2f} s)')
    fig.tight_layout()

    hb_history = (data['Jsat'] @ w_history.T).T
    plt.figure(figsize=(7, 4.5))
    plt.plot(t_history, hb_history[:, 0], label=r'$h_x$')
    plt.plot(t_history, hb_history[:, 1], label=r'$h_y$')
    plt.plot(t_history, hb_history[:, 2], label=r'$h_z$')
    plt.plot(t_nodes, x_value[4:7, :].T, 'o', ms=3, alpha=0.45)
    plt.xlabel('Time [s]'); plt.ylabel('Body angular momentum [N m s]')
    plt.title('Body Angular Momentum: Open-Loop vs SCvx')
    plt.legend(); plt.grid(True); plt.tight_layout()

    plt.figure(figsize=(7, 4))
    iters = [h['iteration'] for h in result['history']]
    tfs = [h['tf'] for h in result['history']]
    plt.plot(iters, tfs, 'o-')
    plt.axhline(data['tf_init'], linestyle=':', alpha=0.7, label=f"initial guess {data['tf_init']} s")
    plt.axhline(tf_opt, linestyle='--', alpha=0.7, label=f'optimal {tf_opt:.2f} s')
    plt.xlabel('SCvx iteration'); plt.ylabel(r'$t_f$ [s]')
    plt.title('Free Final Time Convergence')
    plt.legend(); plt.grid(True); plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    run_birds_vat_dilation()
