"""
The SCvx loop: successive convexification with a trust region and the
classic rho acceptance test.

Each iteration
    1. unscales the reference trajectory (x_bar, u_bar, p_bar),
    2. discretizes the linearized dynamics around it (RK4, ZOH or FOH,
       with or without time dilation),
    3. solves the convex subproblem, obtaining the linearized cost L,
    4. evaluates the true nonlinear augmented cost J at the new iterate,
    5. compares the predicted and actual improvements,
           rho = (J_bar - J) / (J_bar - L)
       and grows, shrinks, or rejects the trust region accordingly.
"""

import time

import numpy as np

from .cost_dsl import coefficient_map, node_count, resolve_coefficient, resolve_term
from .discretization import (ScalingBundle, discretize_constraints,
                             discretize_interval, rk4_flow_map)
from .modes import FOH
from .scaling import scale_p, scale_u, scale_x, unscale_p, unscale_u, unscale_x


# ----------------------------------------------------------------------
# Nonlinear augmented cost J
# ----------------------------------------------------------------------

def eval_user_cost(cost_terms, x, u, T, u_nodes, tau, lamb):
    """
    Numpy mirror of build_convex_cost -- same DSL, same numbers.

    `x` and `u` are scaled trajectories.
    """
    sqrt_tau = float(tau) ** 0.5
    coeff_map = coefficient_map(sqrt_tau, tau, float(tau) * float(lamb))
    var_map = {'x': x, 'u': u}

    cost = 0.0
    for term in cost_terms:
        rt = resolve_term(term)
        arr = var_map[rt.var]
        coeff = resolve_coefficient(rt, coeff_map)
        K = node_count(rt.k_range, T, u_nodes)
        for k in range(K):
            vec = arr[rt.slc, k:k + 1]
            expr = coeff * vec if rt.offset is None else coeff * (vec - rt.offset)
            if rt.kind == 'sumsq':
                cost += rt.weight * float(np.linalg.norm(expr)) ** 2
            elif rt.kind == 'norm1':
                cost += rt.weight * float(np.linalg.norm(expr, ord=1))
            elif rt.kind == 'norm2':
                cost += rt.weight * float(np.linalg.norm(expr, ord=2))
    return cost


def compute_J(x, u, p, T, u_nodes, tau, lamb, dyn_par, ng, funcs, sc,
              cost_terms, hold, dilation, alpha_p, defect_substeps=1,
              defect_penalty_space='unscaled'):
    """
    Augmented nonlinear cost: user cost + defect penalty + constraint
    violation penalty + terminal cost on the final time.

    `x`, `u`, `p` are scaled.

    `defect_penalty_space` picks the units the defect penalty is measured in:

        'unscaled'  physical units.  This is what the reference notebooks do,
                    and it is the default so that results reproduce.  Note
                    that the convex subproblem penalizes the slack vc in
                    *scaled* units, so with an aggressive S_x the two are not
                    on the same footing and rho can come out negative even
                    though the step improves J.
        'scaled'    S_x^-1 applied to the defect, which makes this penalty
                    numerically identical to the vc penalty in L.  Use this
                    if rho is behaving strangely.
    """
    cost = eval_user_cost(cost_terms, x, u, T, u_nodes, tau, lamb)

    scaled_defect = (defect_penalty_space == 'scaled')

    f_func = funcs['f']
    for k in range(T):
        u_next = u[:, k + 1:k + 2] if (hold == FOH and k + 1 < u_nodes) else None
        flow = rk4_flow_map(x[:, k:k + 1], u[:, k:k + 1], u_next, p,
                            k * tau, tau, dyn_par, f_func, sc,
                            hold=hold, dilation=dilation,
                            substeps=defect_substeps)
        defect = unscale_x(x[:, k + 1:k + 2], sc.S_x, sc.c_x) - flow
        if scaled_defect:
            defect = sc.inv_S_x @ defect
        cost += tau * float(np.linalg.norm(lamb * defect, ord=1))

    if ng > 0 and funcs.get('g') is not None:
        hp = unscale_p(p, sc.S_p, sc.c_p) if dilation else np.zeros((0, 1))
        for k in range(T + 1):
            uk = min(k, u_nodes - 1)
            hx = unscale_x(x[:, k:k + 1], sc.S_x, sc.c_x)
            hu = unscale_u(u[:, uk:uk + 1], sc.S_u, sc.c_u)
            g_val = funcs['g'](hx, hu, hp, k * tau, dyn_par)
            for i in range(g_val.shape[0]):
                cost += tau * abs(lamb * max(float(g_val[i, 0]), 0.0))

    if dilation and float(alpha_p) != 0.0:
        tf_val = float(unscale_p(p, sc.S_p, sc.c_p)[0, 0])
        cost += float(alpha_p) * tf_val

    return cost


# ----------------------------------------------------------------------
# SCvx
# ----------------------------------------------------------------------

def solve_scvx(prob_dict, funcs, config):
    """
    Run SCvx on a problem built by build_problem().

    `config` keys
        tf, size_N, dynamic_parameters, S_x, c_x, S_u, c_u, [S_p, c_p]
        start_pos, end_pos, warm_start_x, warm_start_u, [warm_start_tf]
        rho0, rho1, rho2, eta0, eta1, eta_init, beta_sh, beta_gr
        lamb, e_tol, epsilon_stop, max_iter, solver, verbose
        [alpha_p], defect_substeps, solver_options

    Returns a dict with the scaled and unscaled optimal trajectories, the
    optimal final time, convergence status and the per-iteration history.
    """
    nx = prob_dict['nx']
    nu = prob_dict['nu']
    ng = prob_dict['ng']
    T = prob_dict['T']
    u_nodes = prob_dict['u_nodes']
    hold = prob_dict['hold']
    dilation = prob_dict['dilation']
    cost_terms = prob_dict['cost_terms']

    size_N = config['size_N']
    dyn_par = config['dynamic_parameters']
    lamb_val = float(config['lamb'])
    alpha_p = float(config.get('alpha_p', 0.0))
    defect_substeps = int(config.get('defect_substeps', 1))
    defect_space = config.get('defect_penalty_space', 'unscaled')
    max_stall = int(config.get('max_stall', 3))

    sc = ScalingBundle(config['S_x'], config['c_x'],
                       config['S_u'], config['c_u'],
                       config.get('S_p'), config.get('c_p'))

    # Step size: normalized under dilation, absolute otherwise.
    tau_val = 1.0 / T if dilation else float(config['tf']) / T

    rho0, rho1, rho2 = config['rho0'], config['rho1'], config['rho2']
    eta0, eta1 = config['eta0'], config['eta1']
    eta_val = config['eta_init']
    beta_sh, beta_gr = config['beta_sh'], config['beta_gr']
    e_tol = config['e_tol']
    epsilon_stop = config['epsilon_stop']
    max_iter = config['max_iter']
    solver_name = config.get('solver', 'ECOS')
    solver_options = dict(config.get('solver_options') or {})
    verbose = config.get('verbose', True)

    problem = prob_dict['problem']
    variables = prob_dict['variables']

    # ------------------------------------------------------------------
    # Reference trajectory (scaled)
    # ------------------------------------------------------------------
    warm_x = np.asarray(config['warm_start_x'], dtype=float)
    warm_u = np.asarray(config['warm_start_u'], dtype=float)
    if warm_u.shape[1] < u_nodes:            # FOH needs one extra control node
        pad = np.repeat(warm_u[:, -1:], u_nodes - warm_u.shape[1], axis=1)
        warm_u = np.hstack([warm_u, pad])

    ox = np.zeros((nx, T + 1))
    ou = np.zeros((nu, u_nodes))
    for k in range(T + 1):
        ox[:, k:k + 1] = scale_x(warm_x[:, k:k + 1], sc.inv_S_x, sc.c_x)
    for k in range(u_nodes):
        ou[:, k:k + 1] = scale_u(warm_u[:, k:k + 1], sc.inv_S_u, sc.c_u)

    if dilation:
        tf_guess = np.array([[float(config['warm_start_tf'])]])
        op = scale_p(tf_guess, sc.inv_S_p, sc.c_p)
    else:
        op = np.zeros((0, 1))

    # ------------------------------------------------------------------
    # Fixed parameter values
    # ------------------------------------------------------------------
    prob_dict['eta'].value = eta_val
    prob_dict['start_pos'].value = scale_x(config['start_pos'], sc.inv_S_x, sc.c_x)
    prob_dict['end_pos'].value = scale_x(config['end_pos'], sc.inv_S_x, sc.c_x)

    def push_reference():
        for k in range(T + 1):
            prob_dict['ox_params'][k].value = ox[:, k:k + 1]
        for k in range(u_nodes):
            prob_dict['ou_params'][k].value = ou[:, k:k + 1]
        if dilation:
            prob_dict['op_param'].value = op

    push_reference()

    # ------------------------------------------------------------------
    # Loop
    # ------------------------------------------------------------------
    history = []
    converged = False
    past_startup = False
    rho_i = None
    final_iter = 0
    stalled_at_floor = 0
    status = 'max_iter'
    t_start = time.time()

    for iteration in range(1, max_iter + 1):
        final_iter = iteration

        # 1. Unscale the reference
        ohx = np.zeros((nx, T + 1))
        ohu = np.zeros((nu, u_nodes))
        for k in range(T + 1):
            ohx[:, k:k + 1] = unscale_x(ox[:, k:k + 1], sc.S_x, sc.c_x)
        for k in range(u_nodes):
            ohu[:, k:k + 1] = unscale_u(ou[:, k:k + 1], sc.S_u, sc.c_u)
        ohp = unscale_p(op, sc.S_p, sc.c_p) if dilation else None

        # 2. Discretize
        t_disc = time.time()
        for k in range(T):
            k_next = min(k + 1, u_nodes - 1)
            A_k, Bm_k, Bp_k, F_k, y_k = discretize_interval(
                funcs,
                ohx[:, k:k + 1], ohx[:, k + 1:k + 2],
                ohu[:, k:k + 1], ohu[:, k_next:k_next + 1],
                ohp, k * tau_val, tau_val, size_N, dyn_par, sc,
                hold=hold, dilation=dilation)

            prob_dict['A_params'][k].value = A_k
            prob_dict['Bm_params'][k].value = Bm_k
            prob_dict['y_params'][k].value = y_k
            if hold == FOH:
                prob_dict['Bp_params'][k].value = Bp_k
            if dilation:
                prob_dict['F_params'][k].value = F_k

        if ng > 0:
            for k in range(T + 1):
                uk = min(k, u_nodes - 1)
                C_k, D_k, G_k, z_k = discretize_constraints(
                    funcs, ohx[:, k:k + 1], ohu[:, uk:uk + 1], ohp,
                    k * tau_val, dyn_par, sc, dilation)
                prob_dict['C_params'][k].value = C_k
                prob_dict['D_params'][k].value = D_k
                prob_dict['z_params'][k].value = z_k
                if dilation:
                    prob_dict['G_params'][k].value = G_k
        disc_ms = (time.time() - t_disc) * 1000.0

        # 3. Solve the convex subproblem
        t_solve = time.time()
        L_opt = problem.solve(solver=solver_name, ignore_dpp=True, **solver_options)
        solve_ms = (time.time() - t_solve) * 1000.0

        if variables['x'].value is None:
            # Losing a long run to an exception helps nobody: keep the last
            # accepted iterate and report why the loop stopped.
            status = 'subproblem_' + str(problem.status)
            if verbose:
                print(f"SCvx stopped at iteration {iteration}: the convex "
                      f"subproblem returned no solution (status: "
                      f"{problem.status}).")
                if 'infeasible' in str(problem.status).lower():
                    print("  A shrinking trust region usually goes infeasible "
                          "because some term in it is a constant the trust "
                          "region cannot influence. Check that the warm start "
                          "agrees with `start`/`end`, and keep "
                          "pin_aware_trust_region=True.")
            final_iter = iteration - 1
            break

        x_opt = variables['x'].value.copy()
        u_opt = variables['u'].value.copy()
        p_opt = variables['p'].value.copy() if dilation else np.zeros((0, 1))

        # The convex objective drops the constant alpha_p*c_p; add it back so
        # that L and J are directly comparable.
        if dilation and alpha_p != 0.0:
            L_opt = L_opt + alpha_p * float(np.asarray(sc.c_p).reshape(-1)[0])

        # 4. Nonlinear cost at the new iterate and at the reference
        J_opt = compute_J(x_opt, u_opt, p_opt, T, u_nodes, tau_val, lamb_val,
                          dyn_par, ng, funcs, sc, cost_terms, hold, dilation,
                          alpha_p, defect_substeps, defect_space)
        J_ref = compute_J(ox, ou, op, T, u_nodes, tau_val, lamb_val,
                          dyn_par, ng, funcs, sc, cost_terms, hold, dilation,
                          alpha_p, defect_substeps, defect_space)

        norm_diff_2 = float(np.max(np.linalg.norm(x_opt - ox, ord=2, axis=0)))
        norm_diff_1 = float(np.max(np.linalg.norm(x_opt - ox, ord=1, axis=0)))
        p_diff = float(np.linalg.norm(p_opt - op, ord=1)) if dilation else 0.0
        tf_now = (float(unscale_p(p_opt, sc.S_p, sc.c_p)[0, 0]) if dilation
                  else float(config['tf']))

        if verbose:
            rho_txt = '   None' if rho_i is None else f'{rho_i:7.4f}'
            tf_txt = f' | tf={tf_now:8.3f}s' if dilation else ''
            print(f"Iter {iteration:3d} | eta={eta_val:.2e} | rho={rho_txt}"
                  f" | L={L_opt:12.4f} | J={J_opt:12.4f} | Jref={J_ref:12.4f}"
                  f" | dx={norm_diff_2:.3e}{tf_txt}"
                  f" | disc={disc_ms:6.1f}ms solve={solve_ms:6.1f}ms")

        history.append({
            'iteration': iteration, 'eta': eta_val, 'rho': rho_i,
            'L_SCVx': L_opt, 'J_SCVx': J_opt, 'oJ_SCVx': J_ref,
            'norm_diff': norm_diff_2, 'tf': tf_now,
            'discretization_ms': disc_ms, 'solve_ms': solve_ms,
        })

        # 5. Convergence
        delta_J = J_ref - J_opt
        delta_L = J_ref - L_opt

        # Convergence uses |delta_L|, not delta_L.  A rejected step has
        # delta_L << 0, which would satisfy a signed test trivially and
        # declare convergence precisely when the linearization is at its
        # least trustworthy.
        step_is_tiny = norm_diff_1 < epsilon_stop and p_diff < epsilon_stop
        if past_startup and (abs(delta_L) < e_tol * abs(J_ref) or step_is_tiny):
            converged = True
            status = 'converged'
            if verbose:
                print("SCvx converged.")
            break

        # 6. Trust region update
        rho_i = delta_J / delta_L if delta_L != 0 else 0.0

        if rho_i < rho0:
            # Step rejected.  If the trust region is already on its floor,
            # shrinking it again changes nothing and the next iteration will
            # reproduce this one exactly -- bail out instead of spinning.
            if eta_val <= eta0 * (1.0 + 1e-12):
                stalled_at_floor += 1
                if stalled_at_floor >= max_stall:
                    status = 'stalled'
                    if verbose:
                        print(f"SCvx stalled: the step is still rejected with "
                              f"the trust region at its floor (eta0={eta0:g}). "
                              f"Stopping at iteration {iteration}.")
                    break
            eta_val = max(eta0, eta_val / beta_sh)
        else:
            stalled_at_floor = 0
            if rho_i < rho1:
                eta_val = max(eta0, eta_val / beta_sh)
            elif rho_i >= rho2:
                eta_val = min(eta1, beta_gr * eta_val)
            ox = x_opt.copy()
            ou = u_opt.copy()
            if dilation:
                op = p_opt.copy()

        prob_dict['eta'].value = eta_val
        push_reference()

        if iteration >= 3:
            past_startup = True

    # ------------------------------------------------------------------
    # Unscale the accepted solution
    # ------------------------------------------------------------------
    x_final = np.zeros((nx, T + 1))
    u_final = np.zeros((nu, u_nodes))
    for k in range(T + 1):
        x_final[:, k:k + 1] = unscale_x(ox[:, k:k + 1], sc.S_x, sc.c_x)
    for k in range(u_nodes):
        u_final[:, k:k + 1] = unscale_u(ou[:, k:k + 1], sc.S_u, sc.c_u)

    tf_opt = (float(unscale_p(op, sc.S_p, sc.c_p)[0, 0]) if dilation
              else float(config['tf']))

    return {
        'x_opt': ox, 'u_opt': ou, 'p_opt': op,
        'x_opt_unscaled': x_final, 'u_opt_unscaled': u_final,
        'tf': tf_opt,
        'tau': tf_opt / T,                       # absolute step of the solution
        'tau_normalized': tau_val,
        't_nodes': np.linspace(0.0, tf_opt, T + 1),
        'converged': converged, 'iterations': final_iter,
        'status': status,
        'history': history,
        'time_mode': prob_dict['time_mode'],
        'hold': hold, 'dilation': dilation,
        'wall_time': time.time() - t_start,
    }
