import numpy as np
from .scaling import (scale_x, scale_u, inv_scale_x, inv_scale_u,
                      scale_A, scale_B, scale_y, scale_C, scale_D, scale_z)
from .cost_dsl import resolve_term


# ==================================================
# DISCRETIZACION RK4 DE PHI (Cell 40 del template)
# ==================================================

def discretize_ABy(ohx, ohu, t0, dyn_par, tau, size_N,
                   f_A, f_B, f_y,
                   S_x, inv_S_x, c_x, S_u, inv_S_u, c_u):
    """
    RK4 integration of the state transition matrix Phi to get discrete A, B, y.
    ZOH assumption: linearization point (ohx, ohu) held constant over the interval.

    ohx, ohu: unscaled reference state/control at node k (column vectors)
    t0: time at start of interval
    tau: step time between nodes
    size_N: number of RK4 sub-steps
    f_A, f_B, f_y: code-gen functions (operate in unscaled space)
    S_x, inv_S_x, c_x, S_u, inv_S_u, c_u: scaling matrices/vectors

    Returns: (A_disc, B_disc, y_disc) — scaled and discrete
    """

    n = ohx.shape[0]
    m = ohu.shape[0]

    Phi = np.eye(n)
    G = np.zeros((n, m))
    zz = np.zeros((n, 1))

    dt = tau / (size_N - 1)
    t = t0

    for _ in range(size_N - 1):

        # ===== k1 =====
        A1_no = f_A(ohx, ohu, t, dyn_par)
        B1_no = f_B(ohx, ohu, t, dyn_par)
        y1_no = f_y(ohx, ohu, t, dyn_par)

        A1 = scale_A(A1_no, inv_S_x, S_x)
        B1 = scale_B(B1_no, inv_S_x, S_u)
        y1 = scale_y(A1_no, B1_no, y1_no, inv_S_x, c_x, c_u)

        k1_Phi = A1 @ Phi
        k1_G = A1 @ G + B1
        k1_z = A1 @ zz + y1

        # ===== k2 =====
        t2 = t + 0.5 * dt
        Phi2 = Phi + 0.5 * dt * k1_Phi
        G2 = G + 0.5 * dt * k1_G
        z2 = zz + 0.5 * dt * k1_z

        A2_no = f_A(ohx, ohu, t2, dyn_par)
        B2_no = f_B(ohx, ohu, t2, dyn_par)
        y2_no = f_y(ohx, ohu, t2, dyn_par)

        A2 = scale_A(A2_no, inv_S_x, S_x)
        B2 = scale_B(B2_no, inv_S_x, S_u)
        y2 = scale_y(A2_no, B2_no, y2_no, inv_S_x, c_x, c_u)

        k2_Phi = A2 @ Phi2
        k2_G = A2 @ G2 + B2
        k2_z = A2 @ z2 + y2

        # ===== k3 =====
        Phi3 = Phi + 0.5 * dt * k2_Phi
        G3 = G + 0.5 * dt * k2_G
        z3 = zz + 0.5 * dt * k2_z

        A3_no = f_A(ohx, ohu, t2, dyn_par)
        B3_no = f_B(ohx, ohu, t2, dyn_par)
        y3_no = f_y(ohx, ohu, t2, dyn_par)

        A3 = scale_A(A3_no, inv_S_x, S_x)
        B3 = scale_B(B3_no, inv_S_x, S_u)
        y3 = scale_y(A3_no, B3_no, y3_no, inv_S_x, c_x, c_u)

        k3_Phi = A3 @ Phi3
        k3_G = A3 @ G3 + B3
        k3_z = A3 @ z3 + y3

        # ===== k4 =====
        t4 = t + dt
        Phi4 = Phi + dt * k3_Phi
        G4 = G + dt * k3_G
        z4 = zz + dt * k3_z

        A4_no = f_A(ohx, ohu, t4, dyn_par)
        B4_no = f_B(ohx, ohu, t4, dyn_par)
        y4_no = f_y(ohx, ohu, t4, dyn_par)

        A4 = scale_A(A4_no, inv_S_x, S_x)
        B4 = scale_B(B4_no, inv_S_x, S_u)
        y4 = scale_y(A4_no, B4_no, y4_no, inv_S_x, c_x, c_u)

        k4_Phi = A4 @ Phi4
        k4_G = A4 @ G4 + B4
        k4_z = A4 @ z4 + y4

        # ===== Update =====
        Phi += (dt / 6) * (k1_Phi + 2 * k2_Phi + 2 * k3_Phi + k4_Phi)
        G += (dt / 6) * (k1_G + 2 * k2_G + 2 * k3_G + k4_G)
        zz += (dt / 6) * (k1_z + 2 * k2_z + 2 * k3_z + k4_z)

        t += dt

    return Phi, G, zz


# ==================================================
# DISCRETIZACION C, D, z (bypass ZOH)
# ==================================================

def discretize_CDz(ohx, ohu, t0, dyn_par,
                   f_C, f_D, f_z,
                   S_x, c_x, S_u, c_u):
    """
    Bypass — constraints don't need RK4 discretization.
    Just evaluate and scale C, D, z at the linearization point.
    """

    C_no = f_C(ohx, ohu, t0, dyn_par)
    D_no = f_D(ohx, ohu, t0, dyn_par)
    z_no = f_z(ohx, ohu, t0, dyn_par)

    return (scale_C(C_no, S_x),
            scale_D(D_no, S_u),
            scale_z(C_no, D_no, z_no, c_x, c_u))


# ==================================================
# RK4 STEP (para evaluar J_SCVx real)
# ==================================================

def rk4_step(xk_scaled, uk_scaled, t, dt, dyn_par, f_func,
             S_x, c_x, S_u, c_u):
    """
    Single RK4 step of the actual dynamics.
    Takes SCALED x, u. Returns UNSCALED next state.
    """

    xk = inv_scale_x(xk_scaled, S_x, c_x)
    uk = inv_scale_u(uk_scaled, S_u, c_u)

    k1 = f_func(xk, uk, t, dyn_par)
    k2 = f_func(xk + 0.5 * dt * k1, uk, t + 0.5 * dt, dyn_par)
    k3 = f_func(xk + 0.5 * dt * k2, uk, t + 0.5 * dt, dyn_par)
    k4 = f_func(xk + dt * k3, uk, t + dt, dyn_par)

    return xk + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


# ==================================================
# COSTO REAL J_SCVx (Cell 42 del template)
# ==================================================

def eval_user_cost(cost_terms, x, u, T, sqrt_tau, tau_val, tau_lamb):
    """
    Evalua el costo del usuario en numpy — espejo numerico de
    _build_cost_from_terms (problem.py), que lo construye en CVXPY.
    x, u: trayectorias ESCALADAS (numpy). Misma fuente de verdad (cost_terms).
    """
    coeff_map = {'sqrt_tau': sqrt_tau, 'tau': float(tau_val), 'tau_lamb': tau_lamb}
    var_map = {'u': u, 'x': x}
    cost = 0.0
    for term in cost_terms:
        rt = resolve_term(term)
        v_full = var_map[rt.var]
        coeff = coeff_map[rt.coeff] if isinstance(rt.coeff, str) else float(rt.coeff)
        K = T + 1 if rt.k_range == 'T+1' else T
        for k in range(K):
            vec = v_full[rt.slc, k:k+1]
            expr = coeff * vec if rt.offset is None else coeff * (vec - rt.offset)
            if rt.kind == 'sumsq':
                cost += rt.weight * np.linalg.norm(expr) ** 2
            elif rt.kind == 'norm1':
                cost += rt.weight * np.linalg.norm(expr, ord=1)
            elif rt.kind == 'norm2':
                cost += rt.weight * np.linalg.norm(expr, ord=2)
    return cost


def compute_J(x, u, T, tau, lamb, dyn_par, ng,
              f_func, g_func, S_x, c_x, S_u, c_u, tf, cost_terms=None):
    """
    Real SCVx cost evaluated with actual dynamics (not linearization).
    x, u: SCALED trajectories (numpy arrays).
    cost_terms: DSL del costo del usuario (mismo que build_problem).
    Returns scalar cost.
    """
    sqrt_tau = tau ** 0.5
    tau_lamb = tau * lamb

    # Costo del usuario (control, tracking, etc.) desde el DSL
    cost = eval_user_cost(cost_terms, x, u, T, sqrt_tau, tau, tau_lamb)

    for k in range(T):
        # Defect: difference between propagated and actual next state
        flow_map = rk4_step(x[:, k:k+1], u[:, k:k+1],
                            k / T * tf, tau, dyn_par, f_func,
                            S_x, c_x, S_u, c_u)
        x_next_unscaled = inv_scale_x(x[:, k+1:k+2], S_x, c_x)
        defect = x_next_unscaled - flow_map
        cost += tau * np.linalg.norm(lamb * defect, ord=1)

    if ng > 0 and g_func is not None:
        for k in range(T + 1):
            x_unscaled = inv_scale_x(x[:, k:k+1], S_x, c_x)
            uk_idx = min(k, T - 1)
            u_unscaled = inv_scale_u(u[:, uk_idx:uk_idx+1], S_u, c_u)
            g_val = g_func(x_unscaled, u_unscaled, k / T * tf, dyn_par)
            for i in range(g_val.shape[0]):
                cost += tau * np.abs(lamb * max(g_val[i, 0], 0.0))

    return cost


# ==================================================
# LOOP SCVx COMPLETO (Cell 43 del template)
# ==================================================

def solve_scvx(prob_dict, funcs, config):
    """
    Full SCVx loop with rho logic.

    Parametros:
        prob_dict: output de build_problem()
        funcs: dict con funciones code-gen:
            'f': f(hx, hu, t, dyn_par) -> dinamica
            'A': A(hx, hu, t, dyn_par) -> df/dx
            'B': B(hx, hu, t, dyn_par) -> df/du
            'y': y(hx, hu, t, dyn_par) -> f - Ax - Bu
            'g', 'C', 'D', 'z': igual (requeridas si ng > 0)
        config: dict con parametros numericos:
            'tf': tiempo final
            'size_N': precision discretizacion RK4
            'dynamic_parameters': array numpy
            'S_x', 'c_x', 'S_u', 'c_u': matrices/vectores de escalado
            'start_pos', 'end_pos': condiciones frontera (unscaled, nx x 1)
            'warm_start_x': trayectoria inicial x (unscaled, nx x (T+1))
            'warm_start_u': trayectoria inicial u (unscaled, nu x T)
            'rho0', 'rho1', 'rho2': umbrales rho
            'etta0', 'etta1': limites etta
            'etta_init': etta inicial
            'beta_sh', 'beta_gr': factores shrink/grow
            'lamb': peso penalizacion virtual
            'e_tol': tolerancia convergencia
            'epsilon_stop': criterio norma parada
            'max_iter': iteraciones maximas
            'solver': nombre solver CVXPY (default 'ECOS')
            'verbose': imprimir info por iteracion (default True)
            'extra_param_values': dict {nombre: valor} para cp.Parameters extra

    Retorna:
        dict con:
            'x_opt': estados optimos (scaled)
            'u_opt': controles optimos (scaled)
            'x_opt_unscaled': estados optimos (unscaled)
            'u_opt_unscaled': controles optimos (unscaled)
            'converged': bool
            'iterations': int
            'history': lista de info por iteracion
    """

    # --- Dimensiones ---
    nx = prob_dict['nx']
    nu = prob_dict['nu']
    ng = prob_dict['ng']
    T = prob_dict['T']
    cost_terms = prob_dict.get('cost_terms')

    # --- Funciones code-gen ---
    f_func = funcs['f']
    f_A = funcs['A']
    f_B = funcs['B']
    f_y = funcs['y']
    f_g = funcs.get('g')
    f_C = funcs.get('C')
    f_D = funcs.get('D')
    f_z = funcs.get('z')

    # --- Config numerica ---
    tf = config['tf']
    size_N = config['size_N']
    dyn_par = config['dynamic_parameters']
    S_x = config['S_x']
    c_x = config['c_x']
    S_u = config['S_u']
    c_u = config['c_u']
    inv_S_x = np.linalg.inv(S_x)
    inv_S_u = np.linalg.inv(S_u)

    rho0 = config['rho0']
    rho1 = config['rho1']
    rho2 = config['rho2']
    etta0 = config['etta0']
    etta1 = config['etta1']
    etta_val = config['etta_init']
    beta_sh = config['beta_sh']
    beta_gr = config['beta_gr']
    lamb_val = config['lamb']
    tau_val = tf / T
    sqrt_tau_val = tau_val ** 0.5

    e_tol = config['e_tol']
    epsilon_stop = config['epsilon_stop']
    max_iter = config['max_iter']
    solver_name = config.get('solver', 'ECOS')
    verbose = config.get('verbose', True)
    extra_vals = config.get('extra_param_values', {})

    # --- Referencias CVXPY ---
    problem = prob_dict['problem']
    vars_ = prob_dict['variables']

    # --- Buffers de discretizacion ---
    aux_A = np.zeros((nx, nx * T))
    aux_B = np.zeros((nx, nu * T))
    aux_y = np.zeros((nx, T))
    if ng > 0:
        aux_C = np.zeros((ng, nx * (T + 1)))
        aux_D = np.zeros((ng, nu * (T + 1)))
        aux_z = np.zeros((ng, T + 1))

    # --- Escalar warm start ---
    warm_x = config['warm_start_x']    # unscaled (nx, T+1)
    warm_u = config['warm_start_u']    # unscaled (nu, T)

    ox = np.zeros((nx, T + 1))
    ou = np.zeros((nu, T))
    for k in range(T + 1):
        ox[:, k:k+1] = scale_x(warm_x[:, k:k+1], inv_S_x, c_x)
    for k in range(T):
        ou[:, k:k+1] = scale_u(warm_u[:, k:k+1], inv_S_u, c_u)

    # --- Asignar parametros fijos ---
    prob_dict['etta'].value = etta_val
    prob_dict['start_pos'].value = scale_x(config['start_pos'], inv_S_x, c_x)
    prob_dict['end_pos'].value = scale_x(config['end_pos'], inv_S_x, c_x)
    for k in range(T + 1):
        prob_dict['ox_params'][k].value = ox[:, k:k+1]
        prob_dict['ou_params'][k].value = ou[:, k:k+1] if k < T else ou[:, T-1:T]

    # --- Loop SCVx ---
    history = []
    converged = False
    no_first_iterations = False
    rho_i = None
    final_iter = 0

    for iteration in range(1, max_iter + 1):
        final_iter = iteration

        # 1. Desescalar trayectoria de referencia
        ohx = np.zeros((nx, T + 1))
        ohu = np.zeros((nu, T + 1))      # T+1 para eval C,D,z en k=T
        for k in range(T + 1):
            ohx[:, k:k+1] = inv_scale_x(ox[:, k:k+1], S_x, c_x)
        for k in range(T):
            ohu[:, k:k+1] = inv_scale_u(ou[:, k:k+1], S_u, c_u)
        ohu[:, T:T+1] = ohu[:, T-1:T]    # repetir ultimo control

        # 2. Discretizar en cada nodo
        for k in range(T):
            t_k = k / T * tf
            A_k, B_k, y_k = discretize_ABy(
                ohx[:, k:k+1], ohu[:, k:k+1], t_k, dyn_par, tau_val, size_N,
                f_A, f_B, f_y,
                S_x, inv_S_x, c_x, S_u, inv_S_u, c_u)
            aux_A[:, nx*k:nx*(k+1)] = A_k
            aux_B[:, nu*k:nu*(k+1)] = B_k
            aux_y[:, k:k+1] = y_k

        if ng > 0:
            for k in range(T + 1):
                t_k = k / T * tf
                C_k, D_k, z_k = discretize_CDz(
                    ohx[:, k:k+1], ohu[:, k:k+1], t_k, dyn_par,
                    f_C, f_D, f_z,
                    S_x, c_x, S_u, c_u)
                aux_C[:, nx*k:nx*(k+1)] = C_k
                aux_D[:, nu*k:nu*(k+1)] = D_k
                aux_z[:, k:k+1] = z_k

        # 3. Asignar valores a parametros CVXPY (por paso)
        for k in range(T):
            prob_dict['A_params'][k].value = aux_A[:, nx*k:nx*(k+1)]
            prob_dict['B_params'][k].value = aux_B[:, nu*k:nu*(k+1)]
            prob_dict['y_params'][k].value = aux_y[:, k:k+1]
        if ng > 0:
            for k in range(T + 1):
                prob_dict['C_params'][k].value = aux_C[:, nx*k:nx*(k+1)]
                prob_dict['D_params'][k].value = aux_D[:, nu*k:nu*(k+1)]
                prob_dict['z_params'][k].value = aux_z[:, k:k+1]

        # 4. Resolver
        val = problem.solve(solver=solver_name, ignore_dpp=True)
        x_opt = vars_['x'].value.copy()
        u_opt = vars_['u'].value.copy()

        # 5. Evaluar costos
        L_SCVx_opt = val
        J_SCVx_opt = compute_J(x_opt, u_opt, T, tau_val, lamb_val, dyn_par, ng,
                               f_func, f_g, S_x, c_x, S_u, c_u, tf, cost_terms)
        oJ_SCVx = compute_J(ox, ou, T, tau_val, lamb_val, dyn_par, ng,
                            f_func, f_g, S_x, c_x, S_u, c_u, tf, cost_terms)

        norm_diff = np.max(np.linalg.norm(x_opt - ox, ord=2, axis=0))
        norm_diff_1 = np.max(np.linalg.norm(x_opt - ox, ord=1, axis=0))

        if verbose:
            print(f"Iter {iteration:3d} | etta={etta_val:.2e} | rho={rho_i}"
                  f" | L={L_SCVx_opt:.4f} | J={J_SCVx_opt:.4f}"
                  f" | oJ={oJ_SCVx:.4f} | dx={norm_diff:.6f}")

        history.append({
            'iteration': iteration,
            'etta': etta_val,
            'rho': rho_i,
            'L_SCVx': L_SCVx_opt,
            'J_SCVx': J_SCVx_opt,
            'oJ_SCVx': oJ_SCVx,
            'norm_diff': norm_diff,
        })

        # 6. Verificar convergencia
        Delta_J = oJ_SCVx - J_SCVx_opt
        Delta_L = oJ_SCVx - L_SCVx_opt

        if no_first_iterations and (
                Delta_L < e_tol * np.abs(oJ_SCVx)
                or norm_diff_1 < epsilon_stop):
            converged = True
            if verbose:
                print("SCVx converged!")
            break

        # 7. Logica rho
        rho_i = Delta_J / Delta_L if Delta_L != 0 else 0.0

        if rho_i < rho0:
            etta_val = max(etta0, etta_val / beta_sh)
            # Paso rechazado: ox, ou no cambian
        elif rho_i < rho1:
            etta_val = max(etta0, etta_val / beta_sh)
            ox = x_opt.copy()
            ou = u_opt.copy()
        elif rho_i < rho2:
            ox = x_opt.copy()
            ou = u_opt.copy()
        else:    # rho_i >= rho2
            etta_val = min(etta1, beta_gr * etta_val)
            ox = x_opt.copy()
            ou = u_opt.copy()

        # 8. Actualizar parametros CVXPY
        prob_dict['etta'].value = etta_val
        for k in range(T + 1):
            prob_dict['ox_params'][k].value = ox[:, k:k+1]
            prob_dict['ou_params'][k].value = ou[:, k:k+1] if k < T else ou[:, T-1:T]

        if iteration >= 3:
            no_first_iterations = True

    # --- Desescalar resultado final ---
    x_final = np.zeros((nx, T + 1))
    u_final = np.zeros((nu, T))
    for k in range(T + 1):
        x_final[:, k:k+1] = inv_scale_x(ox[:, k:k+1], S_x, c_x)
    for k in range(T):
        u_final[:, k:k+1] = inv_scale_u(ou[:, k:k+1], S_u, c_u)

    return {
        'x_opt': ox,
        'u_opt': ou,
        'x_opt_unscaled': x_final,
        'u_opt_unscaled': u_final,
        'converged': converged,
        'iterations': final_iter,
        'history': history,
    }
