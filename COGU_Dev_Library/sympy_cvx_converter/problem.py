import cvxpy as cp
import numpy as np


def build_problem(nx, nu, T, ng=0,
                  tau_val=1.0, lamb_val=1.0,
                  S_x=None, c_x=None, S_u=None, c_u=None,
                  state_bounds=None, control_bounds=None):
    """
    Construye un problema CVXPY para SCVx compatible con CVXPYgen + CVXPY 1.8.

    Dos cambios respecto al notebook de Mariel (CVXPY 1.6):
    1. tau/lamb como floats en el costo (no cp.Parameter) — evita bug mul_elem
       del backend COO de CVXPY 1.8 con parametros escalares x slices 2D.
    2. Parametros por paso A_0..A_{T-1} en vez de A_discrete apilado — la
       matriz parametrica DPP es dispersa (~17K filas vs ~76M con stacked).

    S_x, c_x, S_u, c_u y limites de bounds son constantes numpy.

    Parametros que cambian entre iteraciones SCVx (cp.Parameter):
        A_{k}, B_{k}, y_{k}   — dinamica linealizada (por paso)
        C_{k}, D_{k}, z_{k}   — restricciones no-convexas (por paso)
        ox_{k}, ou_{k}         — trayectoria de referencia (por paso)
        etta, start_pos, end_pos
    """

    sqrt_tau = float(tau_val ** 0.5)
    tau_lamb = float(tau_val * lamb_val)

    if S_x is None: S_x = np.eye(nx)
    if c_x is None: c_x = np.zeros((nx, 1))
    if S_u is None: S_u = np.eye(nu)
    if c_u is None: c_u = np.zeros((nu, 1))

    norm_map = {'norm1': 1, 'norm2': 2, 'norminf': 'inf'}

    # ==================================================
    # VARIABLES
    # ==================================================
    x  = cp.Variable((nx, T+1), name='x')
    u  = cp.Variable((nu, T),   name='u')
    vc = cp.Variable((nx, T),   name='vc')
    vi = cp.Variable((ng, T+1), name='vi') if ng > 0 else None

    # ==================================================
    # PARAMETROS por paso — locales a un unico constraint
    # → matriz parametrica DPP dispersa → cabe en RAM para CVXPYgen
    # ==================================================

    A_params = [cp.Parameter((nx, nx), name=f'A_{k}') for k in range(T)]
    B_params = [cp.Parameter((nx, nu), name=f'B_{k}') for k in range(T)]
    y_params = [cp.Parameter((nx, 1),  name=f'y_{k}') for k in range(T)]

    if ng > 0:
        C_params = [cp.Parameter((ng, nx), name=f'C_{k}') for k in range(T+1)]
        D_params = [cp.Parameter((ng, nu), name=f'D_{k}') for k in range(T+1)]
        z_params = [cp.Parameter((ng, 1),  name=f'z_{k}') for k in range(T+1)]
    else:
        C_params = D_params = z_params = None

    ox_params = [cp.Parameter((nx, 1), name=f'ox_{k}') for k in range(T+1)]
    ou_params = [cp.Parameter((nu, 1), name=f'ou_{k}') for k in range(T+1)]

    etta      = cp.Parameter(name='etta')
    start_pos = cp.Parameter((nx, 1),  name='start_pos')
    end_pos   = cp.Parameter((nx, 1),  name='end_pos')

    # ==================================================
    # COSTO — sqrt_tau y tau_lamb son floats, sin param x var
    # ==================================================
    cost = 0
    for k in range(T):
        cost += cp.sum_squares(sqrt_tau * u[:, k:k+1])
        cost += cp.norm(tau_lamb * vc[:, k:k+1], 1)

    if ng > 0:
        for k in range(T+1):
            cost += cp.norm(tau_lamb * vi[:, k:k+1], 1)
        for k in range(T+1):
            cost += cp.norm(tau_lamb * vi[:, k:k+1], 1)

    # ==================================================
    # RESTRICCIONES
    # ==================================================
    constraints = []

    constraints += [x[:, 0] == start_pos[:, 0]]
    constraints += [x[:, T] == end_pos[:, 0]]

    # Dinamica linealizada (por paso)
    for k in range(T):
        constraints += [
            x[:, k+1:k+2] == A_params[k] @ x[:, k:k+1]
                            + B_params[k] @ u[:, k:k+1]
                            + y_params[k]
                            + vc[:, k:k+1]
        ]

    # Trust region
    for k in range(T):
        constraints += [
            cp.norm(x[:, k:k+1] - ox_params[k], 'inf')
            + cp.norm(u[:, k:k+1] - ou_params[k], 'inf')
            <= etta
        ]
    constraints += [
        cp.norm(x[:, T:T+1]   - ox_params[T],   'inf')
        + cp.norm(u[:, T-1:T] - ou_params[T-1], 'inf')
        <= etta
    ]

    # Restricciones no-convexas linealizadas (por paso)
    if ng > 0:
        for k in range(T):
            constraints += [
                C_params[k] @ x[:, k:k+1]
                + D_params[k] @ u[:, k:k+1]
                + z_params[k]
                <= vi[:, k:k+1]
            ]
        constraints += [
            C_params[T] @ x[:, T:T+1]
            + D_params[T] @ u[:, T-1:T]
            + z_params[T]
            <= vi[:, T:T+1]
        ]

    # Bounds de estado (constantes numpy)
    if state_bounds:
        for slc, norm_type, limit in state_bounds:
            p_val = norm_map[norm_type]
            for k in range(T+1):
                constraints += [
                    cp.norm(S_x[slc, slc] @ x[slc, k:k+1]
                            + c_x[slc, 0:1], p_val) <= float(limit)
                ]

    # Bounds de control (constantes numpy)
    if control_bounds:
        for slc, norm_type, limit in control_bounds:
            p_val = norm_map[norm_type]
            for k in range(T):
                constraints += [
                    cp.norm(S_u[slc, slc] @ u[slc, k:k+1]
                            + c_u[slc, 0:1], p_val) <= float(limit)
                ]

    # ==================================================
    # PROBLEMA
    # ==================================================
    problem = cp.Problem(cp.Minimize(cost), constraints)

    return {
        'problem':   problem,
        'nx': nx, 'nu': nu, 'ng': ng, 'T': T,
        'variables': {'x': x, 'u': u, 'vc': vc, 'vi': vi},
        'A_params': A_params, 'B_params': B_params, 'y_params': y_params,
        'C_params': C_params, 'D_params': D_params, 'z_params': z_params,
        'ox_params': ox_params, 'ou_params': ou_params,
        'etta':      etta,
        'start_pos': start_pos,
        'end_pos':   end_pos,
        'S_x': S_x, 'c_x': c_x, 'S_u': S_u, 'c_u': c_u,
        'tau_val': tau_val, 'lamb_val': lamb_val,
    }
