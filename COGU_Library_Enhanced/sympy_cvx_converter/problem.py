import cvxpy as cp
import numpy as np

from .cost_dsl import resolve_term


def _build_cost_from_terms(cost_terms, vars_, T, sqrt_tau, tau_val, tau_lamb):
    coeff_map = {'sqrt_tau': sqrt_tau, 'tau': float(tau_val), 'tau_lamb': tau_lamb}
    cost = 0
    for term in cost_terms:
        rt = resolve_term(term)
        v_full = vars_[rt.var]
        coeff = coeff_map[rt.coeff] if isinstance(rt.coeff, str) else float(rt.coeff)
        K = T + 1 if rt.k_range == 'T+1' else T
        for k in range(K):
            v = v_full[rt.slc, k:k+1]
            expr = coeff * v if rt.offset is None else coeff * (v - rt.offset)
            if rt.kind == 'sumsq':
                cost += rt.weight * cp.sum_squares(expr)
            elif rt.kind == 'norm1':
                cost += rt.weight * cp.norm(expr, 1)
            elif rt.kind == 'norm2':
                cost += rt.weight * cp.norm(expr, 2)
    return cost


def _build_bound_constraints(bounds, var, S, c, K):
    norm_map = {'norm1': 1, 'norm2': 2, 'norminf': 'inf'}
    result = []
    for spec in bounds:
        if isinstance(spec, tuple):
            slc, norm_type, limit = spec
            spec = {'kind': 'norm', 'slice': slc, 'norm_type': norm_type, 'limit': limit}
        kind = spec['kind']
        slc  = spec['slice']
        if kind == 'norm':
            p_val = norm_map[spec['norm_type']]
            # S/c custom (opcional): matriz de transformacion propia en vez del sub-bloque de S global
            S_local = spec.get('S')
            c_local = spec.get('c')
            for k in range(K):
                if S_local is not None:
                    affine = S_local @ var[slc, k:k+1] + c_local
                else:
                    affine = S[slc, slc] @ var[slc, k:k+1] + c[slc, 0:1]
                result.append(cp.norm(affine, p_val) <= float(spec['limit']))
        elif kind == 'box':
            lower = spec.get('lower')
            upper = spec.get('upper')
            for k in range(K):
                affine = S[slc, slc] @ var[slc, k:k+1] + c[slc, 0:1]
                if lower is not None:
                    result.append(affine >= float(lower))
                if upper is not None:
                    result.append(affine <= float(upper))
        else:
            raise ValueError(f"_build_bound_constraints: unknown kind={kind!r}")
    return result


def build_problem(nx, nu, T, ng=0,
                  tau_val=1.0, lamb_val=1.0,
                  S_x=None, c_x=None, S_u=None, c_u=None,
                  state_bounds=None, control_bounds=None,
                  cost_terms=None):
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
    if cost_terms is None:
        raise ValueError("cost_terms es obligatorio: declara el costo del problema "
                         "(ej. [{'kind':'sumsq','var':'u','coeff':'sqrt_tau'}]).")
    cost = _build_cost_from_terms(cost_terms, {'u': u, 'x': x}, T, sqrt_tau, tau_val, tau_lamb)
    for k in range(T):
        cost += cp.norm(tau_lamb * vc[:, k:k+1], 1)   # SCVx penalty — always fixed

    if ng > 0:
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

    # Bounds de estado y control (soporta tuple format y dict format con kind='norm'|'box')
    if state_bounds:
        constraints += _build_bound_constraints(state_bounds, x, S_x, c_x, T+1)
    if control_bounds:
        constraints += _build_bound_constraints(control_bounds, u, S_u, c_u, T)

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
        'cost_terms': cost_terms,

        # ── Fase 3 — metadatos para generacion C generica (Fases 5 y 6) ──
        'constant_params': {
            'tau_val':  float(tau_val),
            'sqrt_tau': sqrt_tau,
            'tau_lamb': tau_lamb,
            'S_x': S_x, 'c_x': c_x,
            'S_u': S_u, 'c_u': c_u,
        },
        'iteration_params': {
            'per_step': {
                'A':  {'count': T,   'shape': (nx, nx)},
                'B':  {'count': T,   'shape': (nx, nu)},
                'y':  {'count': T,   'shape': (nx, 1)},
                **({'C': {'count': T+1, 'shape': (ng, nx)},
                    'D': {'count': T+1, 'shape': (ng, nu)},
                    'z': {'count': T+1, 'shape': (ng, 1)}} if ng > 0 else {}),
                'ox': {'count': T+1, 'shape': (nx, 1)},
                'ou': {'count': T+1, 'shape': (nu, 1)},
            },
            'scalar':   ['etta'],
            'boundary': ['start_pos', 'end_pos'],
        },
    }
