"""
Construction of the SCvx convex subproblem in CVXPY.

The subproblem solved at every SCvx iteration is

    min   L_lambda(x, u, p, vc, vi) + alpha_p (S_p p + c_p)
    s.t.  x_k+1 = A_k x_k + Bm_k u_k [+ Bp_k u_k+1] [+ F_k p] + y_k + vc_k
          C_k x_k + D_k u_k [+ G_k p] + z_k <= vi_k
          ||x_k - x_bar_k||_inf + ||u_k - u_bar_k||_inf [+ ||p - p_bar||_inf] <= eta
          [tf_min <= S_p p + c_p <= tf_max]
          x_k in X_k,  u_k in U_k

The terms in brackets are present only in the relevant time mode.

Two deliberate choices carried over from the reference implementation, both
needed for CVXPYgen:

1. tau and lambda enter the cost as plain Python floats rather than
   cp.Parameter.  A scalar parameter multiplying a 2-D slice trips a bug in
   the COO backend of recent CVXPY versions.

2. The linearized dynamics are exposed as per-step parameters A_0..A_T-1
   instead of one stacked A matrix.  The resulting DPP parameter matrix is
   sparse, which is the difference between a canonicalization that fits in
   RAM and one that does not.
"""

import cvxpy as cp
import numpy as np

from .cost_dsl import coefficient_map, node_count, resolve_coefficient, resolve_term
from .modes import FOH, control_node_count, parse_time_mode


# ----------------------------------------------------------------------
# Cost
# ----------------------------------------------------------------------

def build_convex_cost(cost_terms, variables, T, u_nodes, tau, lamb):
    """User running cost, expressed in CVXPY."""
    sqrt_tau = float(tau) ** 0.5
    tau_lamb = float(tau) * float(lamb)
    coeff_map = coefficient_map(sqrt_tau, tau, tau_lamb)

    cost = 0
    for term in cost_terms:
        rt = resolve_term(term)
        var = variables[rt.var]
        coeff = resolve_coefficient(rt, coeff_map)
        K = node_count(rt.k_range, T, u_nodes)
        available = var.shape[1]
        if K > available:
            raise ValueError(
                f"cost term {term!r} spans {K} nodes but variable "
                f"'{rt.var}' only has {available} columns"
            )
        for k in range(K):
            v = var[rt.slc, k:k + 1]
            expr = coeff * v if rt.offset is None else coeff * (v - rt.offset)
            if rt.kind == 'sumsq':
                cost += rt.weight * cp.sum_squares(expr)
            elif rt.kind == 'norm1':
                cost += rt.weight * cp.norm(expr, 1)
            elif rt.kind == 'norm2':
                cost += rt.weight * cp.norm(expr, 2)
    return cost


# ----------------------------------------------------------------------
# Convex state / control bounds
# ----------------------------------------------------------------------

def build_bound_constraints(bounds, var, S, c, K):
    """
    Translate the bound DSL into CVXPY constraints on the *unscaled* quantity
    S[slice] var[slice, k] + c[slice].

    Supported entries:
        (slice, 'norm1'|'norm2'|'norminf', limit)         -- shorthand tuple
        {'kind': 'norm', 'slice':, 'norm_type':, 'limit':,
         'S': optional, 'c': optional}
        {'kind': 'box', 'slice':, 'lower':, 'upper':}
        {'kind': 'box_affine', 'c_affine': (1,n), 'lower':, 'upper':}

    An explicit 'S'/'c' pair replaces the global scaling sub-block, which is
    how a bound on a derived quantity is expressed -- e.g. bounding the
    angular velocity inv(J) h when the state is the angular momentum h.
    """
    norm_map = {'norm1': 1, 'norm2': 2, 'norminf': 'inf'}
    result = []

    for spec in bounds:
        if isinstance(spec, tuple):
            slc, norm_type, limit = spec
            spec = {'kind': 'norm', 'slice': slc,
                    'norm_type': norm_type, 'limit': limit}
        kind = spec['kind']

        if kind == 'norm':
            slc = spec['slice']
            p_val = norm_map[spec['norm_type']]
            S_local = spec.get('S')
            c_local = spec.get('c')
            for k in range(K):
                if S_local is not None:
                    affine = S_local @ var[slc, k:k + 1] + c_local
                else:
                    affine = S[slc, slc] @ var[slc, k:k + 1] + c[slc, 0:1]
                result.append(cp.norm(affine, p_val) <= float(spec['limit']))

        elif kind == 'box':
            slc = spec['slice']
            lower = spec.get('lower')
            upper = spec.get('upper')
            for k in range(K):
                affine = S[slc, slc] @ var[slc, k:k + 1] + c[slc, 0:1]
                if lower is not None:
                    result.append(affine >= float(lower))
                if upper is not None:
                    result.append(affine <= float(upper))

        elif kind == 'box_affine':
            lower = spec.get('lower')
            upper = spec.get('upper')
            c_affine = spec['c_affine']
            for k in range(K):
                affine = c_affine @ (S @ var[:, k:k + 1] + c[:, 0:1])
                if lower is not None:
                    result.append(affine >= float(lower))
                if upper is not None:
                    result.append(affine <= float(upper))

        else:
            raise ValueError(f"unknown bound kind: {kind!r}")

    return result


# ----------------------------------------------------------------------
# Problem
# ----------------------------------------------------------------------

def build_problem(nx, nu, T, ng=0,
                  time_mode='ZOH',
                  tau_val=1.0, lamb_val=1.0,
                  S_x=None, c_x=None, S_u=None, c_u=None,
                  S_p=None, c_p=None,
                  alpha_p=0.0, tf_min=None, tf_max=None,
                  state_bounds=None, control_bounds=None,
                  cost_terms=None, pin_aware_trust_region=True):
    """
    Build the CVXPY convex subproblem.

    All arguments are plain numbers / numpy arrays / slices, which keeps the
    whole call serializable -- the C generation worker rebuilds an identical
    problem in a clean subprocess from exactly these kwargs.

    Parameters that change between SCvx iterations are cp.Parameter:
        A_k, Bm_k, [Bp_k], [F_k], y_k     linearized dynamics, per step
        C_k, D_k, [G_k], z_k              linearized nonconvex constraints
        ox_k, ou_k, [op]                  reference trajectory
        eta, start_pos, end_pos

    pin_aware_trust_region
        Leave the boundary-pinned states x_0 and x_T out of the trust region.
        Set False only to reproduce the legacy formulation bit for bit.
    """
    mode_name, hold, dilation = parse_time_mode(time_mode)
    u_nodes = control_node_count(hold, T)

    if cost_terms is None:
        raise ValueError(
            "cost_terms is required: declare the running cost, e.g. "
            "[{'kind': 'sumsq', 'var': 'u', 'coeff': 'sqrt_tau'}]"
        )

    if S_x is None:
        S_x = np.eye(nx)
    if c_x is None:
        c_x = np.zeros((nx, 1))
    if S_u is None:
        S_u = np.eye(nu)
    if c_u is None:
        c_u = np.zeros((nu, 1))
    if dilation and (S_p is None or c_p is None):
        raise ValueError("time dilation requires S_p and c_p")

    tau_val = float(tau_val)
    lamb_val = float(lamb_val)
    tau_lamb = tau_val * lamb_val

    # ------------------------------------------------------------------
    # Variables
    # ------------------------------------------------------------------
    x = cp.Variable((nx, T + 1), name='x')
    u = cp.Variable((nu, u_nodes), name='u')
    vc = cp.Variable((nx, T), name='vc')
    vi = cp.Variable((ng, T + 1), name='vi') if ng > 0 else None
    p = cp.Variable((1, 1), name='p') if dilation else None

    # ------------------------------------------------------------------
    # Parameters (per step -> sparse DPP parameter matrix)
    # ------------------------------------------------------------------
    A_params = [cp.Parameter((nx, nx), name=f'A_{k}') for k in range(T)]
    Bm_params = [cp.Parameter((nx, nu), name=f'Bm_{k}') for k in range(T)]
    Bp_params = ([cp.Parameter((nx, nu), name=f'Bp_{k}') for k in range(T)]
                 if hold == FOH else None)
    F_params = ([cp.Parameter((nx, 1), name=f'F_{k}') for k in range(T)]
                if dilation else None)
    y_params = [cp.Parameter((nx, 1), name=f'y_{k}') for k in range(T)]

    if ng > 0:
        C_params = [cp.Parameter((ng, nx), name=f'C_{k}') for k in range(T + 1)]
        D_params = [cp.Parameter((ng, nu), name=f'D_{k}') for k in range(T + 1)]
        G_params = ([cp.Parameter((ng, 1), name=f'G_{k}') for k in range(T + 1)]
                    if dilation else None)
        z_params = [cp.Parameter((ng, 1), name=f'z_{k}') for k in range(T + 1)]
    else:
        C_params = D_params = G_params = z_params = None

    ox_params = [cp.Parameter((nx, 1), name=f'ox_{k}') for k in range(T + 1)]
    ou_params = [cp.Parameter((nu, 1), name=f'ou_{k}') for k in range(u_nodes)]
    op_param = cp.Parameter((1, 1), name='op') if dilation else None

    eta = cp.Parameter(name='eta')
    start_pos = cp.Parameter((nx, 1), name='start_pos')
    end_pos = cp.Parameter((nx, 1), name='end_pos')

    # ------------------------------------------------------------------
    # Cost
    # ------------------------------------------------------------------
    cost = build_convex_cost(cost_terms, {'x': x, 'u': u}, T, u_nodes,
                             tau_val, lamb_val)

    for k in range(T):
        cost += cp.norm(tau_lamb * vc[:, k:k + 1], 1)

    if ng > 0:
        for k in range(T + 1):
            cost += cp.norm(tau_lamb * vi[:, k:k + 1], 1)

    # Convex terminal cost on the final time: alpha_p (S_p p + c_p).  The
    # constant alpha_p c_p does not change the argmin, so only the linear part
    # is built here; the reported objective adds the offset back.
    if dilation and float(alpha_p) != 0.0:
        cost += float(alpha_p) * float(np.asarray(S_p).reshape(-1)[0]) * p[0, 0]

    # ------------------------------------------------------------------
    # Constraints
    # ------------------------------------------------------------------
    constraints = [
        x[:, 0:1] == start_pos,
        x[:, T:T + 1] == end_pos,
    ]

    # Discretized (possibly dilated) dynamics
    for k in range(T):
        rhs = (A_params[k] @ x[:, k:k + 1]
               + Bm_params[k] @ u[:, k:k + 1]
               + y_params[k]
               + vc[:, k:k + 1])
        if hold == FOH:
            rhs = rhs + Bp_params[k] @ u[:, k + 1:k + 2]
        if dilation:
            rhs = rhs + F_params[k] @ p
        constraints += [x[:, k + 1:k + 2] == rhs]

    # Trust region, including the dilation parameter when it is free.
    #
    # x_0 and x_T are pinned by the boundary equality constraints above, so
    # ||x_0 - ox_0||_inf and ||x_T - ox_T||_inf are *constants*, not decision
    # quantities.  Charging them against eta does nothing useful and is
    # actively harmful: if the warm start disagrees with `start` or `end` --
    # which is the norm, e.g. any warm start whose angular velocity is built
    # by finite differences starts at zero -- that constant becomes a hard
    # floor on eta, and the subproblem turns infeasible the moment the trust
    # region shrinks below it.  Skip the state term at the pinned nodes.
    for k in range(T + 1):
        uk = min(k, u_nodes - 1)
        tr = cp.norm(u[:, uk:uk + 1] - ou_params[uk], 'inf')
        if not (pin_aware_trust_region and k in (0, T)):
            tr = tr + cp.norm(x[:, k:k + 1] - ox_params[k], 'inf')
        if dilation:
            tr = tr + cp.norm(p - op_param, 'inf')
        constraints += [tr <= eta]

    # Linearized nonconvex path constraints
    if ng > 0:
        for k in range(T + 1):
            uk = min(k, u_nodes - 1)
            lhs = (C_params[k] @ x[:, k:k + 1]
                   + D_params[k] @ u[:, uk:uk + 1]
                   + z_params[k])
            if dilation:
                lhs = lhs + G_params[k] @ p
            constraints += [lhs <= vi[:, k:k + 1]]

    # Convex bounds
    if state_bounds:
        constraints += build_bound_constraints(state_bounds, x, S_x, c_x, T + 1)
    if control_bounds:
        constraints += build_bound_constraints(control_bounds, u, S_u, c_u, u_nodes)

    # Free final time bounds
    if dilation:
        S_p_f = float(np.asarray(S_p).reshape(-1)[0])
        c_p_f = float(np.asarray(c_p).reshape(-1)[0])
        constraints += [S_p_f * p[0, 0] + c_p_f >= float(tf_min)]
        constraints += [S_p_f * p[0, 0] + c_p_f <= float(tf_max)]

    problem = cp.Problem(cp.Minimize(cost), constraints)

    return {
        'problem': problem,
        'time_mode': mode_name, 'hold': hold, 'dilation': dilation,
        'nx': nx, 'nu': nu, 'ng': ng, 'T': T, 'u_nodes': u_nodes,
        'variables': {'x': x, 'u': u, 'vc': vc, 'vi': vi, 'p': p},
        'A_params': A_params, 'Bm_params': Bm_params, 'Bp_params': Bp_params,
        'F_params': F_params, 'y_params': y_params,
        'C_params': C_params, 'D_params': D_params,
        'G_params': G_params, 'z_params': z_params,
        'ox_params': ox_params, 'ou_params': ou_params, 'op_param': op_param,
        'eta': eta, 'start_pos': start_pos, 'end_pos': end_pos,
        'S_x': S_x, 'c_x': c_x, 'S_u': S_u, 'c_u': c_u,
        'S_p': S_p, 'c_p': c_p,
        'alpha_p': float(alpha_p), 'tf_min': tf_min, 'tf_max': tf_max,
        'tau_val': tau_val, 'lamb_val': lamb_val,
        'pin_aware_trust_region': pin_aware_trust_region,
        'cost_terms': cost_terms,
    }
