"""
Discretization of the linearized, scaled LTV system.

Over one interval [t_k, t_k+1] of length tau, the scaled linearized dynamics
are

    xdot = A(t) x + B(t) u(t) + F(t) p + y(t)

and the exact solution at t_k+1 is obtained by propagating the state
transition matrix together with the convolution integrals.  Writing
Phi(t) := Phi(t, t_k), all of the required quantities satisfy linear ODEs
that can be integrated in one shot with RK4:

    Phi'   = A Phi                         Phi(t_k) = I
    Gm'    = A Gm + lam_minus(t) B         Gm(t_k)  = 0
    Gp'    = A Gp + lam_plus(t)  B         Gp(t_k)  = 0     (FOH only)
    H'     = A H  + F                      H(t_k)   = 0     (dilation only)
    zz'    = A zz + y                      zz(t_k)  = 0

which gives the discrete update

    ZOH:  x_k+1 = A_k x_k + Bm_k u_k                + F_k p + y_k
    FOH:  x_k+1 = A_k x_k + Bm_k u_k + Bp_k u_k+1   + F_k p + y_k

with A_k = Phi(t_k+1), Bm_k = Gm(t_k+1), Bp_k = Gp(t_k+1), F_k = H(t_k+1),
y_k = zz(t_k+1).  Under ZOH lam_minus == 1 and lam_plus == 0; under FOH

    lam_minus(t) = (t_k+1 - t)/tau,   lam_plus(t) = (t - t_k)/tau.

Linearization reference inside the interval
-------------------------------------------
Under ZOH the Jacobians are evaluated at the node reference (x_bar_k, u_bar_k)
held constant across the interval -- the same convention as the reference
notebooks.  Under FOH they are evaluated along the linearly interpolated
reference (x_bar, u_bar), which is what makes the first-order hold worth its
extra cost.
"""

import numpy as np

from .modes import FOH
from .scaling import (scale_A, scale_B, scale_C, scale_D, scale_F, scale_G,
                      scale_y, scale_z)


class ScalingBundle:
    """The scaling matrices needed by the discretizer, precomputed once."""

    def __init__(self, S_x, c_x, S_u, c_u, S_p=None, c_p=None):
        self.S_x = S_x
        self.c_x = c_x
        self.S_u = S_u
        self.c_u = c_u
        self.S_p = S_p
        self.c_p = c_p
        self.inv_S_x = np.linalg.inv(S_x)
        self.inv_S_u = np.linalg.inv(S_u)
        self.inv_S_p = np.linalg.inv(S_p) if S_p is not None else None


def _scaled_dynamics_matrices(funcs, hx, hu, hp, t, dyn_par, sc, dilation):
    """Evaluate and scale A, B, F, y at one point."""
    A_hat = funcs['A'](hx, hu, hp, t, dyn_par)
    B_hat = funcs['B'](hx, hu, hp, t, dyn_par)
    y_hat = funcs['y'](hx, hu, hp, t, dyn_par)
    F_hat = funcs['F'](hx, hu, hp, t, dyn_par) if dilation else None

    A = scale_A(A_hat, sc.inv_S_x, sc.S_x)
    B = scale_B(B_hat, sc.inv_S_x, sc.S_u)
    F = scale_F(F_hat, sc.inv_S_x, sc.S_p) if dilation else None
    y = scale_y(A_hat, B_hat, F_hat, y_hat,
                sc.inv_S_x, sc.c_x, sc.c_u, sc.c_p if dilation else None)
    return A, B, F, y


def discretize_interval(funcs, ohx_k, ohx_k1, ohu_k, ohu_k1, ohp,
                        t0, tau, size_N, dyn_par, sc,
                        hold='ZOH', dilation=False):
    """
    RK4 discretization of one interval.

    Parameters
    ----------
    funcs : dict
        Evaluation functions 'A', 'B', 'F', 'y' with the common signature
        (hx, hu, hp, t, dyn_par).  All operate in unscaled coordinates.
    ohx_k, ohx_k1 : (nx,1)
        Unscaled reference state at the interval endpoints.  `ohx_k1` is only
        read under FOH.
    ohu_k, ohu_k1 : (nu,1)
        Unscaled reference control at the interval endpoints.
    ohp : (1,1) or None
        Unscaled reference dilation parameter (tf).
    t0, tau : float
        Interval start and length, in whichever time the formulation uses
        (absolute without dilation, normalized with it).
    size_N : int
        Number of RK4 sub-steps + 1 (matches the reference implementations:
        the interval is split into size_N-1 sub-steps).

    Returns
    -------
    (A_k, Bm_k, Bp_k, F_k, y_k)
        `Bp_k` is None under ZOH, `F_k` is None without dilation.
    """
    nx = ohx_k.shape[0]
    nu = ohu_k.shape[0]
    foh = (hold == FOH)

    Phi = np.eye(nx)
    Gm = np.zeros((nx, nu))
    Gp = np.zeros((nx, nu)) if foh else None
    H = np.zeros((nx, 1)) if dilation else None
    zz = np.zeros((nx, 1))

    n_sub = max(int(size_N) - 1, 1)
    dt = tau / n_sub
    t = t0

    def stage(t_eval):
        """Scaled A, B, F, y and the hold weights at a sub-time."""
        if foh:
            alpha = (t_eval - t0) / tau
            alpha = min(max(alpha, 0.0), 1.0)
            hx = (1.0 - alpha) * ohx_k + alpha * ohx_k1
            hu = (1.0 - alpha) * ohu_k + alpha * ohu_k1
            lam_m, lam_p = 1.0 - alpha, alpha
        else:
            hx, hu = ohx_k, ohu_k
            lam_m, lam_p = 1.0, 0.0
        A, B, F, y = _scaled_dynamics_matrices(
            funcs, hx, hu, ohp, t_eval, dyn_par, sc, dilation)
        return A, B, F, y, lam_m, lam_p

    def derivatives(t_eval, Phi_s, Gm_s, Gp_s, H_s, zz_s):
        A, B, F, y, lam_m, lam_p = stage(t_eval)
        d_Phi = A @ Phi_s
        d_Gm = A @ Gm_s + lam_m * B
        d_Gp = (A @ Gp_s + lam_p * B) if foh else None
        d_H = (A @ H_s + F) if dilation else None
        d_zz = A @ zz_s + y
        return d_Phi, d_Gm, d_Gp, d_H, d_zz

    for _ in range(n_sub):
        t_half = t + 0.5 * dt
        t_full = t + dt

        k1 = derivatives(t, Phi, Gm, Gp, H, zz)

        s2 = _advance(Phi, Gm, Gp, H, zz, k1, 0.5 * dt, foh, dilation)
        k2 = derivatives(t_half, *s2)

        s3 = _advance(Phi, Gm, Gp, H, zz, k2, 0.5 * dt, foh, dilation)
        k3 = derivatives(t_half, *s3)

        s4 = _advance(Phi, Gm, Gp, H, zz, k3, dt, foh, dilation)
        k4 = derivatives(t_full, *s4)

        w = dt / 6.0
        Phi = Phi + w * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Gm = Gm + w * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if foh:
            Gp = Gp + w * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        if dilation:
            H = H + w * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        zz = zz + w * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4])

        t = t_full

    return Phi, Gm, Gp, H, zz


def _advance(Phi, Gm, Gp, H, zz, k, step, foh, dilation):
    """Euler predictor used to build the intermediate RK4 stages."""
    return (
        Phi + step * k[0],
        Gm + step * k[1],
        (Gp + step * k[2]) if foh else None,
        (H + step * k[3]) if dilation else None,
        zz + step * k[4],
    )


def discretize_constraints(funcs, hx, hu, hp, t0, dyn_par, sc, dilation):
    """
    Linearized nonconvex path constraints at one node.

    These are pointwise, so no RK4 propagation is needed -- evaluate the
    Jacobians at the linearization point and scale them.
    """
    C_hat = funcs['C'](hx, hu, hp, t0, dyn_par)
    D_hat = funcs['D'](hx, hu, hp, t0, dyn_par)
    z_hat = funcs['z'](hx, hu, hp, t0, dyn_par)
    G_hat = funcs['G'](hx, hu, hp, t0, dyn_par) if dilation else None

    C = scale_C(C_hat, sc.S_x)
    D = scale_D(D_hat, sc.S_u)
    G = scale_G(G_hat, sc.S_p) if dilation else None
    z = scale_z(C_hat, D_hat, G_hat, z_hat,
                sc.c_x, sc.c_u, sc.c_p if dilation else None)
    return C, D, G, z


def rk4_flow_map(x_k_scaled, u_k_scaled, u_k1_scaled, p_scaled, t0, dt,
                 dyn_par, f_func, sc, hold='ZOH', dilation=False,
                 substeps=1):
    """
    Propagate the *true* nonlinear dynamics across one interval.

    Takes scaled quantities and returns the unscaled state at t0 + dt.  Under
    dilation the generated `f` already carries the factor p, so `dt` is the
    normalized step and the result is still the physical state.
    """
    from .scaling import unscale_p, unscale_u, unscale_x

    x = unscale_x(x_k_scaled, sc.S_x, sc.c_x)
    u_k = unscale_u(u_k_scaled, sc.S_u, sc.c_u)
    hp = unscale_p(p_scaled, sc.S_p, sc.c_p) if dilation else np.zeros((0, 1))

    if hold == FOH and u_k1_scaled is not None:
        u_k1 = unscale_u(u_k1_scaled, sc.S_u, sc.c_u)
    else:
        u_k1 = u_k

    def control_at(t_eval):
        if hold != FOH:
            return u_k
        alpha = (t_eval - t0) / dt
        alpha = min(max(alpha, 0.0), 1.0)
        return (1.0 - alpha) * u_k + alpha * u_k1

    h = dt / max(int(substeps), 1)
    t = t0
    for _ in range(max(int(substeps), 1)):
        k1 = f_func(x, control_at(t), hp, t, dyn_par)
        k2 = f_func(x + 0.5 * h * k1, control_at(t + 0.5 * h), hp, t + 0.5 * h, dyn_par)
        k3 = f_func(x + 0.5 * h * k2, control_at(t + 0.5 * h), hp, t + 0.5 * h, dyn_par)
        k4 = f_func(x + h * k3, control_at(t + h), hp, t + h, dyn_par)
        x = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
    return x
