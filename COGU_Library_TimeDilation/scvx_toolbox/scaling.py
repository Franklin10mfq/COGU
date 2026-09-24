"""
Affine scaling of states, controls and the dilation parameter.

Convention (hat = physical / unscaled, plain = scaled):

    x_hat = S_x x + c_x
    u_hat = S_u u + c_u
    p_hat = S_p p + c_p        (p_hat = tf, only under time dilation)

The linearized dynamics are scaled accordingly:

    A = S_x^-1 A_hat S_x
    B = S_x^-1 B_hat S_u
    F = S_x^-1 F_hat S_p
    y = S_x^-1 (A_hat c_x + B_hat c_u + F_hat c_p + y_hat)

and the linearized nonconvex constraints:

    C = C_hat S_x,  D = D_hat S_u,  G = G_hat S_p
    z = C_hat c_x + D_hat c_u + G_hat c_p + z_hat
"""

import numpy as np


# ----------------------------------------------------------------------
# States, controls, dilation parameter
# ----------------------------------------------------------------------

def scale_x(x_unscaled, inv_S_x, c_x):
    """x = S_x^-1 (x_hat - c_x)"""
    return inv_S_x @ (x_unscaled - c_x)


def scale_u(u_unscaled, inv_S_u, c_u):
    """u = S_u^-1 (u_hat - c_u)"""
    return inv_S_u @ (u_unscaled - c_u)


def scale_p(p_unscaled, inv_S_p, c_p):
    """p = S_p^-1 (p_hat - c_p)"""
    return inv_S_p @ (p_unscaled - c_p)


def unscale_x(x_scaled, S_x, c_x):
    """x_hat = S_x x + c_x"""
    return S_x @ x_scaled + c_x


def unscale_u(u_scaled, S_u, c_u):
    """u_hat = S_u u + c_u"""
    return S_u @ u_scaled + c_u


def unscale_p(p_scaled, S_p, c_p):
    """p_hat = S_p p + c_p"""
    return S_p @ p_scaled + c_p


# ----------------------------------------------------------------------
# Dynamics Jacobians
# ----------------------------------------------------------------------

def scale_A(A_hat, inv_S_x, S_x):
    return inv_S_x @ A_hat @ S_x


def scale_B(B_hat, inv_S_x, S_u):
    return inv_S_x @ B_hat @ S_u


def scale_F(F_hat, inv_S_x, S_p):
    return inv_S_x @ F_hat @ S_p


def scale_y(A_hat, B_hat, F_hat, y_hat, inv_S_x, c_x, c_u, c_p):
    """
    y = S_x^-1 (A_hat c_x + B_hat c_u + F_hat c_p + y_hat).

    F_hat / c_p may be None when time dilation is disabled.
    """
    acc = A_hat @ c_x + B_hat @ c_u + y_hat
    if F_hat is not None and c_p is not None:
        acc = acc + F_hat @ c_p
    return inv_S_x @ acc


# ----------------------------------------------------------------------
# Nonconvex constraint Jacobians
# ----------------------------------------------------------------------

def scale_C(C_hat, S_x):
    return C_hat @ S_x


def scale_D(D_hat, S_u):
    return D_hat @ S_u


def scale_G(G_hat, S_p):
    return G_hat @ S_p


def scale_z(C_hat, D_hat, G_hat, z_hat, c_x, c_u, c_p):
    acc = C_hat @ c_x + D_hat @ c_u + z_hat
    if G_hat is not None and c_p is not None:
        acc = acc + G_hat @ c_p
    return acc


# ----------------------------------------------------------------------
# Defaults
# ----------------------------------------------------------------------

def default_scaling(nx, nu, S_x=None, c_x=None, S_u=None, c_u=None):
    """Fill in identity/zero scaling for whatever the user left out."""
    S_x = np.eye(nx) if S_x is None else np.asarray(S_x, dtype=float)
    c_x = np.zeros((nx, 1)) if c_x is None else np.asarray(c_x, dtype=float).reshape(nx, 1)
    S_u = np.eye(nu) if S_u is None else np.asarray(S_u, dtype=float)
    c_u = np.zeros((nu, 1)) if c_u is None else np.asarray(c_u, dtype=float).reshape(nu, 1)
    return S_x, c_x, S_u, c_u


def dilation_scaling(tf_min, tf_max):
    """
    Scaling of the dilation parameter so that the scaled p lives in [0, 1]
    when tf ranges over [tf_min, tf_max]:

        p_hat = (tf_max - tf_min) p + tf_min
    """
    span = float(tf_max) - float(tf_min)
    if span <= 0.0:
        raise ValueError(f"tf_max ({tf_max}) must be strictly greater than tf_min ({tf_min})")
    return np.array([[span]]), np.array([[float(tf_min)]])
