import numpy as np


# ==========================================
# ESCALADO DE ESTADOS Y CONTROLES
# ==========================================

def scale_x(x_no_scaled, inv_S_x, c_x):
    """x_scaled = inv(S_x) @ (x - c_x)"""
    return inv_S_x @ (x_no_scaled - c_x)


def scale_u(u_no_scaled, inv_S_u, c_u):
    """u_scaled = inv(S_u) @ (u - c_u)"""
    return inv_S_u @ (u_no_scaled - c_u)


def inv_scale_x(x_scaled, S_x, c_x):
    """x = S_x @ x_scaled + c_x"""
    return S_x @ x_scaled + c_x


def inv_scale_u(u_scaled, S_u, c_u):
    """u = S_u @ u_scaled + c_u"""
    return S_u @ u_scaled + c_u


# ==========================================
# ESCALADO DE MATRICES DE DINAMICA (A, B, y)
# ==========================================

def scale_A(A_no_scaled, inv_S_x, S_x):
    """A_scaled = inv(S_x) @ A @ S_x"""
    return inv_S_x @ A_no_scaled @ S_x


def scale_B(B_no_scaled, inv_S_x, S_u):
    """B_scaled = inv(S_x) @ B @ S_u"""
    return inv_S_x @ B_no_scaled @ S_u


def scale_y(A_no_scaled, B_no_scaled, y_no_scaled, inv_S_x, c_x, c_u):
    """y_scaled = inv(S_x) @ (A @ c_x + B @ c_u + y)"""
    return inv_S_x @ (A_no_scaled @ c_x + B_no_scaled @ c_u + y_no_scaled)


# ==========================================
# ESCALADO DE MATRICES DE RESTRICCIONES (C, D, z)
# ==========================================

def scale_C(C_no_scaled, S_x):
    """C_scaled = C @ S_x"""
    return C_no_scaled @ S_x


def scale_D(D_no_scaled, S_u):
    """D_scaled = D @ S_u"""
    return D_no_scaled @ S_u


def scale_z(C_no_scaled, D_no_scaled, z_no_scaled, c_x, c_u):
    """z_scaled = C @ c_x + D @ c_u + z"""
    return C_no_scaled @ c_x + D_no_scaled @ c_u + z_no_scaled
