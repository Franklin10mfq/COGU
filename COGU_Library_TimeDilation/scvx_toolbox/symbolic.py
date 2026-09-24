"""
Symbolic model handling: SymPy in, numpy (or C) evaluation functions out.

The user supplies the dynamics in *absolute* time,

    d(x_hat)/dt = f(x_hat, u_hat, t)

and, when time dilation is enabled, the library forms the dilated dynamics on
the normalized interval s in [0, 1] itself:

    d(x_hat)/ds = p_hat * f(x_hat, u_hat)   with p_hat = tf

so the same model definition can be reused across all four time modes.

Linearization (first-order Taylor around the reference x_bar, u_bar, p_bar):

    A_hat = df/dx,  B_hat = df/du,  F_hat = df/dp
    y_hat = f - A_hat x_bar - B_hat u_bar - F_hat p_bar

and identically for the nonconvex inequalities g(x, u, p) <= 0:

    C_hat = dg/dx,  D_hat = dg/du,  G_hat = dg/dp
    z_hat = g - C_hat x_bar - D_hat u_bar - G_hat p_bar

Every generated function shares one signature,

    fn(hx, hu, hp, t, dyn_par) -> numpy array

where `hx`, `hu`, `hp` are unscaled column vectors.  `hp` is a (1,1) array
holding tf under dilation, and an empty array otherwise.
"""

import re

import numpy as np
import sympy as sp

DILATION_SYMBOL_NAME = 'p_dil'


# ----------------------------------------------------------------------
# Symbol -> array index substitution
# ----------------------------------------------------------------------

def _substitute_symbols(expr_str, state_names, control_names,
                        dilation_names, dyn_param_names,
                        fmt_state, fmt_control, fmt_dilation, fmt_param):
    """Replace symbol names by indexed array accesses, longest name first.

    Longest-first matters: without it a state named `q1` could be clobbered
    while rewriting a state named `q10`.
    """
    replacements = []
    for i, name in enumerate(state_names):
        replacements.append((name, fmt_state(i)))
    for i, name in enumerate(control_names):
        replacements.append((name, fmt_control(i)))
    for i, name in enumerate(dilation_names):
        replacements.append((name, fmt_dilation(i)))
    for i, name in enumerate(dyn_param_names):
        replacements.append((name, fmt_param(i)))

    for name, target in sorted(replacements, key=lambda kv: -len(kv[0])):
        expr_str = re.sub(rf"\b{re.escape(name)}\b", target, expr_str)
    return expr_str


def matrix_to_numpy_string(matrix, state_names, control_names,
                           dilation_names, dyn_param_names):
    """Render a SymPy matrix as the source of a numpy 2-D array literal."""
    matrix = sp.Matrix(matrix)
    n_rows, n_cols = matrix.shape
    if n_rows == 0:
        return "np.zeros((0, 1))"

    rows = []
    for i in range(n_rows):
        rows.append(",".join(sp.pycode(matrix[i, j]) for j in range(n_cols)))
    array_str = "np.array([[" + "],[".join(rows) + "]])"

    return _substitute_symbols(
        array_str, state_names, control_names, dilation_names, dyn_param_names,
        fmt_state=lambda i: f"hx[{i},0]",
        fmt_control=lambda i: f"hu[{i},0]",
        fmt_dilation=lambda i: f"hp[{i},0]",
        fmt_param=lambda i: f"dyn_par[{i}]",
    )


def matrix_to_c_assignments(matrix, state_names, control_names,
                            dilation_names, dyn_param_names, indent="    "):
    """Render a SymPy matrix as a block of C `out[idx] = expr;` statements."""
    matrix = sp.Matrix(matrix)
    n_rows, n_cols = matrix.shape
    if n_rows == 0:
        return f"{indent}/* empty */"

    lines = []
    for i in range(n_rows):
        for j in range(n_cols):
            expr_str = sp.ccode(matrix[i, j])
            expr_str = _substitute_symbols(
                expr_str, state_names, control_names, dilation_names, dyn_param_names,
                fmt_state=lambda k: f"hx[{k}]",
                fmt_control=lambda k: f"hu[{k}]",
                fmt_dilation=lambda k: f"hp[{k}]",
                fmt_param=lambda k: f"dyn_par[{k}]",
            )
            lines.append(f"{indent}out[{i * n_cols + j}] = {expr_str};")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Symbolic model
# ----------------------------------------------------------------------

class SymbolicModel:
    """Container for the symbolic dynamics and all of its Jacobians."""

    def __init__(self, states, controls, dynamics,
                 nonconvex_constraints=None, dynamic_parameters=None,
                 dilation=False, dilation_symbol=None):

        self.states = list(states)
        self.controls = list(controls)
        self.dynamic_parameters = list(dynamic_parameters) if dynamic_parameters else []
        self.dilation = bool(dilation)

        x_vec = sp.Matrix(self.states)
        u_vec = sp.Matrix(self.controls)
        f_user = sp.Matrix(list(dynamics))

        if f_user.shape[0] != len(self.states):
            raise ValueError(
                f"dynamics has {f_user.shape[0]} equations but there are "
                f"{len(self.states)} states"
            )

        # -- dilation parameter -------------------------------------------
        if self.dilation:
            if dilation_symbol is None:
                dilation_symbol = sp.Symbol(DILATION_SYMBOL_NAME)
            self.dilation_params = [dilation_symbol]
            p_vec = sp.Matrix(self.dilation_params)

            # If the user already wrote the dilation factor into the dynamics
            # (as the reference notebook does), do not multiply it a second
            # time.  Otherwise apply f -> p*f, i.e. dx/ds = p f(x,u).
            already_dilated = dilation_symbol in f_user.free_symbols
            self.auto_dilated = not already_dilated
            f_sym = f_user if already_dilated else dilation_symbol * f_user
        else:
            self.dilation_params = []
            self.auto_dilated = False
            p_vec = sp.Matrix([])
            f_sym = f_user

        self.f = f_sym
        self.n_p = len(self.dilation_params)

        # -- dynamics Jacobians -------------------------------------------
        self.A = f_sym.jacobian(x_vec)
        self.B = f_sym.jacobian(u_vec)
        if self.dilation:
            self.F = f_sym.jacobian(p_vec)
            self.y = f_sym - self.A @ x_vec - self.B @ u_vec - self.F @ p_vec
        else:
            self.F = sp.zeros(len(self.states), 0)
            self.y = f_sym - self.A @ x_vec - self.B @ u_vec

        # -- nonconvex constraint Jacobians -------------------------------
        if nonconvex_constraints:
            g_sym = sp.Matrix(list(nonconvex_constraints))
            self.g = g_sym
            self.C = g_sym.jacobian(x_vec)
            self.D = g_sym.jacobian(u_vec)
            if self.dilation:
                self.G = g_sym.jacobian(p_vec)
                self.z = g_sym - self.C @ x_vec - self.D @ u_vec - self.G @ p_vec
            else:
                self.G = sp.zeros(g_sym.shape[0], 0)
                self.z = g_sym - self.C @ x_vec - self.D @ u_vec
        else:
            self.g = sp.zeros(0, 1)
            self.C = sp.zeros(0, len(self.states))
            self.D = sp.zeros(0, len(self.controls))
            self.G = sp.zeros(0, max(self.n_p, 0))
            self.z = sp.zeros(0, 1)

        self.nx = len(self.states)
        self.nu = len(self.controls)
        self.ng = self.g.shape[0]

    # -- name lists used by the code emitters ------------------------------

    @property
    def state_names(self):
        return [str(s) for s in self.states]

    @property
    def control_names(self):
        return [str(s) for s in self.controls]

    @property
    def dilation_names(self):
        return [str(s) for s in self.dilation_params]

    @property
    def dyn_param_names(self):
        return [str(s) for s in self.dynamic_parameters]

    def _numpy(self, matrix):
        return matrix_to_numpy_string(matrix, self.state_names, self.control_names,
                                      self.dilation_names, self.dyn_param_names)

    def numpy_sources(self):
        """Dict of numpy array-literal sources, one per generated function."""
        return {
            'f': self._numpy(self.f),
            'A': self._numpy(self.A),
            'B': self._numpy(self.B),
            'F': self._numpy(self.F) if self.dilation else "np.zeros((%d, 0))" % self.nx,
            'y': self._numpy(self.y),
            'g': self._numpy(self.g),
            'C': self._numpy(self.C),
            'D': self._numpy(self.D),
            'G': self._numpy(self.G) if self.dilation else "np.zeros((%d, 0))" % self.ng,
            'z': self._numpy(self.z),
        }


# ----------------------------------------------------------------------
# numpy module generation
# ----------------------------------------------------------------------

_MODULE_TEMPLATE = '''\
"""Auto-generated by scvx_toolbox -- do not edit by hand."""

import math

import numpy as np


def f_dynamics(hx, hu, hp, t, dyn_par):
    return {f}


def A_jacobian(hx, hu, hp, t, dyn_par):
    return {A}


def B_jacobian(hx, hu, hp, t, dyn_par):
    return {B}


def F_jacobian(hx, hu, hp, t, dyn_par):
    return {F}


def y_affine(hx, hu, hp, t, dyn_par):
    return {y}


def g_constraints(hx, hu, hp, t, dyn_par):
    return {g}


def C_jacobian(hx, hu, hp, t, dyn_par):
    return {C}


def D_jacobian(hx, hu, hp, t, dyn_par):
    return {D}


def G_jacobian(hx, hu, hp, t, dyn_par):
    return {G}


def z_affine(hx, hu, hp, t, dyn_par):
    return {z}
'''


def generate_numpy_module(model, output_file=None):
    """
    Emit the numpy evaluation module for a SymbolicModel.

    Returns the module source; also writes it to `output_file` when given.
    """
    sources = model.numpy_sources()
    content = _MODULE_TEMPLATE.format(**sources)
    if output_file is not None:
        with open(output_file, 'w', encoding='utf-8') as fp:
            fp.write(content)
    return content
