import os
import re
import sympy as sp
import numpy as np

from .cost_dsl import resolve_term


# ==========================================
# REEMPLAZO DE VARIABLES EN STRINGS
# ==========================================

def replace_variables_in_string(f_str,
                                state_parameters_string,
                                control_parameters_string,
                                dynamic_parameters_string):
    """
    Toma un string con nombres de variables simbolicas y los reemplaza
    por indices de array numpy.

    Ejemplo: 'rx' → 'hx[0,0]', 'ax' → 'hu[0,0]', 'Ixx' → 'dyn_par[0]'
    """

    # Reemplazar estados
    for i, var in enumerate(state_parameters_string):
        f_str = re.sub(rf"\b{var}\b", f"hx[{i},0]", f_str)

    # Reemplazar controles
    for i, var in enumerate(control_parameters_string):
        f_str = re.sub(rf"\b{var}\b", f"hu[{i},0]", f_str)

    # Reemplazar parametros dinamicos
    for i, var in enumerate(dynamic_parameters_string):
        f_str = re.sub(rf"\b{var}\b", f"dyn_par[{i}]", f_str)

    return f_str


# ==========================================
# CONVERSION DE MATRIZ SYMPY A STRING NUMPY
# ==========================================

def matrix_to_numpy_string(matrix, state_params, control_params, dynamic_params):
    """
    Convierte una matriz simbolica SymPy a un string que representa un np.array.

    Ejemplo:
        sp.Matrix([[x + y], [x*y]])  →  "np.array([[hx[0,0] + hx[1,0]],[hx[0,0]*hx[1,0]]])"

    Usa sp.pycode() para convertir cada elemento a Python,
    luego replace_variables_in_string para cambiar nombres por indices.
    """

    n_rows, *rest = np.shape(matrix)
    n_cols = rest[0] if rest else 1

    matrix = sp.Matrix(n_rows, n_cols, lambda i, j: matrix[i * n_cols + j])

    array_str = "np.array([["

    for i in range(n_rows):
        for j in range(n_cols):
            array_str += sp.pycode(matrix[i, j])
            if j < n_cols - 1:
                array_str += ","
        if i < n_rows - 1:
            array_str += "],["
        else:
            array_str += "]])"

    array_str = replace_variables_in_string(array_str,
                                            state_params,
                                            control_params,
                                            dynamic_params)
    return array_str


# ==========================================
# GENERACION COMPLETA DE FUNCIONES
# ==========================================

def generate_functions(states, controls, dynamics, nonconvex_constraints=None,
                       dynamic_parameters=None, output_file=None):
    """
    Recibe el modelo simbolico y genera funciones numpy para f, A, B, y, (C, D, z).

    Parametros:
        states:    lista de simbolos SymPy de estados [rx, ry, vx, vy, ...]
        controls:  lista de simbolos SymPy de controles [ax, ay, ...]
        dynamics:  lista/matriz SymPy con las ecuaciones dx/dt = f(x,u)
        nonconvex_constraints: lista/matriz SymPy con g(x,u) <= 0 (opcional)
        dynamic_parameters: lista de simbolos SymPy de parametros dinamicos (opcional)
        output_file: ruta del archivo .py a generar (opcional, si None no escribe archivo)

    Retorna:
        dict con los strings generados: 'f', 'A', 'B', 'y', ('g', 'C', 'D', 'z' si hay constraints)
    """

    # Convertir a matrices SymPy
    x_states = sp.Matrix(states)
    u_input = sp.Matrix(controls)
    f = sp.Matrix(dynamics)

    # Listas de nombres como strings (para replace_variables_in_string)
    state_str = [str(s) for s in states]
    control_str = [str(s) for s in controls]
    dyn_par_str = [str(p) for p in dynamic_parameters] if dynamic_parameters else []

    # --- Dinamica: f ---
    string_f = matrix_to_numpy_string(f, state_str, control_str, dyn_par_str)

    # --- Jacobianos de la dinamica: A = df/dx, B = df/du ---
    A_symbolic = f.jacobian(x_states)
    B_symbolic = f.jacobian(u_input)

    string_A = matrix_to_numpy_string(A_symbolic, state_str, control_str, dyn_par_str)
    string_B = matrix_to_numpy_string(B_symbolic, state_str, control_str, dyn_par_str)

    # --- Termino afin: y = f - A*x - B*u ---
    y_symbolic = f - A_symbolic @ x_states - B_symbolic @ u_input

    string_y = matrix_to_numpy_string(y_symbolic, state_str, control_str, dyn_par_str)

    result = {
        'f': string_f,
        'A': string_A,
        'B': string_B,
        'y': string_y,
    }

    # --- Restricciones no convexas (si existen) ---
    string_g = "0"
    string_C = "0"
    string_D = "0"
    string_z = "0"

    if nonconvex_constraints is not None and len(nonconvex_constraints) > 0:
        g = sp.Matrix(nonconvex_constraints)

        string_g = matrix_to_numpy_string(g, state_str, control_str, dyn_par_str)

        C_symbolic = g.jacobian(x_states)
        D_symbolic = g.jacobian(u_input)
        z_symbolic = g - C_symbolic @ x_states - D_symbolic @ u_input

        string_C = matrix_to_numpy_string(C_symbolic, state_str, control_str, dyn_par_str)
        string_D = matrix_to_numpy_string(D_symbolic, state_str, control_str, dyn_par_str)
        string_z = matrix_to_numpy_string(z_symbolic, state_str, control_str, dyn_par_str)

    result['g'] = string_g
    result['C'] = string_C
    result['D'] = string_D
    result['z'] = string_z

    # --- Escribir archivo .py si se pide ---
    if output_file is not None:
        file_content = f"""\
import numpy as np
import math

def f_SCVx(hx, hu, t, dyn_par):
    return {string_f}

def g_SCVx(hx, hu, t, dyn_par):
    return {string_g}

def A_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_A}

def B_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_B}

def y_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_y}

def C_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_C}

def D_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_D}

def z_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return {string_z}
"""
        with open(output_file, "w", encoding="utf-8") as fp:
            fp.write(file_content)

    return result


# ==========================================
# HELPERS PARA GENERACION EN C
# ==========================================

def _replace_variables_c(expr_str, state_params, control_params, dynamic_params):
    """Reemplaza nombres de simbolos SymPy por indexado plano C (hx[i], hu[i], dyn_par[i])."""
    for i, var in enumerate(state_params):
        expr_str = re.sub(rf"\b{var}\b", f"hx[{i}]", expr_str)
    for i, var in enumerate(control_params):
        expr_str = re.sub(rf"\b{var}\b", f"hu[{i}]", expr_str)
    for i, var in enumerate(dynamic_params):
        expr_str = re.sub(rf"\b{var}\b", f"dyn_par[{i}]", expr_str)
    return expr_str


def _matrix_to_c_assignments(matrix, state_params, control_params, dynamic_params):
    """Convierte una matriz SymPy a sentencias C 'out[idx] = expr;' con sp.ccode()."""
    n_rows, *rest = np.shape(matrix)
    n_cols = rest[0] if rest else 1
    matrix = sp.Matrix(n_rows, n_cols, lambda i, j: matrix[i * n_cols + j])

    lines = []
    for i in range(n_rows):
        for j in range(n_cols):
            idx = i * n_cols + j
            expr_str = sp.ccode(matrix[i, j])
            expr_str = _replace_variables_c(expr_str, state_params, control_params, dynamic_params)
            lines.append(f"    out[{idx}] = {expr_str};")
    return "\n".join(lines)


def generate_c_user_cost(cost_terms, nx, nu):
    """
    Genera el bloque C del costo del usuario a partir de cost_terms.

    Tercera materializacion del mismo DSL: espejo de _build_cost_from_terms
    (CVXPY, problem.py) y eval_user_cost (numpy, solver.py).

    El codigo generado asume en scope las variables C:
        x  (estados, stride ns=STATES_SIZE)
        u  (controles, stride mi=INPUTS_SIZE)
        T, tau, lam, cost (acumulador)
    y declara sqrt_tau / tau_lamb solo si algun termino los usa.

    Retorna un string C (a insertar en J_cost vía el placeholder {USER_COST}).
    """
    coeff_c = {'sqrt_tau': 'sqrt_tau', 'tau': 'tau', 'tau_lamb': 'tau_lamb'}
    resolved = [resolve_term(t) for t in cost_terms]

    # Declarar solo los coeficientes string que se usan (evita warnings unused)
    used = {rt.coeff for rt in resolved if isinstance(rt.coeff, str)}
    lines = []
    if 'sqrt_tau' in used:
        lines.append('    double sqrt_tau=sqrt(tau);')
    if 'tau_lamb' in used:
        lines.append('    double tau_lamb=tau*lam;')

    off_idx = 0
    for rt in resolved:
        arr, stride, dim = ('u', 'mi', nu) if rt.var == 'u' else ('x', 'ns', nx)
        lo = 0 if rt.slc.start is None else rt.slc.start
        hi = dim if rt.slc.stop is None else rt.slc.stop
        coeff = coeff_c[rt.coeff] if isinstance(rt.coeff, str) else repr(float(rt.coeff))
        weight = repr(rt.weight)
        kcond = 'k<=T' if rt.k_range == 'T+1' else 'k<T'

        if rt.offset is not None:
            off = np.asarray(rt.offset, dtype=float).reshape(-1)
            off_name = f'_off{off_idx}'
            off_idx += 1
            off_vals = ','.join(repr(float(v)) for v in off)
            lines.append(f'    static const double {off_name}[]={{{off_vals}}};')
            e_expr = f'{coeff}*({arr}k[i]-{off_name}[i-{lo}])'
        else:
            e_expr = f'{coeff}*{arr}k[i]'

        accum = 's+=e*e;' if rt.kind in ('sumsq', 'norm2') else 's+=fabs(e);'
        lines.append(f'    for(int k=0;{kcond};k++){{')
        lines.append(f'        const double*{arr}k={arr}+k*{stride}; double s=0;')
        lines.append(f'        for(int i={lo};i<{hi};i++){{ double e={e_expr}; {accum} }}')
        if rt.kind == 'norm2':
            lines.append(f'        cost+={weight}*sqrt(s);')
        else:   # sumsq, norm1
            lines.append(f'        cost+={weight}*s;')
        lines.append('    }')

    return '\n'.join(lines)


def _c_array(name, flat_values):
    """Declara un array C const: 'static const double name[]={...};'"""
    body = ','.join(repr(float(v)) for v in flat_values)
    return f'    static const double {name}[]={{{body}}};'


def generate_c_embedded_data(embed, nx, nu, T):
    """
    Genera el bloque C de datos del problema embebidos (Fase 6).

    `embed` es un dict con (todo SIN escalar, valores que el usuario pasa en Python):
        tf, lam (escalares); start, end (nx,1); dp (np,);
        S_x (nx,nx), c_x (nx,1), S_u (nu,nu), c_u (nu,1);
        warm_x (nx,T+1), warm_u (nu,T).

    Layouts (deben coincidir con el uso del template):
        S_x/S_u row-major (Sx[i*ns+j]); warm time-major (ohx[i+t*ns]).

    Retorna un string C (declaraciones a insertar vía {EMBEDDED_DATA}).
    """
    start = np.asarray(embed['start'], dtype=float).reshape(-1)
    end   = np.asarray(embed['end'],   dtype=float).reshape(-1)
    dp    = np.asarray(embed['dp'],    dtype=float).reshape(-1)
    Sx    = np.asarray(embed['S_x'],   dtype=float)
    cx    = np.asarray(embed['c_x'],   dtype=float).reshape(-1)
    Su    = np.asarray(embed['S_u'],   dtype=float)
    cu    = np.asarray(embed['c_u'],   dtype=float).reshape(-1)
    wx    = np.asarray(embed['warm_x'], dtype=float)   # (nx, T+1)
    wu    = np.asarray(embed['warm_u'], dtype=float)   # (nu, T)

    # warm en layout time-major: idx = i + t*ns  ->  recorrer t, luego i
    wx_tm = [wx[i, t] for t in range(T + 1) for i in range(nx)]
    wu_tm = [wu[i, t] for t in range(T)     for i in range(nu)]

    lines = [
        f'    double tf={repr(float(embed["tf"]))};',
        f'    double lam={repr(float(embed["lam"]))};',
        _c_array('_start', start),
        _c_array('_end',   end),
        _c_array('_dp',    dp),
        _c_array('_Sx',    Sx.reshape(-1)),   # row-major
        _c_array('_cx',    cx),
        _c_array('_Su',    Su.reshape(-1)),   # row-major
        _c_array('_cu',    cu),
        _c_array('_warm_x', wx_tm),
        _c_array('_warm_u', wu_tm),
    ]
    return '\n'.join(lines)


# ==========================================
# GENERACION DE CODIGO C
# ==========================================

def generate_c_functions(states, controls, dynamics, nonconvex_constraints=None,
                         dynamic_parameters=None, output_dir=None):
    """
    Genera dynamics.c y dynamics.h a partir del modelo simbolico SymPy.

    Funciones generadas (todas con firma identica):
        void name(double *hx, double *hu, double t, double *dyn_par, double *out)

        f_dynamics    — f(x,u):        out size nx
        A_jacobian    — df/dx:         out size nx*nx (row-major)
        B_jacobian    — df/du:         out size nx*nu (row-major)
        y_affine      — f - A*x - B*u: out size nx
        g_constraints — g(x,u):        out size ng  (0 si sin restricciones)
        C_jac         — dg/dx:         out size ng*nx
        D_jac         — dg/du:         out size ng*nu
        z_affine      — g - C*x - D*u: out size ng

    Parametros:
        states:                lista de simbolos SymPy de estados
        controls:              lista de simbolos SymPy de controles
        dynamics:              lista/matriz SymPy con dx/dt = f(x,u)
        nonconvex_constraints: lista/matriz SymPy con g(x,u) <= 0 (opcional)
        dynamic_parameters:    lista de simbolos SymPy de parametros (opcional)
        output_dir:            directorio donde escribir dynamics.c y dynamics.h
                               (si None, no escribe archivos)

    Retorna:
        dict con claves 'c_source' y 'h_source' (strings del contenido generado)
    """

    x_states = sp.Matrix(states)
    u_input = sp.Matrix(controls)
    f = sp.Matrix(dynamics)

    state_str = [str(s) for s in states]
    control_str = [str(s) for s in controls]
    dyn_par_str = [str(p) for p in dynamic_parameters] if dynamic_parameters else []

    nx = len(states)
    nu = len(controls)

    # --- Dinamica f ---
    body_f = _matrix_to_c_assignments(f, state_str, control_str, dyn_par_str)

    # --- Jacobianos A = df/dx, B = df/du ---
    A_sym = f.jacobian(x_states)
    B_sym = f.jacobian(u_input)
    body_A = _matrix_to_c_assignments(A_sym, state_str, control_str, dyn_par_str)
    body_B = _matrix_to_c_assignments(B_sym, state_str, control_str, dyn_par_str)

    # --- Termino afin y = f - A*x - B*u ---
    y_sym = f - A_sym @ x_states - B_sym @ u_input
    body_y = _matrix_to_c_assignments(y_sym, state_str, control_str, dyn_par_str)

    # --- Restricciones no convexas (opcionales) ---
    ng = 0
    if nonconvex_constraints is not None and len(nonconvex_constraints) > 0:
        g = sp.Matrix(nonconvex_constraints)
        ng = g.shape[0]
        C_sym = g.jacobian(x_states)
        D_sym = g.jacobian(u_input)
        z_sym = g - C_sym @ x_states - D_sym @ u_input
        body_g = _matrix_to_c_assignments(g, state_str, control_str, dyn_par_str)
        body_C = _matrix_to_c_assignments(C_sym, state_str, control_str, dyn_par_str)
        body_D = _matrix_to_c_assignments(D_sym, state_str, control_str, dyn_par_str)
        body_z = _matrix_to_c_assignments(z_sym, state_str, control_str, dyn_par_str)
    else:
        body_g = "    /* no nonconvex constraints */"
        body_C = "    /* no nonconvex constraints */"
        body_D = "    /* no nonconvex constraints */"
        body_z = "    /* no nonconvex constraints */"

    sig = "double *hx, double *hu, double t, double *dyn_par, double *out"

    c_content = f"""\
#include <math.h>
#include "dynamics.h"

void f_dynamics({sig}) {{
{body_f}
}}

void A_jacobian({sig}) {{
{body_A}
}}

void B_jacobian({sig}) {{
{body_B}
}}

void y_affine({sig}) {{
{body_y}
}}

void g_constraints({sig}) {{
{body_g}
}}

void C_jac({sig}) {{
{body_C}
}}

void D_jac({sig}) {{
{body_D}
}}

void z_affine({sig}) {{
{body_z}
}}
"""

    h_content = f"""\
#ifndef DYNAMICS_H
#define DYNAMICS_H

/* nx={nx}, nu={nu}, ng={ng} */

void f_dynamics(double *hx, double *hu, double t, double *dyn_par, double *out);
void A_jacobian(double *hx, double *hu, double t, double *dyn_par, double *out);
void B_jacobian(double *hx, double *hu, double t, double *dyn_par, double *out);
void y_affine(double *hx, double *hu, double t, double *dyn_par, double *out);
void g_constraints(double *hx, double *hu, double t, double *dyn_par, double *out);
void C_jac(double *hx, double *hu, double t, double *dyn_par, double *out);
void D_jac(double *hx, double *hu, double t, double *dyn_par, double *out);
void z_affine(double *hx, double *hu, double t, double *dyn_par, double *out);

#endif /* DYNAMICS_H */
"""

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, "dynamics.c"), "w", encoding="utf-8") as fp:
            fp.write(c_content)
        with open(os.path.join(output_dir, "dynamics.h"), "w", encoding="utf-8") as fp:
            fp.write(h_content)

    return {
        'c_source': c_content,
        'h_source': h_content,
    }
