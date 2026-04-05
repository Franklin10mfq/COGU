import re
import sympy as sp
import numpy as np


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
