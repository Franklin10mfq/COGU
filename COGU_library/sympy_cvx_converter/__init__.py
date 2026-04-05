# sympy_cvx_converter
# Traductor automatico SymPy → CVXPY + pipeline SCVx completo

from .translator import (
    sympy_to_cvx,
    convert_constraint,
    create_variables,
    detect_quad_form,
    cvx_norm,
    cvx_sum_squares,
    cvx_huber,
)

from .codegen import (
    generate_functions,
    replace_variables_in_string,
    matrix_to_numpy_string,
)

from .scaling import (
    scale_x, scale_u,
    inv_scale_x, inv_scale_u,
    scale_A, scale_B, scale_y,
    scale_C, scale_D, scale_z,
)

from .problem import build_problem

from .solver import (
    discretize_ABy,
    discretize_CDz,
    rk4_step,
    compute_J,
    solve_scvx,
)
