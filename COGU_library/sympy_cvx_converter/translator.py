import sympy as sp
from sympy.matrices.expressions.slice import MatrixSlice
from sympy.matrices.expressions.matexpr import MatrixElement
import cvxpy as cp
import numpy as np


# ==========================================
# FUNCIONES SYMPY PERSONALIZADAS
# ==========================================
# Estas funciones no existen nativamente en SymPy.
# El usuario las usa en expresiones SymPy y el converter
# las traduce a las funciones CVXPY correspondientes.

cvx_norm = sp.Function('norm')
cvx_sum_squares = sp.Function('sum_squares')
cvx_huber = sp.Function('huber')


# ==========================================
# CREACION DE VARIABLES
# ==========================================

def create_variables(states, controls, T, parameters=None, constants=None):
    """
    Crea variables CVXPY a partir de la definicion del problema.
    Los tamaños se deducen automaticamente de states, controls y T.

    - states: lista de simbolos SymPy de estados
    - controls: lista de simbolos SymPy de controles
    - T: numero de pasos de tiempo
    - parameters: lista de simbolos → cp.Parameter (datos que cambian entre iteraciones)
    - constants: dict {simbolo: valor} → float (valores fijos)
    """

    var_map = {}

    nx = len(states)
    nu = len(controls)

    # 1. Constantes (valores fijos)
    if constants:
        for v, val in constants.items():
            var_map[v] = float(val)

    # 2. Parametros (cambian entre iteraciones, no se optimizan)
    if parameters:
        for v in parameters:
            var_map[v] = cp.Parameter(name=str(v))

    # 3. Variables de optimizacion (tamaños deducidos)
    var_map['x'] = cp.Variable((nx, T + 1), name='x')
    var_map['u'] = cp.Variable((nu, T), name='u')
    var_map['vc'] = cp.Variable((nx, T), name='vc')

    # 4. Parametros matriciales (tamaños deducidos)
    var_map['A_discrete'] = cp.Parameter((nx, nx * T), name='A_discrete')
    var_map['B_discrete'] = cp.Parameter((nx, nu * T), name='B_discrete')
    var_map['y_discrete'] = cp.Parameter((nx, T), name='y_discrete')

    return var_map


# ==========================================
# DETECCION x^T Q x
# ==========================================

def detect_quad_form(expr, var_map):

    if not isinstance(expr, sp.MatMul):
        return None

    args = expr.args

    if len(args) != 3:
        return None

    a, b, c = args

    if isinstance(a, sp.Transpose) and isinstance(c, sp.MatrixSymbol):

        x = a.args[0]

        if x == c:

            if isinstance(b, sp.MatrixSymbol):

                return cp.quad_form(var_map[x], np.eye(var_map[x].shape[0]))

    return None


# ==========================================
# CONVERSION PRINCIPAL
# ==========================================

def sympy_to_cvx(expr, var_map):

    if expr.is_Number:
        return float(expr)

    if expr.is_Symbol:
        return var_map[expr]

    # -------------------
    # SUMA
    # -------------------

    if isinstance(expr, sp.Add):

        return sum(sympy_to_cvx(a, var_map) for a in expr.args)

    # -------------------
    # SUMA MATRICIAL (MatAdd)
    # -------------------

    if isinstance(expr, sp.MatAdd):

        return sum(sympy_to_cvx(a, var_map) for a in expr.args)

    # -------------------
    # PATRONES DE MULTIPLICACION (antes de Mul generico)
    # -------------------

    # Patron: -x * log(x) → cp.entr(x)
    if isinstance(expr, sp.Mul):
        args = expr.args
        if len(args) == 3:
            has_neg1 = sp.S.NegativeOne in args
            log_args = [a for a in args if isinstance(a, sp.log)]
            other = [a for a in args if a != sp.S.NegativeOne and not isinstance(a, sp.log)]
            if has_neg1 and len(log_args) == 1 and len(other) == 1:
                x_sym = other[0]
                log_inner = log_args[0].args[0]
                if x_sym == log_inner:
                    return cp.entr(sympy_to_cvx(x_sym, var_map))

    # Patron: x * exp(x) → cp.xexp(x)
    if isinstance(expr, sp.Mul):
        args = expr.args
        if len(args) == 2:
            a, b = args
            if isinstance(b, sp.exp) and a == b.args[0]:
                return cp.xexp(sympy_to_cvx(a, var_map))
            if isinstance(a, sp.exp) and b == a.args[0]:
                return cp.xexp(sympy_to_cvx(b, var_map))

    # -------------------
    # MULTIPLICACION GENERICA
    # -------------------

    if isinstance(expr, sp.Mul):

        result = sympy_to_cvx(expr.args[0], var_map)

        for a in expr.args[1:]:

            result = result * sympy_to_cvx(a, var_map)

        return result

    # -------------------
    # POTENCIAS
    # -------------------

    if isinstance(expr, sp.Pow):

        base, exponent = expr.args

        if exponent == 2:
            return cp.square(sympy_to_cvx(base, var_map))

        if exponent == sp.Rational(1, 2):
            return cp.sqrt(sympy_to_cvx(base, var_map))

        # sp.Pow(x, -1) → cp.inv_pos(x)
        if exponent == -1:
            return cp.inv_pos(sympy_to_cvx(base, var_map))

        # sp.Pow(x, p) general → cp.power(x, p)
        return cp.power(sympy_to_cvx(base, var_map), float(exponent))

    # -------------------
    # FUNCIONES
    # -------------------

    if isinstance(expr, sp.Abs):

        return cp.abs(sympy_to_cvx(expr.args[0], var_map))

    # --- Patron: log(1 + exp(x)) → cp.logistic(x) [ANTES de log generico] ---

    if expr.func == sp.log:
        inner = expr.args[0]
        if isinstance(inner, sp.Add) and len(inner.args) == 2:
            a, b = inner.args
            if a == 1 and isinstance(b, sp.exp):
                return cp.logistic(sympy_to_cvx(b.args[0], var_map))
            if b == 1 and isinstance(a, sp.exp):
                return cp.logistic(sympy_to_cvx(a.args[0], var_map))

    # --- Patron: log(1 + x) → cp.log1p(x) [ANTES de log generico] ---

    if expr.func == sp.log:
        inner = expr.args[0]
        if isinstance(inner, sp.Add) and len(inner.args) == 2:
            a, b = inner.args
            if a == 1:
                return cp.log1p(sympy_to_cvx(b, var_map))
            if b == 1:
                return cp.log1p(sympy_to_cvx(a, var_map))

    if expr.func == sp.exp:

        return cp.exp(sympy_to_cvx(expr.args[0], var_map))

    if expr.func == sp.log:

        return cp.log(sympy_to_cvx(expr.args[0], var_map))

    if expr.func == sp.sqrt:

        return cp.sqrt(sympy_to_cvx(expr.args[0], var_map))

    # -------------------
    # MAX / MIN (con deteccion de pos/neg)
    # -------------------

    if isinstance(expr, sp.Max):
        args = expr.args
        if len(args) == 2:
            a, b = args
            if b == 0 or b is sp.S.Zero:
                return cp.pos(sympy_to_cvx(a, var_map))
            if a == 0 or a is sp.S.Zero:
                return cp.pos(sympy_to_cvx(b, var_map))
        cvx_args = [sympy_to_cvx(a, var_map) for a in args]
        result = cvx_args[0]
        for ca in cvx_args[1:]:
            result = cp.maximum(result, ca)
        return result

    if isinstance(expr, sp.Min):
        args = expr.args
        cvx_args = [sympy_to_cvx(a, var_map) for a in args]
        result = cvx_args[0]
        for ca in cvx_args[1:]:
            result = cp.minimum(result, ca)
        return result

    # -------------------
    # TRACE
    # -------------------

    if isinstance(expr, sp.Trace):
        return cp.trace(sympy_to_cvx(expr.args[0], var_map))

    # -------------------
    # MATMUL
    # -------------------

    if isinstance(expr, sp.MatMul):
        q = detect_quad_form(expr, var_map)
        if q is not None:
            return q
        result = sympy_to_cvx(expr.args[0], var_map)
        for a in expr.args[1:]:
            result = result @ sympy_to_cvx(a, var_map)
        return result

    # -------------------
    # SLICING 2D: X[3:6, k:k+1]
    # -------------------

    if isinstance(expr, MatrixSlice):
        parent_cvx = sympy_to_cvx(expr.parent, var_map)
        r_start, r_stop, r_step = expr.rowslice
        c_start, c_stop, c_step = expr.colslice
        return parent_cvx[r_start:r_stop:r_step, c_start:c_stop:c_step]

    # -------------------
    # ELEMENTO INDIVIDUAL: X[i, j]
    # -------------------

    if isinstance(expr, MatrixElement):
        parent_cvx = sympy_to_cvx(expr.parent, var_map)
        i = int(expr.i)
        j = int(expr.j)
        return parent_cvx[i, j]

    # -------------------
    # MATRICES
    # -------------------

    if isinstance(expr, sp.MatrixSymbol):
        if expr in var_map:
            return var_map[expr]
        raise ValueError(
            f"MatrixSymbol '{expr}' no encontrado en var_map."
        )

    # -------------------
    # NORMA (patron sqrt(x^2 + y^2 + ...))
    # -------------------

    if expr.func == sp.sqrt and isinstance(expr.args[0], sp.Add):

        inside = expr.args[0]

        terms = []

        for t in inside.args:

            if isinstance(t, sp.Pow) and t.args[1] == 2:
                terms.append(sympy_to_cvx(t.args[0], var_map))

        if len(terms) > 1:

            vec = cp.hstack(terms)

            return cp.norm(vec)

    # -------------------
    # FUNCIONES PERSONALIZADAS
    # -------------------

    # norm(x, p) → cp.norm(x, p)
    if hasattr(expr, 'func') and hasattr(expr.func, 'name') and expr.func.name == 'norm':
        args = expr.args
        x_cvx = sympy_to_cvx(args[0], var_map)
        if len(args) >= 2:
            p_val = args[1]
            if p_val == sp.oo:
                return cp.norm(x_cvx, 'inf')
            return cp.norm(x_cvx, int(p_val))
        return cp.norm(x_cvx)

    # sum_squares(x) → cp.sum_squares(x)
    if hasattr(expr, 'func') and hasattr(expr.func, 'name') and expr.func.name == 'sum_squares':
        return cp.sum_squares(sympy_to_cvx(expr.args[0], var_map))

    # huber(x, M) → cp.huber(x, M)
    if hasattr(expr, 'func') and hasattr(expr.func, 'name') and expr.func.name == 'huber':
        args = expr.args
        x_cvx = sympy_to_cvx(args[0], var_map)
        if len(args) >= 2:
            return cp.huber(x_cvx, float(args[1]))
        return cp.huber(x_cvx)

    # -------------------
    # QUADRATIC FORM
    # -------------------

    q = detect_quad_form(expr, var_map)

    if q is not None:
        return q

    raise NotImplementedError(f"Expresion no soportada: {expr}")


# ==========================================
# CONVERSION DE RESTRICCIONES
# ==========================================

def convert_constraint(constr, var_map):

    if isinstance(constr, sp.LessThan):

        return sympy_to_cvx(constr.lhs, var_map) <= sympy_to_cvx(constr.rhs, var_map)

    if isinstance(constr, sp.GreaterThan):

        return sympy_to_cvx(constr.lhs, var_map) >= sympy_to_cvx(constr.rhs, var_map)

    if isinstance(constr, sp.Equality):

        return sympy_to_cvx(constr.lhs, var_map) == sympy_to_cvx(constr.rhs, var_map)

    raise NotImplementedError("Tipo de restriccion no soportado")
