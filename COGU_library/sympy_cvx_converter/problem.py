import cvxpy as cp


def build_problem(nx, nu, T, ng=0, convex_constraints_fn=None, extra_params=None):
    """
    Construye un problema CVXPY generico para SCVx, DPP-compliant.

    Generaliza Cell 36 del template Astrobee: las dimensiones (nx, nu, ng)
    se pasan como argumentos, el costo es estandar SCVx, y las restricciones
    convexas del usuario se inyectan via callback.

    Parametros:
        nx: int — numero de estados
        nu: int — numero de controles
        T:  int — numero de pasos de tiempo (discretizacion)
        ng: int — numero de restricciones no-convexas (0 si no hay)
        convex_constraints_fn: callable(prob_dict) -> list[cp.Constraint]
            Recibe el dict del problema y retorna restricciones convexas adicionales.
            Ejemplo Astrobee: norm(vel) <= vel_max, norm(accel) <= acc_max
        extra_params: dict {nombre: shape} — cp.Parameters adicionales que el
            usuario necesita para sus restricciones convexas.
            shape puede ser () para escalar, o (n,) / (n,m) para vectores/matrices.
            Ejemplo: {'vel_max': (), 'acc_max': (), 'torq_max': ()}

    Retorna:
        dict con:
            'problem':     cp.Problem listo para resolver
            'variables':   dict de cp.Variable {x, u, vc, vi}
            'params':      dict de cp.Parameter {A_discrete, B_discrete, ...}
            'extra_params': dict de cp.Parameter adicionales del usuario
            'nx', 'nu', 'ng', 'T': dimensiones
    """

    # ==================================================
    # VARIABLES DE OPTIMIZACION
    # ==================================================

    x = cp.Variable((nx, T + 1), name='x')
    u = cp.Variable((nu, T), name='u')
    vc = cp.Variable((nx, T), name='vc')           # virtual control (dinamica)
    vi = cp.Variable((ng, T + 1), name='vi') if ng > 0 else None

    # ==================================================
    # PARAMETROS (datos que cambian entre iteraciones SCVx)
    # ==================================================

    # Matrices de linealizacion (concatenadas horizontalmente por nodo k)
    A_discrete = cp.Parameter((nx, nx * T), name='A_discrete')
    B_discrete = cp.Parameter((nx, nu * T), name='B_discrete')
    y_discrete = cp.Parameter((nx, T), name='y_discrete')

    # Restricciones no-convexas linealizadas
    if ng > 0:
        C_discrete = cp.Parameter((ng, nx * (T + 1)), name='C_discrete')
        D_discrete = cp.Parameter((ng, nu * (T + 1)), name='D_discrete')
        z_discrete = cp.Parameter((ng, T + 1), name='z_discrete')
    else:
        C_discrete = None
        D_discrete = None
        z_discrete = None

    # Tiempo
    tau = cp.Parameter(name='tau')
    sqrt_tau = cp.Parameter(name='sqrt_tau')

    # Condiciones iniciales y finales
    start_pos = cp.Parameter((nx, 1), name='start_pos')
    end_pos = cp.Parameter((nx, 1), name='end_pos')

    # Punto de referencia (iteracion anterior)
    ox_cvxpy = cp.Parameter((nx, T + 1), name='ox_cvxpy')
    ou_cvxpy = cp.Parameter((nu, T), name='ou_cvxpy')

    # Escalado
    S_x_scaling = cp.Parameter((nx, nx), name='S_x_scaling')
    S_u_scaling = cp.Parameter((nu, nu), name='S_u_scaling')
    c_x_scaling = cp.Parameter((nx, 1), name='c_x_scaling')
    c_u_scaling = cp.Parameter((nu, 1), name='c_u_scaling')

    # Pesos SCVx
    lamb = cp.Parameter(name='lamb')
    tau_lamb = cp.Parameter(name='tau_lamb')
    etta = cp.Parameter(name='etta')

    # ==================================================
    # PARAMETROS EXTRA DEL USUARIO
    # ==================================================

    user_extra = {}
    if extra_params:
        for name, shape in extra_params.items():
            if shape == () or shape is None:
                user_extra[name] = cp.Parameter(name=name)
            else:
                user_extra[name] = cp.Parameter(shape, name=name)

    # ==================================================
    # COSTO
    # ==================================================
    # Estandar SCVx: esfuerzo de control + penalizacion virtual

    cost = 0

    for k in range(T):
        # Esfuerzo de control: ||sqrt(tau) * u_k||^2
        cost += cp.sum_squares(sqrt_tau * u[:, k:k+1])
        # Penalizacion virtual control (dinamica): tau*lamb * ||vc_k||_1
        cost += cp.norm(tau_lamb * vc[:, k:k+1], 1)

    # Penalizacion virtual inequalities (restricciones no-convexas)
    # Doble penalizacion: 2*lamb para vi vs 1*lamb para vc
    # Prioriza satisfacer restricciones no-convexas sobre cerrar gap de dinamica
    if ng > 0:
        for k in range(T + 1):
            cost += cp.norm(tau_lamb * vi[:, k:k+1], 1)
    for k in range(T + 1):
        cost += cp.norm(tau_lamb * vi[:, k:k+1], 1)

    # ==================================================
    # RESTRICCIONES
    # ==================================================

    constraints = []

    # --- Condiciones iniciales y finales ---
    constraints += [x[:, 0] == start_pos[:, 0]]
    constraints += [x[:, T] == end_pos[:, 0]]

    # --- Dinamica linealizada: x_{k+1} = A_k x_k + B_k u_k + y_k + vc_k ---
    for k in range(T):
        constraints += [
            x[:, k+1:k+2] == A_discrete[:, nx*k:nx*(k+1)] @ x[:, k:k+1]
            + B_discrete[:, nu*k:nu*(k+1)] @ u[:, k:k+1]
            + y_discrete[:, k:k+1]
            + vc[:, k:k+1]
        ]

    # --- Trust region: ||x_k - ox_k||_inf + ||u_k - ou_k||_inf <= etta ---
    for k in range(T):
        constraints += [
            cp.norm(x[:, k:k+1] - ox_cvxpy[:, k:k+1], 'inf')
            + cp.norm(u[:, k:k+1] - ou_cvxpy[:, k:k+1], 'inf')
            <= etta
        ]
    # Ultimo nodo (k=T): x usa columna T, u repite columna T-1
    constraints += [
        cp.norm(x[:, T:T+1] - ox_cvxpy[:, T:T+1], 'inf')
        + cp.norm(u[:, T-1:T] - ou_cvxpy[:, T-1:T], 'inf')
        <= etta
    ]

    # --- Restricciones no-convexas linealizadas: C_k x_k + D_k u_k + z_k <= vi_k ---
    if ng > 0:
        for k in range(T):
            constraints += [
                C_discrete[:, nx*k:nx*(k+1)] @ x[:, k:k+1]
                + D_discrete[:, nu*k:nu*(k+1)] @ u[:, k:k+1]
                + z_discrete[:, k:k+1]
                <= vi[:, k:k+1]
            ]
        # Ultimo nodo (k=T): u repite columna T-1
        constraints += [
            C_discrete[:, nx*T:nx*(T+1)] @ x[:, T:T+1]
            + D_discrete[:, nu*T:nu*(T+1)] @ u[:, T-1:T]
            + z_discrete[:, T:T+1]
            <= vi[:, T:T+1]
        ]

    # --- Restricciones convexas del usuario ---
    prob_dict = _build_prob_dict(
        nx=nx, nu=nu, ng=ng, T=T,
        x=x, u=u, vc=vc, vi=vi,
        A_discrete=A_discrete, B_discrete=B_discrete, y_discrete=y_discrete,
        C_discrete=C_discrete, D_discrete=D_discrete, z_discrete=z_discrete,
        tau=tau, sqrt_tau=sqrt_tau,
        start_pos=start_pos, end_pos=end_pos,
        ox_cvxpy=ox_cvxpy, ou_cvxpy=ou_cvxpy,
        S_x_scaling=S_x_scaling, S_u_scaling=S_u_scaling,
        c_x_scaling=c_x_scaling, c_u_scaling=c_u_scaling,
        lamb=lamb, tau_lamb=tau_lamb, etta=etta,
        extra_params=user_extra,
    )

    if convex_constraints_fn is not None:
        user_constraints = convex_constraints_fn(prob_dict)
        constraints += user_constraints

    # ==================================================
    # PROBLEMA
    # ==================================================

    objective = cp.Minimize(cost)
    problem = cp.Problem(objective, constraints)

    prob_dict['problem'] = problem
    return prob_dict


def _build_prob_dict(*, nx, nu, ng, T,
                     x, u, vc, vi,
                     A_discrete, B_discrete, y_discrete,
                     C_discrete, D_discrete, z_discrete,
                     tau, sqrt_tau,
                     start_pos, end_pos,
                     ox_cvxpy, ou_cvxpy,
                     S_x_scaling, S_u_scaling, c_x_scaling, c_u_scaling,
                     lamb, tau_lamb, etta,
                     extra_params):
    """Organiza todo en un dict accesible por el solver y el usuario."""

    return {
        # Dimensiones
        'nx': nx, 'nu': nu, 'ng': ng, 'T': T,

        # Variables
        'variables': {
            'x': x, 'u': u, 'vc': vc, 'vi': vi,
        },

        # Parametros de linealizacion
        'params': {
            'A_discrete': A_discrete,
            'B_discrete': B_discrete,
            'y_discrete': y_discrete,
            'C_discrete': C_discrete,
            'D_discrete': D_discrete,
            'z_discrete': z_discrete,
        },

        # Parametros SCVx
        'tau': tau,
        'sqrt_tau': sqrt_tau,
        'start_pos': start_pos,
        'end_pos': end_pos,
        'ox_cvxpy': ox_cvxpy,
        'ou_cvxpy': ou_cvxpy,

        # Escalado
        'S_x_scaling': S_x_scaling,
        'S_u_scaling': S_u_scaling,
        'c_x_scaling': c_x_scaling,
        'c_u_scaling': c_u_scaling,

        # Pesos
        'lamb': lamb,
        'tau_lamb': tau_lamb,
        'etta': etta,

        # Extra del usuario
        'extra_params': extra_params,
    }
