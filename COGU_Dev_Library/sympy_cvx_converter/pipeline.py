import numpy as np
from cvxpygen import cpg
import importlib.util
import tempfile
import os
from .codegen import generate_functions
from .problem import build_problem
from .solver import solve_scvx


def solve_trajectory(states, controls, dynamics, start, end, T, tf,
                     nonconvex_constraints=None, dynamic_parameters_sym=None,
                     dynamic_parameters_val=None,
                     S_x=None, c_x=None, S_u=None, c_u=None,
                     warm_start_x=None, warm_start_u=None,
                     state_bounds=None, control_bounds=None,
                     size_N=20, lamb=1000.0, max_iter=25,
                     rho0=0.0, rho1=0.1, rho2=0.7,
                     etta0=1e-8, etta1=10.0, etta_init=1.0,
                     beta_sh=2.0, beta_gr=2.0,
                     e_tol=0.005, epsilon_stop=1e-5,
                     solver='ECOS', verbose=True,
                     generate_c=False, c_output_dir="c_output"):
# AGREGADOS GENERATE C & C OUTPUT DIR 
    """
    Resuelve un problema de trayectoria optima via SCVx.
    El usuario define todo en SymPy y numpy — la libreria hace el resto.

    Parametros requeridos:
        states:    lista de simbolos SymPy de estados [rx, vx, ...]
        controls:  lista de simbolos SymPy de controles [ax, ...]
        dynamics:  lista SymPy dx/dt = f(x, u) [vx, ax, ...]
        start:     condicion inicial (nx, 1) numpy array
        end:       condicion final (nx, 1) numpy array
        T:         numero de pasos de tiempo
        tf:        tiempo final

    Parametros opcionales:
        nonconvex_constraints:  lista SymPy g(x,u) <= 0
        dynamic_parameters_sym: lista de simbolos SymPy de parametros dinamicos
        dynamic_parameters_val: numpy array con valores de los parametros
        S_x, c_x, S_u, c_u:   matrices de escalado (default: identidad)
        warm_start_x:          trayectoria inicial estados (nx, T+1)
        warm_start_u:          trayectoria inicial controles (nu, T)
        state_bounds:          lista de (slice, 'norm1'|'norm2'|'norminf', limit)
        control_bounds:        lista de (slice, 'norm1'|'norm2'|'norminf', limit)
        size_N, lamb, max_iter, rho0-2, etta0-1, beta_sh/gr, e_tol, epsilon_stop:
                               parametros SCVx (defaults del template)
        solver:                solver CVXPY (default 'ECOS')
        verbose:               imprimir progreso (default True)

    Retorna:
        dict con x_opt, u_opt (unscaled), converged, iterations, history
    """

    nx = len(states)
    nu = len(controls)
    ng = len(nonconvex_constraints) if nonconvex_constraints else 0

    # ── [A] Code generation ──
    codegen_file = tempfile.NamedTemporaryFile(
        suffix='.py', delete=False, dir='.', prefix='scvx_codegen')
    codegen_file.close()
    codegen_path = codegen_file.name

    try:
        generate_functions(
            states=states, controls=controls, dynamics=dynamics,
            nonconvex_constraints=nonconvex_constraints,
            dynamic_parameters=dynamic_parameters_sym,
            output_file=codegen_path,
        )

        spec = importlib.util.spec_from_file_location("_scvx_codegen", codegen_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        funcs = {
            'f': mod.f_SCVx,
            'A': mod.A_no_discrete_no_scaled_SCVx,
            'B': mod.B_no_discrete_no_scaled_SCVx,
            'y': mod.y_no_discrete_no_scaled_SCVx,
        }
        if ng > 0:
            funcs['g'] = mod.g_SCVx
            funcs['C'] = mod.C_no_discrete_no_scaled_SCVx
            funcs['D'] = mod.D_no_discrete_no_scaled_SCVx
            funcs['z'] = mod.z_no_discrete_no_scaled_SCVx

        # ── [B] Auto-scaling ──
        if S_x is None:
            S_x = np.eye(nx)
        if c_x is None:
            c_x = np.zeros((nx, 1))
        if S_u is None:
            S_u = np.eye(nu)
        if c_u is None:
            c_u = np.zeros((nu, 1))

        # ── [C] Auto warm-start ──
        if warm_start_x is None:
            warm_start_x = np.zeros((nx, T + 1))
            for k in range(T + 1):
                warm_start_x[:, k:k+1] = start + (end - start) * k / T
        if warm_start_u is None:
            warm_start_u = np.zeros((nu, T))

        # ── [D] Build problem ──
        tau_val = tf / T

        prob_dict = build_problem(
            nx=nx, nu=nu, T=T, ng=ng,
            tau_val=tau_val, lamb_val=lamb,
            S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
            state_bounds=state_bounds,
            control_bounds=control_bounds,
        )
        # ── [D.5] Generar solver C con CVXPYgen (opcional) ──
        if generate_c:
            solver_dir = os.path.join(c_output_dir, 'solver')
            os.makedirs(solver_dir, exist_ok=True)
            try:
                cpg.generate_code(prob_dict['problem'], code_dir=solver_dir, solver='ECOS')
            except MemoryError as e:
                raise MemoryError(
                    f"CVXPYgen: problema demasiado grande para CVXPY 1.8 DPP expansion "
                    f"(T={T}, nx={nx}, nu={nu}). Usar T<=10 con nx<=6 para generar C. "
                    f"Original: {e}"
                ) from e
            except (ModuleNotFoundError, Exception) as e:
                # C files generated OK; Python wrapper compilation failed (CMake missing, etc.)
                c_src = os.path.join(solver_dir, 'c', 'src')
                if os.path.isdir(c_src):
                    print(f"[CVXPYgen] C code generated at {solver_dir}/c/")
                    print(f"[CVXPYgen] Python wrapper compilation failed: {e}")
                    print("[CVXPYgen] Para compilar: instalar CMake y correr setup.py en solver/")
                else:
                    raise

        # ── [E] Solve SCVx ──
        dyn_par = dynamic_parameters_val if dynamic_parameters_val is not None else np.array([])

        result = solve_scvx(prob_dict, funcs, config={
            'tf': tf,
            'size_N': size_N,
            'dynamic_parameters': dyn_par,
            'S_x': S_x, 'c_x': c_x,
            'S_u': S_u, 'c_u': c_u,
            'start_pos': start,
            'end_pos': end,
            'warm_start_x': warm_start_x,
            'warm_start_u': warm_start_u,
            'rho0': rho0, 'rho1': rho1, 'rho2': rho2,
            'etta0': etta0, 'etta1': etta1, 'etta_init': etta_init,
            'beta_sh': beta_sh, 'beta_gr': beta_gr,
            'lamb': lamb,
            'e_tol': e_tol,
            'epsilon_stop': epsilon_stop,
            'max_iter': max_iter,
            'solver': solver,
            'verbose': verbose,
            'extra_param_values': {},
        })

        return result

    finally:
        # Limpiar archivo temporal
        try:
            os.unlink(codegen_path)
        except OSError:
            pass