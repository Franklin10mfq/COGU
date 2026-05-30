import numpy as np
import shutil
import importlib.util
import tempfile
import os
import sys
import pickle
import subprocess
from .codegen import generate_functions, generate_c_functions
from .template_filler import fill_scvx_template
from .problem import build_problem
from .solver import solve_scvx


def _generate_solver_subprocess(build_kwargs, solver_dir):
    """
    Genera el solver C llamando a CVXPYgen en un SUBPROCESO limpio.

    Por que subproceso:
        CVXPYgen necesita un pico de RAM contiguo de varios cientos de MB en la
        canonicalizacion DPP. En el proceso principal, SymPy y el solve SCVx ya
        ocuparon y fragmentaron la memoria -> MemoryError. Un proceso hijo
        arranca limpio y dispone del bloque contiguo completo.

    Parametros:
        build_kwargs: dict con los argumentos de build_problem() — todos
                      serializables (ints, floats, numpy arrays, slices).
                      Se transfiere al hijo via pickle (pickle maneja slice y
                      numpy de forma nativa, a diferencia de JSON).
        solver_dir:   carpeta destino del codigo C generado.

    Lanza RuntimeError si el hijo no logra escribir cpg_solve.c.
    """
    # Raiz a importar en el hijo: el padre de la carpeta del paquete.
    # __file__ = .../sympy_cvx_converter/pipeline.py
    pkg_dir = os.path.dirname(os.path.abspath(__file__))      # .../sympy_cvx_converter
    pkg_root = os.path.dirname(pkg_dir)                        # carpeta que contiene el paquete
    worker = os.path.join(pkg_dir, "_cpg_worker.py")

    # Empaquetar inputs en un pickle temporal
    payload = {"pkg_root": pkg_root, "build_kwargs": build_kwargs, "solver_dir": solver_dir}
    with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False, mode="wb") as fp:
        pickle.dump(payload, fp)
        payload_path = fp.name

    try:
        proc = subprocess.run(
            [sys.executable, worker, payload_path],
            capture_output=True, text=True,
        )
        # Criterio de exito: el archivo C existe (el wrapper Python es opcional)
        cpg_solve_c = os.path.join(solver_dir, "c", "src", "cpg_solve.c")
        if not os.path.isfile(cpg_solve_c):
            raise RuntimeError(
                "CVXPYgen (subproceso) no genero el solver C.\n"
                f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
            )
    finally:
        try:
            os.unlink(payload_path)
        except OSError:
            pass


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
                     generate_c=False, c_output_dir="c_output",
                     solver_source_dir=None):
    """
    Resuelve un problema de trayectoria optima via SCVx.
    El usuario define todo en SymPy y numpy — la libreria hace el resto.

    Flujo (version Enhanced):
        1. Resuelve en Python (solve_scvx).
        2. SOLO si el solve corrio sin lanzar excepcion y generate_c=True,
           genera el codigo C. La generacion del solver con CVXPYgen corre en
           un subproceso limpio para evitar el MemoryError de canonicalizacion.

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
        generate_c:            si True, genera C tras un solve Python exitoso
        c_output_dir:          carpeta destino del codigo C
        solver_source_dir:     si se da, copia un solver preexistente en vez de
                               generarlo (atajo opcional; ya NO es necesario)

    Retorna:
        dict con x_opt, u_opt (unscaled), converged, iterations, history
    """

    nx = len(states)
    nu = len(controls)
    ng = len(nonconvex_constraints) if nonconvex_constraints else 0

    # ── [A] Code generation (SymPy → funciones numpy para el solve Python) ──
    codegen_file = tempfile.NamedTemporaryFile(
        suffix='.py', delete=False, prefix='scvx_codegen')
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

        # Argumentos serializables de build_problem: se reutilizan tal cual en el
        # subproceso de generacion C (donde se reconstruye el cp.Problem en limpio).
        build_kwargs = dict(
            nx=nx, nu=nu, T=T, ng=ng,
            tau_val=tau_val, lamb_val=lamb,
            S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
            state_bounds=state_bounds,
            control_bounds=control_bounds,
        )
        prob_dict = build_problem(**build_kwargs)

        # ── [E] Resolver en Python PRIMERO ──
        # Si solve_scvx lanza una excepcion, se propaga aqui y NO se llega al
        # bloque de generacion C: no se genera codigo de un solve fallido.
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

        # ── [F] Generar C SOLO si el solve Python corrio sin problemas ──
        if generate_c:
            os.makedirs(c_output_dir, exist_ok=True)

            # 1. Solver C
            solver_dir = os.path.join(c_output_dir, 'solver')
            cpg_solve_c = os.path.join(solver_dir, 'c', 'src', 'cpg_solve.c')

            if solver_source_dir is not None:
                # Atajo opcional: copiar un solver preexistente
                if os.path.isdir(solver_dir):
                    shutil.rmtree(solver_dir)
                shutil.copytree(solver_source_dir, solver_dir)
                print(f"[CVXPYgen] Solver copiado desde {solver_source_dir}")
            elif os.path.isfile(cpg_solve_c):
                # Ya existe en c_output_dir — reutilizar
                print(f"[CVXPYgen] Solver ya existe en {solver_dir}/c/ — omitiendo regeneracion")
            else:
                # Generar desde cero en subproceso limpio (evita MemoryError)
                os.makedirs(solver_dir, exist_ok=True)
                print("[CVXPYgen] Python OK -> generando solver C en subproceso (RAM limpia)...")
                _generate_solver_subprocess(build_kwargs, solver_dir)
                print(f"[CVXPYgen] Solver C generado en {solver_dir}/c/")

            # 2. Dinamica en C (dynamics.c + dynamics.h)
            generate_c_functions(
                states=states,
                controls=controls,
                dynamics=dynamics,
                nonconvex_constraints=nonconvex_constraints,
                dynamic_parameters=dynamic_parameters_sym,
                output_dir=c_output_dir,
            )

            # 3. Loop SCVx C (scvx_main.c + cpg_compat.h)
            np_ = len(dynamic_parameters_sym) if dynamic_parameters_sym else 0
            fill_scvx_template(nx, nu, ng, np_, T, c_output_dir)

        return result

    finally:
        # Limpiar archivo temporal de codegen
        try:
            os.unlink(codegen_path)
        except OSError:
            pass
