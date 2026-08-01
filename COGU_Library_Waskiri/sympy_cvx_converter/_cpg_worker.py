"""
Worker de subproceso para generar el solver C con CVXPYgen en RAM limpia.

Lo invoca pipeline.py (_generate_solver_subprocess) como proceso hijo:
    python _cpg_worker.py <payload.pkl>

El payload (pickle) contiene:
    pkg_root     : ruta a agregar a sys.path para importar el paquete
    build_kwargs : argumentos de build_problem() (incluye slices y numpy arrays)
    solver_dir   : carpeta destino del codigo C generado

Por que un proceso aparte:
    CVXPYgen necesita un bloque de RAM contiguo de varios cientos de MB durante
    la canonicalizacion DPP. En el proceso principal, SymPy (generate_functions)
    y el solve SCVx ya consumieron y fragmentaron la memoria -> MemoryError.
    Un proceso hijo arranca limpio y dispone del bloque contiguo completo.
"""
import sys
import os
import pickle


def main():
    # argv[1] = ruta al pickle con los datos que mando el proceso padre
    payload_path = sys.argv[1]
    with open(payload_path, "rb") as fp:
        payload = pickle.load(fp)

    # cmake (MSYS2) en PATH: CVXPYgen intenta compilar un wrapper Python opcional.
    # Si cmake/VS no estan, ese wrapper falla pero el codigo C (lo que importa)
    # ya queda escrito en disco. Por eso lo agregamos y, mas abajo, toleramos el fallo.
    os.environ["PATH"] = r"C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + os.environ.get("PATH", "")

    # Importar la libreria desde la raiz que paso el padre (no asumimos cwd)
    sys.path.insert(0, payload["pkg_root"])
    from sympy_cvx_converter.problem import build_problem
    from cvxpygen import cpg

    # Reconstruir el problema CVXPY en limpio (cp.Problem no es serializable,
    # pero sus inputs si lo son -> lo regeneramos identico aqui)
    prob_dict = build_problem(**payload["build_kwargs"])
    solver_dir = payload["solver_dir"]

    try:
        cpg.generate_code(prob_dict["problem"], code_dir=solver_dir, solver="ECOS")
    except Exception as e:
        # Distinguir "fallo del wrapper Python" (tolerable) de "no se genero C" (fatal)
        c_src = os.path.join(solver_dir, "c", "src", "cpg_solve.c")
        if os.path.isfile(c_src):
            print(f"[worker] Codigo C generado; wrapper Python omitido ({e})")
        else:
            # No hay C -> fallo real (p.ej. MemoryError antes de escribir). Propagar.
            raise


if __name__ == "__main__":
    main()
