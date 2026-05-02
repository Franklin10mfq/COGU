"""
Test rapido: puede CVXPYgen generar codigo C para T=30, nx=13, nu=6, ng=3 (Astrobee)?

Prueba solo build_problem() + CVXPYgen, sin correr el loop SCVx (que tarda 40s).
Si esto pasa, astrobee_simple.py con T=30 y generate_c=True funcionara.
"""
import sys, os, time, tracemalloc
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'COGU_library'))
from sympy_cvx_converter.problem import build_problem
from cvxpygen import cpg

# --- Dimensiones del problema Astrobee ---
nx, nu, ng = 13, 6, 3
T  = 30
tf = 200.0
tau_val = tf / T

acc_max   = 20 / 7.2 * 0.001
torq_max  = 100 * 0.000001
vel_max   = 5.0
omega_max = 5.0 * np.pi / 180

S_x = np.diag([60,60,60, 2*vel_max,2*vel_max,2*vel_max,
               2,2,2,2, 2*0.1*omega_max,2*0.1*omega_max,2*0.1*omega_max])
c_x = np.array([[-30],[-30],[-30],
                [-vel_max],[-vel_max],[-vel_max],
                [-1],[-1],[-1],[-1],
                [-(0.1*omega_max)],[-(0.1*omega_max)],[-(0.1*omega_max)]])
S_u = np.diag([2*10*acc_max]*3 + [2*0.1*torq_max]*3)
c_u = np.array([[-(10*acc_max)]]*3 + [[-(0.1*torq_max)]]*3)

state_bounds = [
    (slice(3,6),  'norm2', vel_max),
    (slice(10,13),'norm2', omega_max),
]
control_bounds = [
    (slice(0,3), 'norm1', acc_max),
    (slice(3,6), 'norm1', torq_max),
]

print(f"Construyendo problema CVXPY: nx={nx}, nu={nu}, ng={ng}, T={T} ...")
t0 = time.time()
prob_dict = build_problem(
    nx=nx, nu=nu, T=T, ng=ng,
    tau_val=tau_val, lamb_val=1000.0,
    S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
    state_bounds=state_bounds,
    control_bounds=control_bounds,
)
print(f"  -> build_problem OK ({time.time()-t0:.1f}s)")

output_dir = "solver_t30/solver"
os.makedirs(output_dir, exist_ok=True)

print(f"\nCorriendo CVXPYgen (T={T}) ...")
tracemalloc.start()
t1 = time.time()
try:
    cpg.generate_code(prob_dict['problem'], code_dir=output_dir, solver='ECOS')
    elapsed = time.time() - t1
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"\n=== EXITO ===")
    print(f"Tiempo CVXPYgen: {elapsed:.1f}s")
    print(f"Memoria pico:    {peak/1e6:.0f} MB")
    print(f"Archivos en:     {output_dir}/c/")
    # Verificar que los archivos C existen
    c_src = os.path.join(output_dir, 'c', 'src')
    if os.path.isdir(c_src):
        for f in os.listdir(c_src):
            fpath = os.path.join(c_src, f)
            print(f"  {f}  ({os.path.getsize(fpath)/1024:.0f} KB)")
except MemoryError as e:
    elapsed = time.time() - t1
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"\n=== FALLO: MemoryError ===")
    print(f"Tiempo antes de fallar: {elapsed:.1f}s")
    print(f"Memoria pico:           {peak/1e6:.0f} MB")
    print(f"Necesitas CVXPY 1.6.4 con entorno virtual.")
    print(f"Error: {e}")
except Exception as e:
    elapsed = time.time() - t1
    tracemalloc.stop()
    # Puede ser que los archivos C se generaron pero fallo la compilacion Python
    c_src = os.path.join(output_dir, 'c', 'src')
    if os.path.isdir(c_src) and any(f.endswith('.c') for f in os.listdir(c_src)):
        print(f"\n=== C GENERADO (wrapper Python fallo — normal sin CMake) ===")
        for f in os.listdir(c_src):
            fpath = os.path.join(c_src, f)
            print(f"  {f}  ({os.path.getsize(fpath)/1024:.0f} KB)")
    else:
        print(f"\n=== FALLO: {type(e).__name__} ===")
        traceback.print_exc()
