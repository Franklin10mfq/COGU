"""
Verificacion final del pipeline COGU.

Genera C en astrobee_c_output/, compila y muestra resultados
para comparar con SCVx_official_template_cvxpy_cvxpygen_ZOH.ipynb.

Run desde c:/Users/luiss/COGU/:
    python COGU_Dev_Library/verify_final.py
"""

import sys, os, time
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from sympy_cvx_converter import solve_trajectory
from sympy_cvx_converter.utils import slerp, angular_vel

# =============================================================================
# 1. MODELO SYMPY (del Cell 21 del notebook de referencia)
# =============================================================================
ax, ay, az = sp.symbols('ax ay az')
taux, tauy, tauz = sp.symbols('taux tauy tauz')
rx, ry, rz = sp.symbols('rx ry rz')
vx, vy, vz = sp.symbols('vx vy vz')
q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
wx, wy, wz = sp.symbols('wx wy wz')
Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')
d_obs1, d_obs2, d_obs3 = sp.symbols('d_obs1 d_obs2 d_obs3')
c_obs1_x, c_obs1_y, c_obs1_z = sp.symbols('c_obs1_x c_obs1_y c_obs1_z')
c_obs2_x, c_obs2_y, c_obs2_z = sp.symbols('c_obs2_x c_obs2_y c_obs2_z')
c_obs3_x, c_obs3_y, c_obs3_z = sp.symbols('c_obs3_x c_obs3_y c_obs3_z')

I = sp.Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
w_sym = sp.Matrix([wx, wy, wz])
tau_sym = sp.Matrix([taux, tauy, tauz])
f_rb = sp.simplify(I.inv() * (-w_sym.cross(I * w_sym) + tau_sym))

dynamics = [
    vx, vy, vz, ax, ay, az,
    0.5*(-wx*q1 - wy*q2 - wz*q3),
    0.5*(+wx*q0 + wz*q2 - wy*q3),
    0.5*(+wy*q0 - wz*q1 + wx*q3),
    0.5*(+wz*q0 + wy*q1 - wx*q2),
    f_rb[0], f_rb[1], f_rb[2],
]

diff1 = sp.Matrix([rx,ry,rz]) - sp.Matrix([c_obs1_x, c_obs1_y, c_obs1_z])
diff2 = sp.Matrix([rx,ry,rz]) - sp.Matrix([c_obs2_x, c_obs2_y, c_obs2_z])
diff3 = sp.Matrix([rx,ry,rz]) - sp.Matrix([c_obs3_x, c_obs3_y, c_obs3_z])
nonconvex = [
    -sp.sqrt(diff1.dot(diff1)) + d_obs1,
    -sp.sqrt(diff2.dot(diff2)) + d_obs2,
    -sp.sqrt(diff3.dot(diff3)) + d_obs3,
]

states   = [rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz]
controls = [ax, ay, az, taux, tauy, tauz]
dyn_params_sym = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz,
                  c_obs1_x, c_obs1_y, c_obs1_z,
                  c_obs2_x, c_obs2_y, c_obs2_z,
                  c_obs3_x, c_obs3_y, c_obs3_z,
                  d_obs1, d_obs2, d_obs3]

# =============================================================================
# 2. PARAMETROS NUMERICOS (del Cell 29 del notebook)
# =============================================================================
T = 30
tf = 200.0
acc_max  = 20 / 7.2 * 0.001
torq_max = 100 * 0.000001
vel_max  = 5.0
omega_max = 5.0 * np.pi / 180

start_pos = np.array([[0.0],[-1.0],[1.0],
                       [0.001],[-0.002],[0.0],
                       [0],[1],[0],[0],
                       [0.0],[0.0],[0.0]], dtype=float)
end_pos = np.array([[5.0],[2.0],[1.6],
                     [0.0],[0.0],[0.0],
                     [0],[0],[1],[0],
                     [0.0],[0.0],[0.0]], dtype=float)

dyn_params_val = np.array([
    0.1083, 0.1083, 0.1083, 0.01, 0.02, 0.01,
    1.1, -0.5, 1.0,   2.6, 0.5, 1.1,   4.0, 1.6, 1.2,
    0.8, 0.8, 0.8,
])

S_x = np.diag([60,60,60, 2*vel_max,2*vel_max,2*vel_max,
               2,2,2,2, 2*0.1*omega_max,2*0.1*omega_max,2*0.1*omega_max])
c_x = np.array([[-30],[-30],[-30],
                [-vel_max],[-vel_max],[-vel_max],
                [-1],[-1],[-1],[-1],
                [-(0.1*omega_max)],[-(0.1*omega_max)],[-(0.1*omega_max)]])
S_u = np.diag([2*10*acc_max]*3 + [2*0.1*torq_max]*3)
c_u = np.array([[-(10*acc_max)]]*3 + [[-(0.1*torq_max)]]*3)

# =============================================================================
# 3. WARM START (del Cell 33)
# =============================================================================
tau_val = tf / T
pos = np.column_stack([np.linspace(start_pos[i,0], end_pos[i,0], T+1) for i in range(3)]).T
vel = np.zeros_like(pos)
for t in range(1, T): vel[:, t] = (pos[:, t] - pos[:, t-1]) / tau_val
vel[:, 0:1], vel[:, T:T+1] = start_pos[3:6], end_pos[3:6]
quat = slerp(start_pos[6:10,0], end_pos[6:10,0], T+1).T
angvel = angular_vel(quat.T, tau_val).T

warm_x = np.zeros((13, T+1))
warm_x[0:3,:] = pos; warm_x[3:6,:] = vel; warm_x[6:10,:] = quat; warm_x[10:13,:] = angvel
warm_u = np.zeros((6, T))

# =============================================================================
# 4. RESOLVER con generate_c=True en astrobee_c_output/
# =============================================================================
C_OUTPUT = "COGU_Dev_Library/astrobee_c_output"
# Solver preexistente generado en sesion anterior (evita MemoryError de CVXPYgen en T=30)
SOLVER_SOURCE = "COGU_Dev_Library/references/solver_t30/solver"

print("=" * 60)
print("COGU pipeline -- Verificacion final")
print(f"Generando C en: {C_OUTPUT}/")
print(f"Solver copiado desde: {SOLVER_SOURCE}/")
print("=" * 60)
t0 = time.time()

result = solve_trajectory(
    states=states,
    controls=controls,
    dynamics=dynamics,
    start=start_pos,
    end=end_pos,
    T=T, tf=tf,
    nonconvex_constraints=nonconvex,
    dynamic_parameters_sym=dyn_params_sym,
    dynamic_parameters_val=dyn_params_val,
    S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
    warm_start_x=warm_x,
    warm_start_u=warm_u,
    state_bounds=[
        (slice(3,6), 'norm2', vel_max),
        (slice(10,13), 'norm2', omega_max),
    ],
    control_bounds=[
        (slice(0,3), 'norm1', acc_max),
        (slice(3,6), 'norm1', torq_max),
    ],
    size_N=20,
    max_iter=25,
    generate_c=True,
    c_output_dir=C_OUTPUT,
    solver_source_dir=SOLVER_SOURCE,
)

total = time.time() - t0

# =============================================================================
# 5. RESULTADOS PYTHON — comparar con notebook Cell 43
# =============================================================================
x = result['x_opt_unscaled']
h = result['history']

print("\n" + "=" * 60)
print("RESULTADOS PYTHON")
print("=" * 60)
print(f"Tiempo total:  {total:.1f}s")
print(f"Iteraciones:   {result['iterations']}")
print(f"Convergio:     {result['converged']}")
print(f"Start pos:     {x[0:3, 0].round(4)}")
print(f"End pos:       {x[0:3,-1].round(4)}")
print(f"Err start:     {np.linalg.norm(x[0:3,0:1]-start_pos[0:3]):.2e}")
print(f"Err end:       {np.linalg.norm(x[0:3,-1:]-end_pos[0:3]):.2e}")
print(f"Max vel:       {np.max(np.linalg.norm(x[3:6,:],axis=0)):.4f}  (limit {vel_max})")
print(f"Max omega:     {np.max(np.linalg.norm(x[10:13,:],axis=0)):.6f}  (limit {omega_max:.6f})")

print("\n--- Comparacion con notebook Cell 43 ---")
print(f"Iter 1: L={h[0]['L_SCVx']:.4f}  J={h[0]['J_SCVx']:.4f}   (ref: L=3267.2262  J=8899.0704)")
print(f"Iter 5: L={h[4]['L_SCVx']:.4f}  J={h[4]['J_SCVx']:.4f}   (ref: L=723.5744   J=733.6812)")
print(f"Iter 8: L={h[7]['L_SCVx']:.4f}  J={h[7]['J_SCVx']:.4f}   (ref: L=723.5393   J=729.4670)")

print("\n--- Archivos C generados ---")
for f in ['scvx_main.c', 'dynamics.c', 'dynamics.h', 'cpg_compat.h']:
    path = os.path.join(C_OUTPUT, f)
    size = os.path.getsize(path) if os.path.exists(path) else 0
    status = "OK" if size > 0 else "FALTA"
    print(f"  {f:<20} {size:>8} bytes  [{status}]")
solver_c = os.path.join(C_OUTPUT, 'solver', 'c', 'src', 'cpg_solve.c')
size = os.path.getsize(solver_c) if os.path.exists(solver_c) else 0
print(f"  {'solver/c/src/cpg_solve.c':<20} {size:>8} bytes  [{'OK' if size>0 else 'FALTA'}]")

print(f"\nC output en: {os.path.abspath(C_OUTPUT)}/")
print("Listo para compilar con build_astrobee_c_output.ps1")
