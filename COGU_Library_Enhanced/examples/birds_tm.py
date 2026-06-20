"""
BIRDS-TM (Reaction Wheels + Thrusters) via solve_trajectory() — Fase 4.

Porta el notebook BIRDS_TM_RW_Th_cvxpy_cvxpygen_l2_q.ipynb a la libreria Enhanced.
Control de actitud puro: nx=10 (quaternion + momento cuerpo + momento ruedas),
nu=7 (3 torques de rueda + 4 magnitudes de thruster), ng=0, np=30.

Es Python-only (generate_c=False); el C llega en Fases 5-7.

Usa dos extensiones DSL de Fase 4:
  - cost_terms con 'offset'  -> seguimiento  sqrt_tau*(x - x_goal)
  - state_bounds con 'S'/'c' -> transformacion inv(J)@S_x  (momento -> velocidad angular)

Referencia del notebook (loop ECOS, mismo problema con escalares como float):
  Iter 1: L=2318.04  J=6425.60
  Iter 5: L=1557.53  J=1643.41
  Iter 8: L=1557.13  J=1614.23
  Converge: L->~1557.1  J->~1614   quat_final=[0,0,1,0]  hb,hw final=0

Run desde c:/Users/luiss/COGU/:
    python COGU_Library_Enhanced/examples/birds_tm.py
"""

import sys, os, time
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sympy_cvx_converter import solve_trajectory
from sympy_cvx_converter.utils import slerp, angular_vel

# =============================================================================
# 1. MODELO SYMPY (copiado EXACTO del Cell 9 del notebook)
# =============================================================================
tauWx, tauWy, tauWz = sp.symbols('tauWx tauWy tauWz')
b1th, b2th, b3th, b4th = sp.symbols('b1th b2th b3th b4th')

q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
hbx, hby, hbz = sp.symbols('hbx hby hbz')
hw1, hw2, hw3 = sp.symbols('hw1 hw2 hw3')

Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')

rth1x, rth1y, rth1z = sp.symbols('rth1x rth1y rth1z')
rth2x, rth2y, rth2z = sp.symbols('rth2x rth2y rth2z')
rth3x, rth3y, rth3z = sp.symbols('rth3x rth3y rth3z')
rth4x, rth4y, rth4z = sp.symbols('rth4x rth4y rth4z')
rth1 = sp.Matrix([rth1x, rth1y, rth1z])
rth2 = sp.Matrix([rth2x, rth2y, rth2z])
rth3 = sp.Matrix([rth3x, rth3y, rth3z])
rth4 = sp.Matrix([rth4x, rth4y, rth4z])

fth1x, fth1y, fth1z = sp.symbols('fth1x fth1y fth1z')
fth2x, fth2y, fth2z = sp.symbols('fth2x fth2y fth2z')
fth3x, fth3y, fth3z = sp.symbols('fth3x fth3y fth3z')
fth4x, fth4y, fth4z = sp.symbols('fth4x fth4y fth4z')
fth1 = sp.Matrix([fth1x, fth1y, fth1z])
fth2 = sp.Matrix([fth2x, fth2y, fth2z])
fth3 = sp.Matrix([fth3x, fth3y, fth3z])
fth4 = sp.Matrix([fth4x, fth4y, fth4z])

J = sp.Matrix([[Ixx, Ixy, Ixz],
               [Ixy, Iyy, Iyz],
               [Ixz, Iyz, Izz]])

hb_sympy   = sp.Matrix([hbx, hby, hbz])
hw_sympy   = sp.Matrix([hw1, hw2, hw3])
tauW_sympy = sp.Matrix([tauWx, tauWy, tauWz])

f_quaternion = 0.5 * sp.Matrix([[-q1, -q2, -q3],
                                [ q0, -q3,  q2],
                                [ q3,  q0, -q1],
                                [-q2,  q1,  q0]]) * (J.inv() * hb_sympy)

torque_th = (b1th*(rth1.cross(fth1)) + b2th*(rth2.cross(fth2))
             + b3th*(rth3.cross(fth3)) + b4th*(rth4.cross(fth4)))

f_ridig_body  = (-(J.inv()*hb_sympy).cross(hb_sympy + hw_sympy) - tauW_sympy + torque_th)
f_ridig_body  = sp.simplify(f_ridig_body)
f_rigid_wheel = tauW_sympy

dynamics = [
    f_quaternion[0], f_quaternion[1], f_quaternion[2], f_quaternion[3],
    f_ridig_body[0], f_ridig_body[1], f_ridig_body[2],
    f_rigid_wheel[0], f_rigid_wheel[1], f_rigid_wheel[2],
]

states   = [q0, q1, q2, q3, hbx, hby, hbz, hw1, hw2, hw3]
controls = [tauWx, tauWy, tauWz, b1th, b2th, b3th, b4th]
dyn_params_sym = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz,
                  rth1x, rth1y, rth1z, rth2x, rth2y, rth2z,
                  rth3x, rth3y, rth3z, rth4x, rth4y, rth4z,
                  fth1x, fth1y, fth1z, fth2x, fth2y, fth2z,
                  fth3x, fth3y, fth3z, fth4x, fth4y, fth4z]

# =============================================================================
# 2. PARAMETROS NUMERICOS (copiado del Cell 17 y Cell 19)
# =============================================================================
T = 50
tf = 100.0
tau_val = tf / T

J1, J2, J3 = 0.0333, 0.0333, 0.0067
Jxy = Jxz = Jyz = 0.0001
IW1 = IW2 = IW3 = 4.774648292756861e-06

Jsat = np.array([[J1, Jxy, Jxz], [Jxy, J2, Jyz], [Jxz, Jyz, J3]])
IWmatrix = np.array([[IW1, 0, 0], [0, IW2, 0], [0, 0, IW3]])

rth1_d = np.array([[-26.43], [-26.61], [-147.785]]) * 0.001
rth2_d = np.array([[ 26.43], [-26.61], [-147.785]]) * 0.001
rth3_d = np.array([[ 26.43], [ 26.61], [-147.785]]) * 0.001
rth4_d = np.array([[-26.43], [ 26.61], [-147.785]]) * 0.001

fth1_d = np.array([[0], [ 0.643], [0.766]]) * 0.001
fth2_d = np.array([[0], [ 0.643], [0.766]]) * 0.001
fth3_d = np.array([[0], [-0.643], [0.766]]) * 0.001
fth4_d = np.array([[0], [-0.643], [0.766]]) * 0.001
fth1_d = fth1_d / np.linalg.norm(fth1_d)
fth2_d = fth2_d / np.linalg.norm(fth2_d)
fth3_d = fth3_d / np.linalg.norm(fth3_d)
fth4_d = fth4_d / np.linalg.norm(fth4_d)

dyn_params_val = np.array([
    J1, J2, J3, Jxy, Jxz, Jyz,
    rth1_d[0, 0], rth1_d[1, 0], rth1_d[2, 0],
    rth2_d[0, 0], rth2_d[1, 0], rth2_d[2, 0],
    rth3_d[0, 0], rth3_d[1, 0], rth3_d[2, 0],
    rth4_d[0, 0], rth4_d[1, 0], rth4_d[2, 0],
    fth1_d[0, 0], fth1_d[1, 0], fth1_d[2, 0],
    fth2_d[0, 0], fth2_d[1, 0], fth2_d[2, 0],
    fth3_d[0, 0], fth3_d[1, 0], fth3_d[2, 0],
    fth4_d[0, 0], fth4_d[1, 0], fth4_d[2, 0],
])

# Limites fisicos
torqW_max  = 0.001                 # N*m
omegab_max = 10.0 * np.pi / 180    # rad/s
speedW_max = 8000 * 2*np.pi / 60   # rad/s
Fth_max    = 0.632 * 0.001         # N

# =============================================================================
# 3. BOUNDARY CONDITIONS (Cell 17)
# =============================================================================
start_quat = np.array([[1], [0], [0], [0]])
end_quat   = np.array([[0], [0], [1], [0]])
start_hb = Jsat @ np.array([[0.1*np.pi/180], [0.2*np.pi/180], [-0.1*np.pi/180]])
start_hW = IWmatrix @ np.array([[500*2*np.pi/60], [-200*2*np.pi/60], [700*2*np.pi/60]])
end_hb   = Jsat @ np.zeros((3, 1))
end_hW   = IWmatrix @ np.zeros((3, 1))

start_pos = np.block([[start_quat], [start_hb], [start_hW]]).astype(float)
end_pos   = np.block([[end_quat],   [end_hb],   [end_hW]]).astype(float)

# =============================================================================
# 4. ESCALADO (Cell 17)
# =============================================================================
S_x = np.diag([2*1, 2*1, 2*1, 2*1,
               2*J1*omegab_max, 2*J2*omegab_max, 2*J3*omegab_max,
               2*IW1*speedW_max, 2*IW2*speedW_max, 2*IW3*speedW_max])
c_x = np.array([[-1], [-1], [-1], [-1],
                [-J1*omegab_max], [-J2*omegab_max], [-J3*omegab_max],
                [-IW1*speedW_max], [-IW2*speedW_max], [-IW3*speedW_max]])
S_u = np.diag([2*0.1*torqW_max]*3 + [2*10*Fth_max]*4)
c_u = np.array([[-(0.1*torqW_max)]]*3 + [[-(10*Fth_max)]]*4)

# =============================================================================
# 5. WARM START (Cell 21) — slerp quaternion + Jsat@omega para hb, ceros resto
# =============================================================================
quat = slerp(start_pos[0:4, 0], end_pos[0:4, 0], T+1).T   # (4, T+1)
omega = angular_vel(quat.T, tau_val).T                     # (3, T+1)
hb_ws = Jsat @ omega

warm_x = np.zeros((10, T+1))
warm_x[0:4, :] = quat
warm_x[4:7, :] = hb_ws
warm_x[7:10, :] = 0.0
warm_u = np.zeros((7, T))

# =============================================================================
# 6. COST TERMS (usa offset de Fase 4) — x_goal es el estado objetivo ESCALADO
# =============================================================================
inv_S_x = np.linalg.inv(S_x)
x_goal = (inv_S_x @ (end_pos - c_x))[0:10]   # end_pos escalado (igual que el notebook)

cost_terms = [
    {'kind': 'sumsq', 'var': 'u', 'slice': slice(0, 3), 'coeff': 'sqrt_tau'},
    {'kind': 'norm1', 'var': 'u', 'slice': slice(3, 7), 'coeff': 'tau', 'weight': 1.0},
    {'kind': 'sumsq', 'var': 'x', 'slice': slice(0, 10), 'coeff': 'sqrt_tau',
     'weight': 100, 'k_range': 'T+1', 'offset': x_goal},
]

# =============================================================================
# 7. CONSTRAINTS (usa S/c custom de Fase 4 para las de estado)
# =============================================================================
inv_Jb_S = np.linalg.inv(Jsat) @ S_x[4:7, 4:7]
inv_Jb_c = np.linalg.inv(Jsat) @ c_x[4:7, 0:1]
inv_JW_S = np.linalg.inv(IWmatrix) @ S_x[7:10, 7:10]
inv_JW_c = np.linalg.inv(IWmatrix) @ c_x[7:10, 0:1]

state_bounds = [
    {'kind': 'norm', 'slice': slice(4, 7),  'norm_type': 'norm2',
     'limit': omegab_max, 'S': inv_Jb_S, 'c': inv_Jb_c},
    {'kind': 'norm', 'slice': slice(7, 10), 'norm_type': 'norminf',
     'limit': speedW_max, 'S': inv_JW_S, 'c': inv_JW_c},
]
control_bounds = [
    (slice(0, 3), 'norminf', torqW_max),
    {'kind': 'box',  'slice': slice(3, 7), 'lower': 0.0, 'upper': Fth_max},
    {'kind': 'norm', 'slice': slice(3, 5), 'norm_type': 'norm1', 'limit': Fth_max},
    {'kind': 'norm', 'slice': slice(5, 7), 'norm_type': 'norm1', 'limit': Fth_max},
]

# =============================================================================
# 8. RESOLVER (Python-only)
# =============================================================================
print("=" * 60)
print("BIRDS-TM (RW + Th) -- Fase 4 Python-only")
print("=" * 60)
t0 = time.time()

result = solve_trajectory(
    states=states,
    controls=controls,
    dynamics=dynamics,
    start=start_pos,
    end=end_pos,
    T=T, tf=tf,
    dynamic_parameters_sym=dyn_params_sym,
    dynamic_parameters_val=dyn_params_val,
    S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
    warm_start_x=warm_x,
    warm_start_u=warm_u,
    cost_terms=cost_terms,
    state_bounds=state_bounds,
    control_bounds=control_bounds,
    size_N=20,
    lamb=5000.0,
    max_iter=25,
    generate_c=False,
)

total = time.time() - t0

# =============================================================================
# 9. RESULTADOS — comparar con referencia del notebook
# =============================================================================
x = result['x_opt_unscaled']
h = result['history']

print("\n" + "=" * 60)
print("RESULTADOS BIRDS-TM")
print("=" * 60)
print(f"Tiempo total:  {total:.1f}s")
print(f"Iteraciones:   {result['iterations']}")
print(f"Convergio:     {result['converged']}")
print(f"Quat final:    {x[0:4, -1].round(5)}   (objetivo [0 0 1 0])")
print(f"Quat start:    {x[0:4, 0].round(5)}   (objetivo [1 0 0 0])")
print(f"hb  final:     {x[4:7, -1].round(6)}   (objetivo [0 0 0])")
print(f"hw  final:     {x[7:10, -1].round(6)}  (objetivo [0 0 0])")

# Verificacion de constraints fisicas
omega_traj = np.linalg.inv(Jsat) @ x[4:7, :]
wheel_traj = np.linalg.inv(IWmatrix) @ x[7:10, :]
print(f"\nMax ||omega||:  {np.max(np.linalg.norm(omega_traj, axis=0)):.5f}  (limit {omegab_max:.5f})")
print(f"Max |wheel|:    {np.max(np.abs(wheel_traj)):.3f}  (limit {speedW_max:.3f})")

print("\n--- Comparacion con notebook (referencia ECOS) ---")
print(f"Iter 1: L={h[0]['L_SCVx']:.4f}  J={h[0]['J_SCVx']:.4f}   (ref: L=2318.04  J=6425.60)")
print(f"Iter 5: L={h[4]['L_SCVx']:.4f}  J={h[4]['J_SCVx']:.4f}   (ref: L=1557.53  J=1643.41)")
print(f"Iter 8: L={h[7]['L_SCVx']:.4f}  J={h[7]['J_SCVx']:.4f}   (ref: L=1557.13  J=1614.23)")
