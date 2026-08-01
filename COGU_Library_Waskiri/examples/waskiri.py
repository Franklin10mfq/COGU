"""
waskiri.py -- Waskiri-Sat (1U CubeSat, solo magnetorquers).

Modelo de referencia: COGU_algorithms/Waskiri-Sat_trials_working/
SCVx_Waskiri_Sat_cvxpy_cvxpygen_l2_q_dm_ZOH_opt.ipynb (confirmado por Franklin
como el modelo correcto: discretizacion ZOH, campo magnetico terrestre via
dipolo inclinado analitico -- NO usa WMM/propagador orbital).

Estado (7): cuaternion q0..q3 + momento angular en cuerpo hbx,hby,hbz.
Control (3): dipolo magnetico mx,my,mz.
Parametros dinamicos (11): inv(I) (6 componentes), n (movimiento medio orbital),
Omega (RAAN), inc (inclinacion), B0 (escala del campo magnetico), t0 (tiempo inicial).

Costo: energia del dipolo + tasa de cambio del dipolo (evita conmutaciones
bruscas del magnetorquer) + seguimiento de actitud objetivo durante toda la
maniobra. El termino de tasa de cambio requirio extender el DSL de costo con
un nuevo 'kind': 'diffsumsq' (ver sympy_cvx_converter/cost_dsl.py).

Restricciones: velocidad angular (norma de h_B transformado por inv(J)) y
dipolo maximo (norma infinito).
"""

import sys, os
import numpy as np
import sympy as sp
from sympy.physics.vector import ReferenceFrame

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sympy_cvx_converter import solve_trajectory
from sympy_cvx_converter.utils import slerp, angular_vel

# =============================================================================
# 1. MODELO SIMBOLICO (identico a la notebook de referencia)
# =============================================================================
mx, my, mz = sp.symbols('mx my mz')
t = sp.symbols('t')
q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
hbx, hby, hbz = sp.symbols('hbx hby hbz')
invIxx, invIyy, invIzz = sp.symbols('invIxx invIyy invIzz')
invIxy, invIxz, invIyz = sp.symbols('invIxy invIxz invIyz')
n_cubesat, Omega, inc, B0, t0 = sp.symbols('n_cubesat Omega inc B0 t0')

# Cinematica de cuerpo rigido: DCM cuerpo<-inercial (BN)
BN = sp.Matrix([
    [q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 + q0*q3), 2*(-q0*q2 + q1*q3)],
    [2*(q1*q2 - q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q0*q1 + q2*q3)],
    [2*(q0*q2 + q1*q3), 2*(-q0*q1 + q2*q3), q0**2 - q1**2 - q2**2 + q3**2],
])

invI = sp.Matrix([
    [invIxx, invIxy, invIxz],
    [invIxy, invIyy, invIyz],
    [invIxz, invIyz, invIzz],
])

# Campo magnetico terrestre: dipolo inclinado analitico en el frame orbital (O)
O_B = sp.Matrix([
    -2*sp.sin(inc)*sp.sin(n_cubesat*(t + t0)),
    sp.sin(inc)*sp.cos(n_cubesat*(t + t0)),
    sp.cos(inc),
]) * B0

Nframe = ReferenceFrame('Nframe')
Oframe = Nframe.orientnew('Oframe', 'Body', [Omega, inc, n_cubesat*(t + t0)], 'ZXZ')
ON = Oframe.dcm(Nframe)
BO = BN * ON.T

hb = sp.Matrix([hbx, hby, hbz])
m_sym = sp.Matrix([mx, my, mz])

f_rigid = -(invI * hb).cross(hb) + m_sym.cross(BO * O_B)
f_quat = 0.5 * sp.Matrix([
    [-q1, -q2, -q3],
    [ q0, -q3,  q2],
    [ q3,  q0, -q1],
    [-q2,  q1,  q0],
]) * (invI * hb)

dynamics = [f_quat[0], f_quat[1], f_quat[2], f_quat[3], f_rigid[0], f_rigid[1], f_rigid[2]]
states = [q0, q1, q2, q3, hbx, hby, hbz]
controls = [mx, my, mz]
dynamic_parameters_sym = [invIxx, invIyy, invIzz, invIxy, invIxz, invIyz,
                          n_cubesat, Omega, inc, B0, t0]

# =============================================================================
# 2. PARAMETROS NUMERICOS (identicos a la notebook de referencia, salvo T)
# =============================================================================
# El notebook usa T=200 (alta fidelidad, solo validacion en tierra). CVXPYgen no
# puede generar el solver C con T=200 en esta maquina (canonicalizacion de ECOS
# se queda sin memoria por la cantidad de parametros por paso -- confirmado que
# es un limite de escala, no un error del modelo: con T=50 genera sin problema).
# T=50 (misma escala que BIRDS-TM) se usa para el C embebido; puede subirse si
# se corre en una maquina con mas RAM o se resuelve el limite de CVXPYgen.
T = 50
tf = 200.0
tau = tf / T
size_N = 5

start_quat = np.array([[1], [0], [0], [0]])
end_quat = np.array([[0], [0], [1], [0]])
start_w_xyz = np.array([[0.1*np.pi/180], [0.1*np.pi/180], [-0.1*np.pi/180]])
end_w_xyz = np.array([[0.0], [0.0], [0.0]])

Jsat = np.array([[0.002, 0.00002, 0.00002],
                 [0.00002, 0.002, 0.00002],
                 [0.00002, 0.00002, 0.002]])
invJsat = np.linalg.inv(Jsat)

start_hb = Jsat @ start_w_xyz
end_hb = Jsat @ end_w_xyz
start = np.block([[start_quat], [start_hb]])
end = np.block([[end_quat], [end_hb]])

m_max = 0.2                    # A m^2
omega_max = 5.0 * np.pi / 180  # rad/s

# Orbita circular SSO -- parametros del dipolo inclinado analitico
mu_earth = 3.986004418e14
R_E = 6378e3
h_cubesat = 500e3
mu0_over_4pi = 1e-7
m_E = 7.94e22
a_cubesat = R_E + h_cubesat

n_cubesat_val = np.sqrt(mu_earth / a_cubesat**3)
Omega_val = 45 * np.pi / 180
inc_val = 97.6 * np.pi / 180
B0_val = mu0_over_4pi * m_E / a_cubesat**3
t0_val = 0.0

dynamic_parameters_val = np.array([
    invJsat[0, 0], invJsat[1, 1], invJsat[2, 2],
    invJsat[0, 1], invJsat[0, 2], invJsat[1, 2],
    n_cubesat_val, Omega_val, inc_val, B0_val, t0_val,
])

alpha_m = 1000.0   # peso de la tasa de cambio del dipolo
alpha_q = 1000.0   # peso del seguimiento de actitud objetivo

# Escalado
S_u = np.diag([2*(0.5*m_max)]*3)
c_u = np.array([[-(0.5*m_max)]]*3)

S_x = np.diag([2, 2, 2, 2,
               2*Jsat[0, 0]*omega_max, 2*Jsat[1, 1]*omega_max, 2*Jsat[2, 2]*omega_max]).astype(float)
c_x = np.array([[-1], [-1], [-1], [-1],
                [-Jsat[0, 0]*omega_max], [-Jsat[1, 1]*omega_max], [-Jsat[2, 2]*omega_max]])

# Objetivo de seguimiento, ya escalado (offset del cost_term de tracking)
inv_S_x = np.linalg.inv(S_x)
xq_goal_scaled = inv_S_x @ (end - c_x)

# =============================================================================
# 3. WARM START (slerp + velocidad angular por diferencias finitas)
# =============================================================================
ohx_quat_rows = slerp(start[0:4, 0], end[0:4, 0], T + 1)      # (T+1, 4)
ohx_omega_rows = angular_vel(ohx_quat_rows, tau)               # (T+1, 3)

warm_start_x = np.zeros((7, T + 1))
warm_start_x[0:4, :] = ohx_quat_rows.T
warm_start_x[4:7, :] = Jsat @ ohx_omega_rows.T
warm_start_u = np.zeros((3, T))

# =============================================================================
# 4. COSTO -- DSL (energia + tasa de cambio + seguimiento de actitud)
# =============================================================================
cost_terms = [
    {'kind': 'sumsq', 'var': 'u', 'coeff': 'sqrt_tau'},
    {'kind': 'diffsumsq', 'var': 'u', 'coeff': 'sqrt_tau', 'weight': alpha_m},
    {'kind': 'sumsq', 'var': 'x', 'coeff': 'sqrt_tau', 'weight': alpha_q,
     'offset': xq_goal_scaled, 'k_range': 'T+1'},
]

# =============================================================================
# 5. RESTRICCIONES -- DSL
# =============================================================================
state_bounds = [
    {'kind': 'norm', 'slice': slice(4, 7), 'norm_type': 'norm2', 'limit': omega_max,
     'S': invJsat @ S_x[4:7, 4:7], 'c': invJsat @ c_x[4:7, 0:1]},
]
control_bounds = [
    (slice(0, 3), 'norminf', m_max),
]

# =============================================================================
# 6. RESOLVER
# =============================================================================
if __name__ == '__main__':
    result = solve_trajectory(
        states, controls, dynamics, start, end, T, tf,
        dynamic_parameters_sym=dynamic_parameters_sym,
        dynamic_parameters_val=dynamic_parameters_val,
        S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
        warm_start_x=warm_start_x, warm_start_u=warm_start_u,
        state_bounds=state_bounds, control_bounds=control_bounds,
        cost_terms=cost_terms,
        size_N=size_N, lamb=20000.0, max_iter=25,
        etta0=1e-6,
        solver='ECOS', verbose=True,
        generate_c=True, c_output_dir=os.path.join(os.path.dirname(__file__), '..', 'waskiri_c_output'),
    )

    print("\n--- Resultado ---")
    print("Convergio:", result['converged'])
    print("Iteraciones:", result['iterations'])
    x = result['x_opt_unscaled']
    print("Cuaternion final:", np.round(x[0:4, -1], 4), " (objetivo:", end_quat.ravel(), ")")
    print("h_B final:", np.round(x[4:7, -1], 6))
