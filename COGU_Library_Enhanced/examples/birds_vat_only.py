"""
BIRDS-TM (VAT only) via solve_trajectory() — Fase 4
+ Open-loop RK4 validation
Control de actitud puro:
    nx = 7  -> quaternion + body angular momentum
    nu = 4  -> VAT thrust magnitudes
The SCVx solution is validated independently using:
    - Zero-order hold control
    - RK4 integration
    - Original nonlinear rotational dynamics
"""

import sys
import os
import time
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# =============================================================================
# IMPORT LIBRARY
# =============================================================================
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sympy_cvx_converter import solve_trajectory
from sympy_cvx_converter.utils import slerp, angular_vel

# =============================================================================
# 1. SYMPY MODEL
# =============================================================================
b1th, b2th, b3th, b4th = sp.symbols('b1th b2th b3th b4th')
q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
hbx, hby, hbz = sp.symbols('hbx hby hbz')
Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')

# -------------------------------------------------------------------------
# Thruster positions
# -------------------------------------------------------------------------
rth1x, rth1y, rth1z = sp.symbols('rth1x rth1y rth1z')
rth2x, rth2y, rth2z = sp.symbols('rth2x rth2y rth2z')
rth3x, rth3y, rth3z = sp.symbols('rth3x rth3y rth3z')
rth4x, rth4y, rth4z = sp.symbols('rth4x rth4y rth4z')

rth1 = sp.Matrix([rth1x, rth1y, rth1z])
rth2 = sp.Matrix([rth2x, rth2y, rth2z])
rth3 = sp.Matrix([rth3x, rth3y, rth3z])
rth4 = sp.Matrix([rth4x, rth4y, rth4z])

# -------------------------------------------------------------------------
# Thruster force directions
# -------------------------------------------------------------------------
fth1x, fth1y, fth1z = sp.symbols('fth1x fth1y fth1z')
fth2x, fth2y, fth2z = sp.symbols('fth2x fth2y fth2z')
fth3x, fth3y, fth3z = sp.symbols('fth3x fth3y fth3z')
fth4x, fth4y, fth4z = sp.symbols('fth4x fth4y fth4z')

fth1 = sp.Matrix([fth1x, fth1y, fth1z])
fth2 = sp.Matrix([fth2x, fth2y, fth2z])
fth3 = sp.Matrix([fth3x, fth3y, fth3z])
fth4 = sp.Matrix([fth4x, fth4y, fth4z])

# -------------------------------------------------------------------------
# Inertia
# -------------------------------------------------------------------------
J = sp.Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
hb_sympy = sp.Matrix([hbx, hby, hbz])

# =============================================================================
# Quaternion dynamics
# =============================================================================
f_quaternion = 0.5 * sp.Matrix([[-q1, -q2, -q3], [q0, -q3, q2], [q3, q0, -q1], [-q2, q1, q0]]) * (J.inv() * hb_sympy)

# =============================================================================
# Thruster torque
# =============================================================================
torque_th = b1th * (rth1.cross(fth1)) + b2th * (rth2.cross(fth2)) + b3th * (rth3.cross(fth3)) + b4th * (rth4.cross(fth4))

# =============================================================================
# Rigid body dynamics
# =============================================================================
f_ridig_body = -(J.inv() * hb_sympy).cross(hb_sympy) + torque_th
f_ridig_body = sp.simplify(f_ridig_body)

dynamics = [f_quaternion[0], f_quaternion[1], f_quaternion[2], f_quaternion[3], f_ridig_body[0], f_ridig_body[1], f_ridig_body[2]]
states = [q0, q1, q2, q3, hbx, hby, hbz]
controls = [b1th, b2th, b3th, b4th]

dyn_params_sym = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz, rth1x, rth1y, rth1z, rth2x, rth2y, rth2z, rth3x, rth3y, rth3z, rth4x, rth4y, rth4z, fth1x, fth1y, fth1z, fth2x, fth2y, fth2z, fth3x, fth3y, fth3z, fth4x, fth4y, fth4z]

# =============================================================================
# 2. NUMERICAL PARAMETERS
# =============================================================================
T = 31 - 1
tf = 480.0
tau_val = tf / T

# -------------------------------------------------------------------------
# Inertia
# -------------------------------------------------------------------------
J1 = 0.0333
J2 = 0.0333
J3 = 0.0067
Jxy = 0.0001
Jxz = 0.0001
Jyz = 0.0001

Jsat = np.array([[J1, Jxy, Jxz], [Jxy, J2, Jyz], [Jxz, Jyz, J3]])
Jsat_inv = np.linalg.inv(Jsat)

# -------------------------------------------------------------------------
# Thruster locations
# -------------------------------------------------------------------------
rth1_d = np.array([[-26.43], [-26.61], [-147.785]]) * 0.001
rth2_d = np.array([[26.43], [-26.61], [-147.785]]) * 0.001
rth3_d = np.array([[26.43], [26.61], [-147.785]]) * 0.001
rth4_d = np.array([[-26.43], [26.61], [-147.785]]) * 0.001

# -------------------------------------------------------------------------
# VAT directions
# -------------------------------------------------------------------------
theta_VAT = 40 * np.pi / 180

fth1_d = np.array([[0], [np.sin(theta_VAT)], [np.cos(theta_VAT)]])
fth2_d = np.array([[0], [np.sin(theta_VAT)], [np.cos(theta_VAT)]])
fth3_d = np.array([[0], [-np.sin(theta_VAT)], [np.cos(theta_VAT)]])
fth4_d = np.array([[0], [-np.sin(theta_VAT)], [np.cos(theta_VAT)]])

# Normalize directions
fth1_d /= np.linalg.norm(fth1_d)
fth2_d /= np.linalg.norm(fth2_d)
fth3_d /= np.linalg.norm(fth3_d)
fth4_d /= np.linalg.norm(fth4_d)

# =============================================================================
# Dynamic parameters
# =============================================================================
dyn_params_val = np.array([J1, J2, J3, Jxy, Jxz, Jyz, rth1_d[0, 0], rth1_d[1, 0], rth1_d[2, 0], rth2_d[0, 0], rth2_d[1, 0], rth2_d[2, 0], rth3_d[0, 0], rth3_d[1, 0], rth3_d[2, 0], rth4_d[0, 0], rth4_d[1, 0], rth4_d[2, 0], fth1_d[0, 0], fth1_d[1, 0], fth1_d[2, 0], fth2_d[0, 0], fth2_d[1, 0], fth2_d[2, 0], fth3_d[0, 0], fth3_d[1, 0], fth3_d[2, 0], fth4_d[0, 0], fth4_d[1, 0], fth4_d[2, 0]])

# =============================================================================
# Physical limits
# =============================================================================
omegab_max = 1.5 * np.pi / 180
Fth_max = 25e-6

# =============================================================================
# 3. BOUNDARY CONDITIONS
# =============================================================================

#start_quat = np.array([[0.5],[0.5],[0.5],[0.5]])                       # |q| = 1
#end_quat   = np.array([[0.8660254],[0.2886751],[0.2886751],[0.2886751]])  # |q| = 1

start_quat = np.array([[1], [0], [0], [0]])
end_quat = np.array([[0], [1], [0], [0]])

start_hb = Jsat @ np.array([[0.2 * np.pi / 180], [0.1 * np.pi / 180], [-0.2 * np.pi / 180]])
end_hb = Jsat @ np.zeros((3, 1))

start_pos = np.block([[start_quat], [start_hb]]).astype(float)
end_pos = np.block([[end_quat], [end_hb]]).astype(float)

# =============================================================================
# 4. SCALING
# =============================================================================
S_x = np.diag([2, 2, 2, 2, 2 * J1 * omegab_max, 2 * J2 * omegab_max, 2 * J3 * omegab_max])
c_x = np.array([[-1], [-1], [-1], [-1], [-J1 * omegab_max], [-J2 * omegab_max], [-J3 * omegab_max]])
S_u = np.diag([2 * Fth_max, 2 * Fth_max, 2 * Fth_max, 2 * Fth_max])
c_u = np.array([[-Fth_max], [-Fth_max], [-Fth_max], [-Fth_max]])

# =============================================================================
# 5. WARM START
# =============================================================================
quat = slerp(start_pos[0:4, 0], end_pos[0:4, 0], T + 1).T
omega = angular_vel(quat.T, tau_val).T
hb_ws = Jsat @ omega

warm_x = np.zeros((7, T + 1))
warm_x[0:4, :] = quat
warm_x[4:7, :] = hb_ws
warm_u = np.zeros((4, T))

# =============================================================================
# 6. COST
# =============================================================================
inv_S_x = np.linalg.inv(S_x)
x_goal = (inv_S_x @ (end_pos - c_x))[0:7]

cost_terms = [{'kind': 'norm1', 'var': 'u', 'slice': slice(0, 4), 'coeff': 'tau', 'weight': 0.5}]

# =============================================================================
# 7. CONSTRAINTS
# =============================================================================
inv_Jb_S = np.linalg.inv(Jsat) @ S_x[4:7, 4:7]
inv_Jb_c = np.linalg.inv(Jsat) @ c_x[4:7, 0:1]

state_bounds = [{'kind': 'norm', 'slice': slice(4, 7), 'norm_type': 'norm2', 'limit': omegab_max, 'S': inv_Jb_S, 'c': inv_Jb_c}]

control_bounds = [
    {'kind': 'box', 'slice': slice(0, 4), 'lower': 0.0, 'upper': Fth_max},
    {'kind': 'box_affine', 'c_affine': np.array([[1.0, 1.0, 0.0, 0.0]]), 'lower': 0.0, 'upper': Fth_max},
    {'kind': 'box_affine', 'c_affine': np.array([[0.0, 0.0, 1.0, 1.0]]), 'lower': 0.0, 'upper': Fth_max}
]

# =============================================================================
# 8. SOLVE SCVX
# =============================================================================
print("=" * 60)
print("BIRDS-TM (VAT only) -- Fase 4")
print("=" * 60)

t0 = time.time()

result = solve_trajectory(states=states, controls=controls, dynamics=dynamics, start=start_pos, end=end_pos, T=T, tf=tf, dynamic_parameters_sym=dyn_params_sym, dynamic_parameters_val=dyn_params_val, S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u, warm_start_x=warm_x, warm_start_u=warm_u, cost_terms=cost_terms, state_bounds=state_bounds, control_bounds=control_bounds, size_N=20, lamb=5000.0, max_iter = 20, rho0=0.0, rho1=0.1, rho2=0.7, etta0=0.001, etta1=10.0, etta_init=1.0, beta_sh=2.0, beta_gr=2.0, e_tol=0.005, epsilon_stop=0.00001, solver='ECOS', verbose=True, generate_c=True, c_output_dir="COGU_Library_Enhanced/birds_vat_c_output")

total = time.time() - t0

# =============================================================================
# 9. SOLVER RESULTS
# =============================================================================
x_value = result['x_opt_unscaled']
u_value = result['u_opt_unscaled']
history = result['history']

print("\n" + "=" * 60)
print("RESULTADOS BIRDS-TM")
print("=" * 60)
print(f"Tiempo total:  {total:.1f} s")
print(f"Iteraciones:   {result['iterations']}")
print(f"Convergio:     {result['converged']}")
print(f"Quat final:    {x_value[0:4, -1].round(6)}")
print(f"Quat start:    {x_value[0:4, 0].round(6)}")
print(f"hb final:      {x_value[4:7, -1].round(9)}")

# =============================================================================
# 10. GUIDANCE ANGULAR VELOCITY
# =============================================================================
omega_guidance = Jsat_inv @ x_value[4:7, :]
omega_guidance_norm = np.linalg.norm(omega_guidance, axis=0)

print(f"\nMax ||omega|| guidance: {np.max(omega_guidance_norm):.6f} rad/s ({np.max(omega_guidance_norm) * 180/np.pi:.6f} deg/s)")
print(f"Omega limit:            {omegab_max:.6f} rad/s ({omegab_max * 180/np.pi:.6f} deg/s)")

# =============================================================================
# 11. THRUSTER TORQUE MATRIX
# =============================================================================
B = np.column_stack([np.cross(rth1_d[:, 0], fth1_d[:, 0]), np.cross(rth2_d[:, 0], fth2_d[:, 0]), np.cross(rth3_d[:, 0], fth3_d[:, 0]), np.cross(rth4_d[:, 0], fth4_d[:, 0])])

print("\nThruster torque matrix B [m]:")
print(B)

# =============================================================================
# 12. OPEN-LOOP VALIDATION
# =============================================================================

def u_interpol_zoh(t, u_value, dt_control):
    """Zero-order hold of the SCVx control."""
    T_control = u_value.shape[1]
    index = int(np.floor(t / dt_control))
    index = min(index, T_control - 1)
    return u_value[:, index]

def attitude_dynamics_open_loop(t, x, u, params):
    q = x[0:4]
    w = x[4:7]
    I = params["I"]
    I_inv = params["I_inv"]
    B = params["B"]

    torque = B @ u

    w_dot = I_inv @ (torque - np.cross(w, I @ w))

    wx, wy, wz = w
    q0, q1, q2, q3 = q
    q_dot = 0.5 * np.array([-wx*q1 - wy*q2 - wz*q3, wx*q0 + wz*q2 - wy*q3, wy*q0 - wz*q1 + wx*q3, wz*q0 + wy*q1 - wx*q2])
    return np.concatenate([q_dot, w_dot])

# =============================================================================
# RK4
# =============================================================================
def rk4_step(f, t, x, u, dt, params):
    k1 = f(t, x, u, params)
    k2 = f(t + dt/2, x + dt/2 * k1, u, params)
    k3 = f(t + dt/2, x + dt/2 * k2, u, params)
    k4 = f(t + dt, x + dt * k3, u, params)
    return x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

# =============================================================================
# OPEN-LOOP SIMULATION SETUP
# =============================================================================
params = {"I": Jsat, "I_inv": Jsat_inv, "B": B}

x_sim = np.zeros(7)
x_sim[0:4] = x_value[0:4, 0]
x_sim[4:7] = Jsat_inv @ x_value[4:7, 0]

# -------------------------------------------------------------------------
# Integration step
# -------------------------------------------------------------------------
dt_state = 0.1
dt_control = tau_val
sim_time = tf

N = int(round(sim_time / dt_state)) + 1

if abs(dt_control / dt_state - round(dt_control / dt_state)) > 1e-9:
    print(f"[warning] dt_state={dt_state}s does not divide dt_control={dt_control}s")

# =============================================================================
# HISTORY ARRAYS
# =============================================================================
q_history = np.zeros((N, 4))
w_history = np.zeros((N, 3))
u_history = np.zeros((N, 4))
t_history = np.zeros(N)

# =============================================================================
# OPEN-LOOP SIMULATION
# =============================================================================
t = 0.0

for i in range(N):
    # -------------------------------------------------------------
    # Store
    # -------------------------------------------------------------
    q_history[i] = x_sim[0:4]
    w_history[i] = x_sim[4:7]
    t_history[i] = t

    # -------------------------------------------------------------
    # ZOH control
    # -------------------------------------------------------------
    u_current = u_interpol_zoh(t, u_value, dt_control)
    u_history[i] = u_current

    # -------------------------------------------------------------
    # RK4
    # -------------------------------------------------------------
    if i < N - 1:
        # Do not integrate beyond tf
        dt = min(dt_state, tf - t)
        if dt > 0:
            x_sim = rk4_step(attitude_dynamics_open_loop, t, x_sim, u_current, dt, params)

            # -----------------------------------------------------
            # Quaternion normalization
            # -----------------------------------------------------
            q_norm = np.linalg.norm(x_sim[0:4])
            if q_norm > 0:
                x_sim[0:4] /= q_norm
            t += dt

# =============================================================================
# 13. OPEN-LOOP FINAL RESULTS
# =============================================================================
q_final_ol = q_history[-1]
w_final_ol = w_history[-1]
hb_final_ol = Jsat @ w_final_ol

print("\n" + "=" * 60)
print("OPEN-LOOP VALIDATION")
print("=" * 60)
print(f"Final quaternion OL: {q_final_ol.round(8)}")
print(f"Target quaternion:   {end_quat[:,0]}")
print(f"Final omega OL:      {w_final_ol.round(10)} rad/s")
print(f"Final omega OL:      {(w_final_ol * 180/np.pi).round(8)} deg/s")
print(f"Final hb OL:         {hb_final_ol.round(10)}")
print(f"Max ||omega|| OL:    {np.max(np.linalg.norm(w_history, axis=1)):.8f} rad/s")
print(f"Omega limit:         {omegab_max:.8f} rad/s")

# =============================================================================
# 14. FINAL ERRORS
# =============================================================================
q_target = end_quat[:, 0]
q_error = q_final_ol - q_target
q_error_norm = np.linalg.norm(q_error)
w_error_norm = np.linalg.norm(w_final_ol)

print("\n--- OPEN-LOOP ERRORS ---")
print(f"||q_final - q_target||: {q_error_norm:.8e}")
print(f"||omega_final||:        {w_error_norm:.8e} rad/s")

# =============================================================================
# 15. CONTROL STATISTICS
# =============================================================================
print("\n--- CONTROL STATISTICS ---")

for i in range(4):
    print(f"VAT{i+1}: max={np.max(u_value[i,:])*1e6:.4f} uN, mean={np.mean(u_value[i,:])*1e6:.4f} uN, active={np.sum(u_value[i,:] > 1e-12)}/{T}")

print(f"Max b1+b2: {np.max(u_value[0,:] + u_value[1,:])*1e6:.4f} uN")
print(f"Max b3+b4: {np.max(u_value[2,:] + u_value[3,:])*1e6:.4f} uN")

# =============================================================================
# 16. PLOT 1 — QUATERNION
# =============================================================================
t_nodes = np.arange(T + 1) * tau_val

plt.figure(figsize=(7, 4.5))
plt.plot(t_history, q_history[:, 0], label='q0')
plt.plot(t_history, q_history[:, 1], label='q1')
plt.plot(t_history, q_history[:, 2], label='q2')
plt.plot(t_history, q_history[:, 3], label='q3')
plt.plot(t_history, np.linalg.norm(q_history, axis=1), '--', label=r'$\|q\|$')
plt.plot(t_nodes, x_value[0:4, :].T, 'o', ms=3, alpha=0.45)
plt.xlabel('Time [s]')
plt.ylabel('Quaternion')
plt.title('Quaternion Evolution: Open-Loop vs SCVx')
plt.legend()
plt.grid(True)
plt.tight_layout()

# =============================================================================
# 17. PLOT 2 — ANGULAR VELOCITY
# =============================================================================
plt.figure(figsize=(7, 4.5))
plt.plot(t_history, w_history[:, 0] * 180 / np.pi, label=r'$\omega_x$')
plt.plot(t_history, w_history[:, 1] * 180 / np.pi, label=r'$\omega_y$')
plt.plot(t_history, w_history[:, 2] * 180 / np.pi, label=r'$\omega_z$')
plt.plot(t_history, np.linalg.norm(w_history, axis=1) * 180 / np.pi, '--', label=r'$\|\omega\|$')
plt.plot(t_nodes, omega_guidance.T * 180 / np.pi, 'o', ms=3, alpha=0.45)
plt.axhline(omegab_max * 180 / np.pi, linestyle='--', alpha=0.7, label=r'$\omega_{\max}$')
plt.axhline(-omegab_max * 180 / np.pi, linestyle='--', alpha=0.7)
plt.xlabel('Time [s]')
plt.ylabel(r'Angular velocity [deg/s]')
plt.title('Angular Velocity: Open-Loop vs SCVx')
plt.legend()
plt.grid(True)
plt.tight_layout()

# =============================================================================
# 18. PLOT 3 — VAT CONTROL
# =============================================================================
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(7, 7))

for i in range(4):
    axs[i].step(t_history, u_history[:, i] * 1e6, where='post', label=f'VAT{i+1}')
    axs[i].set_ylabel(f'VAT{i+1} [uN]')
    axs[i].grid(True)
    axs[i].legend(loc='upper right')

axs[-1].set_xlabel('Time [s]')
fig.suptitle('VAT Thrust Sequence — SCVx ZOH')
fig.tight_layout()

# =============================================================================
# 19. PLOT 4 — OPEN-LOOP VS GUIDANCE MOMENTUM
# =============================================================================
hb_history = (Jsat @ w_history.T).T

plt.figure(figsize=(7, 4.5))
plt.plot(t_history, hb_history[:, 0], label=r'$h_x$')
plt.plot(t_history, hb_history[:, 1], label=r'$h_y$')
plt.plot(t_history, hb_history[:, 2], label=r'$h_z$')
plt.plot(t_nodes, x_value[4:7, :].T, 'o', ms=3, alpha=0.45)
plt.xlabel('Time [s]')
plt.ylabel(r'Body angular momentum [N m s]')
plt.title('Body Angular Momentum: Open-Loop vs SCVx')
plt.legend()
plt.grid(True)
plt.tight_layout()

# =============================================================================
# 20. SHOW
# =============================================================================
plt.show()