"""
Astrobee example — full SCVx trajectory optimization using sympy_cvx_converter.
Reproduces the SCVx official template (ZOH) results.

Run from project root:
    python examples/astrobee.py
"""

import sys
import os
import time
import importlib

import numpy as np
import sympy as sp
import cvxpy as cp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy_cvx_converter import generate_functions, build_problem, solve_scvx


# ==============================================================================
# 1. SymPy MODEL (Cell 21 del template)
# ==============================================================================

def define_model():
    """Define Astrobee dynamics and obstacle constraints in SymPy."""

    # Controls
    ax, ay, az = sp.symbols('ax ay az')
    taux, tauy, tauz = sp.symbols('taux tauy tauz')

    # States
    rx, ry, rz = sp.symbols('rx ry rz')
    vx, vy, vz = sp.symbols('vx vy vz')
    q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
    wx, wy, wz = sp.symbols('wx wy wz')

    # Dynamic parameters (inertia + obstacles)
    Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
    Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')
    d_obs1, d_obs2, d_obs3 = sp.symbols('d_obs1 d_obs2 d_obs3')
    c_obs1_x, c_obs1_y, c_obs1_z = sp.symbols('c_obs1_x c_obs1_y c_obs1_z')
    c_obs2_x, c_obs2_y, c_obs2_z = sp.symbols('c_obs2_x c_obs2_y c_obs2_z')
    c_obs3_x, c_obs3_y, c_obs3_z = sp.symbols('c_obs3_x c_obs3_y c_obs3_z')

    # Inertia matrix and rigid body kinematics
    I = sp.Matrix([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz]
    ])
    w_sympy = sp.Matrix([wx, wy, wz])
    tau_sympy = sp.Matrix([taux, tauy, tauz])
    f_rigid_body = sp.simplify(I.inv() * (-w_sympy.cross(I * w_sympy) + tau_sympy))

    # Dynamics: dx/dt = f(x, u)
    dynamics = [
        vx, vy, vz,                                     # dr/dt = v
        ax, ay, az,                                      # dv/dt = a (control)
        0.5 * (-wx*q1 - wy*q2 - wz*q3),                 # dq0/dt
        0.5 * (+wx*q0 + wz*q2 - wy*q3),                 # dq1/dt
        0.5 * (+wy*q0 - wz*q1 + wx*q3),                 # dq2/dt
        0.5 * (+wz*q0 + wy*q1 - wx*q2),                 # dq3/dt
        f_rigid_body[0], f_rigid_body[1], f_rigid_body[2] # dw/dt (rigid body)
    ]

    # Non-convex constraints: obstacle avoidance g(x) <= 0
    diff1 = sp.Matrix([rx, ry, rz]) - sp.Matrix([c_obs1_x, c_obs1_y, c_obs1_z])
    diff2 = sp.Matrix([rx, ry, rz]) - sp.Matrix([c_obs2_x, c_obs2_y, c_obs2_z])
    diff3 = sp.Matrix([rx, ry, rz]) - sp.Matrix([c_obs3_x, c_obs3_y, c_obs3_z])

    constraints_nc = [
        -sp.sqrt(diff1.dot(diff1)) + d_obs1,
        -sp.sqrt(diff2.dot(diff2)) + d_obs2,
        -sp.sqrt(diff3.dot(diff3)) + d_obs3,
    ]

    states = [rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz]
    controls = [ax, ay, az, taux, tauy, tauz]
    dyn_params = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz,
                  c_obs1_x, c_obs1_y, c_obs1_z,
                  c_obs2_x, c_obs2_y, c_obs2_z,
                  c_obs3_x, c_obs3_y, c_obs3_z,
                  d_obs1, d_obs2, d_obs3]

    return states, controls, dynamics, constraints_nc, dyn_params


# ==============================================================================
# 2. CODE GENERATION
# ==============================================================================

def run_codegen(states, controls, dynamics, constraints_nc, dyn_params):
    """Generate codegen file (slow — ~2-5 min for Astrobee). Skips if file exists."""

    codegen_path = os.path.join(os.path.dirname(__file__), 'astrobee_codegen.py')

    if os.path.exists(codegen_path):
        print("Code-gen file already exists, skipping generation.")
    else:
        print("Generating code-gen functions (this takes a few minutes)...")
        t0 = time.time()
        generate_functions(
            states=states,
            controls=controls,
            dynamics=dynamics,
            nonconvex_constraints=constraints_nc,
            dynamic_parameters=dyn_params,
            output_file=codegen_path,
        )
        print(f"Code-gen done in {time.time() - t0:.1f}s")

    # Import the generated module
    spec = importlib.util.spec_from_file_location("astrobee_codegen", codegen_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    return {
        'f': mod.f_SCVx,
        'A': mod.A_no_discrete_no_scaled_SCVx,
        'B': mod.B_no_discrete_no_scaled_SCVx,
        'y': mod.y_no_discrete_no_scaled_SCVx,
        'g': mod.g_SCVx,
        'C': mod.C_no_discrete_no_scaled_SCVx,
        'D': mod.D_no_discrete_no_scaled_SCVx,
        'z': mod.z_no_discrete_no_scaled_SCVx,
    }


# ==============================================================================
# 3. NUMERICAL PARAMETERS (Cell 29 del template)
# ==============================================================================

def get_parameters():
    """Return all numerical parameters for the Astrobee problem."""

    T = 50
    tf = 200.0
    tau = tf / T
    size_N = 20

    # Boundary conditions
    start_r = np.array([[0.0], [-1.0], [1.0]])
    end_r = np.array([[5.0], [2.0], [1.6]])
    start_v = np.array([[0.001], [-0.002], [0.0]])
    end_v = np.array([[0.0], [0.0], [0.0]])
    start_quat = np.array([[0], [1], [0], [0]], dtype=float)
    end_quat = np.array([[0], [0], [1], [0]], dtype=float)
    start_w = np.array([[0.0], [0.0], [0.0]])
    end_w = np.array([[0.0], [0.0], [0.0]])

    start_pos = np.block([[start_r], [start_v], [start_quat], [start_w]])
    end_pos = np.block([[end_r], [end_v], [end_quat], [end_w]])

    # Obstacles
    c_obs1 = np.array([[1.1], [-0.5], [1.0]])
    c_obs2 = np.array([[2.6], [0.5], [1.1]])
    c_obs3 = np.array([[4.0], [1.6], [1.2]])
    d_obs1 = 0.8
    d_obs2 = 0.8
    d_obs3 = 0.8

    # Inertia
    J1, J2, J3 = 0.1083, 0.1083, 0.1083

    # Limits
    acc_max = 20 / 7.2 * 0.001
    torq_max = 100 * 0.000001
    vel_max = 5.0
    omega_max = 5.0 * np.pi / 180

    # Dynamic parameters vector (must match dyn_params order)
    dynamic_parameters = np.array([
        J1, J2, J3,           # Ixx, Iyy, Izz
        0.01, 0.02, 0.01,     # Ixy, Ixz, Iyz
        c_obs1[0, 0], c_obs1[1, 0], c_obs1[2, 0],
        c_obs2[0, 0], c_obs2[1, 0], c_obs2[2, 0],
        c_obs3[0, 0], c_obs3[1, 0], c_obs3[2, 0],
        d_obs1, d_obs2, d_obs3,
    ])

    # Scaling matrices
    S_u = np.diag([
        2 * 10 * acc_max, 2 * 10 * acc_max, 2 * 10 * acc_max,
        2 * 0.1 * torq_max, 2 * 0.1 * torq_max, 2 * 0.1 * torq_max,
    ])
    c_u = np.array([
        [-(10 * acc_max)], [-(10 * acc_max)], [-(10 * acc_max)],
        [-(0.1 * torq_max)], [-(0.1 * torq_max)], [-(0.1 * torq_max)],
    ])

    S_x = np.diag([
        2*30, 2*30, 2*30,
        2*vel_max, 2*vel_max, 2*vel_max,
        2*1, 2*1, 2*1, 2*1,
        2*0.1*omega_max, 2*0.1*omega_max, 2*0.1*omega_max,
    ])
    c_x = np.array([
        [-30], [-30], [-30],
        [-vel_max], [-vel_max], [-vel_max],
        [-1], [-1], [-1], [-1],
        [-(0.1 * omega_max)], [-(0.1 * omega_max)], [-(0.1 * omega_max)],
    ])

    return {
        'T': T, 'tf': tf, 'tau': tau, 'size_N': size_N,
        'start_pos': start_pos, 'end_pos': end_pos,
        'dynamic_parameters': dynamic_parameters,
        'S_x': S_x, 'c_x': c_x, 'S_u': S_u, 'c_u': c_u,
        'acc_max': acc_max, 'torq_max': torq_max,
        'vel_max': vel_max, 'omega_max': omega_max,
    }


# ==============================================================================
# 4. WARM START (Cell 33 del template)
# ==============================================================================

def slerp(q1, q2, num_samples):
    """Spherical linear interpolation between two quaternions."""
    dot = np.dot(q1, q2)
    if dot < 0.0:
        q2 = -q2
        dot = -dot
    dot = np.clip(dot, -1.0, 1.0)
    theta_0 = np.arccos(dot)

    if np.abs(theta_0) < 1e-6:
        return np.linspace(q1, q2, num_samples)

    sin_theta_0 = np.sin(theta_0)
    quaternions = []
    for i in range(num_samples):
        t = i / (num_samples - 1)
        theta = theta_0 * t
        sin_theta = np.sin(theta)
        s0 = np.cos(theta) - dot * sin_theta / sin_theta_0
        s1 = sin_theta / sin_theta_0
        quaternions.append(s0 * q1 + s1 * q2)

    return np.array(quaternions)


def compute_angular_velocity(quaternions, dt):
    """Compute angular velocity from quaternion trajectory."""
    from scipy.spatial.transform import Rotation as R

    rotations = R.from_quat(quaternions)
    angular_velocities = [[0, 0, 0]]

    for i in range(len(rotations) - 1):
        delta_rot = rotations[i + 1] * rotations[i].inv()
        log_rot = delta_rot.as_rotvec() / dt
        angular_velocities.append(log_rot)

    return np.array(angular_velocities)


def build_warm_start(params):
    """Build initial trajectory guess (unscaled)."""

    T = params['T']
    tau = params['tau']
    start_pos = params['start_pos']
    end_pos = params['end_pos']
    nx = 13
    nu = 6

    # Position: linear interpolation
    ohx_pos = np.column_stack([
        np.linspace(start_pos[0, 0], end_pos[0, 0], T + 1),
        np.linspace(start_pos[1, 0], end_pos[1, 0], T + 1),
        np.linspace(start_pos[2, 0], end_pos[2, 0], T + 1),
    ]).T  # (3, T+1)

    # Velocity: finite differences
    ohx_vel = np.zeros_like(ohx_pos)
    for t in range(1, T):
        ohx_vel[:, t] = (ohx_pos[:, t] - ohx_pos[:, t - 1]) / tau
    ohx_vel[:, 0:1] = start_pos[3:6, 0:1]
    ohx_vel[:, T:T+1] = end_pos[3:6, 0:1]

    # Quaternion: slerp
    ohx_quat = slerp(start_pos[6:10, 0], end_pos[6:10, 0], T + 1).T  # (4, T+1)

    # Angular velocity from quaternion trajectory
    ohx_angvel = compute_angular_velocity(ohx_quat.T, tau).T  # (3, T+1)

    # Assemble
    warm_x = np.zeros((nx, T + 1))
    warm_x[0:3, :] = ohx_pos
    warm_x[3:6, :] = ohx_vel
    warm_x[6:10, :] = ohx_quat
    warm_x[10:13, :] = ohx_angvel

    warm_u = np.zeros((nu, T))

    return warm_x, warm_u


# ==============================================================================
# 5. CONVEX CONSTRAINTS CALLBACK (Cell 36 del template)
# ==============================================================================

def astrobee_convex_constraints(p):
    """Astrobee-specific convex constraints: velocity, angular velocity, accel, torque."""

    x = p['variables']['x']
    u = p['variables']['u']
    S_x = p['S_x_scaling']
    c_x = p['c_x_scaling']
    S_u = p['S_u_scaling']
    c_u = p['c_u_scaling']
    T = p['T']
    ep = p['extra_params']

    constraints = []

    for k in range(T + 1):
        # Velocity limit: ||v|| <= vel_max  (states 3:6)
        constraints += [cp.norm(S_x[3:6, 3:6] @ x[3:6, k:k+1]
                                + c_x[3:6, 0:1], 2) <= ep['vel_max']]
        # Angular velocity limit: ||w|| <= omega_max  (states 10:13)
        constraints += [cp.norm(S_x[10:13, 10:13] @ x[10:13, k:k+1]
                                + c_x[10:13, 0:1], 2) <= ep['omega_max']]

    for k in range(T):
        # Acceleration limit: ||a||_1 <= acc_max  (controls 0:3)
        constraints += [cp.norm(S_u[0:3, 0:3] @ u[0:3, k:k+1]
                                + c_u[0:3, 0:1], 1) <= ep['acc_max']]
        # Torque limit: ||tau||_1 <= torq_max  (controls 3:6)
        constraints += [cp.norm(S_u[3:6, 3:6] @ u[3:6, k:k+1]
                                + c_u[3:6, 0:1], 1) <= ep['torq_max']]

    return constraints


# ==============================================================================
# 6. MAIN
# ==============================================================================

def main():
    print("=" * 60)
    print("  ASTROBEE SCVx — sympy_cvx_converter example")
    print("=" * 60)

    # 1. Model
    print("\n[1/6] Defining SymPy model...")
    states, controls, dynamics, constraints_nc, dyn_params = define_model()
    print(f"  nx={len(states)}, nu={len(controls)}, ng={len(constraints_nc)}")

    # 2. Code generation
    print("\n[2/6] Code generation...")
    funcs = run_codegen(states, controls, dynamics, constraints_nc, dyn_params)

    # 3. Parameters
    print("\n[3/6] Loading parameters...")
    params = get_parameters()
    T = params['T']
    print(f"  T={T}, tf={params['tf']}, tau={params['tau']:.1f}")

    # 4. Warm start
    print("\n[4/6] Building warm start...")
    warm_x, warm_u = build_warm_start(params)
    print(f"  warm_x shape: {warm_x.shape}, warm_u shape: {warm_u.shape}")

    # 5. Build CVXPY problem
    print("\n[5/6] Building CVXPY problem...")
    prob_dict = build_problem(
        nx=13, nu=6, T=T, ng=3,
        convex_constraints_fn=astrobee_convex_constraints,
        extra_params={
            'vel_max': (), 'omega_max': (),
            'acc_max': (), 'torq_max': (),
        },
    )
    print(f"  DPP: {prob_dict['problem'].is_dpp()}")
    print(f"  Constraints: {len(prob_dict['problem'].constraints)}")

    # 6. Solve SCVx
    print("\n[6/6] Running SCVx loop...")
    t0 = time.time()

    result = solve_scvx(prob_dict, funcs, config={
        'tf': params['tf'],
        'size_N': params['size_N'],
        'dynamic_parameters': params['dynamic_parameters'],
        'S_x': params['S_x'], 'c_x': params['c_x'],
        'S_u': params['S_u'], 'c_u': params['c_u'],
        'start_pos': params['start_pos'],
        'end_pos': params['end_pos'],
        'warm_start_x': warm_x,
        'warm_start_u': warm_u,
        'rho0': 0.0,
        'rho1': 0.1,
        'rho2': 0.7,
        'etta0': 1e-8,
        'etta1': 10.0,
        'etta_init': 1.0,
        'beta_sh': 2.0,
        'beta_gr': 2.0,
        'lamb': 1000.0,
        'e_tol': 0.005,
        'epsilon_stop': 1e-5,
        'max_iter': 25,
        'solver': 'ECOS',
        'verbose': True,
        'extra_param_values': {
            'vel_max': params['vel_max'],
            'omega_max': params['omega_max'],
            'acc_max': params['acc_max'],
            'torq_max': params['torq_max'],
        },
    })

    total_time = time.time() - t0

    # 7. Results
    print("\n" + "=" * 60)
    print("  RESULTS")
    print("=" * 60)
    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Total time: {total_time:.2f}s")

    x_opt = result['x_opt_unscaled']
    u_opt = result['u_opt_unscaled']

    print(f"\n  Start position: {x_opt[0:3, 0].round(4)}")
    print(f"  End position:   {x_opt[0:3, -1].round(4)}")
    print(f"  Expected start: {params['start_pos'][0:3, 0].round(4)}")
    print(f"  Expected end:   {params['end_pos'][0:3, 0].round(4)}")

    err_start = np.linalg.norm(x_opt[0:3, 0:1] - params['start_pos'][0:3])
    err_end = np.linalg.norm(x_opt[0:3, -1:] - params['end_pos'][0:3])
    print(f"\n  Position error start: {err_start:.2e}")
    print(f"  Position error end:   {err_end:.2e}")

    # Velocity profile
    max_vel = np.max(np.linalg.norm(x_opt[3:6, :], axis=0))
    max_omega = np.max(np.linalg.norm(x_opt[10:13, :], axis=0))
    print(f"\n  Max velocity: {max_vel:.4f} m/s (limit: {params['vel_max']})")
    print(f"  Max omega:    {max_omega:.6f} rad/s (limit: {params['omega_max']:.6f})")

    return result


if __name__ == '__main__':
    main()
