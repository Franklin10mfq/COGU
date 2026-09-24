"""
Astrobee free-flyer model, shared by every Astrobee example and test.

The numbers come from `COGU_algorithms/Astrobee_template/
SCVx_official_template_cvxpy_cvxpygen_ZOH.ipynb` (cells 21, 29 and 33), so a
run built from `astrobee_problem()` is directly comparable with the notebook.

13 states  r(3) v(3) q(4) w(3)
 6 controls a(3) tau(3)
 3 nonconvex constraints: keep-out spheres around three obstacles
"""

import numpy as np
import sympy as sp


# ----------------------------------------------------------------------
# Symbolic model (notebook cell 21)
# ----------------------------------------------------------------------

def astrobee_symbols():
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
        0.5 * (-wx * q1 - wy * q2 - wz * q3),
        0.5 * (+wx * q0 + wz * q2 - wy * q3),
        0.5 * (+wy * q0 - wz * q1 + wx * q3),
        0.5 * (+wz * q0 + wy * q1 - wx * q2),
        f_rb[0], f_rb[1], f_rb[2],
    ]

    r = sp.Matrix([rx, ry, rz])
    diff1 = r - sp.Matrix([c_obs1_x, c_obs1_y, c_obs1_z])
    diff2 = r - sp.Matrix([c_obs2_x, c_obs2_y, c_obs2_z])
    diff3 = r - sp.Matrix([c_obs3_x, c_obs3_y, c_obs3_z])

    nonconvex = [
        -sp.sqrt(diff1.dot(diff1)) + d_obs1,
        -sp.sqrt(diff2.dot(diff2)) + d_obs2,
        -sp.sqrt(diff3.dot(diff3)) + d_obs3,
    ]

    states = [rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz]
    controls = [ax, ay, az, taux, tauy, tauz]
    dyn_params_sym = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz,
                      c_obs1_x, c_obs1_y, c_obs1_z,
                      c_obs2_x, c_obs2_y, c_obs2_z,
                      c_obs3_x, c_obs3_y, c_obs3_z,
                      d_obs1, d_obs2, d_obs3]

    return states, controls, dynamics, nonconvex, dyn_params_sym


# ----------------------------------------------------------------------
# Numerical data (notebook cell 29)
# ----------------------------------------------------------------------

ACC_MAX = 20 / 7.2 * 0.001
TORQ_MAX = 100 * 1e-6
VEL_MAX = 5.0
OMEGA_MAX = 5.0 * np.pi / 180

START_POS = np.array([[0.0], [-1.0], [1.0],
                      [0.001], [-0.002], [0.0],
                      [0.0], [1.0], [0.0], [0.0],
                      [0.0], [0.0], [0.0]])
END_POS = np.array([[5.0], [2.0], [1.6],
                    [0.0], [0.0], [0.0],
                    [0.0], [0.0], [1.0], [0.0],
                    [0.0], [0.0], [0.0]])

DYN_PARAMS_VAL = np.array([
    0.1083, 0.1083, 0.1083, 0.01, 0.02, 0.01,
    1.1, -0.5, 1.0,   2.6, 0.5, 1.1,   4.0, 1.6, 1.2,
    0.8, 0.8, 0.8,
])

OBSTACLES = [                      # (centre, radius), for the plots
    (np.array([1.1, -0.5, 1.0]), 0.8),
    (np.array([2.6, 0.5, 1.1]), 0.8),
    (np.array([4.0, 1.6, 1.2]), 0.8),
]


def astrobee_scaling():
    S_x = np.diag([60, 60, 60,
                   2 * VEL_MAX, 2 * VEL_MAX, 2 * VEL_MAX,
                   2, 2, 2, 2,
                   2 * 0.1 * OMEGA_MAX, 2 * 0.1 * OMEGA_MAX, 2 * 0.1 * OMEGA_MAX])
    c_x = np.array([[-30], [-30], [-30],
                    [-VEL_MAX], [-VEL_MAX], [-VEL_MAX],
                    [-1], [-1], [-1], [-1],
                    [-(0.1 * OMEGA_MAX)], [-(0.1 * OMEGA_MAX)], [-(0.1 * OMEGA_MAX)]],
                   dtype=float)
    S_u = np.diag([2 * 10 * ACC_MAX] * 3 + [2 * 0.1 * TORQ_MAX] * 3)
    c_u = np.array([[-(10 * ACC_MAX)]] * 3 + [[-(0.1 * TORQ_MAX)]] * 3)
    return S_x, c_x, S_u, c_u


# ----------------------------------------------------------------------
# Warm start (notebook cell 33)
# ----------------------------------------------------------------------

def astrobee_warm_start(T, tf, u_nodes=None, legacy_quaternion=True):
    """
    Straight line in position, slerp in attitude, finite differences for the
    rest -- exactly the notebook's construction.

    legacy_quaternion
        True reproduces the notebook bit for bit: it feeds scalar-first
        quaternions to scipy's scalar-last `from_quat`, so the warm-start
        angular velocity is not the angular velocity of the slerp.  False
        uses the corrected convention.
    """
    from scvx_toolbox.utils import angular_velocity, slerp

    if u_nodes is None:
        u_nodes = T
    tau = tf / T

    pos = np.column_stack([
        np.linspace(START_POS[i, 0], END_POS[i, 0], T + 1) for i in range(3)
    ]).T
    vel = np.zeros_like(pos)
    for k in range(1, T):
        vel[:, k] = (pos[:, k] - pos[:, k - 1]) / tau
    vel[:, 0:1] = START_POS[3:6]
    vel[:, T:T + 1] = END_POS[3:6]

    quat = slerp(START_POS[6:10, 0], END_POS[6:10, 0], T + 1).T
    angvel = angular_velocity(quat.T, tau,
                              scalar_first=not legacy_quaternion).T

    warm_x = np.zeros((13, T + 1))
    warm_x[0:3, :] = pos
    warm_x[3:6, :] = vel
    warm_x[6:10, :] = quat
    warm_x[10:13, :] = angvel
    warm_u = np.zeros((6, u_nodes))
    return warm_x, warm_u


# ----------------------------------------------------------------------
# One call that assembles everything solve_trajectory() needs
# ----------------------------------------------------------------------

def astrobee_problem(T=30, tf=200.0, time_mode='ZOH', legacy_quaternion=True):
    """kwargs for `solve_trajectory`, minus the SCvx tuning."""
    from scvx_toolbox.modes import control_node_count, parse_time_mode

    states, controls, dynamics, nonconvex, dyn_sym = astrobee_symbols()
    S_x, c_x, S_u, c_u = astrobee_scaling()
    _, hold, _ = parse_time_mode(time_mode)
    u_nodes = control_node_count(hold, T)
    warm_x, warm_u = astrobee_warm_start(T, tf, u_nodes, legacy_quaternion)

    return dict(
        states=states, controls=controls, dynamics=dynamics,
        start=START_POS, end=END_POS, T=T, tf=tf,
        time_mode=time_mode,
        nonconvex_constraints=nonconvex,
        dynamic_parameters_sym=dyn_sym,
        dynamic_parameters_val=DYN_PARAMS_VAL,
        S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
        warm_start_x=warm_x, warm_start_u=warm_u,
        cost_terms=[
            {'kind': 'sumsq', 'var': 'u', 'slice': slice(None),
             'coeff': 'sqrt_tau', 'weight': 1.0, 'k_range': 'U'},
        ],
        state_bounds=[
            (slice(3, 6), 'norm2', VEL_MAX),
            (slice(10, 13), 'norm2', OMEGA_MAX),
        ],
        control_bounds=[
            (slice(0, 3), 'norm1', ACC_MAX),
            (slice(3, 6), 'norm1', TORQ_MAX),
        ],
    )
