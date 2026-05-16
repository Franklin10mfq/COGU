"""
plot_trajectory.py -- Graficador de trayectoria SCVx para Astrobee.

Resuelve el problema de guiado usando la COGU_Dev_Library (sympy_cvx_converter)
y grafica los resultados en el mismo formato de referencia que el notebook
SCVx_official_template_cvxpy_cvxpygen_ZOH.ipynb (6 plots estado/control +
animacion 3D Plotly). Los datos graficados provienen del solver de la libreria,
no del notebook.

"""

import sys, os, time
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from scipy.spatial.transform import Rotation

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from sympy_cvx_converter import solve_trajectory

# =============================================================================
# 1. MODELO SYMPY (identico a verify_final.py)
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

I_mat = sp.Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
w_sym = sp.Matrix([wx, wy, wz])
tau_sym = sp.Matrix([taux, tauy, tauz])
f_rb = sp.simplify(I_mat.inv() * (-w_sym.cross(I_mat * w_sym) + tau_sym))

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
# 2. PARAMETROS NUMERICOS
# =============================================================================
T = 30
tf = 200.0
acc_max   = 20 / 7.2 * 0.001
torq_max  = 100 * 0.000001
vel_max   = 5.0
omega_max = 5.0 * np.pi / 180

start_pos = np.array([[0.0],[-1.0],[1.0],
                       [0.001],[-0.002],[0.0],
                       [0],[1],[0],[0],
                       [0.0],[0.0],[0.0]], dtype=float)
end_pos   = np.array([[5.0],[2.0],[1.6],
                       [0.0],[0.0],[0.0],
                       [0],[0],[1],[0],
                       [0.0],[0.0],[0.0]], dtype=float)

dyn_params_val = np.array([
    0.1083, 0.1083, 0.1083, 0.01, 0.02, 0.01,
    1.1, -0.5, 1.0,   2.6, 0.5, 1.1,   4.0, 1.6, 1.2,
    0.8, 0.8, 0.8,
])

# Obstaculos: centro [x,y,z] y radio desde dyn_params_val
obstacles = [
    {'h': dyn_params_val[6:9],  'rc': dyn_params_val[15]},
    {'h': dyn_params_val[9:12], 'rc': dyn_params_val[16]},
    {'h': dyn_params_val[12:15],'rc': dyn_params_val[17]},
]

S_x = np.diag([60,60,60, 2*vel_max,2*vel_max,2*vel_max,
               2,2,2,2, 2*0.1*omega_max,2*0.1*omega_max,2*0.1*omega_max])
c_x = np.array([[-30],[-30],[-30],
                [-vel_max],[-vel_max],[-vel_max],
                [-1],[-1],[-1],[-1],
                [-(0.1*omega_max)],[-(0.1*omega_max)],[-(0.1*omega_max)]])
S_u = np.diag([2*10*acc_max]*3 + [2*0.1*torq_max]*3)
c_u = np.array([[-(10*acc_max)]]*3 + [[-(0.1*torq_max)]]*3)

# =============================================================================
# 3. WARM START
# =============================================================================
def _slerp(qa, qb, n):
    dot = np.clip(np.dot(qa, qb), -1.0, 1.0)
    if dot < 0: qb, dot = -qb, -dot
    theta = np.arccos(dot)
    if abs(theta) < 1e-6: return np.linspace(qa, qb, n)
    return np.array([(np.cos(theta*t/(n-1)) - dot*np.sin(theta*t/(n-1))/np.sin(theta))*qa
                     + np.sin(theta*t/(n-1))/np.sin(theta)*qb for t in range(n)])

def _angular_vel(quats, dt):
    rots = Rotation.from_quat(quats)
    w = [[0,0,0]]
    for i in range(len(rots)-1):
        w.append((rots[i+1]*rots[i].inv()).as_rotvec()/dt)
    return np.array(w)

tau_val = tf / T
pos = np.column_stack([np.linspace(start_pos[i,0], end_pos[i,0], T+1) for i in range(3)]).T
vel = np.zeros_like(pos)
for t in range(1, T): vel[:, t] = (pos[:, t] - pos[:, t-1]) / tau_val
vel[:, 0:1], vel[:, T:T+1] = start_pos[3:6], end_pos[3:6]
quat = _slerp(start_pos[6:10,0], end_pos[6:10,0], T+1).T
angvel = _angular_vel(quat.T, tau_val).T

warm_x = np.zeros((13, T+1))
warm_x[0:3,:] = pos; warm_x[3:6,:] = vel; warm_x[6:10,:] = quat; warm_x[10:13,:] = angvel
warm_u = np.zeros((6, T))

# =============================================================================
# 4. RESOLVER (sin generar C -- solo Python)
# =============================================================================
print("=" * 60)
print("COGU plot_trajectory.py -- resolviendo...")
print("=" * 60)
t0 = time.time()

result = solve_trajectory(
    states=states, controls=controls, dynamics=dynamics,
    start=start_pos, end=end_pos,
    T=T, tf=tf,
    nonconvex_constraints=nonconvex,
    dynamic_parameters_sym=dyn_params_sym,
    dynamic_parameters_val=dyn_params_val,
    S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
    warm_start_x=warm_x, warm_start_u=warm_u,
    state_bounds=[
        (slice(3,6), 'norm2', vel_max),
        (slice(10,13), 'norm2', omega_max),
    ],
    control_bounds=[
        (slice(0,3), 'norm1', acc_max),
        (slice(3,6), 'norm1', torq_max),
    ],
    size_N=20, max_iter=25,
    generate_c=False,
)

print(f"Resuelto en {time.time()-t0:.1f}s | iters={result['iterations']} | convergio={result['converged']}")

x = result['x_opt_unscaled']   # (13, T+1)
u = result['u_opt_unscaled']   # (6,  T)
time_x = np.linspace(0, tf, T+1)
time_u = np.linspace(0, tf, T)

# =============================================================================
# 5. PLOTS MATPLOTLIB (6 figuras de estado/control)
# =============================================================================
plt.figure(1)
plt.plot(time_x, x[6,:], label='q0')
plt.plot(time_x, x[7,:], label='q1')
plt.plot(time_x, x[8,:], label='q2')
plt.plot(time_x, x[9,:], label='q3')
plt.plot(time_x, np.linalg.norm(x[6:10,:], axis=0), label='||q||', linestyle='--')
plt.ylabel('Quaternion sequence')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

plt.figure(2)
plt.plot(time_x, x[0,:], label='x')
plt.plot(time_x, x[1,:], label='y')
plt.plot(time_x, x[2,:], label='z')
plt.ylabel('Position sequence [m]')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

plt.figure(3)
plt.plot(time_x, x[3,:], label='vx')
plt.plot(time_x, x[4,:], label='vy')
plt.plot(time_x, x[5,:], label='vz')
plt.ylabel('Velocity sequence [m/s]')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

plt.figure(4)
plt.plot(time_x, x[10,:], label='wx')
plt.plot(time_x, x[11,:], label='wy')
plt.plot(time_x, x[12,:], label='wz')
plt.ylabel('Angular velocity sequence [rad/s]')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

plt.figure(5)
plt.step(time_u, u[0,:], label='ax')
plt.step(time_u, u[1,:], label='ay')
plt.step(time_u, u[2,:], label='az')
plt.ylabel('Control acceleration sequence [m/s^2]')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

plt.figure(6)
plt.step(time_u, u[3,:], label='taux')
plt.step(time_u, u[4,:], label='tauy')
plt.step(time_u, u[5,:], label='tauz')
plt.ylabel('Control torque sequence [Nm]')
plt.xlabel('Time [s]')
plt.legend(); plt.grid(True); plt.tight_layout()

# =============================================================================
# 6. ANIMACION 3D PLOTLY (Astrobee + obstaculos)
# =============================================================================
def _create_cube(center, size, quaternion):
    half = size / 2.0
    verts = np.array([
        [-half,-half,-half],[half,-half,-half],[half,half,-half],[-half,half,-half],
        [-half,-half, half],[half,-half, half],[half,half, half],[-half,half, half]
    ])
    # quaternion del modelo: [q0,q1,q2,q3] con q0=escalar
    # scipy from_quat espera [x,y,z,w] -> reordenar
    q_scipy = np.array([quaternion[1], quaternion[2], quaternion[3], quaternion[0]])
    r = Rotation.from_quat(q_scipy)
    rv = r.apply(verts) + center

    faces = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[0,3,7,4],[1,2,6,5]]
    i_f = [f[0] for f in faces] + [f[0] for f in faces]
    j_f = [f[1] for f in faces] + [f[2] for f in faces]
    k_f = [f[2] for f in faces] + [f[3] for f in faces]

    cube = go.Mesh3d(x=rv[:,0], y=rv[:,1], z=rv[:,2],
                     i=i_f, j=j_f, k=k_f, opacity=1.0, color='white')

    eye_positions = np.array([[half,0.13,0],[half,-0.13,0]])
    re = r.apply(eye_positions) + center
    eye_traces = []
    for pos in re:
        ug, vg = np.mgrid[0:2*np.pi:10j, 0:np.pi:5j]
        er = 0.08
        xe = er*np.cos(ug)*np.sin(vg)+pos[0]
        ye = er*np.sin(ug)*np.sin(vg)+pos[1]
        ze = er*np.cos(vg)+pos[2]
        c_arr = np.full_like(xe, 0.5)
        eye_traces.append(go.Surface(x=xe, y=ye, z=ze, surfacecolor=c_arr,
                                     colorscale=[[0,'deepskyblue'],[1,'deepskyblue']],
                                     opacity=1.0, showscale=False))
    return cube, eye_traces

def _create_spheres(obs_list):
    traces = []
    for obs in obs_list:
        h, rc = obs['h'], obs['rc']
        ug, vg = np.mgrid[0:2*np.pi:16j, 0:np.pi:8j]
        xs = rc*np.cos(ug)*np.sin(vg)+h[0]
        ys = rc*np.sin(ug)*np.sin(vg)+h[1]
        zs = rc*np.cos(vg)+h[2]
        traces.append(go.Surface(x=xs, y=ys, z=zs, colorscale='Reds',
                                 opacity=0.3, showscale=False))
    return traces

x_traj, y_traj, z_traj = x[0,:], x[1,:], x[2,:]
quaternions = x[6:10,:]

traj_line = go.Scatter3d(x=x_traj, y=y_traj, z=z_traj,
                         mode='lines', line=dict(color='blue', width=2),
                         name='Trayectoria')

sphere_traces = _create_spheres(obstacles)
cube0, eyes0 = _create_cube([x_traj[0],y_traj[0],z_traj[0]], 0.5, quaternions[:,0])

# Limites de ejes con escala uniforme
pad = 0.5
max_r = max(x_traj.max()-x_traj.min(), y_traj.max()-y_traj.min(), z_traj.max()-z_traj.min()) / 2
mx, my, mz = (x_traj.max()+x_traj.min())/2, (y_traj.max()+y_traj.min())/2, (z_traj.max()+z_traj.min())/2
x_rng = [mx-max_r-pad, mx+max_r+pad]
y_rng = [my-max_r-pad, my+max_r+pad]
z_rng = [mz-max_r-pad, mz+max_r+pad]

frames = []
for k in range(T+1):
    cube_k, eyes_k = _create_cube([x_traj[k],y_traj[k],z_traj[k]], 0.5, quaternions[:,k])
    frames.append(go.Frame(data=[
        go.Scatter3d(x=x_traj, y=y_traj, z=z_traj, mode='lines',
                     line=dict(color='blue', width=2)),
        cube_k,
    ] + _create_spheres(obstacles) + eyes_k))

layout = go.Layout(
    title='Astrobee guidance COGU (SCVx, ECOS) -- COGU_Dev_Library',
    scene=dict(
        xaxis_title='x [m]', yaxis_title='y [m]', zaxis_title='z [m]',
        xaxis=dict(range=x_rng, autorange=False),
        yaxis=dict(range=y_rng, autorange=False),
        zaxis=dict(range=z_rng, autorange=False),
        aspectmode='cube',
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
    ),
    updatemenus=[dict(type='buttons', buttons=[
        dict(label='Play', method='animate',
             args=[None, dict(frame=dict(duration=80, redraw=True), fromcurrent=True)]),
        dict(label='Pause', method='animate',
             args=[[None], dict(frame=dict(duration=0, redraw=False), mode='immediate')])
    ])]
)

fig = go.Figure(data=[traj_line, cube0] + sphere_traces + eyes0, layout=layout, frames=frames)
fig.show()

plt.show()
