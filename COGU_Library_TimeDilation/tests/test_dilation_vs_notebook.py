"""
Replicate cell 45 of BIRDS_TM_Th_cvxpy_cvxpygen_l2_time_dilation.ipynb exactly.

Notebook output for reference:
    Iter 1: L=34.9272  J=302.1056  oJ=707.6715  tf=298.6339
    Iter 2: L=49.3521  J=251.4258  oJ=302.1056  tf=585.6905  rho=0.60285
    Iter 3: L=34.0228  J=155.0715  oJ=251.4258  tf=410.7513  rho=0.20051
    Iter 4: L=35.3089  J= 43.7126  oJ=155.0715  tf=445.8738  rho=0.44321

Note: the notebook's L omits the constant alpha_p*c_p = 0.05*10 = 0.5 that its
J includes.  This library adds the constant back to L so that L and J are on
the same footing (rho is otherwise biased), hence the expected +0.5 offset.
"""
import sys, os
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from scvx_toolbox.pipeline import solve_trajectory
from scvx_toolbox.utils import slerp, angular_velocity

b1th, b2th, b3th, b4th = sp.symbols('b1th b2th b3th b4th')
q0, q1, q2, q3 = sp.symbols('q0 q1 q2 q3')
hbx, hby, hbz = sp.symbols('hbx hby hbz')
Ixx, Iyy, Izz = sp.symbols('Ixx Iyy Izz')
Ixy, Ixz, Iyz = sp.symbols('Ixy Ixz Iyz')
rth1x, rth1y, rth1z = sp.symbols('rth1x rth1y rth1z')
rth2x, rth2y, rth2z = sp.symbols('rth2x rth2y rth2z')
rth3x, rth3y, rth3z = sp.symbols('rth3x rth3y rth3z')
rth4x, rth4y, rth4z = sp.symbols('rth4x rth4y rth4z')
fth1x, fth1y, fth1z = sp.symbols('fth1x fth1y fth1z')
fth2x, fth2y, fth2z = sp.symbols('fth2x fth2y fth2z')
fth3x, fth3y, fth3z = sp.symbols('fth3x fth3y fth3z')
fth4x, fth4y, fth4z = sp.symbols('fth4x fth4y fth4z')

rth1 = sp.Matrix([rth1x, rth1y, rth1z]); fth1 = sp.Matrix([fth1x, fth1y, fth1z])
rth2 = sp.Matrix([rth2x, rth2y, rth2z]); fth2 = sp.Matrix([fth2x, fth2y, fth2z])
rth3 = sp.Matrix([rth3x, rth3y, rth3z]); fth3 = sp.Matrix([fth3x, fth3y, fth3z])
rth4 = sp.Matrix([rth4x, rth4y, rth4z]); fth4 = sp.Matrix([fth4x, fth4y, fth4z])

J = sp.Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
hb = sp.Matrix([hbx, hby, hbz])
f_quat = 0.5*sp.Matrix([[-q1,-q2,-q3],[q0,-q3,q2],[q3,q0,-q1],[-q2,q1,q0]])*(J.inv()*hb)
torque = (b1th*rth1.cross(fth1) + b2th*rth2.cross(fth2)
          + b3th*rth3.cross(fth3) + b4th*rth4.cross(fth4))
f_body = sp.simplify(-(J.inv()*hb).cross(hb) + torque)

dynamics = [f_quat[0], f_quat[1], f_quat[2], f_quat[3], f_body[0], f_body[1], f_body[2]]
states = [q0,q1,q2,q3,hbx,hby,hbz]
controls = [b1th,b2th,b3th,b4th]
dps = [Ixx,Iyy,Izz,Ixy,Ixz,Iyz,
       rth1x,rth1y,rth1z, rth2x,rth2y,rth2z, rth3x,rth3y,rth3z, rth4x,rth4y,rth4z,
       fth1x,fth1y,fth1z, fth2x,fth2y,fth2z, fth3x,fth3y,fth3z, fth4x,fth4y,fth4z]

T = 30
tf_min, tf_max, tf_init = 10.0, 720.0, 600.0
tau_abs = tf_init / T

J1, J2, J3 = 0.0333, 0.0333, 0.0067
Jxy = Jxz = Jyz = 0.0001
Jsat = np.array([[J1,Jxy,Jxz],[Jxy,J2,Jyz],[Jxz,Jyz,J3]])

sat_CG = np.array([[0],[0],[-5]])*0.001
rth1_d = np.array([[-26.43],[-26.61],[-147.785]])*0.001 - sat_CG
rth2_d = np.array([[ 26.43],[-26.61],[-147.785]])*0.001 - sat_CG
rth3_d = np.array([[ 26.43],[ 26.61],[-147.785]])*0.001 - sat_CG
rth4_d = np.array([[-26.43],[ 26.61],[-147.785]])*0.001 - sat_CG

th = 40*np.pi/180
fth1_d = np.array([[0],[ np.sin(th)],[np.cos(th)]]); fth1_d /= np.linalg.norm(fth1_d)
fth2_d = np.array([[0],[ np.sin(th)],[np.cos(th)]]); fth2_d /= np.linalg.norm(fth2_d)
fth3_d = np.array([[0],[-np.sin(th)],[np.cos(th)]]); fth3_d /= np.linalg.norm(fth3_d)
fth4_d = np.array([[0],[-np.sin(th)],[np.cos(th)]]); fth4_d /= np.linalg.norm(fth4_d)

dpv = np.array([J1,J2,J3,Jxy,Jxz,Jyz,
    *rth1_d[:,0], *rth2_d[:,0], *rth3_d[:,0], *rth4_d[:,0],
    *fth1_d[:,0], *fth2_d[:,0], *fth3_d[:,0], *fth4_d[:,0]])

omegab_max = 1.5*np.pi/180
Fth_max = 25e-6

start_quat = np.array([[1.0],[0.0],[0.0],[0.0]])
end_quat   = np.array([[0.0],[1.0],[0.0],[0.0]])
start_hb = Jsat @ np.array([[0.2*np.pi/180],[-0.1*np.pi/180],[-0.2*np.pi/180]])
end_hb   = Jsat @ np.zeros((3,1))
start_pos = np.block([[start_quat],[start_hb]])
end_pos   = np.block([[end_quat],[end_hb]])

# Notebook scaling: note c_u is POSITIVE there.
S_u = np.diag([2*Fth_max]*4)
c_u = np.array([[Fth_max]]*4)
S_x = np.diag([2,2,2,2, 2*J1*omegab_max, 2*J2*omegab_max, 2*J3*omegab_max])
c_x = np.array([[-1],[-1],[-1],[-1],
                [-J1*omegab_max],[-J2*omegab_max],[-J3*omegab_max]])

quat = slerp(start_pos[0:4,0], end_pos[0:4,0], T+1).T
omega = angular_velocity(quat.T, tau_abs, scalar_first=False).T   # notebook convention
warm_x = np.zeros((7, T+1))
warm_x[0:4,:] = quat
warm_x[4:7,:] = Jsat @ omega
warm_u = np.zeros((4, T))

inv_Jb_S = np.linalg.inv(Jsat) @ S_x[4:7,4:7]
inv_Jb_c = np.linalg.inv(Jsat) @ c_x[4:7,0:1]

res = solve_trajectory(
    states=states, controls=controls, dynamics=dynamics,
    start=start_pos, end=end_pos, T=T,
    time_mode='ZOH+dilation',
    tf_init=tf_init, tf_min=tf_min, tf_max=tf_max, alpha_p=0.05,
    dynamic_parameters_sym=dps, dynamic_parameters_val=dpv,
    S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
    warm_start_x=warm_x, warm_start_u=warm_u,
    cost_terms=[{'kind':'norm1','var':'u','slice':slice(0,4),
                 'coeff':'tau','weight':10.0}],
    state_bounds=[{'kind':'norm','slice':slice(4,7),'norm_type':'norm2',
                   'limit':omegab_max,'S':inv_Jb_S,'c':inv_Jb_c}],
    control_bounds=[
        {'kind':'box','slice':slice(0,4),'lower':0.0,'upper':Fth_max},
        {'kind':'box_affine','c_affine':np.array([[1.,1.,0.,0.]]),'upper':Fth_max},
        {'kind':'box_affine','c_affine':np.array([[0.,0.,1.,1.]]),'upper':Fth_max},
    ],
    size_N=10, lamb=5000.0, max_iter=15,
    rho0=0.0, rho1=0.1, rho2=0.7, eta0=0.001, eta1=10.0, eta_init=1.0,
    beta_sh=2.0, beta_gr=2.0, e_tol=0.005, epsilon_stop=1e-5,
    solver='ECOS', solver_options={'max_iters': 200},
    # Bit-exact reproduction needs the legacy trust region, which also
    # charges eta for the boundary-pinned states x_0 and x_T.
    pin_aware_trust_region=False,
)

ref = [(34.9272, 302.1056, 707.6715, 298.6339),
       (49.3521, 251.4258, 302.1056, 585.6905),
       (34.0228, 155.0715, 251.4258, 410.7513),
       (35.3089,  43.7126, 155.0715, 445.8738)]
print("\n  it |     L (mine)   L(nb)+0.5 |     J (mine)      J (nb) |    tf (mine)    tf (nb)")
for i, (L, Jv, oJ, tfv) in enumerate(ref):
    h = res['history'][i]
    print(f"  {i+1:2d} | {h['L_SCVx']:12.4f} {L+0.5:11.4f} |"
          f" {h['J_SCVx']:12.4f} {Jv:11.4f} | {h['tf']:12.4f} {tfv:11.4f}")
print(f"\nfinal tf = {res['tf']:.4f} s   (notebook converges to ~443.9 s)")
