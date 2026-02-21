import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time
import ecos
import math
from scipy.spatial.transform import Rotation as R

print("numpy version",np.__version__)
print("cvxpy version",cp.__version__)

T = 31-1 # 101-1 means 101 discretization points
tf = 200.0
tau = tf/(T)
size_N=20
print("Step: ",tau," [s]")

start_quat=np.array([[0],[1],[0],[0]])
#end_quat=np.array([[0.3**0.5],[0.4**0.5],[0.1**0.5],[0.2**0.5]])
end_quat=np.array([[0],[0],[0],[1]])

#r = R.from_euler('xyz', [45*np.pi/180, 25*np.pi/180, 115*np.pi/180])
#end_quat=np.reshape(np.array(r.as_quat()),(4,1))

print(end_quat)

startpos=np.block([[start_quat],[0],[0],[0]])
endpos=np.block([[end_quat],[0],[0],[0]])

J1_double=0.0333
J2_double=0.0333
J3_double=0.0067

u_max_torq_double = 0.001 #Aspina RW

omega_max_double = 1.0*np.pi/180 #deg/s

#Guidance parameters
rho0=0.0
rho1=0.1
rho2=0.7
etta0=0.001
etta1=5
beta_sh=2
beta_gr=2

lamb_double=2000
etta_double=0.1

e_tol=0.05
epsilon_stop_norm=0.04

u_qw_scaling=np.array([[1/(1*u_max_torq_double),0,0],[0,1/(1*u_max_torq_double),0],[0,0,1/(1*u_max_torq_double)]])

nx = cp.Variable((7, T + 1), name='nx')
u = cp.Variable((3, T), name='u')
vc = cp.Variable((7, T), name='vc')

startpos_cvxpy = cp.Parameter((7,1), name='start_pos')
endpos_cvxpy = cp.Parameter((7,1), name='end_pos')

ox_cvxpy = cp.Parameter((7,T + 1), name='ox_cvxpy')
ou_cvxpy = cp.Parameter((3,T), name='ou_cvxpy')

A_discrete_qw = cp.Parameter((7,7*T), name='A_discrete_qw')
B_discrete_qw_scaled = cp.Parameter((7,3*T), name='B_discrete_qw_scaled')
w_discrete_qw = cp.Parameter((7,T), name='w_discrete_qw')

lamb = cp.Parameter(name='lamb')
etta = cp.Parameter(name='etta')

vel_max = cp.Parameter(name='vel_max')
omega_max = cp.Parameter(name='omega_max')

constraints = [
    nx[:, 0] == startpos_cvxpy[:,0],
    nx[:, T] == endpos_cvxpy[:,0],
    #u[:,0]==0,u[:,-1]==0
]
cost = 0

for k in range(0, T): # from 0 to T-1

    constraints += [nx[0:7, k+1] == A_discrete_qw[:,7*k:7*k+7] @ nx[0:7, k] + B_discrete_qw_scaled[:,3*k:3*k+3] @ (u[0:3, k]) + w_discrete_qw[:,k] + vc[0:7, k]] #q1_k+1
    
    constraints += [cp.norm(u[0:3, k], 2) <= 1]

    cost += tau*cp.sum_squares(u[0:3,k])
    cost += tau*cp.norm(lamb*vc[:, k], 1)
    constraints  += [cp.norm(nx[:, k]-ox_cvxpy[:,k],'inf')+cp.norm(u[:, k]-ou_cvxpy[:,k],'inf')<=etta]

for k in range(0, T+1):
    constraints += [cp.norm(nx[4:7,k], 2)<=omega_max]

objective = cp.Minimize(cost)
problem = cp.Problem(objective, constraints)

print("Is DPP? ", problem.is_dcp(dpp=True))

aux_A_discrete_qw = np.zeros((7,7*T))
aux_B_discrete_qw_scaled = np.zeros((7,3*T))
aux_w_discrete_qw = np.zeros((7,T))

def scaling_begin(u,u_scaling,T):
    for k in range(0, T):
        u[:,k:k+1]=u_scaling@u[:,k:k+1]
    return u
def scaling_end(u,u_scaling,T):
    for k in range(0, T):
        u[:,k:k+1]=np.linalg.inv(u_scaling)@u[:,k:k+1]
    return u
    
def slerp(q1, q2, num_samples):
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

    rotations = R.from_quat(quaternions)
    angular_velocities = [[0,0,0]]

    for i in range(len(rotations) - 1):
        delta_rot = rotations[i + 1] * rotations[i].inv()
        log_rot = delta_rot.as_rotvec() / dt
        angular_velocities.append(log_rot)

    return np.array(angular_velocities)

def exp_matrix_taylor_A(A,h,n):
    size_n,size_aux=np.shape(A)
    sum=np.eye(size_n)+A*h
    for i in range(2,n+2,1):
        sum=sum+1/math.factorial(i)*np.linalg.matrix_power(h*A,i)
    return sum

def exp_matrix_taylor_B(A,B,h,n):
    size_n,size_aux=np.shape(A)
    sum=h*np.eye(size_n)
    for i in range(2,n+2,1):
        sum=sum+1/math.factorial(i)*np.linalg.matrix_power(A,i-1)*h**i
    return sum@B

def f_qw(x,u):
    aux_f=np.zeros((7,1))

    oq1=x[0,0]
    oq2=x[1,0]
    oq3=x[2,0]
    oq4=x[3,0]
    ow1=x[4,0]
    ow2=x[5,0]
    ow3=x[6,0]

    ou1=u[0,0]
    ou2=u[1,0]
    ou3=u[2,0]

    aux_f[0,0]=0.5*(oq4*ow1-oq3*ow2+oq2*ow3)
    aux_f[1,0]=0.5*(oq3*ow1+oq4*ow2-oq1*ow3)
    aux_f[2,0]=0.5*(-oq2*ow1+oq1*ow2+oq4*ow3)
    aux_f[3,0]=0.5*(-oq1*ow1-oq2*ow2-oq3*ow3)

    aux_f[4,0]=(1/J1_double)*(-(J3_double-J2_double)*ow2*ow3+ou1)
    aux_f[5,0]=(1/J2_double)*(-(J1_double-J3_double)*ow3*ow1+ou2)
    aux_f[6,0]=(1/J3_double)*(-(J2_double-J1_double)*ow1*ow2+ou3)
    return aux_f

def A_qw(oxqw):
  oq1=oxqw[0,0]
  oq2=oxqw[1,0]
  oq3=oxqw[2,0]
  oq4=oxqw[3,0]
  ow1=oxqw[4,0]
  ow2=oxqw[5,0]
  ow3=oxqw[6,0]
  aux_A_qw=np.zeros((7,7))
  aux_A_qw[0:4,:]=0.5*np.array([[0,ow3,-ow2,ow1,oq4,-oq3,oq2],
                            [-ow3,0,ow1,ow2,oq3,oq4,-oq1],
                            [ow2,-ow1,0,ow3,-oq2,oq1,oq4],
                            [-ow1,-ow2,-ow3,0,-oq1,-oq2,-oq3]])
  aux_A_qw[4:7,:] = np.array([
                     [0,0,0,0,0,1/J1_double*(-(J3_double-J2_double)*ow3),1/J1_double*(-(J3_double-J2_double)*ow2)],
                     [0,0,0,0,1/J2_double*(-(J1_double-J3_double)*ow3),0,1/J2_double*(-(J1_double-J3_double)*ow1)],
                     [0,0,0,0,1/J3_double*(-(J2_double-J1_double)*ow2),1/J3_double*(-(J2_double-J1_double)*ow1),0]])
  return aux_A_qw
def B_qw(oxqw):
  aux_B_qw = np.block([[np.zeros((4,3))],[1/J1_double,0,0],[0,1/J2_double,0],[0,0,1/J3_double]])
  return aux_B_qw
def w_qw(oxqw,ou):
  aux_w_qw=f_qw(oxqw,ou)-A_qw(oxqw)@oxqw
  return aux_w_qw

ox_quat = slerp(startpos[0:4,0], endpos[0:4,0],T+1).T
ox_angvel = compute_angular_velocity(ox_quat.T, tau).T

ox=np.zeros((7,T+1))
ox[0:4,:]=ox_quat
ox[4:7,:]=ox_angvel

ou=np.zeros((3,T))

omega_max.value=omega_max_double

def f_qw_rk4_step(xk, uk, dt):
    k1 = f_qw(xk, uk)
    k2 = f_qw(xk + 0.5 * dt * k1, uk)
    k3 = f_qw(xk + 0.5 * dt * k2, uk)
    k4 = f_qw(xk + dt * k3, uk)
    return xk + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def f_qw_rk5_step(xk, uk, dt):
    k1 = f_qw(xk, uk)
    k2 = f_qw(xk + (1/4) * dt * k1, uk)
    k3 = f_qw(xk + (3/8) * dt * k2, uk)
    k4 = f_qw(xk + (12/13) * dt * k3, uk)
    k5 = f_qw(xk + dt * k4, uk)
    k6 = f_qw(xk + (1/2) * dt * k5, uk)
    
    return xk + dt * (16/135 * k1 + 6656/12825 * k3 + 28561/56430 * k4 - 9/50 * k5 + 2/55 * k6)

def f_qw_Euler(xk, uk, dt):
    return xk + dt * f_qw(xk, uk)

def f_SCVx(xk,uk):
    xk1=np.zeros((7,1))
    xk_aux=np.copy(xk[0:7,0:1])
    for N in range(0,size_N,1):
        xk_aux=f_qw_rk5_step(xk_aux, uk[0:3,0:1], tau/size_N)
    
    xk1[0:7,0:1]=np.copy(xk_aux)

    return xk1

def J_SCVx(x,u,T):
    cost = 0
    for k in range(0, T): # from 0 to T-1
        cost += tau*np.linalg.norm(u[0:3,k], ord=2)**2
    print("Cost 1: ",cost)
    u[0:3,:] = scaling_end(u[0:3,:],u_qw_scaling,T) #scaling end

    for k in range(0, T): # from 0 to T-1
        cost += tau*np.linalg.norm(lamb_double*(x[:, k+1:k+2]-f_SCVx(x[:, k:k+1],u[:, k:k+1])), ord=1)

    u[0:3,:]=scaling_begin(u[0:3,:],u_qw_scaling,T) #scaling begin
    print("Cost 2: ",cost)
    return cost

def L_SCVx(x,u,vc,T):
    cost = 0
    for k in range(0, T): # from 0 to T-1
        cost += tau*np.linalg.norm(u[0:3,k], ord=2)**2
        cost += tau*np.linalg.norm(lamb_double*(vc[:, k]), ord=1)

    return cost

startpos_cvxpy.value=startpos
endpos_cvxpy.value=endpos

t0 = time.time()

lamb.value = lamb_double
etta.value = etta_double

ou[0:3,:]=scaling_begin(ou[0:3,:],u_qw_scaling,T) #scaling begin

ox_cvxpy.value=np.copy(ox) #trajectory initialization (solver, scaled)
ou_cvxpy.value=np.copy(ou) #trajectory initialization (solver, scaled)

i=1
no_first_iterations = False
while True:

    ou[0:3,:] = scaling_end(ou[0:3,:],u_qw_scaling,T) #scaling end

    for k in range(0, T):
        aux_A_discrete_qw[0:7,7*k:7*k+7] = exp_matrix_taylor_A(A_qw(ox[0:7,k:k+1]),tau,7)
        #aux_matrix_Bqw=np.block([[B_qw(ox[6:13,k:k+1])[0:4,:]],[B_qw(ox[6:13,k:k+1])[4:7,:]@np.linalg.inv(u_qw_scaling)]])
        aux_B_discrete_qw_scaled[0:7,3*k:3*k+3] = exp_matrix_taylor_B(A_qw(ox[0:7,k:k+1]),B_qw(ox[0:7,k:k+1])@np.linalg.inv(u_qw_scaling),tau,7)
        aux_w_discrete_qw[0:7,k:k+1] = exp_matrix_taylor_B(A_qw(ox[0:7,k:k+1]),w_qw(ox[0:7,k:k+1],np.zeros((3,1))),tau,7)

    A_discrete_qw.value = np.copy(aux_A_discrete_qw)
    B_discrete_qw_scaled.value = np.copy(aux_B_discrete_qw_scaled)
    w_discrete_qw.value = np.copy(aux_w_discrete_qw)

    val=problem.solve(solver="ECOS",ignore_dpp=True) #IF YOU USE CVXPYGEN: ignore_dpp=False


    vc_opt=np.copy(vc.value)
    x_opt=np.copy(nx.value)
    u_opt=np.copy(u.value)

    ou[0:3,:]=scaling_begin(ou[0:3,:],u_qw_scaling,T) #scaling begin

    oJ_SCVx=J_SCVx(ox,ou,T)
    J_SCVx_opt=J_SCVx(x_opt,u_opt,T)
    L_SCVx_opt=L_SCVx(x_opt,u_opt,vc_opt,T)
    

    Delta_J_SCVx=oJ_SCVx-J_SCVx_opt
    Delta_L_SCVx=oJ_SCVx-L_SCVx_opt
    print("oJ_SCVx: ",oJ_SCVx,"J_SCVx_opt",J_SCVx_opt,"L_SCVx_opt",L_SCVx_opt,"cvxpy_L",val,"Norm_x_diff: ",np.max(np.linalg.norm((x_opt-ox), ord=2,axis=0)))
    
    if (Delta_L_SCVx<e_tol*np.abs(oJ_SCVx) or np.max(np.linalg.norm((x_opt-ox), ord=1,axis=0))<epsilon_stop_norm) and no_first_iterations:
      
      ou[0:3,:] = scaling_end(ou[0:3,:],u_qw_scaling,T) #scaling end

      x_global=np.copy(ox)
      u_global=np.copy(ou)
      break;
    
    else:
      rho_i=Delta_J_SCVx/Delta_L_SCVx
      if rho_i<rho0:
        etta.value=max([etta0,etta.value/beta_sh])
        ox=np.copy(ox)
        ou=np.copy(ou)
      if rho_i>=rho0 and rho_i<rho1:
        etta.value=max([etta0,etta.value/beta_sh])
        ox=np.copy(x_opt)
        ou=np.copy(u_opt)
      if rho_i>=rho1 and rho_i<rho2:
        ox=np.copy(x_opt)
        ou=np.copy(u_opt)
      if rho_i>=rho2:
        etta.value=min([etta1,beta_gr*etta.value])
        ox=np.copy(x_opt)
        ou=np.copy(u_opt)

      print(" Iteration number: ",i," Cost function: ", val," Etta: ",etta.value, " Rho: ",rho_i)

    ox_cvxpy.value=np.copy(ox)
    ou_cvxpy.value=np.copy(ou)

    if i==3:
        no_first_iterations = True
    i=i+1
    
t1 = time.time()
print('\nCVXPY\nSolve time: %.3f ms' % (1000 * (t1 - t0)))

u_value=np.copy(u_global)
nx_value=np.copy(x_global)
#u_value = scaling_begin(u_value,u_qw_scaling,T)
#print(u_value)

print(np.linalg.norm(vc.value,2))