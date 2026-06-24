import numpy as np

def dcm_ZXZ(Omega, inc, theta):
    cO = np.cos(Omega)
    sO = np.sin(Omega)
    ci = np.cos(inc)
    si = np.sin(inc)
    ct = np.cos(theta)
    st = np.sin(theta)

    Rz1 = np.array([[cO, -sO, 0],
                    [sO,  cO, 0],
                    [0,   0,  1]])

    Rx  = np.array([[1, 0, 0],
                    [0, ci, -si],
                    [0, si,  ci]])

    Rz2 = np.array([[ct, -st, 0],
                    [st,  ct, 0],
                    [0,   0,  1]])

    return (Rz1 @ Rx @ Rz2).T

def ECI_B_function(t, inc, Omega, B0, n, t0):

    O_B = np.array([
            [-2*np.sin(inc)*np.cos(n*(t+t0))],
            [-np.sin(inc)*np.sin(n*(t+t0))],
            [np.cos(inc)]
        ]) * B0

    ON = dcm_ZXZ(Omega, inc, n*(t+t0))
    N_B = ON.T@O_B
    return N_B

def mRMM_function(t,mRMx,mRMy,mRMz):
    return np.array([[mRMx],[mRMy],[mRMz]])
    
    
