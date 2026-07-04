import numpy as np
import ctypes
import time

def mRMM_function(t,mRMx,mRMy,mRMz):
    return np.array([[mRMx],[mRMy],[mRMz]])
# Global precomputed matrix (must exist)
# ECI_B_precomputed: shape (3, N)

def ECI_B_function(t, N_B_precomputed, size_N_B_precomputed, tf):
    """
    Linear interpolation of ECI_B (3x1) from precomputed table.

    Inputs:
        t  : query time (scalar)
        N_B_precomputed : array (3, N)
        size_N_B_precomputed : N
        tf : final time

    Output:
        ECI_B : list of size 3
    """

    # Precompute dt (uniform grid)
    dt = tf / (size_N_B_precomputed - 1)
    inv_dt = 1.0 / dt

    # Compute index (O(1))
    idx = int(t * inv_dt)

    # Clamp
    if idx < 0:
        return N_B_precomputed[0:3, 0:1]
    if idx >= size_N_B_precomputed - 1:
        return N_B_precomputed[0:3, size_N_B_precomputed - 1:size_N_B_precomputed]

    # Interpolation factor
    t_idx = idx * dt
    alpha = (t - t_idx) * inv_dt

    # Manual interpolation (faster for 3 elements)
    return np.array([
        [(1.0 - alpha) * N_B_precomputed[0, idx] + alpha * N_B_precomputed[0, idx + 1]],
        [(1.0 - alpha) * N_B_precomputed[1, idx] + alpha * N_B_precomputed[1, idx + 1]],
        [(1.0 - alpha) * N_B_precomputed[2, idx] + alpha * N_B_precomputed[2, idx + 1]]])

def precompute_B_wmm(v0_double,r0_double,JD_double,dt_wmm_double,tf):
    lib = ctypes.CDLL("./WMM_orb_prop_C/WMM_orb_prop_C_caller.dll")

    # firma de la función
    lib.WMM_ECI_B.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # v (3)
        ctypes.POINTER(ctypes.c_double),  # r (3)
        ctypes.c_double,                  # JD
        ctypes.c_double,                  # dt
        ctypes.c_int,                     # T
        ctypes.POINTER(ctypes.c_double)   # salida (flatten)
    ]
    lib.WMM_ECI_B.restype = None

    B = np.zeros(3 * 2000, dtype=np.float64)

    # llamada
    t0_auxaux_=time.time()
    lib.WMM_ECI_B(
        v0_double.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        r0_double.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_double(JD_double),
        ctypes.c_double(dt_wmm_double),
        ctypes.c_int(int(tf / dt_wmm_double) + 1),
        B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    tf_auxaux_=time.time()
    print("Precomputed WMM, inference time [ms]",(tf_auxaux_-t0_auxaux_)*1000)
    # reshape a [3 x 500]
    B = B.reshape((3, 2000))
    B_final = B[:,0:int(tf / dt_wmm_double) + 1]
    return B_final   
