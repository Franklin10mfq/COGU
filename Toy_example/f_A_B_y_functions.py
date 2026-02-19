
import numpy as np
import math

def f_SCVx(hx, hu, t, dyn_par):
    return np.array([[-1 + hx[1,0]],[hu[0,0] + 1]])

def A_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return np.array([[0,1],[0,0]])

def B_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return np.array([[0],[1]])

def y_no_discrete_no_scaled_SCVx(hx, hu, t, dyn_par):
    return np.array([[-1],[1]])
