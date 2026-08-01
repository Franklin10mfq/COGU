"""
Helpers de cuaterniones para construir warm starts en problemas
SCVx con dinamica de actitud (cuerpo rigido, SO(3)).

slerp:        interpolacion esferica entre dos cuaterniones unitarios.
angular_vel:  velocidad angular omega via diferencias finitas sobre SO(3).

Convencion: cuaterniones (q0, q1, q2, q3) con q0 = parte escalar.
"""

import numpy as np
from scipy.spatial.transform import Rotation as R


def slerp(q1, q2, n):
    """
    Interpolacion lineal esferica entre dos cuaterniones unitarios.

    Mantiene la norma unitaria a lo largo del trayecto (lo que la
    interpolacion lineal componente a componente NO garantiza).

    Parameters
    ----------
    q1 : array_like, shape (4,)
        Cuaternion inicial (q0, q1, q2, q3), norma 1.
    q2 : array_like, shape (4,)
        Cuaternion final, norma 1.
    n : int
        Numero de puntos a generar (incluyendo extremos).

    Returns
    -------
    ndarray, shape (n, 4)
        Trayectoria de cuaterniones unitarios entre q1 y q2.
    """
    dot = np.clip(np.dot(q1, q2), -1.0, 1.0)
    if dot < 0:
        q2, dot = -q2, -dot
    theta = np.arccos(dot)
    if abs(theta) < 1e-6:
        return np.linspace(q1, q2, n)
    return np.array([
        (np.cos(theta * t / (n - 1)) - dot * np.sin(theta * t / (n - 1)) / np.sin(theta)) * q1
        + np.sin(theta * t / (n - 1)) / np.sin(theta) * q2
        for t in range(n)
    ])


def angular_vel(quats, dt):
    """
    Velocidad angular omega de una secuencia de cuaterniones.

    Se calcula via la rotacion relativa entre instantes consecutivos:
        omega_k = log(q_k+1 * q_k^{-1}) / dt

    Parameters
    ----------
    quats : array_like, shape (N, 4)
        Secuencia de cuaterniones unitarios.
    dt : float
        Paso de tiempo entre cuaterniones consecutivos.

    Returns
    -------
    ndarray, shape (N, 3)
        Velocidades angulares. La primera fila es cero (sin
        cuaternion previo del cual derivar).
    """
    rots = R.from_quat(quats)
    w = [[0, 0, 0]]
    for i in range(len(rots) - 1):
        w.append((rots[i + 1] * rots[i].inv()).as_rotvec() / dt)
    return np.array(w)
