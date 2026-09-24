"""
Quaternion helpers for building warm starts on attitude problems.

Convention: quaternions are (q0, q1, q2, q3) with q0 the scalar part.
"""

import numpy as np
from scipy.spatial.transform import Rotation


def _to_scipy(quats):
    """Scalar-first (q0,q1,q2,q3) -> scipy's scalar-last (x,y,z,w)."""
    q = np.asarray(quats, dtype=float)
    return np.column_stack([q[:, 1], q[:, 2], q[:, 3], q[:, 0]])


def slerp(q1, q2, n):
    """
    Spherical linear interpolation between two unit quaternions.

    Unlike component-wise linear interpolation this keeps the norm at 1 all
    the way along, which matters when the result seeds a solver that has a
    unit-norm quaternion in its dynamics.

    Parameters
    ----------
    q1, q2 : array_like, shape (4,)
        Unit quaternions (q0, q1, q2, q3).
    n : int
        Number of samples, endpoints included.

    Returns
    -------
    ndarray, shape (n, 4)
    """
    q1 = np.asarray(q1, dtype=float)
    q2 = np.asarray(q2, dtype=float)

    dot = float(np.clip(np.dot(q1, q2), -1.0, 1.0))
    if dot < 0.0:                       # take the short way round
        q2, dot = -q2, -dot

    theta_0 = np.arccos(dot)
    if abs(theta_0) < 1e-6:
        return np.linspace(q1, q2, n)

    sin_theta_0 = np.sin(theta_0)
    out = []
    for i in range(n):
        t = i / (n - 1)
        theta = theta_0 * t
        s0 = np.cos(theta) - dot * np.sin(theta) / sin_theta_0
        s1 = np.sin(theta) / sin_theta_0
        out.append(s0 * q1 + s1 * q2)
    return np.array(out)


def angular_velocity(quats, dt, scalar_first=True):
    """
    Angular velocity of a quaternion sequence, by finite differences on SO(3):

        omega_k = log(q_k+1 q_k^-1) / dt

    Parameters
    ----------
    quats : array_like, shape (N, 4)
    dt : float
        Time step between consecutive quaternions.  For a time-dilation
        warm start use the *absolute* step tf_init/T -- angular velocity is a
        physical quantity and does not live in normalized time.
    scalar_first : bool
        True if `quats` is (q0, q1, q2, q3) with q0 scalar, which is this
        library's convention.

    Returns
    -------
    ndarray, shape (N, 3).  The first row is zero: there is no preceding
    quaternion to difference against.
    """
    q = np.asarray(quats, dtype=float)
    rotations = Rotation.from_quat(_to_scipy(q) if scalar_first else q)

    omega = [[0.0, 0.0, 0.0]]
    for i in range(len(rotations) - 1):
        delta = rotations[i + 1] * rotations[i].inv()
        omega.append(delta.as_rotvec() / dt)
    return np.array(omega)


# Backwards-compatible alias for the older name.
angular_vel = angular_velocity


def quat_from_axis_angle(axis, angle_deg):
    """Unit quaternion (q0,q1,q2,q3) for a rotation about `axis` by `angle_deg`."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    half = np.deg2rad(angle_deg) / 2.0
    q = np.array([[np.cos(half)],
                  [axis[0] * np.sin(half)],
                  [axis[1] * np.sin(half)],
                  [axis[2] * np.sin(half)]])
    return q / np.linalg.norm(q)
