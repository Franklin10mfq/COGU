"""MPC-RTI para control magnético de actitud de un CubeSat.

ECOS obligatorio; horizonte ZOH de 30 intervalos de 10 s; actualización
receding-horizon cada 1 s; RK4 y sensibilidades variacionales consistentes.
"""
import os
import time
import numpy as np
import cvxpy as cp
import sympy as sp
import matplotlib

if os.name == "nt" or os.environ.get("DISPLAY") or os.environ.get("MPLBACKEND"):
    try:
        matplotlib.use(os.environ.get("MPLBACKEND", "TkAgg"))
    except Exception:
        matplotlib.use("Agg")
else:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==============================================================================
# 0. PARÁMETROS
# ==============================================================================
TS_CTRL = 1.0
TS_NODE = 10.0
NP = 30
SIZE_N = 20
SIM_TIME = 300.0
T_START = 0.0
M_MAX = 0.2                         # [A m^2]
OMEGA_MAX = np.deg2rad(5.0)         # [rad/s]

Q_X = 10.0 * np.diag([10., 10., 10., 10., 5., 5., 5.])

Q_X  = 20*np.diag([20, 20, 20, 20, 1, 1, 1])
R_U = 0.07 * np.eye(3)
R_DU = 15 * np.eye(3)
Q_F = 10.0 * Q_X

W_OMEGA_SLACK = 1.0e3
TERMINAL_EQUALITY = False
ALIGN_QUAT_SIGN = True

ECOS_OPTS = dict(max_iters=200, abstol=1e-5, reltol=1e-5, feastol=1e-5,
                 abstol_inacc=1e-3, reltol_inacc=1e-3, feastol_inacc=1e-3)
STATES, INPUTS = 7, 3
assert TS_CTRL > 0 and TS_NODE > 0 and SIZE_N > 0 and NP > 0
assert TS_CTRL <= TS_NODE


def sqrt_psd(matrix):
    eigval, eigvec = np.linalg.eigh(0.5 * (matrix + matrix.T))
    return (eigvec * np.sqrt(np.clip(eigval, 0.0, None))) @ eigvec.T


Q_H, R_H, RD_H, QF_H = map(sqrt_psd, (Q_X, R_U, R_DU, Q_F))

# ==============================================================================
# 1. MODELO FÍSICO Y DINÁMICA SIMBÓLICA
# ==============================================================================
J = np.array([[0.0022, 0.00002, 0.00002],
              [0.00002, 0.0021, 0.00002],
              [0.00002, 0.00002, 0.0019]], dtype=float)
INV_J = np.linalg.inv(J)
IJ_SYM = (INV_J[0, 0], INV_J[1, 1], INV_J[2, 2],
          INV_J[0, 1], INV_J[0, 2], INV_J[1, 2])

R_E = 6378e3
MU_EARTH = 3.986004418e14
A_ORB = R_E + 550e3
N_ORB = np.sqrt(MU_EARTH / A_ORB**3)
INC = np.deg2rad(97.4)
M_EARTH = 7.94e22
MU0_OVER_4PI = 1e-7
M_E_VEC = M_EARTH * np.array([0., 0., 1.])


def build_symbolic_model():
    mx, my, mz = sp.symbols("mx my mz")
    q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
    hbx, hby, hbz = sp.symbols("hbx hby hbz")
    ixx, iyy, izz, ixy, ixz, iyz = sp.symbols(
        "invIxx invIyy invIzz invIxy invIxz invIyz")
    bx, by, bz = sp.symbols("Bx By Bz")
    xs, us = [q0, q1, q2, q3, hbx, hby, hbz], [mx, my, mz]
    ps = [bx, by, bz, ixx, iyy, izz, ixy, ixz, iyz]
    inv_i = sp.Matrix([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])
    hb, moment = sp.Matrix([hbx, hby, hbz]), sp.Matrix([mx, my, mz])
    rot = sp.Matrix([
        [1-2*(q2**2+q3**2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1**2+q3**2), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1**2+q2**2)]])
    b_body = rot.T @ sp.Matrix([bx, by, bz])
    omega = inv_i @ hb
    om = sp.Matrix([
        [0, -omega[0], -omega[1], -omega[2]],
        [omega[0], 0, -omega[2], omega[1]],
        [omega[1], omega[2], 0, -omega[0]],
        [omega[2], -omega[1], omega[0], 0]
    ])

    qdot = sp.Rational(1, 2) * om @ sp.Matrix([q0, q1, q2, q3])

    hdot = -omega.cross(hb) + moment.cross(b_body)
    f = sp.Matrix(list(qdot) + list(hdot))
    return sp.lambdify(xs + us + ps,
                       (f, f.jacobian(sp.Matrix(xs)), f.jacobian(sp.Matrix(us))),
                       modules="numpy", cse=True)


print("Construyendo modelo simbólico...")
MODEL = build_symbolic_model()


def dynamics_jacobians(x, u, b_eci):
    args = (list(np.asarray(x, float)) + list(np.asarray(u, float)) +
            list(np.asarray(b_eci, float)) + list(IJ_SYM))
    f, a, b = MODEL(*args)
    return (np.asarray(f, float).reshape(7),
            np.asarray(a, float).reshape(7, 7),
            np.asarray(b, float).reshape(7, 3))

# ==============================================================================
# 2. CAMPO MAGNÉTICO, CUATERNIONES Y ESCALADO
# ==============================================================================
def magnetic_field_eci(t):
    t = np.asarray(t, dtype=float)
    nu = N_ORB * t
    r = A_ORB * np.stack([np.cos(nu), np.sin(nu)*np.cos(INC),
                          np.sin(nu)*np.sin(INC)], axis=-1)
    rn = np.linalg.norm(r, axis=-1, keepdims=True)
    rh = r / rn
    dot = rh @ M_E_VEC
    return MU0_OVER_4PI * (3.0*dot[..., None]*rh - M_E_VEC) / rn**3


def quat_to_rotation(q):
    q0, q1, q2, q3 = np.asarray(q, float)
    return np.array([[1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
                     [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
                     [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]])


def normalize_state_quaternion(x):
    x = np.asarray(x, dtype=float).copy()
    if x.ndim == 1:
        n = np.linalg.norm(x[:4])
        if n < 1e-12:
            raise ValueError("Norma de cuaternión prácticamente nula")
        x[:4] /= n
    else:
        n = np.linalg.norm(x[:4, ...], axis=0, keepdims=True)
        if np.any(n < 1e-12):
            raise ValueError("Norma de cuaternión prácticamente nula")
        x[:4, ...] /= n
    return x


def attitude_error_deg(q, qref):
    return 2.0 * np.degrees(np.arccos(np.clip(abs(float(q @ qref)), 0., 1.)))


START_QUAT = np.array([0.0, 0.0, 0., 1.])
END_QUAT = np.array([0., 0., 1., 0.])
START_OMEGA = np.deg2rad(np.array([-0.2, 0.2, 0.2]))
END_OMEGA = np.zeros(3)
S_U = M_MAX * np.eye(3)
INV_S_U = np.linalg.inv(S_U)
H_MAX = np.diag(J) * OMEGA_MAX
S_X = np.diag([1., 1., 1., 1., H_MAX[0], H_MAX[1], H_MAX[2]])
INV_S_X = np.linalg.inv(S_X)
W_OMEGA = INV_J @ S_X[4:7, 4:7]


def scale_x(x):
    return INV_S_X @ np.asarray(x, dtype=float)


def unscale_u(u):
    return S_U @ np.asarray(u, dtype=float)

# ==============================================================================
# 3. INTEGRACIÓN RK4 Y SENSIBILIDADES
# ==============================================================================
def rk4_state_step(x, u, t, dt):
    b0, bm, b1 = magnetic_field_eci(np.array([t, t+0.5*dt, t+dt]))
    k1 = dynamics_jacobians(x, u, b0)[0]
    k2 = dynamics_jacobians(x+0.5*dt*k1, u, bm)[0]
    k3 = dynamics_jacobians(x+0.5*dt*k2, u, bm)[0]
    k4 = dynamics_jacobians(x+dt*k3, u, b1)[0]
    return x + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)


def integrate_zoh(x0, u, t0, duration, substeps):
    x = np.asarray(x0, dtype=float).copy()
    dt = duration / substeps
    for j in range(substeps):
        x = rk4_state_step(x, u, t0+j*dt, dt)
    return x


def augmented_rhs(y, u, t):
    x = y[:7]
    phi = y[7:56].reshape(7, 7)
    gamma = y[56:].reshape(7, 3)
    f, a, b = dynamics_jacobians(x, u, magnetic_field_eci(t))
    return np.concatenate([f, (a @ phi).ravel(), (a @ gamma + b).ravel()])


def integrate_interval_sensitivities(x0, u0, t0):
    y = np.concatenate([x0, np.eye(7).ravel(), np.zeros((7, 3)).ravel()])
    dt = TS_NODE / SIZE_N
    for j in range(SIZE_N):
        tj = t0 + j*dt
        k1 = augmented_rhs(y, u0, tj)
        k2 = augmented_rhs(y+0.5*dt*k1, u0, tj+0.5*dt)
        k3 = augmented_rhs(y+0.5*dt*k2, u0, tj+0.5*dt)
        k4 = augmented_rhs(y+dt*k3, u0, tj+dt)
        y += (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    xnext = y[:7]
    ad = y[7:56].reshape(7, 7)
    bd = y[56:].reshape(7, 3)
    cd = xnext - ad @ x0 - bd @ u0
    return xnext, ad, bd, cd


def rollout_nominal(x0, u_plan_phys, t0):
    xnom = np.zeros((7, NP+1))
    xnom[:, 0] = x0
    for k in range(NP):
        xnom[:, k+1] = integrate_zoh(
            xnom[:, k], u_plan_phys[:, k], t0+k*TS_NODE, TS_NODE, SIZE_N)
    return xnom


def linearize_horizon(xnom, u_plan_phys, t0):
    ad_s = np.zeros((NP, 7, 7))
    bd_s = np.zeros((NP, 7, 3))
    cd_s = np.zeros((NP, 7))
    for k in range(NP):
        _, ad, bd, cd = integrate_interval_sensitivities(
            xnom[:, k], u_plan_phys[:, k], t0+k*TS_NODE)
        ad_s[k] = INV_S_X @ ad @ S_X
        bd_s[k] = INV_S_X @ bd @ S_U
        cd_s[k] = INV_S_X @ cd
    return ad_s, bd_s, cd_s

# ==============================================================================
# 4. OPTIMIZACIÓN CONVEXA (SOCP) Y VALIDACIÓN
# ==============================================================================
def build_mpc_problem():
    # Variables de optimización
    x = cp.Variable((7, NP + 1), name="x_scaled")
    u = cp.Variable((3, NP), name="u_scaled")
    slack = cp.Variable(NP + 1, nonneg=True, name="omega_slack")

    # Parámetros del modelo discreto linealizado
    ad_p = cp.Parameter((7, 7 * NP), name="Ad")
    bd_p = cp.Parameter((7, 3 * NP), name="Bd")
    cd_p = cp.Parameter((7, NP), name="cd")
    x0_p = cp.Parameter(7, name="x0")
    xf_p = cp.Parameter(7, name="xf")
    uprev_p = cp.Parameter(3, name="u_previous_scaled")

    # Función de costo
    cost = 0
    for k in range(NP):
        e_k = x[:, k] - xf_p
        cost += TS_NODE * cp.sum_squares(Q_H @ e_k)
        cost += TS_NODE * cp.sum_squares(R_H @ u[:, k])

        du = u[:, k] - (uprev_p if k == 0 else u[:, k - 1])
        cost += (1.0 / TS_NODE) * cp.sum_squares(RD_H @ du)

    # Costo terminal y holgura
    e_terminal = x[:, NP] - xf_p
    cost += cp.sum_squares(QF_H @ e_terminal)
    cost += W_OMEGA_SLACK * TS_NODE * cp.sum(slack)

    # Restricciones
    cons = [x[:, 0] == x0_p]
    if TERMINAL_EQUALITY:
        cons.append(x[:, NP] == xf_p)

    for k in range(NP):
        ak = ad_p[:, 7*k:7*(k+1)]
        bk = bd_p[:, 3*k:3*(k+1)]
        cons.append(x[:, k+1] == ak @ x[:, k] + bk @ u[:, k] + cd_p[:, k])

    cons += [u <= 1.0, u >= -1.0]

    omega_scaled = W_OMEGA @ x[4:7, :]
    cons.append(cp.norm(omega_scaled, 2, axis=0) <= OMEGA_MAX + slack)

    problem = cp.Problem(cp.Minimize(cost), cons)
    return problem, x, u, slack, ad_p, bd_p, cd_p, x0_p, xf_p, uprev_p


def validate_solution(xsol, usol, ssol, ad, bd, cd, tol=4e-2):
    if xsol is None or usol is None or ssol is None:
        return False, "Variables de solución ausentes"
    if not (np.isfinite(xsol).all() and np.isfinite(usol).all() and np.isfinite(ssol).all()):
        return False, "Solución con NaN/Inf"
    if xsol.shape != (7, NP+1) or usol.shape != (3, NP) or ssol.shape != (NP+1,):
        return False, "Dimensiones inesperadas en la solución"
    if np.max(np.abs(usol)) > 1.0+1e-4:
        return False, f"Límite de control excedido ({np.max(np.abs(usol)):.6f})"
    if np.min(ssol) < -1e-6:
        return False, "Holgura negativa"
    max_res = 0.0
    for k in range(NP):
        pred = ad[k] @ xsol[:, k] + bd[k] @ usol[:, k] + cd[k]
        max_res = max(max_res, float(np.linalg.norm(xsol[:, k+1]-pred)))
    if max_res > tol:
        return False, f"Residuo dinámico {max_res:.3e} > {tol:.3e}"
    omega = np.linalg.norm(W_OMEGA @ xsol[4:7, :], axis=0)
    viol = float(np.max(omega-OMEGA_MAX-ssol))
    if viol > tol:
        return False, f"Restricción de velocidad incumplida ({viol:.3e})"
    return True, f"OK; residuo dinámico={max_res:.2e}"

# ==============================================================================
# 5. INICIALIZACIÓN Y DESPLAZAMIENTO DEL PLAN
# ==============================================================================
def initial_state_and_reference():
    x0 = np.r_[START_QUAT/np.linalg.norm(START_QUAT), J @ START_OMEGA]
    xf = np.r_[END_QUAT/np.linalg.norm(END_QUAT), J @ END_OMEGA]
    return normalize_state_quaternion(x0), normalize_state_quaternion(xf)


def shift_control_plan(uplan, fraction):
    shifted = np.empty_like(uplan)
    shifted[:, :-1] = (1.0-fraction)*uplan[:, :-1] + fraction*uplan[:, 1:]
    shifted[:, -1] = uplan[:, -1]
    return shifted

# ==============================================================================
# 6. SIMULACIÓN EN LAZO CERRADO
# ==============================================================================
def main(sim_time=SIM_TIME, show=True):
    if "ECOS" not in cp.installed_solvers():
        raise RuntimeError("ECOS no está instalado. Ejecuta: pip install ecos")

    print("Iniciando MPC-RTI magnético (ZOH + RK4 + ECOS)...")
    (problem, xvar, uvar, svar, adp, bdp, cdp,
     x0p, xfp, uprevp) = build_mpc_problem()
    xplant, xref = initial_state_and_reference()

    u_plan_scaled = np.zeros((3, NP))
    u_previous_safe = np.zeros(3)

    tlog, xlog, ulog, errlog = [], [], [], []
    omegalog, blog = [], []
    solve_times, lin_times, slack_log, qnorm_log, umag_log = [], [], [], [], []
    statuses, fallback_steps = [], []
    t_sim = T_START
    n_steps = int(round(sim_time/TS_CTRL))
    wall0 = time.perf_counter()

    plant_substeps = max(1, int(round(TS_CTRL/(TS_NODE/SIZE_N))))

    for step in range(n_steps):
        if ALIGN_QUAT_SIGN and np.dot(xplant[:4], xref[:4]) < 0:
            xref[:4] *= -1.0

        u_plan_scaled = shift_control_plan(u_plan_scaled, TS_CTRL/TS_NODE)
        u_plan_phys = S_U @ u_plan_scaled

        solve_t = 0.0
        lin_t = 0.0
        u_phys = u_previous_safe.copy()
        status = "LINEARIZATION_FAILED"
        solve_ok = False
        slack_max = np.nan

        try:
            t0_lin = time.perf_counter()
            xnom = rollout_nominal(xplant, u_plan_phys, t_sim)
            ad_val, bd_val, cd_val = linearize_horizon(xnom, u_plan_phys, t_sim)
            lin_t = time.perf_counter() - t0_lin
            if not (np.isfinite(ad_val).all() and np.isfinite(bd_val).all() and np.isfinite(cd_val).all()):
                raise FloatingPointError("Matrices de linealización no finitas")

            adp.value = np.hstack([ad_val[k] for k in range(NP)])
            bdp.value = np.hstack([bd_val[k] for k in range(NP)])
            cdp.value = cd_val.T
            x0p.value = scale_x(xplant)
            xfp.value = scale_x(xref)
            uprevp.value = INV_S_U @ u_previous_safe

            t0_solve = time.perf_counter()
            problem.solve(solver=cp.ECOS, warm_start=True, **ECOS_OPTS)
            solve_t = time.perf_counter() - t0_solve
            status = str(problem.status)

            if status in ("optimal", "optimal_inaccurate"):
                xsol = None if xvar.value is None else np.asarray(xvar.value, float)
                usol = None if uvar.value is None else np.asarray(uvar.value, float)
                ssol = None if svar.value is None else np.asarray(svar.value, float)
                valid, message = validate_solution(xsol, usol, ssol, ad_val, bd_val, cd_val)
                if valid:
                    u_plan_scaled = usol.copy()
                    u_phys = S_U @ usol[:, 0]
                    u_phys = np.clip(u_phys, -M_MAX, M_MAX)
                    u_previous_safe = u_phys.copy()
                    solve_ok = True
                    slack_max = float(np.max(np.maximum(ssol, 0.0)))
                else:
                    status = "REJECTED: " + message
                    fallback_steps.append(step)
                    print(f"[!] t={t_sim:.1f}s: solución rechazada: {message}")
            else:
                fallback_steps.append(step)
                print(f"[!] t={t_sim:.1f}s: ECOS status={status}; se usa fallback.")

        except Exception as exc:
            status = f"ERROR: {type(exc).__name__}"
            fallback_steps.append(step)
            print(f"[!] t={t_sim:.1f}s: {exc}; se usa control seguro.")

        omega = INV_J @ xplant[4:7]
        b_body = quat_to_rotation(xplant[:4]).T @ magnetic_field_eci(t_sim)
        err = attitude_error_deg(xplant[:4], xref[:4])
        tlog.append(t_sim)
        xlog.append(xplant.copy())
        ulog.append(u_phys.copy())
        errlog.append(err)
        omegalog.append(omega.copy())
        blog.append(b_body.copy())
        solve_times.append(solve_t)
        lin_times.append(lin_t)
        slack_log.append(slack_max)
        qnorm_log.append(float(np.linalg.norm(xplant[:4]) - 1.0))
        umag_log.append(float(np.linalg.norm(u_phys)))
        statuses.append(status)

        if step % 10 == 0:
            wpct = 100.0 * np.linalg.norm(omega) / OMEGA_MAX
            upct = 100.0 * np.max(np.abs(u_phys)) / M_MAX
            q_norm = float(np.linalg.norm(xplant[:4]))
            slack_str = f"slk={slack_max:.2e}" if not np.isnan(slack_max) else "slk=---"
            fb_str = " [FB]" if not solve_ok else ""

            print(f"t={t_sim:7.1f}s | err={err:7.3f}° | q={q_norm:.6f} | "
                  f"|u|={np.linalg.norm(u_phys):.4f} ({upct:5.1f}%) | "
                  f"|w|={np.degrees(np.linalg.norm(omega)):6.3f}°/s ({wpct:5.1f}%) | "
                  f"{status[:18]:18s}{fb_str:5s} | {slack_str:10s} | lin={lin_t:.3f}s ECOS={solve_t:.3f}s")

        xplant = integrate_zoh(xplant, u_phys, t_sim, TS_CTRL, plant_substeps)
        xplant = normalize_state_quaternion(xplant)
        t_sim += TS_CTRL

    # ==============================================================================
    # 7. RESUMEN Y GRÁFICAS
    # ==============================================================================
    xarr, uarr = np.asarray(xlog), np.asarray(ulog)
    omegaarr, barr = np.asarray(omegalog), np.asarray(blog)
    n = len(tlog)
    tail = max(1, n // 5)

    max_q_dev = float(np.max(np.abs(np.linalg.norm(xarr[:, :4], axis=1) - 1.0)))
    max_slack = float(np.nanmax(slack_log)) if not np.all(np.isnan(slack_log)) else 0.0

    print("\n" + "=" * 105)
    print(f"RESUMEN MPC-RTI: {n} pasos de control")
    print(f"Fallos/rechazos (Fallback) : {len(fallback_steps)} ({100.0 * len(fallback_steps) / n:.1f}%)")
    print(f"Linealización [s]          : media={np.mean(lin_times):.4f}, máx={np.max(lin_times):.4f}")
    print(f"ECOS [s]                   : media={np.mean(solve_times):.4f}, máx={np.max(solve_times):.4f}")
    print(f"Desviación máx. norma q    : {max_q_dev:.3e} (Objetivo: < 1e-6)")
    print(f"Holgura (slack) máxima     : {max_slack:.3e} rad/s (>0 indica violación de límite de velocidad)")
    print(f"Error de actitud final     : {errlog[-1]:.3f}°")
    print(f"Error medio (último 20%)   : {np.mean(errlog[-tail:]):.3f}°")
    print(f"Velocidad media (último 20%): {np.degrees(np.mean(np.linalg.norm(omegaarr[-tail:], axis=1))):.4f}°/s")
    print(f"Tiempo total de simulación : {time.perf_counter() - wall0:.1f} s")
    print("=" * 105)

    fig, axs = plt.subplots(7, 1, figsize=(14, 21), sharex=True)
    axs[0].plot(tlog, errlog, "k-", lw=1.8)
    axs[0].set_ylabel("Error [deg]")
    axs[0].set_title("Control magnético de actitud MPC-RTI")
    axs[0].grid(True)

    for i, label in enumerate(["q0", "q1", "q2", "q3"]):
        axs[1].plot(tlog, xarr[:, i], label=label)
    axs[1].set_ylabel("Cuaternión")
    axs[1].legend(ncol=4)
    axs[1].grid(True)

    axs[2].plot(tlog, qnorm_log, "b-", lw=1.2)
    axs[2].set_ylabel("Desviación Norma q\n(||q|| - 1)")
    axs[2].grid(True)

    for i, label in enumerate([r"$\omega_x$", r"$\omega_y$", r"$\omega_z$"]):
        axs[3].plot(tlog, np.degrees(omegaarr[:, i]), label=label)
    axs[3].axhline(np.degrees(OMEGA_MAX), color="k", ls="--", alpha=0.5)
    axs[3].axhline(-np.degrees(OMEGA_MAX), color="k", ls="--", alpha=0.5)
    axs[3].set_ylabel("Velocidad [deg/s]")
    axs[3].legend()
    axs[3].grid(True)

    for i, label in enumerate(["mx", "my", "mz"]):
        axs[4].step(tlog, uarr[:, i], where="post", label=label)
    axs[4].axhline(M_MAX, color="k", ls="--", alpha=0.5)
    axs[4].axhline(-M_MAX, color="k", ls="--", alpha=0.5)
    axs[4].set_ylabel("Momento [A m²]")
    axs[4].legend()
    axs[4].grid(True)

    for i, label in enumerate(["Bx", "By", "Bz"]):
        axs[5].plot(tlog, barr[:, i]*1e6, label=label)
    axs[5].plot(tlog, np.linalg.norm(barr, axis=1)*1e6, "k--", label="|B|")
    axs[5].set_ylabel("Campo [μT]")
    axs[5].legend(ncol=4)
    axs[5].grid(True)

    axs[6].plot(tlog, np.nan_to_num(slack_log, nan=0.0), "m-")
    axs[6].set_ylabel("Holgura [rad/s]")
    axs[6].grid(True)

    plt.tight_layout()

    if matplotlib.get_backend().lower() == "agg":
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "resultados_mpc_rti_corregido.png")
        plt.savefig(out, dpi=130)
        print(f"Gráfica guardada en: {out}")
    elif show:
        plt.show()


if __name__ == "__main__":
    main()