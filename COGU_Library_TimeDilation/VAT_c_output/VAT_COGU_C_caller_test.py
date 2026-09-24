"""
VAT_COGU_C_caller_test.py -- COGU feedforward + real-time disturbance
rejection + terminal-optimal quantized tracking, with the guidance itself
coming from repeated, on-demand calls to COGU_SCVx_main (the compiled C
library) instead of a single pre-solved Python/cvxpy trajectory.

Sigma-Delta scheduler on a 1 s pulse grid, pair-exclusive actuation.

----------------------------------------------------------------------------
THE PROBLEM
-----------
With a 5% inertia error, no purely feedforward scheme can deliver a good
final pointing: the inertia error is a MODEL error, invisible to any
scheduler that only propagates what it believes it is doing. So the strategy
has to close the loop on the state -- and the EKF that estimates the
disturbance already provides it.

FOUR MODES (control_mode)
--------------------------
 "open_loop" : COGU replayed as pulses. Baseline.
 "ff_dist"   : + real-time disturbance rejection (feedforward only).
 "terminal"  : + quantized tracking with TERMINAL-OPTIMAL gains, following a
               single, fixed COGU_SCVx_main solve computed once up front.
 "nmpc"      : same terminal-optimal tracking, but the guidance itself is
               periodically REFRESHED by calling COGU_SCVx_main again mid-
               flight with the current estimated state as the new `start`
               and the remaining time as the new `tf_init` -- a receding-
               horizon / NMPC use of the compiled C solver. This is the
               point of exposing every problem parameter (including the
               boundary conditions) as a COGU_SCVx_main runtime argument:
               nothing stops it from being called again, mid-simulation,
               "whenever desired".

WHY "terminal" / "nmpc"
------------------------
Attitude error is the double integral of torque error. For a double
integrator, the minimum-effort law that nulls BOTH the angle and the rate
error exactly at t_go = 0 is

      a_corr = -(6/t_go^2) dtheta - (4/t_go) domega

i.e. a PD with time-varying gains that tighten as the manoeuvre ends. It
targets the FINAL pointing error directly instead of the instantaneous one.
Feeding tau_ff - tau_d_hat + J*a_corr through the Sigma-Delta quantizer keeps
the pair-exclusivity and 1 s grid constraints exactly.

A NOTE ON STATE CONVENTIONS
----------------------------
COGU_SCVx_main's state is [quaternion(4), angular MOMENTUM h_b = J*omega(3)]
(see build_vat_problem_data / cogu_scvx.h) -- that is what `start`/`end` and
COGU's x_out use. The plant/estimator code below (dynamics, state_estimate,
terminal_correction) works in angular VELOCITY, which is the natural unit
for an EKF's rate channel and for a PD-style correction law. `Guidance`
below converts COGU's h_b output to omega = Jsat_inv @ h_b exactly once,
right after each solve, so every controller in this file sees a plain
[quaternion, angular rate] trajectory; converting back the other way
(omega -> h_b = Jsat_nominal @ omega) is the one line needed before handing
an estimated state to COGU_SCVx_main as a new `start` (see `run_nmpc_resolve`
below).

ACTUATION CONSTRAINT (respected in every mode)
------------------------------------------------
One pulse per second, one thruster per antagonistic pair {VAT1,VAT2} and
{VAT3,VAT4}: 9 admissible actuations per slot. Their torques come in +/-
pairs and rank(B) = 3, so any direction is reachable ON AVERAGE.

<<< EKF HOOKS >>>  disturbance_estimate(t) and state_estimate(t, x)

Build the library first (from this directory):
    cmake -B build -DCMAKE_C_FLAGS="-O3"
    cmake --build build --config Release
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                    # this dir: cogu_bindings.py
sys.path.insert(0, os.path.dirname(_HERE))   # parent dir: birds_vat_only_dilation.py

from cogu_bindings import load_library, cogu_scvx_main, COGU_WarmMode  # noqa: E402
from birds_vat_only_dilation import build_vat_problem_data  # noqa: E402


# =============================================================================
# 1. COGU_SCVx_main CALLER -- the one primitive everything below is built on
# =============================================================================

def solve_cogu(lib, data, start, end, tf_init, warm_mode=COGU_WarmMode.SLERP,
               warm_x=None, warm_u=None, verbose=0, **tuning):
    """
    One call to COGU_SCVx_main. Call it once up front for a fixed guidance
    trajectory, or call it again from inside the simulation loop with a
    fresh `start` (and a fresh `tf_init` = remaining time) for a receding-
    horizon / NMPC guidance refresh -- see `control_mode == 'nmpc'` below.

    tf_init is REQUIRED (no default): this library was generated with time
    dilation enabled, so COGU_SCVx_main has no fallback "good" final time to
    assume -- every caller, including a mid-flight re-solve, must supply one.
    """
    call_kwargs = dict(
        defect_penalty_scaled=1, epsilon_stop=1e-3, e_tol=0.005,
        beta_gr=2.0, beta_sh=2.0, eta_init=1.0, eta1=10.0, eta0=0.001,
        rho2=0.7, rho1=0.1, rho0=0.0, max_stall=3, max_iter=15, size_n=10,
        verbose=verbose,
        tf_fixed=tf_init, lam=data['lamb'],
        Sx=data['S_x'], cx=data['c_x'], Su=data['S_u'], cu=data['c_u'],
        alpha_p=data['alpha_p'], tf_min=data['tf_min'], tf_max=data['tf_max'],
        tf_init=tf_init, Sp=data['S_p'], cp=data['c_p'],
        warm_mode=warm_mode, warm_x=warm_x, warm_u=warm_u,
        quat_offset=(0 if warm_mode == COGU_WarmMode.SLERP else -1),
        omega_offset=(4 if warm_mode == COGU_WarmMode.SLERP else -1),
        omega_transform=(data['Jsat'] if warm_mode == COGU_WarmMode.SLERP else None),
    )
    call_kwargs.update(tuning)
    return cogu_scvx_main(lib, start=start, end=end, dp=data['dyn_params_val'], **call_kwargs)


class Guidance:
    """
    One COGU_SCVx_main solution, anchored at absolute simulation time `t0`
    and re-expressed in [quaternion, angular RATE] (see module docstring).
    Every `*_at(t)` method takes an ABSOLUTE simulation time, so a controller
    never needs to know whether the guidance it is following is the original
    one or came from a mid-flight NMPC re-solve -- it just calls
    `guidance.q_ref_at(t)` and picks up whichever Guidance object is current.
    """

    def __init__(self, cogu_res, t0, Jsat_inv):
        self.res = cogu_res
        self.t0 = float(t0)
        self.tf = float(cogu_res['tf'])
        self.x = cogu_res['x'].copy()
        self.x[4:7, :] = Jsat_inv @ cogu_res['x'][4:7, :]   # h_b -> omega
        self.u = cogu_res['u']
        self.T = self.x.shape[1] - 1
        self.dt_control = self.tf / self.T

    @property
    def t_end(self):
        return self.t0 + self.tf

    def _interp(self, t):
        tr = min(max(t - self.t0, 0.0), self.tf)
        idx = min(int(np.floor(tr / self.dt_control)), self.T - 1)
        a = min(max((tr - idx * self.dt_control) / self.dt_control, 0.0), 1.0)
        return tr, idx, a

    def q_ref_at(self, t):
        tr, idx, a = self._interp(t)
        if tr >= self.tf:
            return self.x[0:4, self.T] / np.linalg.norm(self.x[0:4, self.T])
        q0, q1 = self.x[0:4, idx], self.x[0:4, idx + 1]
        if np.dot(q0, q1) < 0.0:
            q1 = -q1
        qr = (1 - a) * q0 + a * q1
        return qr / np.linalg.norm(qr)

    def w_ref_at(self, t):
        tr, idx, a = self._interp(t)
        if tr >= self.tf:
            return self.x[4:7, self.T]
        return (1 - a) * self.x[4:7, idx] + a * self.x[4:7, idx + 1]

    def u_guid_zoh(self, t):
        tr = min(max(t - self.t0, 0.0), self.tf)
        k = int(np.floor(tr / self.dt_control))
        if k < 0 or k >= self.u.shape[1]:
            return np.zeros(4)
        return self.u[:, k]


# =============================================================================
# 2. QUATERNION HELPERS
# =============================================================================

def quat_mult(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])


def quat_conj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])


# =============================================================================
# 3. PLANT
# =============================================================================

def dynamics(t, x, u, params):
    """x = [quaternion(4), angular RATE(3)]."""
    q = x[0:4]; w = x[4:7]
    I, I_inv = params["I"], params["I_inv"]
    tau_ff = (params["B"] @ u) if params["use_thrusters"] else np.zeros(3)
    tau_d = params["disturbance_torque"](t)
    w_dot = I_inv @ (tau_ff + tau_d - np.cross(w, I @ w))
    wx, wy, wz = w; q0, q1, q2, q3 = q
    q_dot = 0.5 * np.array([
        -wx*q1 - wy*q2 - wz*q3,
         wx*q0 + wz*q2 - wy*q3,
         wy*q0 - wz*q1 + wx*q3,
         wz*q0 + wy*q1 - wx*q2])
    return np.concatenate((q_dot, w_dot))


def rk4_step(f, t, x, u, dt, params):
    k1 = f(t,        x,             u, params)
    k2 = f(t + dt/2, x + dt/2 * k1, u, params)
    k3 = f(t + dt/2, x + dt/2 * k2, u, params)
    k4 = f(t + dt,   x + dt   * k3, u, params)
    return x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)


# =============================================================================
# 4. TERMINAL-OPTIMAL CORRECTION
# =============================================================================

def terminal_correction(t, q_hat, w_hat, guidance, I_nominal, t_go_min, corr_clip, tau_max):
    """a_corr = -(6/tgo^2) dtheta - (4/tgo) domega  (min-effort terminal law)."""
    t_go = max(guidance.t_end - t, t_go_min)
    q_r = guidance.q_ref_at(t)
    q_e = quat_mult(quat_conj(q_r), q_hat)
    if q_e[0] < 0.0:
        q_e = -q_e
    dtheta = 2.0 * q_e[1:4]
    domega = w_hat - guidance.w_ref_at(t)
    a_corr = -(6.0 / t_go**2) * dtheta - (4.0 / t_go) * domega
    tau_corr = I_nominal @ a_corr
    n = np.linalg.norm(tau_corr)
    if n > corr_clip * tau_max:
        tau_corr = tau_corr * (corr_clip * tau_max / n)
    return tau_corr, dtheta, domega


# =============================================================================
# 5. SIGMA-DELTA SCHEDULER
# =============================================================================

class VATScheduler:
    """One .step() per pulse slot. Fed only with information available at
    t_k, plus whichever Guidance object is currently active (fixed for
    'terminal', periodically replaced for 'nmpc')."""

    def __init__(self, M, mode, U_OPT, DH_OPT, NFIRE, W_h, lambda_fuel,
                h_clip, lookahead, dt_pulse, dist_predict_order,
                terminal_correction_fn):
        self.M = M
        self.mode = mode
        self.reject = (mode in ("ff_dist", "terminal", "nmpc"))
        self.closed = (mode in ("terminal", "nmpc"))
        self.U_OPT, self.DH_OPT, self.NFIRE = U_OPT, DH_OPT, NFIRE
        self.W_h, self.lambda_fuel, self.h_clip = W_h, lambda_fuel, h_clip
        self.lookahead, self.dt_pulse = lookahead, dt_pulse
        self.dist_predict_order = dist_predict_order
        self.terminal_correction_fn = terminal_correction_fn
        self.h = np.zeros(3)
        self.tau_prev = None
        self.tau_pred_prev = None
        self.n_windup = 0
        self.h_hist = np.zeros((M + 1, 3))
        self.tau_pred_hist = np.zeros((M, 3))
        self.dth_hist = np.zeros((M, 3))

    def _best_option(self, h, dh_seq, idx, depth):
        best_c, best_j = np.inf, 0
        dh_des = dh_seq[idx]
        for j in range(self.DH_OPT.shape[1]):
            h_next = h + dh_des - self.DH_OPT[:, j]
            c = float(h_next @ self.W_h @ h_next) + self.lambda_fuel * self.NFIRE[j]
            if depth > 1 and idx + 1 < len(dh_seq):
                n = np.linalg.norm(h_next)
                h_cl = h_next * (self.h_clip / n) if n > self.h_clip else h_next
                c += self._best_option(h_cl, dh_seq, idx + 1, depth - 1)[0]
            if c < best_c:
                best_c, best_j = c, j
        return best_c, best_j

    def step(self, k, t_now, tau_hat, q_hat, w_hat, guidance):
        # ---- momentum accumulator ------------------------------------------
        if self.closed:
            # measured momentum debt: absorbs disturbance, inertia error, bias
            self.h = self._I_nominal @ (guidance.w_ref_at(t_now) - w_hat)
        elif self.reject and self.tau_prev is not None:
            tau_mean = 0.5 * (self.tau_prev + tau_hat)
            self.h = self.h + (self.tau_pred_prev - tau_mean) * self.dt_pulse

        # ---- causal prediction of the disturbance over the coming slot ------
        if self.reject:
            if self.dist_predict_order == 1 and self.tau_prev is not None:
                slope = (tau_hat - self.tau_prev) / self.dt_pulse
            else:
                slope = np.zeros(3)
        else:
            tau_hat = np.zeros(3); slope = np.zeros(3)

        # ---- terminal-optimal correction -----------------------------------
        if self.closed:
            tau_corr, dtheta, _ = self.terminal_correction_fn(t_now, q_hat, w_hat, guidance)
            self.dth_hist[k] = dtheta
        else:
            tau_corr = np.zeros(3)

        # ---- desired momentum increments (COGU guidance known ahead) --------
        dh_seq = []
        for i in range(self.lookahead):
            kk = k + i
            if kk >= self.M:
                break
            tm = t_now + i * self.dt_pulse + 0.5 * self.dt_pulse
            tau_cmd = self._B @ guidance.u_guid_zoh(tm)
            if self.reject:
                tau_cmd = tau_cmd - (tau_hat + slope * ((i + 0.5) * self.dt_pulse))
            if self.closed:
                tau_cmd = tau_cmd + tau_corr
            dh_seq.append(tau_cmd * self.dt_pulse)
        if not dh_seq:
            dh_seq = [np.zeros(3)]

        tau_pred_now = tau_hat + slope * (0.5 * self.dt_pulse) if self.reject else np.zeros(3)

        # ---- decide ---------------------------------------------------------
        _, j = self._best_option(self.h, dh_seq, 0, self.lookahead)
        self.h = self.h + dh_seq[0] - self.DH_OPT[:, j]
        n = np.linalg.norm(self.h)
        if n > self.h_clip and not self.closed:
            self.h = self.h * (self.h_clip / n); self.n_windup += 1

        self.tau_prev = tau_hat.copy()
        self.tau_pred_prev = tau_pred_now
        self.h_hist[k + 1] = self.h
        self.tau_pred_hist[k] = tau_pred_now
        return self.U_OPT[j]


# =============================================================================
# 6. FULL CLOSED-LOOP EVALUATION
# =============================================================================

def run_vat_cogu_closed_loop(
    problem_kwargs=None,
    control_mode='terminal', compare_all=True,
    nmpc_resolve_period=60.0,
    lib_path=None,
    # actuator model
    dt_pulse=1.0, pulse_width=1.0, hold_time=0.0,
    # perturbations
    pert_seed=0, thrust_offset=0.0, inertia_var=0.05,
    # EKF hook tuning (identity by default -- override est_scale/est_bias/
    # est_sigma/omega_sigma/quat_sigma to simulate an imperfect EKF)
    est_scale=1.0, est_bias=(0.0, 0.0, 0.0), est_sigma=0.0,
    dist_predict_order=0, omega_sigma=0.0, quat_sigma=0.0,
    # terminal-optimal gains
    t_go_min=15.0, corr_clip=1.0,
    # Sigma-Delta quantizer
    lambda_fuel=0.0, h_clip_bits=3.0, lookahead=3,
    dt_state_target=0.1,
    plot=True, verbose=True,
):
    """
    Evaluate COGU (the compiled C SCvx guidance) + VAT actuation under
    perturbations, in closed loop, exactly reproducing the structure of the
    reference disturbance-rejection / Sigma-Delta script -- except every
    place that script assumed a pre-solved `x_value`/`u_value`/`tf`/`T` now
    gets them from `solve_cogu(...)` (a direct COGU_SCVx_main call), and
    `control_mode='nmpc'` calls it again, periodically, from inside the
    simulation loop.

    Returns (runs, data, guidance0): `runs` maps a mode label to its
    simulation-result dict (t, q, w, u, ucmd, err, plan, final_pt,
    final_rate, ...), `data` is build_vat_problem_data(...)'s dict, and
    `guidance0` is the initial Guidance (COGU_SCVx_main solved once, before
    the loop).
    """
    data = dict(problem_kwargs or {})
    data = build_vat_problem_data(**data)
    lib = load_library(lib_path)

    B = data['B']
    I_nominal = data['Jsat'].copy()
    I_nom_inv = np.linalg.inv(I_nominal)
    Fth_max = data['Fth_max']

    q_initial = np.array([1.0, 0.0, 0.0, 0.0])
    q_final = np.array([0.0, 1.0, 0.0, 0.0])

    #q_initial = np.array([0.5, 0.5, 0.5, 0.5])
    #q_final = np.array([0.8660254, 0.2886751, 0.2886751, 0.2886751])

    #q_initial = np.array([1.0, 0.0, 0.0, 0.0])
    #q_final = np.array([np.cos(np.deg2rad(25)/2), np.sin(np.deg2rad(25)/2), 0.0, 0.0])

    w_initial = np.zeros(3)
    w_final = np.zeros(3)

    q_initial /= np.linalg.norm(q_initial)
    q_final   /= np.linalg.norm(q_final)

    h_initial = I_nominal @ np.asarray(w_initial)
    h_final   = I_nominal @ np.asarray(w_final)

    start0 = np.concatenate([q_initial, h_initial])
    end0   = np.concatenate([q_final, h_final])

    end_quat = q_final.copy()

    end_quat = data['end_quat'][:, 0]

    # ---- 1. Initial guidance: the first of possibly many COGU_SCVx_main calls
    if verbose:
        print("=" * 68)
        print("VAT_COGU_C_caller_test -- COGU+VAT closed-loop evaluation")
        print("=" * 68)
        print("Calling COGU_SCVx_main for the initial guidance "
             "(it can be called again at any time -- see control_mode='nmpc')...")
    res0 = solve_cogu(lib, data, start=start0, end=end0, tf_init=data['tf_init'],
                      warm_mode=COGU_WarmMode.SLERP, verbose=0)
    guidance0 = Guidance(res0, t0=0.0, Jsat_inv=I_nom_inv)
    if verbose:
        print(f"Time dilation solution: tf = {guidance0.tf:.2f} s, "
             f"dt_control = {guidance0.dt_control:.3f} s, "
             f"converged={res0['converged']} iterations={res0['iterations']}")

    # ---- 2. Actuator model: one pulse opportunity per second ----------------
    thrust_nominal = Fth_max
    U_OPT = []
    for a in (None, 0, 1):
        for b in (None, 2, 3):
            u = np.zeros(4)
            if a is not None: u[a] = thrust_nominal
            if b is not None: u[b] = thrust_nominal
            U_OPT.append(u)
    U_OPT = np.array(U_OPT)                          # (9,4)
    DH_OPT = (B @ U_OPT.T) * pulse_width * dt_pulse   # (3,9)
    NFIRE = (U_OPT > 0).sum(axis=1)
    q_ref_norm = np.linalg.norm(DH_OPT, axis=0).max()
    h_clip = h_clip_bits * q_ref_norm

    if verbose:
        print("Torque authority per thruster [uN.m]:",
             np.round(np.linalg.norm(B, axis=0) * thrust_nominal * 1e6, 3))
        print("Momentum quantum per pulse   [uN.m.s]:", np.round(q_ref_norm * 1e6, 4))
        duty_max = max((guidance0.u[0, :] + guidance0.u[1, :]).max(),
                       (guidance0.u[2, :] + guidance0.u[3, :]).max()) / thrust_nominal
        spare = 1.0 - 0.5 * ((guidance0.u[0]+guidance0.u[1]).mean()
                             + (guidance0.u[2]+guidance0.u[3]).mean()) / thrust_nominal
        print(f"Max COGU pair duty = {duty_max:.2f} (pulse_width = {pulse_width:.2f}), "
             f"spare authority ~ {spare*100:.0f} %")
        if duty_max > pulse_width + 1e-9:
            print("  WARNING: pulse_width < required duty -> under-actuated manoeuvre.")

    # ---- 3. Perturbations / uncertainty --------------------------------------
    rng = np.random.default_rng(pert_seed)
    thrust_sigma = thrust_nominal / 10
    thrust_bias_sigma = thrust_nominal / 20
    thrust_bias = thrust_bias_sigma * rng.standard_normal(4)

    def disturbance_torque(t):
        return np.array([
            0.2e-7 * np.sin(2*np.pi*t/5400),
            0.3e-7 * np.sin(2*np.pi*t/4200 + 1.3),
            0.3e-7 * np.sin(2*np.pi*t/6800 + 2.1),
        ])

    _est_rng = np.random.default_rng(pert_seed + 7)
    est_bias_arr = np.asarray(est_bias, dtype=float)

    def disturbance_estimate(t):
        tau = est_scale * disturbance_torque(t) + est_bias_arr
        if est_sigma > 0.0:
            tau = tau + est_sigma * _est_rng.standard_normal(3)
        return tau

    _w_rng = np.random.default_rng(pert_seed + 11)

    def state_estimate(t, x_true):
        q = x_true[0:4].copy(); w = x_true[4:7].copy()
        if omega_sigma > 0.0:
            w = w + omega_sigma * _w_rng.standard_normal(3)
        if quat_sigma > 0.0:
            dth = quat_sigma * _w_rng.standard_normal(3)
            dq = np.concatenate(([1.0], 0.5 * dth))
            q = quat_mult(q, dq); q /= np.linalg.norm(q)
        return q, w

    tau_max = np.linalg.norm(B, axis=0).max() * thrust_nominal

    def terminal_correction_fn(t, q_hat, w_hat, guidance):
        return terminal_correction(t, q_hat, w_hat, guidance, I_nominal,
                                   t_go_min, corr_clip, tau_max)

    while True:  # TRUE perturbed inertia
        P = 1.0 + inertia_var * (2.0 * rng.random((3, 3)) - 1.0)
        P = 0.5 * (P + P.T)
        I_true = I_nominal * P
        if np.all(np.linalg.eigvalsh(I_true) > 0.0):
            break
    if verbose:
        print("Diagonal inertia error (%):",
             np.round(100*(np.diag(I_true)-np.diag(I_nominal))/np.diag(I_nominal), 1))

    plant_params = {"I": I_true, "I_inv": np.linalg.inv(I_true), "B": B,
                    "use_thrusters": True, "disturbance_torque": disturbance_torque}

    steps_per_slot = max(1, int(np.ceil(dt_pulse / dt_state_target)))
    dt_state = dt_pulse / steps_per_slot
    n_on = max(1, int(round(pulse_width * steps_per_slot)))

    # ---- 4. Simulation on the pulse grid: measure, decide, then advance -----
    def run(mode, noise_seed):
        M = int(round((guidance0.tf + hold_time) / dt_pulse))
        sched = VATScheduler(M, mode, U_OPT, DH_OPT, NFIRE,
                             W_h=np.eye(3), lambda_fuel=lambda_fuel, h_clip=h_clip,
                             lookahead=lookahead, dt_pulse=dt_pulse,
                             dist_predict_order=dist_predict_order,
                             terminal_correction_fn=terminal_correction_fn)
        sched._I_nominal = I_nominal
        sched._B = B
        noise = np.random.default_rng(noise_seed).standard_normal((M, 4))

        N_steps = M * steps_per_slot + 1
        q_h = np.zeros((N_steps, 4)); w_h = np.zeros((N_steps, 3))
        u_h = np.zeros((N_steps, 4)); uc_h = np.zeros((N_steps, 4))
        e_h = np.zeros(N_steps);      t_h = np.zeros(N_steps)
        plan = np.zeros((M, 4))
        resolve_log = []

        guidance = guidance0
        last_resolve_t = 0.0

        x = np.concatenate([guidance0.x[0:4, 0], guidance0.x[4:7, 0]])
        x[0:4] /= np.linalg.norm(x[0:4])
        t = 0.0; i = 0
        for k in range(M):
            tau_hat = disturbance_estimate(t)
            q_hat, w_hat = state_estimate(t, x)

            if (mode == "nmpc" and (t - last_resolve_t) >= nmpc_resolve_period
                   and t < guidance0.t_end - t_go_min):
                remaining_tf = max(data['tf_min'], guidance0.t_end - t)
                start_now = np.concatenate([q_hat, I_nominal @ w_hat])  # omega -> h_b
                res_k = solve_cogu(lib, data, start=start_now, end=end0,
                                   tf_init=remaining_tf, warm_mode=COGU_WarmMode.SLERP,
                                   verbose=0)
                guidance = Guidance(res_k, t0=t, Jsat_inv=I_nom_inv)
                last_resolve_t = t
                resolve_log.append(dict(t=t, tf_remaining=guidance.tf,
                                        converged=res_k['converged'],
                                        iterations=res_k['iterations']))
                if verbose:
                    print(f"  [nmpc] re-solved COGU_SCVx_main at t={t:6.1f}s: "
                         f"tf_remaining={guidance.tf:7.2f}s "
                         f"converged={res_k['converged']} iterations={res_k['iterations']}")

            u_cmd = sched.step(k, t, tau_hat, q_hat, w_hat, guidance)
            plan[k] = u_cmd

            fired = u_cmd > 0.0
            u_real = np.where(fired,
                              np.clip(u_cmd + thrust_offset + thrust_bias
                                     + thrust_sigma*noise[k], 0.0, None),
                              0.0)
            for s in range(steps_per_slot):
                on = (s < n_on)
                uu = u_real if on else np.zeros(4)
                uu_c = u_cmd if on else np.zeros(4)
                q_h[i] = x[0:4]; w_h[i] = x[4:7]; t_h[i] = t
                u_h[i] = uu; uc_h[i] = uu_c
                qr = guidance.q_ref_at(t)
                e_h[i] = 2.0*np.degrees(np.arccos(np.clip(abs(np.dot(x[0:4], qr)), 0, 1)))
                x = rk4_step(dynamics, t, x, uu, dt_state, plant_params)
                x[0:4] /= np.linalg.norm(x[0:4])
                t += dt_state; i += 1

        q_h[i] = x[0:4]; w_h[i] = x[4:7]; t_h[i] = t
        e_h[i] = 2.0*np.degrees(np.arccos(np.clip(abs(np.dot(x[0:4], guidance.q_ref_at(t))), 0, 1)))
        return dict(t=t_h, q=q_h, w=w_h, u=u_h, ucmd=uc_h, err=e_h,
                   plan=plan, t_slot=np.arange(M)*dt_pulse, sched=sched,
                   resolve_log=resolve_log, final_guidance=guidance)

    # ---- 5. Run and compare ---------------------------------------------------
    q_end = end_quat
    w_end = guidance0.x[4:7, guidance0.T]

    MODES = {"open_loop": "A: open loop (COGU only)",
             "ff_dist":   "B: + disturbance rejection",
             "terminal":  "C: + terminal-optimal tracking",
             "nmpc":      "D: + NMPC (periodic COGU re-solve)"}
    to_run = list(MODES) if compare_all else [control_mode]

    runs = {}
    for m in to_run:
        out = run(m, noise_seed=pert_seed + 1)
        out["final_pt"] = 2.0*np.degrees(np.arccos(np.clip(abs(np.dot(out["q"][-1], q_end)), 0, 1)))
        out["final_rate"] = np.degrees(np.linalg.norm(out["w"][-1] - w_end))
        runs[MODES[m]] = out
        plan = out["plan"]
        excl = np.all(plan[:, 0]*plan[:, 1] == 0) and np.all(plan[:, 2]*plan[:, 3] == 0)
        if verbose:
            print(f"\n--- {MODES[m]} ---")
            print(f"  pair exclusivity           : {excl}")
            print(f"  pulses per thruster        : {(plan > 0).sum(axis=0)}  (total {(plan > 0).sum()})")
            print(f"  final pointing error       : {out['final_pt']:.3f} deg")
            print(f"  final rate error           : {out['final_rate']:.4f} deg/s")
            print(f"  max deviation vs reference : {out['err'].max():.3f} deg")
            if out["resolve_log"]:
                print(f"  COGU_SCVx_main re-solves   : {len(out['resolve_log'])}")

    if verbose:
        print("\n" + "="*62)
        print(f"{'mode':<34}{'pointing [deg]':>16}{'pulses':>10}")
        for tag, out in runs.items():
            print(f"{tag:<34}{out['final_pt']:>16.3f}{int((out['plan']>0).sum()):>10}")
        print("="*62)

    if plot:
        _plot_vat_cogu_closed_loop(runs, guidance0, MODES, to_run)

    return runs, data, guidance0


def _plot_vat_cogu_closed_loop(runs, guidance0, MODES, to_run):
    t_nodes = np.arange(guidance0.T + 1) * guidance0.dt_control
    ref = runs[MODES[to_run[-1]]]

    plt.figure(figsize=(7, 3))
    plt.plot(ref["t"], ref["q"])
    plt.plot(t_nodes, guidance0.x[0:4, :].T, 'o', ms=3, alpha=0.5)
    plt.axvline(guidance0.tf, color='k', ls=':', lw=1, alpha=0.6)
    plt.legend(["q0", "q1", "q2", "q3"], fontsize=8)
    plt.title("Quaternion  (- achieved, o COGU guidance)")
    plt.xlabel("Time [s]"); plt.grid(); plt.tight_layout()

    plt.figure(figsize=(7, 3))
    plt.plot(ref["t"], ref["w"] * 180/np.pi)
    plt.plot(t_nodes, guidance0.x[4:7, :].T * 180/np.pi, 'o', ms=3, alpha=0.5)
    plt.axvline(guidance0.tf, color='k', ls=':', lw=1, alpha=0.6)
    plt.legend(["wx", "wy", "wz"], fontsize=8)
    plt.title("Angular Velocity [deg/s]")
    plt.xlabel("Time [s]"); plt.grid(); plt.tight_layout()

    plt.figure(figsize=(7, 3.2))
    for tag, out in runs.items():
        plt.plot(out["t"], np.maximum(out["err"], 1e-4), label=tag)
    plt.axvline(guidance0.tf, color='k', ls=':', lw=1, alpha=0.6)
    plt.title("Deviation from the COGU reference (perturbed plant)")
    plt.xlabel("Time [s]"); plt.ylabel("error [deg]")
    plt.legend(fontsize=8); plt.grid(which='both'); plt.tight_layout()

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(6, 6))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    for i in range(4):
        axs[i].step(ref["t"], ref["u"][:, i]*1e6, where='post', color=colors[i], label=f'VAT{i+1}')
        axs[i].step(t_nodes[:-1], guidance0.u[i, :]*1e6, where='post', color='k', ls='--',
                    lw=1.0, alpha=0.6, label='COGU (initial)')
        axs[i].set_ylabel(f'Fth{i+1} [uN]'); axs[i].grid(True); axs[i].legend(loc='upper right', fontsize=7)
    axs[-1].set_xlabel('Time [s]')
    fig.suptitle('VAT pulse train on the 1 s grid')
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode', choices=['open_loop', 'ff_dist', 'terminal', 'nmpc'],
                   default='terminal')
    p.add_argument('--no-compare-all', action='store_true',
                   help='only run --mode instead of all four modes')
    p.add_argument('--nmpc-resolve-period', type=float, default=60.0,
                   help='seconds between COGU_SCVx_main re-solves in nmpc mode')
    p.add_argument('--no-plot', action='store_true')
    p.add_argument('--quiet', action='store_true')
    p.add_argument('--lib', default=None)
    args = p.parse_args()
    run_vat_cogu_closed_loop(
        control_mode=args.mode, compare_all=not args.no_compare_all,
        nmpc_resolve_period=args.nmpc_resolve_period,
        lib_path=args.lib, plot=not args.no_plot, verbose=not args.quiet,
    )
