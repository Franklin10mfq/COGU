"""
C code generation.

What this module emits, given a solved problem:

    scvx_main.c               the SCvx loop, compiled from a template into
                              ONE public function, COGU_SCVx_main() (see
                              cogu_scvx.h) -- no main(). Every number that
                              is specific to a problem *instance* (boundary
                              conditions, physical parameters, scaling, warm
                              start, SCvx tuning, dilation bounds) is a
                              runtime argument of that function; only the
                              problem *shape* (nx, nu, ng, dynamic-parameter
                              count, horizon, ZOH/FOH, dilation on/off) is
                              compiled in, because that shape is what fixes
                              the embedded convex solver's workspace layout.
    cogu_scvx.h               the public C API: COGU_SCVx_main's prototype,
                              the COGU_WarmMode / COGU_Status enums, and the
                              read-only COGU_NX/COGU_NU/... shape constants
    dynamics.c / dynamics.h   the model and every Jacobian, including the
                              dilation Jacobians F = df/dp and G = dg/dp
    solver/                   the convex subproblem solver, produced by
                              CVXPYgen in a clean subprocess
    cpg_compat.h              maps a stacked-index update API onto the
                              per-timestep setters CVXPYgen generates
    embedded_data.h           a REFERENCE snapshot of the problem data
                              (scaling, bounds, warm start, boundary
                              conditions, SCvx tuning) from the Python call
                              that generated this project, as C constants.
                              Not included by scvx_main.c any more -- it is
                              just a convenient set of defaults to pass to
                              COGU_SCVx_main, e.g. from a CLI test harness.
    CMakeLists.txt            builds cogu_scvx as a shared library (.dll/.so)

The SCvx loop itself (`scvx_main.c`) is emitted from a template.  Drop your
template at `scvx_toolbox/templates/scvx_loop_template.c` and it will be
filled in; without it the rest of the project is still generated and a note
is printed.  The placeholders the template may use are listed in
`TEMPLATE_PLACEHOLDERS`.
"""

import os
import pickle
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from .cost_dsl import node_count, resolve_term
from .modes import FOH, parse_time_mode
from .symbolic import matrix_to_c_assignments

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), 'templates',
                             'scvx_loop_template.c')

TEMPLATE_PLACEHOLDERS = (
    '{NX}', '{NU}', '{NG}', '{NP}', '{T}', '{U_NODES}',
    '{DILATION}', '{FOH}', '{TIME_MODE}',
    '{USER_COST}',
)

# SCvx tuning mirrored into the C loop as #defines, so the generated program
# runs the same algorithm as the Python call that produced it.
_SOLVER_DEFINES = ('SIZE_N', 'MAX_ITER', 'MAX_STALL',
                   'RHO0', 'RHO1', 'RHO2',
                   'ETA0', 'ETA1', 'ETA_INIT',
                   'BETA_SH', 'BETA_GR',
                   'E_TOL', 'EPSILON_STOP',
                   'DEFECT_PENALTY_SCALED')

# Fixed compile-time capacity of the hist_tf/hist_J/hist_rho diagnostic
# buffers in COGU_SCVx_main (see generate_cogu_header). A Simulink C Caller
# port needs a fixed array size, but max_iter is a genuine runtime tunable,
# so history recording is capped here independently of max_iter -- the SCvx
# solve itself is never limited by this constant, only what gets recorded.
COGU_HIST_CAP = 64


# ======================================================================
# dynamics.c / dynamics.h
# ======================================================================

_C_SIGNATURE = ("const double *hx, const double *hu, const double *hp, "
                "double t, const double *dyn_par, double *out")

_C_DECL = ("void {name}(const double *hx, const double *hu, const double *hp, "
           "double t, const double *dyn_par, double *out);")


def generate_c_dynamics(model, output_dir=None):
    """
    Emit dynamics.c and dynamics.h for a SymbolicModel.

    Every function has the same signature

        void name(const double *hx, const double *hu, const double *hp,
                  double t, const double *dyn_par, double *out)

    with `out` in row-major order:

        f_dynamics    f(x,u[,p])        nx
        A_jacobian    df/dx             nx*nx
        B_jacobian    df/du             nx*nu
        F_jacobian    df/dp             nx*np   (dilation only)
        y_affine      f - Ax - Bu - Fp  nx
        g_constraints g(x,u)            ng
        C_jacobian    dg/dx             ng*nx
        D_jacobian    dg/du             ng*nu
        G_jacobian    dg/dp             ng*np   (dilation only)
        z_affine      g - Cx - Du - Gp  ng
    """
    names = (model.state_names, model.control_names,
             model.dilation_names, model.dyn_param_names)

    def body(matrix, label):
        if np.prod(np.shape(matrix)) == 0:
            return f"    (void)hx; (void)hu; (void)hp; (void)t; (void)dyn_par; (void)out;  /* {label}: empty */"
        return matrix_to_c_assignments(matrix, *names)

    blocks = [
        ('f_dynamics', body(model.f, 'f')),
        ('A_jacobian', body(model.A, 'A')),
        ('B_jacobian', body(model.B, 'B')),
        ('F_jacobian', body(model.F, 'F')),
        ('y_affine', body(model.y, 'y')),
        ('g_constraints', body(model.g, 'g')),
        ('C_jacobian', body(model.C, 'C')),
        ('D_jacobian', body(model.D, 'D')),
        ('G_jacobian', body(model.G, 'G')),
        ('z_affine', body(model.z, 'z')),
    ]

    c_parts = ['#include <math.h>', '#include "dynamics.h"', '']
    for name, block in blocks:
        c_parts.append(f"void {name}({_C_SIGNATURE}) {{\n{block}\n}}\n")
    c_source = '\n'.join(c_parts)

    h_parts = [
        '#ifndef DYNAMICS_H',
        '#define DYNAMICS_H',
        '',
        f'/* nx={model.nx}, nu={model.nu}, ng={model.ng}, np={model.n_p} */',
        '',
    ]
    for name, _ in blocks:
        h_parts.append(_C_DECL.format(name=name))
    h_parts += ['', '#endif /* DYNAMICS_H */', '']
    h_source = '\n'.join(h_parts)

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, 'dynamics.c'), 'w', encoding='utf-8') as fp:
            fp.write(c_source)
        with open(os.path.join(output_dir, 'dynamics.h'), 'w', encoding='utf-8') as fp:
            fp.write(h_source)

    return {'c_source': c_source, 'h_source': h_source}


# ======================================================================
# User cost in C
# ======================================================================

def generate_c_user_cost(cost_terms, nx, nu, T, u_nodes):
    """
    Emit the body of J_cost from the same cost DSL the Python side uses.

    The generated code assumes these are in scope: `x` (stride ns), `u`
    (stride mi), `T`, `u_nodes`, `tau`, `lam`, and the accumulator `cost`.
    """
    coeff_c = {'sqrt_tau': 'sqrt_tau', 'tau': 'tau', 'tau_lamb': 'tau_lamb'}
    resolved = [resolve_term(t) for t in cost_terms]

    used = {rt.coeff for rt in resolved if isinstance(rt.coeff, str)}
    lines = []
    if 'sqrt_tau' in used:
        lines.append('    const double sqrt_tau = sqrt(tau);')
    if 'tau_lamb' in used:
        lines.append('    const double tau_lamb = tau*lam;')

    offset_idx = 0
    for rt in resolved:
        arr, stride, dim = ('u', 'mi', nu) if rt.var == 'u' else ('x', 'ns', nx)
        lo = 0 if rt.slc.start is None else rt.slc.start
        hi = dim if rt.slc.stop is None else rt.slc.stop
        coeff = coeff_c[rt.coeff] if isinstance(rt.coeff, str) else repr(float(rt.coeff))
        weight = repr(rt.weight)
        n_nodes = node_count(rt.k_range, T, u_nodes)

        if rt.offset is not None:
            offs = np.asarray(rt.offset, dtype=float).reshape(-1)
            name = f'_off{offset_idx}'
            offset_idx += 1
            values = ','.join(repr(float(v)) for v in offs)
            lines.append(f'    static const double {name}[] = {{{values}}};')
            elem = f'{coeff}*({arr}k[i] - {name}[i-{lo}])'
        else:
            elem = f'{coeff}*{arr}k[i]'

        accum = 's += e*e;' if rt.kind in ('sumsq', 'norm2') else 's += fabs(e);'
        lines.append(f'    for (int k = 0; k < {n_nodes}; k++) {{')
        lines.append(f'        const double *{arr}k = {arr} + k*{stride};')
        lines.append('        double s = 0.0;')
        lines.append(f'        for (int i = {lo}; i < {hi}; i++) {{ double e = {elem}; {accum} }}')
        if rt.kind == 'norm2':
            lines.append(f'        cost += {weight}*sqrt(s);')
        else:
            lines.append(f'        cost += {weight}*s;')
        lines.append('    }')

    return '\n'.join(lines)


# ======================================================================
# Embedded problem data
# ======================================================================

def _c_array(name, values):
    body = ','.join(repr(float(v)) for v in np.asarray(values, dtype=float).reshape(-1))
    return f'static const double {name}[] = {{{body}}};'


def generate_c_embedded_data(embed, nx, nu, T, u_nodes, dilation):
    """
    Emit the problem data as C constants.

    Matrices are row-major; the warm-start trajectories are time-major
    (index i + t*stride), matching how the loop walks them.
    """
    wx = np.asarray(embed['warm_x'], dtype=float)
    wu = np.asarray(embed['warm_u'], dtype=float)
    if wu.shape[1] < u_nodes:
        wu = np.hstack([wu, np.repeat(wu[:, -1:], u_nodes - wu.shape[1], axis=1)])

    wx_tm = [wx[i, t] for t in range(T + 1) for i in range(nx)]
    wu_tm = [wu[i, t] for t in range(u_nodes) for i in range(nu)]

    lines = [
        f'static const double tf_fixed = {float(embed["tf"])!r};',
        f'static const double lam = {float(embed["lam"])!r};',
        _c_array('_start', embed['start']),
        _c_array('_end', embed['end']),
        _c_array('_dp', embed['dp']),
        _c_array('_Sx', embed['S_x']),
        _c_array('_cx', embed['c_x']),
        _c_array('_Su', embed['S_u']),
        _c_array('_cu', embed['c_u']),
        _c_array('_warm_x', wx_tm),
        _c_array('_warm_u', wu_tm),
    ]

    # SCvx tuning as #defines, consumed by the loop template
    defines = embed.get('solver_params') or {}
    for name in _SOLVER_DEFINES:
        if name in defines:
            value = defines[name]
            if isinstance(value, bool):
                text = '1' if value else '0'
            elif isinstance(value, int):
                text = str(value)
            else:
                text = repr(float(value))
            lines.insert(0, f'#define {name} {text}')

    if dilation:
        lines += [
            f'static const double alpha_p = {float(embed["alpha_p"])!r};',
            f'static const double tf_min = {float(embed["tf_min"])!r};',
            f'static const double tf_max = {float(embed["tf_max"])!r};',
            f'static const double tf_init = {float(embed["tf_init"])!r};',
            _c_array('_Sp', embed['S_p']),
            _c_array('_cp', embed['c_p']),
        ]

    return '\n'.join(lines)


def write_embedded_header(embed, nx, nu, ng, T, u_nodes, dilation, output_dir):
    """
    Write a REFERENCE snapshot of the problem instance that produced this C
    project, as C constants -- the exact numeric values the Python
    `solve_trajectory(...)` call used for boundary conditions, scaling, warm
    start and SCvx tuning.

    scvx_main.c does NOT include this file: every one of these numbers is a
    runtime argument of COGU_SCVx_main (see cogu_scvx.h), which is what lets
    one compiled library be called again with different values (Monte Carlo
    sampling, NMPC re-solves) without recompiling. This header exists purely
    as a convenient, copy-pasteable set of defaults -- e.g. for a CLI test
    harness that reproduces the original Python solve.
    """
    content = f"""\
/* embedded_data.h -- REFERENCE problem data, generated by scvx_toolbox.
 * NOT included by scvx_main.c -- COGU_SCVx_main takes all of this as
 * runtime arguments (see cogu_scvx.h). This file is only a convenient,
 * copy-pasteable snapshot of the values the Python solve_trajectory(...)
 * call used, for a CLI harness or as a starting point for your own caller.
 */
#ifndef EMBEDDED_DATA_H
#define EMBEDDED_DATA_H

#define STATES_SIZE   {nx}
#define INPUTS_SIZE   {nu}
#define NONCVX_SIZE   {ng}
#define HORIZON_T     {T}
#define CONTROL_NODES {u_nodes}
#define TIME_DILATION {1 if dilation else 0}

{generate_c_embedded_data(embed, nx, nu, T, u_nodes, dilation)}

#endif /* EMBEDDED_DATA_H */
"""
    path = os.path.join(output_dir, 'embedded_data.h')
    with open(path, 'w', encoding='utf-8') as fp:
        fp.write(content)
    return path


# ======================================================================
# cogu_scvx.h -- public API for the generated COGU_SCVx_main
# ======================================================================

def generate_cogu_header(nx, nu, ng, n_dyn_par, T, u_nodes, dilation, foh,
                         output_dir, hist_cap=COGU_HIST_CAP):
    """
    Emit cogu_scvx.h: the COGU_WarmMode / COGU_Status enums, the read-only
    COGU_* shape constants, and the COGU_SCVx_main prototype. Identical for
    every generated project except the enum values, so it is built here in
    Python rather than templated.

    Every COGU_SCVx_main parameter -- scalars included -- is declared with an
    explicit array dimension (`type name[1]` for a scalar), never a bare
    scalar or a `type *name` pointer: see the header's own doc comment for
    why (Simulink's C Caller block needs a literal array dimension on every
    port). `hist_cap` is COGU_HIST_CAP, the fixed capacity of the hist_tf/
    hist_J/hist_rho diagnostic buffers.
    """
    content = f"""\
/*
 * cogu_scvx.h -- public C API for COGU_SCVx_main.
 *
 * COGU_SCVx_main is the SCvx guidance loop compiled ahead of time (O3) for
 * ONE fixed problem shape (state/control/constraint/parameter counts,
 * horizon, ZOH-vs-FOH, dilation on/off). Everything about *this instance* of
 * the problem -- boundary conditions, physical model parameters, scaling,
 * warm start, SCvx tuning, time-dilation bounds -- is a runtime argument, so
 * the same compiled library can be called over and over (Monte Carlo, NMPC,
 * a Simulink C Caller block) without touching the compiler.
 *
 * What CANNOT be a runtime argument: nx, nu, ng, the dynamic-parameter count,
 * the horizon T, the control-node count, and whether dilation/FOH are on.
 * Those fix the size and sparsity of the embedded convex solver's workspace
 * (see cpg_workspace.c) at code-generation time; changing any of them means
 * re-running the Python generator (scvx_toolbox) and recompiling. They are
 * exposed below as read-only enum constants so a caller can size buffers and
 * assert compatibility instead of guessing.
 *
 * EVERY parameter below -- including every scalar -- is declared as an array
 * (`type name[1]` for a scalar, `type name[N]` for a vector/matrix), and
 * there is no bare `type name` or `type *name` anywhere in the prototype.
 * This is deliberate, not stylistic: Simulink's C Caller block parses this
 * header directly, and it only maps ports correctly when every argument has
 * an explicit, fixed array dimension -- a plain scalar or a bare pointer
 * does not parse into a signal the same way. See COGU_SCVx_C_caller.h for a
 * Simulink-oriented copy of this same prototype. Passing a length-1 array by
 * value still means "pass the address" at the ABI level (array parameters
 * decay to pointers in C), so nothing here is actually more expensive than
 * the equivalent pointer would be -- it is purely about what the header
 * text says, for the tools that read it literally.
 *
 * Generated by scvx_toolbox. This file matches this instance's scvx_main.c;
 * do not hand-edit either without keeping them in sync.
 */
#ifndef COGU_SCVX_H
#define COGU_SCVX_H

#ifdef __cplusplus
extern "C" {{
#endif

#if defined(_WIN32) || defined(_WIN64)
#  ifdef COGU_BUILD_DLL
#    define COGU_API __declspec(dllexport)
#  elif defined(COGU_USE_DLL)
#    define COGU_API __declspec(dllimport)
#  else
#    define COGU_API
#  endif
#elif defined(__GNUC__)
#  define COGU_API __attribute__((visibility("default")))
#else
#  define COGU_API
#endif

/* ---- fixed problem shape, baked in at generation time (read-only) ------- */
enum {{
    COGU_NX               = {nx},   /* state dimension                        */
    COGU_NU                = {nu},   /* control dimension                      */
    COGU_NG                = {ng},   /* number of nonconvex/path constraints   */
    COGU_NP_DYN            = {n_dyn_par},  /* length of the dynamic parameter vector */
    COGU_T                 = {T},  /* horizon: T+1 state nodes               */
    COGU_UNODES            = {u_nodes},  /* control nodes: T under ZOH, T+1 under FOH */
    COGU_DILATION_ENABLED  = {1 if dilation else 0},   /* 1: final time p is a decision variable */
    COGU_FOH_ENABLED       = {1 if foh else 0},   /* 1: control is first-order held         */

    /* Fixed capacity of hist_tf/hist_J/hist_rho below. A Simulink C Caller
     * port needs a compile-time array size, but max_iter is a genuine
     * runtime tuning knob -- so the history buffers are sized to this cap
     * instead of to max_iter. If max_iter[0] > COGU_HIST_CAP the SOLVE still
     * runs the full max_iter[0] iterations; only the recorded history is
     * truncated to the first COGU_HIST_CAP of them, and hist_len[0] reports
     * how many entries were actually written (always <= COGU_HIST_CAP). */
    COGU_HIST_CAP          = {hist_cap}
}};

/* ---- warm-start construction ---------------------------------------- */
typedef enum {{
    /* Caller fills warm_x_in / warm_u_in itself and passes them verbatim.
     * Required whenever the state layout is not a plain quaternion+rate
     * block, or whenever a previous COGU solve should seed the next one
     * (receding-horizon / NMPC use). */
    COGU_WARM_USER   = 0,

    /* Elementwise linear interpolation of every state component between
     * `start` and `end` across the T+1 nodes; controls start at zero. */
    COGU_WARM_LINEAR = 1,

    /* Like COGU_WARM_LINEAR, except the 4 components at `quat_offset` are
     * treated as a scalar-first unit quaternion and SLERP-interpolated
     * (so the warm start stays unit-norm), and -- if `omega_offset >= 0` --
     * the 3 components there are filled from the finite-difference angular
     * rate of that SLERP sequence (using the ABSOLUTE step tf_init/T, since
     * a rate is physical and does not live in normalized time), optionally
     * passed through `omega_transform` (e.g. an inertia matrix, to turn a
     * rate into a momentum). Falls back to COGU_WARM_LINEAR wholesale if
     * `quat_offset < 0`. */
    COGU_WARM_SLERP  = 2
}} COGU_WarmMode;

/* ---- solve outcome ---------------------------------------------------- */
typedef enum {{
    COGU_STATUS_CONVERGED = 0,
    COGU_STATUS_MAX_ITER  = 1,
    COGU_STATUS_STALLED   = 2
}} COGU_Status;

/*
 * COGU_SCVx_main -- run the SCvx loop once and return the guidance.
 *
 * All array arguments are dense `double`, row-major for matrices and
 * time-major for trajectories (index = i + t*stride), matching the layout
 * scvx_toolbox uses everywhere else. Every INPUT array may be a pointer to a
 * caller-owned buffer that is only read; nothing is retained after the call
 * returns, so the function is safe to call again immediately (Monte Carlo,
 * NMPC) with no state left over from the previous call.
 */
COGU_API void COGU_SCVx_main(
    /* ---------------- SCvx algorithm tuning ---------------- */
    int    defect_penalty_scaled[1], /* 1: penalize the defect in scaled units (matches how L
                                       * penalizes slack vc); 0: unscaled physical units */
    double epsilon_stop[1],        /* stop when the step size drops below this            */
    double e_tol[1],                /* stop when |Delta_L| < e_tol*|J_ref|                  */
    double beta_gr[1],              /* trust-region growth factor (rho >= rho2)             */
    double beta_sh[1],              /* trust-region shrink factor (rho < rho1, or rejected) */
    double eta_init[1],             /* initial trust-region radius                          */
    double eta1[1],                 /* trust-region ceiling                                 */
    double eta0[1],                 /* trust-region floor                                   */
    double rho2[1],                 /* rho >= rho2  -> grow                                 */
    double rho1[1],                 /* rho1<=rho<rho2 -> hold                               */
    double rho0[1],                 /* rho < rho0   -> reject the step                      */
    int    max_stall[1],           /* consecutive rejections at the eta floor before giving up */
    int    max_iter[1],            /* SCvx iteration cap                                   */
    int    size_n[1],               /* RK4 sub-steps per discretization interval (>=2)      */
    int    verbose[1],              /* 0 silent, 1 one-line summary, 2 full per-iteration log */

    /* ---------------- final time / cost ---------------- */
    double tf_fixed[1],             /* used only if COGU_DILATION_ENABLED==0                */
    double lam[1],                  /* defect / constraint-violation penalty weight in J    */

    /* ---------------- boundary conditions, length COGU_NX ---------------- */
    const double start[COGU_NX],
    const double end[COGU_NX],

    /* ---------------- physical model parameters, length COGU_NP_DYN ---------------- */
    const double dp[COGU_NP_DYN],

    /* ---------------- state scaling: xs = Sx^-1 (x - cx) ---------------- */
    const double Sx[COGU_NX * COGU_NX],  /* row-major */
    const double cx[COGU_NX],

    /* ---------------- control scaling: us = Su^-1 (u - cu) ---------------- */
    const double Su[COGU_NU * COGU_NU],  /* row-major */
    const double cu[COGU_NU],

    /* ---------------- warm start ---------------- */
    int warm_mode[1],               /* a COGU_WarmMode value */
    const double warm_x_in[COGU_NX * (COGU_T + 1)], /* time-major; read iff warm_mode[0]==COGU_WARM_USER,
                                                      * else ignored -- pass zeros */
    const double warm_u_in[COGU_NU * COGU_UNODES],  /* time-major; read iff warm_mode[0]==COGU_WARM_USER,
                                                      * else ignored -- pass zeros */
    int quat_offset[1],             /* index of a scalar-first quaternion block in the state, or -1 */
    int omega_offset[1],            /* index of a 3-vector to receive a finite-difference rate, or -1 */
    const double omega_transform[9],/* 3x3 row-major applied to that rate; pass identity (or zeros,
                                      * if unused) -- ignored whenever quat_offset[0]/omega_offset[0] < 0 */

    /* ---------------- time dilation (ignored if COGU_DILATION_ENABLED==0) ---------------- */
    double alpha_p[1],              /* weight of the convex terminal cost alpha_p*tf        */
    double tf_min[1],
    double tf_max[1],
    double tf_init[1],              /* REQUIRED initial guess for tf; also sets the warm-start time step */
    double Sp[1],                   /* scalar: hp = Sp*p + cp */
    double cp[1],

    /* ---------------- outputs (all caller-allocated) ---------------- */
    double x_out[COGU_NX * (COGU_T + 1)], /* time-major, physical units   [REQUIRED] */
    double u_out[COGU_NU * COGU_UNODES],  /* time-major, physical units  [REQUIRED] */
    double tf_out[1],               /* optimal tf, or tf_fixed echoed back              [REQUIRED] */
    int    converged_out[1],       /* 1 converged, 0 otherwise                         [REQUIRED] */
    int    iterations_out[1],      /* iterations actually run                          [REQUIRED] */
    int    status_out[1],           /* a COGU_Status value */
    double J_out[1],                /* nonlinear augmented cost of the accepted solution */
    double L_out[1],                /* convex objective of the last solve */
    double rho_out[1],              /* rho of the last iteration */
    double disc_ms_out[1],          /* total discretization wall time [ms] */
    double solve_ms_out[1],         /* total convex-solve wall time [ms] */
    double hist_tf[COGU_HIST_CAP],  /* per-iteration tf, first hist_len[0] entries valid */
    double hist_J[COGU_HIST_CAP],   /* per-iteration J,  first hist_len[0] entries valid */
    double hist_rho[COGU_HIST_CAP], /* per-iteration rho (entry 0 unused), first hist_len[0] entries valid */
    int    hist_len[1]              /* number of entries actually written to hist_* (<= COGU_HIST_CAP) */
);

/*
 * All "outputs" above are formally required (the Simulink C Caller /
 * Legacy Code Tool convention: every port is wired to something). The
 * underlying implementation additionally tolerates a literal NULL for
 * status_out, J_out, L_out, rho_out, disc_ms_out, solve_ms_out, hist_tf,
 * hist_J, hist_rho and hist_len (skips writing that one) for callers outside
 * Simulink -- e.g. cogu_bindings.py -- that do not want the extra
 * diagnostics on every call. x_out, u_out, tf_out, converged_out and
 * iterations_out are never NULL-checked and must always point at real
 * caller-owned buffers.
 */

#ifdef __cplusplus
}}
#endif

#endif /* COGU_SCVX_H */
"""
    path = os.path.join(output_dir, 'cogu_scvx.h')
    with open(path, 'w', encoding='utf-8') as fp:
        fp.write(content)
    return path


# ======================================================================
# COGU_SCVx_C_caller.h -- Simulink C Caller restatement of cogu_scvx.h
# ======================================================================

def generate_cogu_c_caller_header(nx, nu, ng, n_dyn_par, T, u_nodes,
                                  hist_cap=COGU_HIST_CAP, output_dir=None):
    """
    Emit COGU_SCVx_C_caller.h: the same COGU_SCVx_main prototype as
    cogu_scvx.h, restated with literal numeric array dimensions instead of
    COGU_NX/COGU_T/...-derived expressions, no `const`, and no COGU_API
    decoration -- the form Simulink's C Caller block (Legacy Code Tool)
    parses most reliably. Array-parameter types decay to pointer types in C,
    so this prototype is ABI-compatible with the real compiled
    COGU_SCVx_main despite being spelled differently; both are legal,
    type-compatible declarations of the same external-linkage function.

    Link the STATIC library (the `{lib_name}_static` CMake target -- see
    generate_cmakelists) against a model that includes this header, not the
    shared library/DLL, which is for the Python ctypes callers instead.
    """
    nx1 = nx * (T + 1)
    nu_nodes = nu * u_nodes
    content = f"""\
/*
 * COGU_SCVx_C_caller.h
 *
 * Generated by scvx_toolbox.cgen (generate_cogu_c_caller_header) -- matches
 * this instance's scvx_main.c / cogu_scvx.h; do not hand-edit either without
 * keeping them in sync (re-run the generator instead).
 *
 * Simulink-oriented restatement of the COGU_SCVx_main() prototype declared
 * in cogu_scvx.h: every argument here uses a literal, fixed array dimension
 * -- `type name[1]` for what is conceptually a scalar, `type name[N]` for
 * everything else -- and nothing is a bare scalar or a `type *name` pointer.
 * Simulink's C Caller block needs an explicit array dimension on every
 * argument to map it to a signal port; this is the same function
 * cogu_scvx.h declares (array-parameter types decay to pointer types in C,
 * so the two prototypes are ABI-compatible), just spelled out with literal
 * numbers instead of COGU_NX/COGU_T/...-derived expressions, since that is
 * the form a header parser reads most reliably. Every parameter -- every
 * problem input, every diagnostic output -- must be wired to something: the
 * function itself keeps no state between calls and reads/writes nothing it
 * is not explicitly given here.
 *
 * Usage in Simulink:
 *   1. Add a "C Caller" block and point it at this header.
 *   2. Link against the STATIC library (cogu_scvx_static.lib on Windows /
 *      libcogu_scvx_static.a elsewhere), built by the "cogu_scvx_static"
 *      CMake target -- NOT the shared library/DLL, which is for the Python
 *      ctypes callers (COGU_C_caller_test.py, VAT_COGU_C_caller_test.py,
 *      cogu_bindings.py) instead.
 *   3. Wire every input/output below. A run that does not use a given
 *      input still needs it driven by something -- e.g. warm_x_in/
 *      warm_u_in are only read when warm_mode[0]==0 (COGU_WARM_USER), but
 *      the ports still need a signal (zeros are fine) the rest of the time.
 *
 * warm_mode[0] values (COGU_WarmMode in cogu_scvx.h):
 *     0 = COGU_WARM_USER   -- caller supplies warm_x_in/warm_u_in verbatim
 *     1 = COGU_WARM_LINEAR -- linear interpolation start -> end, built in C
 *     2 = COGU_WARM_SLERP  -- SLERP + finite-difference rate, built in C
 * status_out[0] values (COGU_Status in cogu_scvx.h):
 *     0 = COGU_STATUS_CONVERGED, 1 = COGU_STATUS_MAX_ITER, 2 = COGU_STATUS_STALLED
 *
 * Fixed dimensions baked into THIS instance (see cogu_scvx.h's COGU_* enum
 * -- regenerate both headers if any of these change):
 *     NX = {nx}, NU = {nu}, NG = {ng}, NP_DYN = {n_dyn_par}, T = {T}, CONTROL_NODES = {u_nodes},
 *     HIST_CAP = {hist_cap}
 *   so: start/end/cx = [{nx}], dp = [{n_dyn_par}], Sx = [{nx * nx}] ({nx}*{nx}), Su = [{nu * nu}] ({nu}*{nu}),
 *       cu = [{nu}], warm_x_in/x_out = [{nx1}] ({nx}*{T + 1}), warm_u_in/u_out = [{nu_nodes}]
 *       ({nu}*{u_nodes}), omega_transform = [9] (3*3), hist_tf/hist_J/hist_rho = [{hist_cap}].
 */

#ifndef COGU_SCVX_C_CALLER_H_
#define COGU_SCVX_C_CALLER_H_

#ifdef __cplusplus
extern "C" {{
#endif

extern void COGU_SCVx_main(
    /* ---- SCvx algorithm tuning ---- */
    int    defect_penalty_scaled[1],
    double epsilon_stop[1],
    double e_tol[1],
    double beta_gr[1],
    double beta_sh[1],
    double eta_init[1],
    double eta1[1],
    double eta0[1],
    double rho2[1],
    double rho1[1],
    double rho0[1],
    int    max_stall[1],
    int    max_iter[1],
    int    size_n[1],
    int    verbose[1],

    /* ---- final time / cost ---- */
    double tf_fixed[1],
    double lam[1],

    /* ---- boundary conditions ---- */
    double start[{nx}],
    double end[{nx}],

    /* ---- physical model parameters ---- */
    double dp[{n_dyn_par}],

    /* ---- state scaling: xs = Sx^-1 (x - cx) ---- */
    double Sx[{nx * nx}],
    double cx[{nx}],

    /* ---- control scaling: us = Su^-1 (u - cu) ---- */
    double Su[{nu * nu}],
    double cu[{nu}],

    /* ---- warm start ---- */
    int    warm_mode[1],
    double warm_x_in[{nx1}],
    double warm_u_in[{nu_nodes}],
    int    quat_offset[1],
    int    omega_offset[1],
    double omega_transform[9],

    /* ---- time dilation ---- */
    double alpha_p[1],
    double tf_min[1],
    double tf_max[1],
    double tf_init[1],
    double Sp[1],
    double cp[1],

    /* Outputs */
    double x_out[{nx1}],
    double u_out[{nu_nodes}],
    double tf_out[1],
    int    converged_out[1],
    int    iterations_out[1],
    int    status_out[1],
    double J_out[1],
    double L_out[1],
    double rho_out[1],
    double disc_ms_out[1],
    double solve_ms_out[1],
    double hist_tf[{hist_cap}],
    double hist_J[{hist_cap}],
    double hist_rho[{hist_cap}],
    int    hist_len[1]);

#ifdef __cplusplus
}}
#endif

#endif /* COGU_SCVX_C_CALLER_H_ */
"""
    path = os.path.join(output_dir, 'COGU_SCVx_C_caller.h')
    with open(path, 'w', encoding='utf-8') as fp:
        fp.write(content)
    return path


# ======================================================================
# cpg_compat.h
# ======================================================================

def generate_cpg_compat(nx, nu, ng, T, u_nodes, hold, dilation, output_path,
                        pin_aware_trust_region=True):
    """
    Bridge a stacked-index update API onto CVXPYgen's per-timestep setters.

    CVXPYgen emits `cpg_update_A_0(i, v) ... cpg_update_A_{T-1}(i, v)`; a loop
    written against a single stacked matrix wants one entry point. This header
    supplies that with a switch on the timestep.

    Two CVXPYgen behaviours have to be matched exactly or the result will not
    compile:

    1. A parameter whose total size is 1 gets the scalar signature
       `cpg_update_NAME(cpg_float val)` -- no index argument. So `ou_k` is
       indexed when nu > 1 but scalar when nu == 1, and `op` (always 1x1) is
       always scalar.

    2. CVXPYgen prunes parameters that appear in no constraint. With
       `pin_aware_trust_region` the states at k=0 and k=T are pinned by the
       boundary equalities and dropped from the trust region, so `ox_0` and
       `ox_T` are never referenced and CVXPYgen does not emit setters for
       them. Those switch cases become no-ops here.
    """
    def switch_block(outer, inner_prefix, stride, n_cases, skip=()):
        """`stride`==1 means the per-step setter takes a scalar, not an index."""
        scalar = (stride == 1)
        lines = [f'static inline void {outer}(int idx, double val) {{']
        if n_cases == 0:
            lines += ['    (void)idx; (void)val;', '}', '']
            return '\n'.join(lines)
        lines.append(f'    int k = idx / {stride};')
        if not scalar:
            lines.append(f'    int inner = idx % {stride};')
        lines.append('    switch (k) {')
        for k in range(n_cases):
            if k in skip:
                lines.append(f'        case {k}: break;  /* pruned: unused parameter */')
            elif scalar:
                lines.append(f'        case {k}: {inner_prefix}_{k}(val); break;')
            else:
                lines.append(f'        case {k}: {inner_prefix}_{k}(inner, val); break;')
        lines += ['        default: break;', '    }', '}', '']
        return '\n'.join(lines)

    ox_skip = (0, T) if pin_aware_trust_region else ()

    blocks = [f"""\
/* cpg_compat.h -- stacked update API over CVXPYgen's per-timestep setters.
 * Generated by scvx_toolbox. Do not edit.
 *
 * nx={nx}, nu={nu}, ng={ng}, T={T}, control nodes={u_nodes},
 * hold={hold}, time dilation={'on' if dilation else 'off'},
 * pin-aware trust region={'on' if pin_aware_trust_region else 'off'}
 */
#ifndef CPG_COMPAT_H
#define CPG_COMPAT_H

#include "cpg_solve.h"
"""]

    blocks.append(switch_block('cpg_update_A_discrete', 'cpg_update_A', nx * nx, T))
    blocks.append(switch_block('cpg_update_Bm_discrete', 'cpg_update_Bm', nx * nu, T))
    if hold == FOH:
        blocks.append(switch_block('cpg_update_Bp_discrete', 'cpg_update_Bp', nx * nu, T))
    if dilation:
        blocks.append(switch_block('cpg_update_F_discrete', 'cpg_update_F', nx, T))
    blocks.append(switch_block('cpg_update_y_discrete', 'cpg_update_y', nx, T))

    if ng > 0:
        blocks.append(switch_block('cpg_update_C_discrete', 'cpg_update_C', ng * nx, T + 1))
        blocks.append(switch_block('cpg_update_D_discrete', 'cpg_update_D', ng * nu, T + 1))
        if dilation:
            blocks.append(switch_block('cpg_update_G_discrete', 'cpg_update_G', ng, T + 1))
        blocks.append(switch_block('cpg_update_z_discrete', 'cpg_update_z', ng, T + 1))

    blocks.append(switch_block('cpg_update_ox_cvxpy', 'cpg_update_ox', nx, T + 1,
                               skip=ox_skip))
    blocks.append(switch_block('cpg_update_ou_cvxpy', 'cpg_update_ou', nu, u_nodes))
    blocks.append('#endif /* CPG_COMPAT_H */\n')

    content = '\n'.join(blocks)
    with open(output_path, 'w', encoding='utf-8') as fp:
        fp.write(content)
    return content


# ======================================================================
# CMake
# ======================================================================


# ======================================================================
# cogu_bindings.py -- ctypes bridge from Python to the compiled library
# ======================================================================

BINDINGS_TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), 'templates',
                                      'cogu_bindings_template.py.in')

# Substituted with str.replace, NOT str.format: the template is real Python
# and contains f-strings ({_HERE}, {n}, {nx}, ...) that .format would try to
# interpret. Only these exact tokens are touched.
BINDINGS_PLACEHOLDERS = (
    '{NX}', '{NU}', '{NG}', '{NP_DYN}', '{T}', '{U_NODES}',
    '{DILATION}', '{FOH}', '{HIST_CAP}', '{OUTPUT_DIR}',
    '{SIZE_N}', '{MAX_ITER}', '{MAX_STALL}',
    '{RHO0}', '{RHO1}', '{RHO2}',
    '{ETA0}', '{ETA1}', '{ETA_INIT}',
    '{BETA_SH}', '{BETA_GR}', '{E_TOL}', '{EPSILON_STOP}',
    '{DEFECT_PENALTY_SCALED}', '{TF_FIXED}', '{LAM}',
    '{ALPHA_P}', '{TF_MIN}', '{TF_MAX}', '{TF_INIT}', '{SP}', '{CP}',
)


def generate_python_bindings(nx, nu, ng, n_dyn_par, T, u_nodes, dilation, foh,
                             embed, output_dir, hist_cap=COGU_HIST_CAP):
    """
    Emit cogu_bindings.py: the ctypes wrapper that lets Python call the
    compiled COGU_SCVx_main.

    This has to be generated rather than hand-written because its constants
    and its ctypes signature must match the cogu_scvx.h emitted alongside it.
    A hand-maintained copy goes out of sync silently, and the failure mode is
    memory corruption at call time rather than an import error.

    The template keeps `load_library`, `cogu_scvx_main` and `COGU_WarmMode`
    with unchanged signatures: VAT_COGU_C_caller_test.py imports those three
    names, and that file is user-owned.

    Returns the path written, or None if the template is missing.
    """
    if not os.path.isfile(BINDINGS_TEMPLATE_PATH):
        return None

    with open(BINDINGS_TEMPLATE_PATH, 'r', encoding='utf-8') as fp:
        content = fp.read()

    sp = embed.get('solver_params', {})

    def _num(value, fallback):
        """Repr of a scalar, falling back when the problem does not define it.

        The dilation knobs (alpha_p, tf_min/max/init, S_p, c_p) are None under
        a fixed-time mode.  COGU_SCVx_main ignores them there, but they still
        need a value the generated Python can carry.
        """
        if value is None:
            return repr(fallback)
        arr = np.asarray(value, dtype=float).reshape(-1)
        return repr(float(arr[0]))

    tf_fixed = float(embed['tf'])
    subs = {
        '{NX}': str(nx), '{NU}': str(nu), '{NG}': str(ng),
        '{NP_DYN}': str(n_dyn_par), '{T}': str(T), '{U_NODES}': str(u_nodes),
        '{DILATION}': '1' if dilation else '0',
        '{FOH}': '1' if foh else '0',
        '{HIST_CAP}': str(int(hist_cap)),
        '{OUTPUT_DIR}': os.path.basename(os.path.normpath(output_dir)),

        '{SIZE_N}': str(int(sp.get('SIZE_N', 20))),
        '{MAX_ITER}': str(int(sp.get('MAX_ITER', 25))),
        '{MAX_STALL}': str(int(sp.get('MAX_STALL', 3))),
        '{RHO0}': repr(float(sp.get('RHO0', 0.0))),
        '{RHO1}': repr(float(sp.get('RHO1', 0.1))),
        '{RHO2}': repr(float(sp.get('RHO2', 0.7))),
        '{ETA0}': repr(float(sp.get('ETA0', 1e-8))),
        '{ETA1}': repr(float(sp.get('ETA1', 10.0))),
        '{ETA_INIT}': repr(float(sp.get('ETA_INIT', 1.0))),
        '{BETA_SH}': repr(float(sp.get('BETA_SH', 2.0))),
        '{BETA_GR}': repr(float(sp.get('BETA_GR', 2.0))),
        '{E_TOL}': repr(float(sp.get('E_TOL', 0.005))),
        '{EPSILON_STOP}': repr(float(sp.get('EPSILON_STOP', 1e-5))),
        '{DEFECT_PENALTY_SCALED}':
            '1' if sp.get('DEFECT_PENALTY_SCALED') else '0',

        '{TF_FIXED}': repr(tf_fixed),
        '{LAM}': repr(float(embed['lam'])),

        '{ALPHA_P}': _num(embed.get('alpha_p'), 0.0),
        '{TF_MIN}': _num(embed.get('tf_min'), tf_fixed),
        '{TF_MAX}': _num(embed.get('tf_max'), tf_fixed),
        '{TF_INIT}': _num(embed.get('tf_init'), tf_fixed),
        '{SP}': _num(embed.get('S_p'), 1.0),
        '{CP}': _num(embed.get('c_p'), 0.0),
    }

    # Longest token first, so replacing '{T}' can never chew into '{TF_MIN}'.
    for key in sorted(subs, key=len, reverse=True):
        content = content.replace(key, subs[key])

    path = os.path.join(output_dir, 'cogu_bindings.py')
    with open(path, 'w', encoding='utf-8') as fp:
        fp.write(content)
    return path

# ======================================================================
# scvx_cli_demo.c -- standalone smoke test for the generated library
# ======================================================================

_CLI_DEMO = r"""/*
 * scvx_cli_demo.c -- optional standalone smoke test for COGU_SCVx_main.
 * Generated by scvx_toolbox; safe to delete.
 *
 * NOT part of the public API (that is scvx_main.c / cogu_scvx.h). This is a
 * small, disposable driver with its own main() so the compiled library can be
 * sanity-checked from the command line, without Python, on the exact problem
 * instance the Python solve_trajectory(...) call generated -- every number
 * comes from embedded_data.h and every dimension from cogu_scvx.h, so this
 * same text serves any model and any time mode.
 *
 * Build (the "scvx_cli_demo" target in CMakeLists.txt):
 *     cmake -B build -DCMAKE_C_FLAGS="-O3"
 *     cmake --build build --config Release
 *     ./build/scvx_cli_demo          (or build/Release/scvx_cli_demo.exe)
 */
#include <stdio.h>
#include "cogu_scvx.h"
#include "embedded_data.h"

int main(void) {
    double x_out[COGU_NX * (COGU_T + 1)];
    double u_out[COGU_NU * COGU_UNODES];
    double tf_out, J_out, L_out, rho_out, disc_ms, solve_ms;
    int converged, iterations, status;

    /* Every COGU_SCVx_main scalar is an array of 1 (see cogu_scvx.h) -- a
     * Simulink C Caller needs a fixed array dimension on every port, so the
     * header never declares a bare scalar or a plain pointer. Length-1 locals
     * here are the plain-C equivalent of that same convention. */
    int    defect_penalty_scaled[1] = { DEFECT_PENALTY_SCALED };
    double epsilon_stop[1] = { EPSILON_STOP };
    double e_tol[1]        = { E_TOL };
    double beta_gr[1]      = { BETA_GR };
    double beta_sh[1]      = { BETA_SH };
    double eta_init[1]     = { ETA_INIT };
    double eta1[1]         = { ETA1 };
    double eta0[1]         = { ETA0 };
    double rho2[1]         = { RHO2 };
    double rho1[1]         = { RHO1 };
    double rho0[1]         = { RHO0 };
    int    max_stall[1]    = { MAX_STALL };
    int    max_iter[1]     = { MAX_ITER };
    int    size_n[1]       = { SIZE_N };
    int    verbose[1]      = { 2 /* full per-iteration log */ };
    double tf_fixed_a[1]   = { tf_fixed };
    double lam_a[1]        = { lam };
    int    warm_mode[1]    = { COGU_WARM_USER };
    int    quat_offset[1]  = { -1 /* unused for COGU_WARM_USER */ };
    int    omega_offset[1] = { -1 /* unused for COGU_WARM_USER */ };
    double omega_transform[9] = { 0 };  /* unused for COGU_WARM_USER */

    /* embedded_data.h defines the dilation constants only when the problem
     * really has a free final time. Under a fixed-time mode COGU_SCVx_main
     * ignores these six arguments, but they still have to exist to compile. */
#if TIME_DILATION
    double alpha_p_a[1] = { alpha_p };
    double tf_min_a[1]  = { tf_min };
    double tf_max_a[1]  = { tf_max };
    double tf_init_a[1] = { tf_init };
    double Sp_a[1]      = { _Sp[0] };
    double cp_a[1]      = { _cp[0] };
#else
    double alpha_p_a[1] = { 0.0 };
    double tf_min_a[1]  = { tf_fixed };
    double tf_max_a[1]  = { tf_fixed };
    double tf_init_a[1] = { tf_fixed };
    double Sp_a[1]      = { 1.0 };
    double cp_a[1]      = { 0.0 };
#endif

    COGU_SCVx_main(
        /* --- SCvx tuning, from embedded_data.h --- */
        defect_penalty_scaled, epsilon_stop, e_tol, beta_gr, beta_sh,
        eta_init, eta1, eta0, rho2, rho1, rho0, max_stall, max_iter, size_n,
        verbose,
        /* --- final time / cost --- */
        tf_fixed_a, lam_a,
        /* --- boundary conditions / physical parameters --- */
        _start, _end, _dp,
        /* --- scaling --- */
        _Sx, _cx, _Su, _cu,
        /* --- warm start: reproduce the original run exactly --- */
        warm_mode, _warm_x, _warm_u, quat_offset, omega_offset, omega_transform,
        /* --- time dilation --- */
        alpha_p_a, tf_min_a, tf_max_a, tf_init_a, Sp_a, cp_a,
        /* --- outputs --- */
        x_out, u_out, &tf_out, &converged, &iterations, &status,
        &J_out, &L_out, &rho_out, &disc_ms, &solve_ms,
        NULL, NULL, NULL, NULL);

    printf("\n--- scvx_cli_demo summary ---\n");
    printf("converged=%d iterations=%d status=%d tf=%.6f J=%.6f L=%.6f\n",
           converged, iterations, status, tf_out, J_out, L_out);
    printf("initial state:");
    for (int i = 0; i < COGU_NX; i++) printf("  %.6f", x_out[i]);
    printf("\nfinal state:  ");
    for (int i = 0; i < COGU_NX; i++) printf("  %.6f", x_out[COGU_T * COGU_NX + i]);
    printf("\n");
    return 0;
}
"""


def generate_cli_demo(output_dir):
    """Write scvx_cli_demo.c -- a main() that reruns the generated problem.

    The driver is problem-independent: every number comes out of
    embedded_data.h and every dimension out of cogu_scvx.h, so the same
    text serves any model and any time mode.
    """
    path = os.path.join(output_dir, 'scvx_cli_demo.c')
    with open(path, 'w', encoding='utf-8') as fp:
        fp.write(_CLI_DEMO)
    return path


def generate_cmakelists(output_dir, lib_name='cogu_scvx', scvx_main_present=True):
    """
    Build `lib_name` as a SHARED library (.dll on Windows, .so on Linux/macOS)
    exporting COGU_SCVx_main -- loadable via ctypes from Python (Monte Carlo,
    NMPC). `dynamics.c` is always compiled in; `scvx_main.c` is added whenever
    the loop template was found (scvx_main_present) -- without it there is no
    COGU_SCVx_main to export, so the target degrades to just the model +
    solver for inspection.

    Also always builds `{lib_name}_static` as a STATIC library (.lib/.a) of
    the exact same sources -- that is what a Simulink C Caller block links
    against (see COGU_SCVx_C_caller.h), since it needs no separate .dll/.so
    sitting next to the built model at runtime. Skipped when scvx_main_present
    is False, same as the shared target, since there is nothing to export
    either way.
    """
    sources = 'scvx_main.c\n    dynamics.c' if scvx_main_present else 'dynamics.c'
    static_block = ''
    demo_block = ''
    if scvx_main_present:
        static_block = f"""

# Static library, for Simulink's C Caller block (see COGU_SCVx_C_caller.h):
# a static .lib/.a is what the C Caller block links against, no separate
# .dll/.so needs to sit next to the built model at runtime. Not used by the
# Python callers -- those load the SHARED library above via ctypes, which
# needs a .dll/.so it can dlopen() at runtime. COGU_API expands to nothing
# here (COGU_BUILD_DLL is not defined for this target), which is correct: a
# static archive has no import/export boundary to decorate.
add_library({lib_name}_static STATIC
    {sources}
    ${{cpg_src}}
)
set_target_properties({lib_name}_static PROPERTIES OUTPUT_NAME {lib_name}_static)
target_include_directories({lib_name}_static PRIVATE
    ${{CMAKE_CURRENT_SOURCE_DIR}}
    ${{cpg_include}}
)
if(NOT MSVC)
    target_link_libraries({lib_name}_static PRIVATE m)
endif()
"""
        demo_block = f"""

# Standalone smoke test: reruns the exact problem instance this project was
# generated from, using the constants in embedded_data.h. Not part of the
# public API -- delete scvx_cli_demo.c and this block if you do not want it.
add_executable(scvx_cli_demo scvx_cli_demo.c)
target_link_libraries(scvx_cli_demo PRIVATE {lib_name}_static)
target_include_directories(scvx_cli_demo PRIVATE
    ${{CMAKE_CURRENT_SOURCE_DIR}}
    ${{cpg_include}}
)
if(NOT MSVC)
    target_link_libraries(scvx_cli_demo PRIVATE m)
endif()
"""
    content = f"""\
cmake_minimum_required(VERSION 3.10)
project({lib_name} C)

# Cross compiling: pass -DCMAKE_TOOLCHAIN_FILE=<your-toolchain>.cmake
# e.g. cmake -B build -DCMAKE_TOOLCHAIN_FILE=arm-none-eabi.cmake

# CVXPYgen solver; exports cpg_include and cpg_src to the parent scope
add_subdirectory(solver/c)

add_library({lib_name} SHARED
    {sources}
    ${{cpg_src}}
)
# COGU_API in cogu_scvx.h expands to __declspec(dllexport) only when this is
# defined, so a *consumer* of the .dll (Python via ctypes, a Simulink S-
# function) does not need to define COGU_USE_DLL to link against it.
target_compile_definitions({lib_name} PRIVATE COGU_BUILD_DLL)

target_include_directories({lib_name} PRIVATE
    ${{CMAKE_CURRENT_SOURCE_DIR}}
    ${{cpg_include}}
)

if(NOT MSVC)
    target_link_libraries({lib_name} PRIVATE m)
endif()
{static_block}{demo_block}"""
    with open(os.path.join(output_dir, 'CMakeLists.txt'), 'w', encoding='utf-8') as fp:
        fp.write(content)


# ======================================================================
# CVXPYgen solver, generated in a clean subprocess
# ======================================================================

def generate_solver_subprocess(build_kwargs, solver_dir, solver='ECOS',
                               extra_path=None):
    """
    Run CVXPYgen in a child process.

    Canonicalizing a DPP problem of this size needs a few hundred MB of
    *contiguous* memory.  In the parent process SymPy and the SCvx solve have
    already allocated and fragmented the heap, which is a reliable way to get
    a MemoryError.  A fresh child starts with the whole block available.
    """
    package_dir = os.path.dirname(os.path.abspath(__file__))
    package_root = os.path.dirname(package_dir)
    worker = os.path.join(package_dir, '_cpg_worker.py')

    payload = {
        'package_root': package_root,
        'build_kwargs': build_kwargs,
        'solver_dir': solver_dir,
        'solver': solver,
        'extra_path': extra_path,
    }
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False, mode='wb') as fp:
        pickle.dump(payload, fp)
        payload_path = fp.name

    try:
        proc = subprocess.run([sys.executable, worker, payload_path],
                              capture_output=True, text=True)
        cpg_solve_c = os.path.join(solver_dir, 'c', 'src', 'cpg_solve.c')
        if not os.path.isfile(cpg_solve_c):
            raise RuntimeError(
                "CVXPYgen did not produce the C solver.\n"
                f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
            )
    finally:
        try:
            os.unlink(payload_path)
        except OSError:
            pass


# ======================================================================
# Project assembly
# ======================================================================

def fill_scvx_template(nx, nu, ng, n_p, T, u_nodes, hold, dilation,
                       output_dir, cost_terms, embed, template_path=None):
    """
    Fill the SCvx loop template and write scvx_main.c.

    Returns True if the template was found and filled, False otherwise.
    """
    template_path = template_path or TEMPLATE_PATH
    if not os.path.isfile(template_path):
        return False

    with open(template_path, 'r', encoding='utf-8') as fp:
        content = fp.read()

    substitutions = {
        '{USER_COST}': generate_c_user_cost(cost_terms, nx, nu, T, u_nodes),
        '{NX}': str(nx),
        '{NU}': str(nu),
        '{NG}': str(ng),
        '{NP}': str(n_p),
        '{T}': str(T),
        '{U_NODES}': str(u_nodes),
        '{DILATION}': '1' if dilation else '0',
        '{FOH}': '1' if hold == FOH else '0',
        '{TIME_MODE}': str(embed.get('time_mode', '')),
    }
    for key, value in substitutions.items():
        content = content.replace(key, value)

    with open(os.path.join(output_dir, 'scvx_main.c'), 'w', encoding='utf-8') as fp:
        fp.write(content)
    return True


def emit_c_project(model, prob_dict, build_kwargs, output_dir, embed,
                   solver_source_dir=None, lib_name='cogu_scvx', verbose=True):
    """Generate the whole C project into `output_dir`."""
    os.makedirs(output_dir, exist_ok=True)

    nx, nu, ng, T = prob_dict['nx'], prob_dict['nu'], prob_dict['ng'], prob_dict['T']
    u_nodes = prob_dict['u_nodes']
    hold, dilation = prob_dict['hold'], prob_dict['dilation']
    _, _, _ = parse_time_mode(prob_dict['time_mode'])

    # 1. Convex solver
    solver_dir = os.path.join(output_dir, 'solver')
    cpg_solve_c = os.path.join(solver_dir, 'c', 'src', 'cpg_solve.c')

    if solver_source_dir is not None:
        if os.path.isdir(solver_dir):
            shutil.rmtree(solver_dir)
        shutil.copytree(solver_source_dir, solver_dir)
        if verbose:
            print(f"[cgen] solver copied from {solver_source_dir}")
    elif os.path.isfile(cpg_solve_c):
        if verbose:
            print(f"[cgen] solver already present in {solver_dir}/c -- reusing")
    else:
        os.makedirs(solver_dir, exist_ok=True)
        if verbose:
            print("[cgen] generating the C solver with CVXPYgen (clean subprocess)...")
        generate_solver_subprocess(build_kwargs, solver_dir)
        if verbose:
            print(f"[cgen] solver written to {solver_dir}/c")

    # 2. Model
    generate_c_dynamics(model, output_dir=output_dir)

    # 3. Headers
    generate_cpg_compat(nx, nu, ng, T, u_nodes, hold, dilation,
                        os.path.join(output_dir, 'cpg_compat.h'),
                        pin_aware_trust_region=prob_dict.get(
                            'pin_aware_trust_region', True))
    write_embedded_header(embed, nx, nu, ng, T, u_nodes, dilation, output_dir)

    # 4. SCvx loop + its public C API
    # {NP} sizes the dynamic-parameter array dp[], so it is the number of
    # model constants -- not the number of dilation parameters.
    n_dyn_par = int(np.asarray(embed['dp'], dtype=float).reshape(-1).size)
    generate_cogu_header(nx, nu, ng, n_dyn_par, T, u_nodes, dilation, hold == FOH,
                         output_dir)
    generate_cogu_c_caller_header(nx, nu, ng, n_dyn_par, T, u_nodes,
                                  output_dir=output_dir)
    scvx_main_present = fill_scvx_template(nx, nu, ng, n_dyn_par, T, u_nodes, hold,
                                           dilation, output_dir,
                                           prob_dict['cost_terms'], embed)
    if verbose and not scvx_main_present:
        print(f"[cgen] no loop template at {TEMPLATE_PATH} -- scvx_main.c not "
              f"generated, so there is no COGU_SCVx_main to export. Everything "
              f"else (solver, dynamics.c/h, cpg_compat.h, cogu_scvx.h, "
              f"COGU_SCVx_C_caller.h, embedded_data.h) is in place; drop your "
              f"template there and re-run to emit the loop.")

    # 5. Build files
    if scvx_main_present:
        generate_cli_demo(output_dir)
        generate_python_bindings(nx, nu, ng, n_dyn_par, T, u_nodes, dilation,
                                 hold == FOH, embed, output_dir)
    generate_cmakelists(output_dir, lib_name=lib_name,
                        scvx_main_present=scvx_main_present)

    if verbose:
        print(f"[cgen] C project written to {os.path.abspath(output_dir)}")
