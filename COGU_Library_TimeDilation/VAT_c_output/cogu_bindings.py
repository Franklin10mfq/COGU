"""
cogu_bindings.py -- ctypes wrapper around the compiled COGU_SCVx_main shared
library (libcogu_scvx.so / cogu_scvx.dll), matching cogu_scvx.h exactly.

Shared by COGU_C_caller_test.py and VAT_COGU_C_caller_test.py so the ctypes
plumbing (array marshalling, struct-free flat argument list, library search)
lives in exactly one place. COGU_SCVx_main is stateless between calls -- every
input array is read-only and nothing is retained after it returns -- so this
module's `cogu_scvx_main(...)` is safe to call repeatedly and rapidly: once
for a single trajectory, thousands of times for a Monte Carlo sweep, or once
per control step in a receding-horizon / NMPC loop.

Build the library first (from VAT_c_output/):
    cmake -B build -DCMAKE_C_FLAGS="-O3"
    cmake --build build --config Release
which produces build/libcogu_scvx.so (Linux/Mac) or
build/Release/cogu_scvx.dll / build/cogu_scvx.dll (Windows, depending on the
generator).
"""

import ctypes
import glob
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))

# ---- fixed problem shape, mirrors the enum in cogu_scvx.h ------------------
COGU_NX = 7
COGU_NU = 4
COGU_NG = 0
COGU_NP_DYN = 30
COGU_T = 30
COGU_UNODES = 30
COGU_DILATION_ENABLED = 1
COGU_FOH_ENABLED = 0
COGU_HIST_CAP = 64  # fixed capacity of hist_tf/hist_J/hist_rho in the C API;
                    # if max_iter > COGU_HIST_CAP the solve still runs in
                    # full, only the recorded history is truncated


class COGU_WarmMode:
    USER = 0    # caller-supplied warm_x/warm_u, passed through verbatim
    LINEAR = 1  # elementwise linear interpolation start -> end
    SLERP = 2   # SLERP on the quaternion block + finite-difference rate


class COGU_Status:
    CONVERGED = 0
    MAX_ITER = 1
    STALLED = 2
    _NAMES = {0: 'CONVERGED', 1: 'MAX_ITER', 2: 'STALLED'}

    @classmethod
    def name(cls, value):
        return cls._NAMES.get(int(value), f'UNKNOWN({value})')


def _candidate_lib_paths():
    names_by_platform = {
        'linux': ['libcogu_scvx.so'],
        'darwin': ['libcogu_scvx.dylib', 'libcogu_scvx.so'],
        'win32': ['cogu_scvx.dll'],
    }
    names = names_by_platform.get(sys.platform, ['libcogu_scvx.so', 'cogu_scvx.dll', 'libcogu_scvx.dylib'])
    search_dirs = [
        _HERE,
        os.path.join(_HERE, 'build'),
        os.path.join(_HERE, 'build', 'Release'),
        os.path.join(_HERE, 'build_final'),
        os.path.join(_HERE, 'build_final', 'Release'),
    ]
    # anything matching build*/**/<name> under VAT_c_output, in case the user
    # used a different build-directory name than "build"
    for name in names:
        for d in search_dirs:
            yield os.path.join(d, name)
        for p in glob.glob(os.path.join(_HERE, 'build*', '**', name), recursive=True):
            yield p


def find_library():
    for path in _candidate_lib_paths():
        if os.path.isfile(path):
            return path
    return None


def load_library(path=None):
    """Load libcogu_scvx and bind COGU_SCVx_main's exact C signature."""
    if path is None:
        path = find_library()
    if path is None or not os.path.isfile(path):
        raise FileNotFoundError(
            "Could not find the compiled COGU shared library "
            "(libcogu_scvx.so / cogu_scvx.dll). Build it first:\n"
            "    cmake -B build -DCMAKE_C_FLAGS=\"-O3\"   (from VAT_c_output/)\n"
            "    cmake --build build --config Release\n"
            f"Searched under: {_HERE}"
        )
    lib = ctypes.CDLL(path)

    D = ctypes.c_double
    I = ctypes.c_int
    PD = ctypes.POINTER(ctypes.c_double)
    PI = ctypes.POINTER(ctypes.c_int)

    # Every parameter of COGU_SCVx_main is declared as an array in
    # cogu_scvx.h -- `type name[1]` for what used to be a plain scalar,
    # `type name[N]` for what was already a pointer -- so that Simulink's
    # C Caller (which needs an explicit array dimension on every port) can
    # parse the header directly. At the ABI level an array parameter is a
    # pointer regardless of its declared length, so every single argument
    # below is POINTER(c_double) or POINTER(c_int) -- D and I (plain-value
    # ctypes) are unused here as a result, kept only as short names PD/PI
    # are built from.
    lib.COGU_SCVx_main.restype = None
    lib.COGU_SCVx_main.argtypes = [
        # ---- SCvx tuning ----
        PI, PD, PD, PD, PD, PD, PD, PD, PD, PD, PD, PI, PI, PI, PI,
        # ---- final time / cost ----
        PD, PD,
        # ---- boundary conditions / physical params ----
        PD, PD, PD,
        # ---- scaling ----
        PD, PD, PD, PD,
        # ---- warm start ----
        PI, PD, PD, PI, PI, PD,
        # ---- time dilation ----
        PD, PD, PD, PD, PD, PD,
        # ---- outputs ----
        PD, PD, PD, PI, PI, PI,
        PD, PD, PD,
        PD, PD,
        PD, PD, PD, PI,
    ]
    return lib


def _dptr(arr, time_major_2d=False):
    """
    numpy array (or None) -> (ctypes double*, keepalive-buffer), or (None, None).

    cogu_scvx.h's convention is: 2-D *matrices* (Sx, Su, omega_transform) are
    flat row-major, flat[row*ncols+col] -- numpy's default C order already
    matches, no transpose needed. 2-D *trajectories* (warm_x/warm_u, x_out/
    u_out) are "time-major", flat[t*n_components+i] -- but this module's
    public API keeps them shaped (n_components, n_time) throughout (matching
    build_vat_problem_data's warm_x/warm_u and scvx_toolbox's x_opt_unscaled
    convention), so those need their *transpose* flattened in C order.
    Set time_major_2d=True for that case.
    """
    if arr is None:
        return None, None
    if time_major_2d and arr.ndim == 2:
        a = np.ascontiguousarray(arr.T, dtype=np.float64)
    else:
        a = np.ascontiguousarray(arr, dtype=np.float64)
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), a  # keep `a` alive


class CoguResult(dict):
    """Dict-like result of a COGU_SCVx_main call; also exposes attribute access."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e


def cogu_scvx_main(
    lib,
    *,
    # ---- SCvx tuning ----
    defect_penalty_scaled=1, epsilon_stop=1e-3, e_tol=0.005,
    beta_gr=2.0, beta_sh=2.0, eta_init=1.0, eta1=10.0, eta0=0.001,
    rho2=0.7, rho1=0.1, rho0=0.0, max_stall=3, max_iter=25, size_n=10,
    verbose=0,
    # ---- final time / cost ----
    tf_fixed=0.0, lam=5000.0,
    # ---- boundary conditions / physical parameters ----
    start, end, dp,
    # ---- scaling ----
    Sx, cx, Su, cu,
    # ---- warm start ----
    warm_mode=COGU_WarmMode.SLERP, warm_x=None, warm_u=None,
    quat_offset=0, omega_offset=4, omega_transform=None,
    # ---- time dilation ----
    alpha_p=0.05, tf_min=60.0, tf_max=900.0, tf_init=480.0, Sp=None, cp=None,
    # ---- shapes (override only for a differently-generated library) ----
    nx=COGU_NX, nu=COGU_NU, T=COGU_T, u_nodes=COGU_UNODES,
    want_history=True,
):
    """
    Call COGU_SCVx_main once and return a CoguResult with every output.

    `start`, `end`: length-nx arrays. `dp`: length-COGU_NP_DYN. `Sx`: nx*nx
    (row-major) or a 1-D array of length nx (interpreted as the diagonal, the
    common case). `cx`: length nx. `Su`/`cu`: analogous, size nu. `Sp`, `cp`:
    scalars for the dilation scaling (required whenever the library was
    generated with COGU_DILATION_ENABLED, i.e. essentially always here).

    warm_mode=COGU_WarmMode.USER requires `warm_x` (nx, T+1) and `warm_u`
    (nu, u_nodes); the other modes build the warm start in C and both may be
    left None. tf_init is REQUIRED whenever dilation is enabled -- it is also
    the absolute time step used by COGU_WARM_SLERP's finite-difference rate.

    Returns a CoguResult with keys: x (nx, T+1), u (nu, u_nodes), tf,
    converged (bool), iterations, status (int, see COGU_Status), status_name,
    J, L, rho, disc_ms, solve_ms, and -- if want_history -- hist_tf, hist_J,
    hist_rho (each trimmed to the actual hist_len) and hist_len.
    """
    def _square_or_diag(mat, n):
        mat = np.asarray(mat, dtype=np.float64)
        if mat.ndim == 1:
            if mat.shape[0] != n:
                raise ValueError(f"expected length-{n} diagonal, got shape {mat.shape}")
            return np.diag(mat)
        if mat.shape != (n, n):
            raise ValueError(f"expected shape ({n},{n}), got {mat.shape}")
        return mat

    start = np.ascontiguousarray(start, dtype=np.float64).reshape(-1)
    end = np.ascontiguousarray(end, dtype=np.float64).reshape(-1)
    dp = np.ascontiguousarray(dp, dtype=np.float64).reshape(-1)
    if start.shape[0] != nx or end.shape[0] != nx:
        raise ValueError(f"start/end must have length nx={nx}")

    Sx_full = _square_or_diag(Sx, nx)
    Su_full = _square_or_diag(Su, nu)
    cx = np.ascontiguousarray(cx, dtype=np.float64).reshape(-1)
    cu = np.ascontiguousarray(cu, dtype=np.float64).reshape(-1)

    if Sp is None or cp is None:
        raise ValueError("Sp and cp (dilation scaling) are required")
    Sp = float(np.asarray(Sp).reshape(-1)[0])
    cp = float(np.asarray(cp).reshape(-1)[0])

    warm_x_arr = None if warm_x is None else np.ascontiguousarray(warm_x, dtype=np.float64)
    warm_u_arr = None if warm_u is None else np.ascontiguousarray(warm_u, dtype=np.float64)
    omega_transform_arr = (None if omega_transform is None
                            else np.ascontiguousarray(omega_transform, dtype=np.float64))

    keepalive = []

    def dptr(arr, time_major_2d=False):
        ptr, kept = _dptr(arr, time_major_2d=time_major_2d)
        if kept is not None:
            keepalive.append(kept)
        return ptr

    def sref(value, ctype=ctypes.c_double):
        """Python scalar -> pointer to a length-1 ctypes buffer (kept alive
        for the call), matching a COGU_SCVx_main `type name[1]` parameter."""
        buf = ctype(value)
        keepalive.append(buf)
        return ctypes.byref(buf)

    # Time-major output buffers, shape (n_time, n_components) so they are
    # already in the library's flat[t*n+i] layout with no transpose needed;
    # transposed back to this module's (n_components, n_time) convention
    # after the call.
    x_buf = np.zeros((T + 1, nx), dtype=np.float64)
    u_buf = np.zeros((u_nodes, nu), dtype=np.float64)
    tf_out = ctypes.c_double(0.0)
    converged_out = ctypes.c_int(0)
    iterations_out = ctypes.c_int(0)
    status_out = ctypes.c_int(0)
    J_out = ctypes.c_double(0.0)
    L_out = ctypes.c_double(0.0)
    rho_out = ctypes.c_double(0.0)
    disc_ms_out = ctypes.c_double(0.0)
    solve_ms_out = ctypes.c_double(0.0)

    if want_history:
        hist_tf = np.zeros(max_iter, dtype=np.float64)
        hist_J = np.zeros(max_iter, dtype=np.float64)
        hist_rho = np.zeros(max_iter, dtype=np.float64)
        hist_len = ctypes.c_int(0)
        hist_tf_ptr, hist_J_ptr, hist_rho_ptr = dptr(hist_tf), dptr(hist_J), dptr(hist_rho)
        hist_len_ptr = ctypes.byref(hist_len)
    else:
        hist_tf = hist_J = hist_rho = None
        hist_tf_ptr = hist_J_ptr = hist_rho_ptr = None
        hist_len_ptr = None

    CI = ctypes.c_int

    lib.COGU_SCVx_main(
        sref(int(defect_penalty_scaled), CI), sref(epsilon_stop), sref(e_tol),
        sref(beta_gr), sref(beta_sh), sref(eta_init), sref(eta1), sref(eta0),
        sref(rho2), sref(rho1), sref(rho0), sref(int(max_stall), CI), sref(int(max_iter), CI),
        sref(int(size_n), CI), sref(int(verbose), CI),
        sref(tf_fixed), sref(lam),
        dptr(start), dptr(end), dptr(dp),
        dptr(Sx_full), dptr(cx), dptr(Su_full), dptr(cu),
        sref(int(warm_mode), CI), dptr(warm_x_arr, time_major_2d=True), dptr(warm_u_arr, time_major_2d=True),
        sref(int(quat_offset), CI), sref(int(omega_offset), CI), dptr(omega_transform_arr),
        sref(alpha_p), sref(tf_min), sref(tf_max), sref(tf_init), sref(Sp), sref(cp),
        dptr(x_buf), dptr(u_buf),
        ctypes.byref(tf_out), ctypes.byref(converged_out), ctypes.byref(iterations_out),
        ctypes.byref(status_out), ctypes.byref(J_out), ctypes.byref(L_out), ctypes.byref(rho_out),
        ctypes.byref(disc_ms_out), ctypes.byref(solve_ms_out),
        hist_tf_ptr, hist_J_ptr, hist_rho_ptr, hist_len_ptr,
    )

    result = CoguResult(
        x=np.ascontiguousarray(x_buf.T), u=np.ascontiguousarray(u_buf.T),
        tf=tf_out.value,
        converged=bool(converged_out.value),
        iterations=iterations_out.value,
        status=status_out.value,
        status_name=COGU_Status.name(status_out.value),
        J=J_out.value, L=L_out.value, rho=rho_out.value,
        disc_ms=disc_ms_out.value, solve_ms=solve_ms_out.value,
    )
    if want_history:
        n = hist_len.value
        result['hist_len'] = n
        result['hist_tf'] = hist_tf[:n].copy()
        result['hist_J'] = hist_J[:n].copy()
        result['hist_rho'] = hist_rho[:n].copy()
    return result
