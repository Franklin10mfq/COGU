/*
 * scvx_cli_demo.c -- optional standalone smoke test for COGU_SCVx_main.
 *
 * NOT part of the public API (that is scvx_main.c / cogu_scvx.h). This is a
 * small, disposable driver with its own main() so the compiled library can
 * be sanity-checked from the command line, without Python, using the exact
 * problem instance the Python solve_trajectory(...) call generated (see
 * embedded_data.h). Reproduces what the old hardcoded main() used to print.
 *
 * Build (added to CMakeLists.txt as the "scvx_cli_demo" target):
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
     * header never declares a bare scalar or a plain pointer. Length-1
     * locals here are the plain-C equivalent of that same convention. */
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
    double alpha_p_a[1]    = { alpha_p };
    double tf_min_a[1]     = { tf_min };
    double tf_max_a[1]     = { tf_max };
    double tf_init_a[1]    = { tf_init };
    double Sp_a[1]         = { _Sp[0] };
    double cp_a[1]         = { _cp[0] };

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
