// cpg_quat_scvx_cpg.c
// SCVx quaternions -> tailored to CVXPYgen/CPG solver
// Column-major updates for all cpg_update_* calls

#include <time.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cpg_workspace.h"   // ADAPT: must provide cpg_update_* and param API
#include "cpg_solve.h"       // ADAPT: must provide cpg_solve() and CPG_Result struct

/**************************** USER TUNABLE MACROS ****************************/
#define T 31-1               // Python T = 31-1 -> 30 (discrete steps)
#define TP1 (T+1)          // number of state columns
#define STATE 7            // 4 quaternion + 3 angular velocities
#define CTL 3              // 3 control torques
#define SIZE_N_SUB 20      // sub-steps for RK5 integration (python size_N)
#define KMAX 20            // maximum SCVx iterations
#define TAYLOR_ORDER 7     // Taylor series order for exp approximation
#define EPS_TINY 1e-12
/******************************************************************************/

/* physical & guidance parameters (copied from python defaults) */
double tf = 200.0;
double tau;                 // tf / T

double J1_double = 0.0333;
double J2_double = 0.0333;
double J3_double = 0.0067;

double u_max_torq_double = 0.001;
double omega_max_double = 1.0 * M_PI / 180.0;

double rho0 = 0.0;
double rho1 = 0.1;
double rho2 = 0.7;
double etta0 = 0.001;
double etta1 = 5.0;
double beta_sh = 2.0;
double beta_gr = 2.0;

double lamb_double = 2000.0;
double etta_double = 0.1;

double e_tol = 0.05;
double epsilon_stop_norm = 0.04;

/* scaling matrix (diagonal) for control */
double u_qw_scaling[3][3];

/* scratch big arrays (A_discrete has 7 rows and 7*T columns, etc.) */
double aux_A_discrete_qw[STATE][STATE * T];
double aux_B_discrete_qw_scaled[STATE][CTL * T];
double aux_w_discrete_qw[STATE][T];

/* state and control trajectories used as current linearization point (ox, ou) */
double ox[STATE][TP1];
double ou[CTL][T]; // control for k=0..T-1

/* solver initial guess copy buffers (flat column-major as solver expects) */
double ox_flat[STATE * TP1];
double ou_flat[CTL * T];

/* helper small functions */
static double sign_double(double x) { return x >= 0.0 ? 1.0 : -1.0; }

void scaling_begin(double u[CTL][T], double S[CTL][CTL])
{
    for (int k = 0; k < T; k++)
        for (int i = 0; i < CTL; i++) {
            double s = 0.0;
            for (int j = 0; j < CTL; j++)
                s += S[i][j] * u[j][k];
            u[i][k] = s;
        }
}

void scaling_end(double u[CTL][T], const double S[CTL][CTL])
{
    double tmp[CTL];

    for (int k = 0; k < T; k++) {
        for (int i = 0; i < CTL; i++) {
            tmp[i] = u[i][k] / S[i][i];   // inv(S) * u
        }
        for (int i = 0; i < CTL; i++)
            u[i][k] = tmp[i];
    }
}


/* quaternion utilities (q = [x,y,z,w]; w is scalar) */
void quat_normalize(double q[4]) {
    double s = 0.0;
    for (int i = 0; i < 4; i++) s += q[i]*q[i];
    s = sqrt(s);
    if (s < EPS_TINY) return;
    for (int i = 0; i < 4; i++) q[i] /= s;
}

/* q = p * r  (Hamilton product, p then r) */
void quat_mul(double p[4], double r[4], double q[4]) {
    // p = [px,py,pz,pw]; r = [rx,ry,rz,rw]
    q[0] = p[3]*r[0] + p[0]*r[3] + p[1]*r[2] - p[2]*r[1];
    q[1] = p[3]*r[1] - p[0]*r[2] + p[1]*r[3] + p[2]*r[0];
    q[2] = p[3]*r[2] + p[0]*r[1] - p[1]*r[0] + p[2]*r[3];
    q[3] = p[3]*r[3] - p[0]*r[0] - p[1]*r[1] - p[2]*r[2];
}

/* quaternion conjugate */
void quat_conj(double q[4], double qc[4]) {
    qc[0] = -q[0]; qc[1] = -q[1]; qc[2] = -q[2]; qc[3] = q[3];
}

/* convert quaternion (unit) to rotation vector (axis * angle) */
void quat_to_rotvec(double q[4], double rotv[3]) {
    double x = q[0], y = q[1], z = q[2], w = q[3];
    double vnorm = sqrt(x*x + y*y + z*z);
    if (vnorm < 1e-12) { rotv[0]=rotv[1]=rotv[2]=0.0; return; }
    double angle = 2.0 * atan2(vnorm, w);
    double invv = 1.0 / vnorm;
    rotv[0] = x * invv * angle;
    rotv[1] = y * invv * angle;
    rotv[2] = z * invv * angle;
}

/* SLERP - produce TP1 samples between q1 and q2; out is q[4][TP1] */
void slerp_quat(double q1_in[4], double q2_in[4], double outq[4][TP1]) {
    double q1[4], q2[4];
    for(int i=0;i<4;i++) { q1[i] = q1_in[i]; q2[i] = q2_in[i]; }

    // Producto punto
    double dot = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];
    if(dot < 0.0) {
        for(int i=0;i<4;i++) q2[i] = -q2[i];
        dot = -dot;
    }
    if(dot > 1.0) dot = 1.0;
    if(dot < -1.0) dot = -1.0;

    double theta_0 = acos(dot);

    if(fabs(theta_0) < 1e-12) {
        for(int i=0;i<TP1;i++) {
            double t = (double)(i) / (TP1 - 1);
            for(int j=0;j<4;j++)
                outq[j][i] = (1.0 - t) * q1[j] + t * q2[j];
        }
        return;
    }

    double sin_theta_0 = sin(theta_0);

    for(int i=0;i<TP1;i++) {
        double t = (double)(i) / (TP1 - 1);
        printf("t: %.4f \n",t);
        double theta = theta_0 * t;
        double sin_theta = sin(theta);

        double s0 = cos(theta) - dot * sin_theta / sin_theta_0;
        double s1 = sin_theta / sin_theta_0;

        for(int j=0;j<4;j++){
            outq[j][i] = s0 * q1[j] + s1 * q2[j];
            printf("outq, %.6f",outq[j][i]);}
    }
}

/* compute angular velocity approx from quaternion sequence q[4][TP1]
   ω_k = rotvec(q_{k+1} * conj(q_k)) / dt for k=0..T-1; last column duplicates previous
   outw[3][TP1] filled.
*/
void compute_angvel_from_quat(double qseq[4][TP1], double dt, double outw[3][TP1]) {
    double qc[4], prod[4], rotv[3];
    for (int k=0;k<T;k++) {
        // prod = q_{k+1} * conj(q_k)
        quat_conj(qseq[0] + 4*k, qc); // NOTE: this usage wrong for flatten — we access differently below
        // Access explicitly:
        double qk[4], qk1[4];
        for (int i=0;i<4;i++) { qk[i] = qseq[i][k]; qk1[i] = qseq[i][k+1]; }
        quat_conj(qk, qc);
        quat_mul(qk1, qc, prod);
        quat_normalize(prod);
        quat_to_rotvec(prod, rotv);
        outw[0][k] = rotv[0] / dt;
        outw[1][k] = rotv[1] / dt;
        outw[2][k] = rotv[2] / dt;
    }
    // duplicate last column
    
    outw[0][T] = outw[0][T-1];
    outw[1][T] = outw[1][T-1];
    outw[2][T] = outw[2][T-1];

    outw[0][0] = 0;
    outw[1][0] = 0;
    outw[2][0] = 0;
}

/* f_qw(x,u)   x: 7x1, u:3x1 -> out: 7x1 */
void f_qw(double x[STATE], double u[CTL], double out[STATE]) {
    double oq1=x[0], oq2=x[1], oq3=x[2], oq4=x[3];
    double ow1=x[4], ow2=x[5], ow3=x[6];
    double ou1=u[0], ou2=u[1], ou3=u[2];

    out[0] = 0.5*(oq4*ow1 - oq3*ow2 + oq2*ow3);
    out[1] = 0.5*(oq3*ow1 + oq4*ow2 - oq1*ow3);
    out[2] = 0.5*(-oq2*ow1 + oq1*ow2 + oq4*ow3);
    out[3] = 0.5*(-oq1*ow1 - oq2*ow2 - oq3*ow3);

    out[4] = (1.0/J1_double) * ( - (J3_double - J2_double) * ow2 * ow3 + ou1 );
    out[5] = (1.0/J2_double) * ( - (J1_double - J3_double) * ow3 * ow1 + ou2 );
    out[6] = (1.0/J3_double) * ( - (J2_double - J1_double) * ow1 * ow2 + ou3 );
}

/* A_qw(x) -> fill A matrix 7x7 */
void A_qw(double ox[STATE], double A[STATE][STATE]) {
    for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) A[i][j]=0.0;
    double oq1=ox[0], oq2=ox[1], oq3=ox[2], oq4=ox[3];
    double ow1=ox[4], ow2=ox[5], ow3=ox[6];

    double top[4][7] = {
        {0,   ow3, -ow2, ow1, oq4, -oq3, oq2},
        {-ow3, 0,  ow1, ow2, oq3,  oq4, -oq1},
        {ow2, -ow1, 0,  ow3, -oq2, oq1, oq4},
        {-ow1, -ow2, -ow3, 0, -oq1, -oq2, -oq3}
    };
    for (int i=0;i<4;i++) for (int j=0;j<7;j++) A[i][j] = 0.5 * top[i][j];

    A[4][5] = (1.0/J1_double) * ( - (J3_double - J2_double) * ow3 );
    A[4][6] = (1.0/J1_double) * ( - (J3_double - J2_double) * ow2 );
    A[5][4] = (1.0/J2_double) * ( - (J1_double - J3_double) * ow3 );
    A[5][6] = (1.0/J2_double) * ( - (J1_double - J3_double) * ow1 );
    A[6][4] = (1.0/J3_double) * ( - (J2_double - J1_double) * ow2 );
    A[6][5] = (1.0/J3_double) * ( - (J2_double - J1_double) * ow1 );
}

/* B_qw(x) -> fill B 7x3 (state independent except zeros on top rows) */
void B_qw(double B[STATE][CTL]) {
    for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) B[i][j]=0.0;
    B[4][0] = 1.0 / J1_double;
    B[5][1] = 1.0 / J2_double;
    B[6][2] = 1.0 / J3_double;
}

/* exp_matrix_taylor_A: approximate exp(hA) by I + hA + ... up to order */
void exp_matrix_taylor_A(double A[STATE][STATE], double h, int order, double out[STATE][STATE]) {
    // out = I + hA + sum_{k=2..order} (h^k / k!) * A^k
    // We'll compute powers iteratively.
    // Initialize out = I + hA
    double hA[STATE][STATE];
    for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) {
        hA[i][j] = h * A[i][j];
        out[i][j] = (i==j) ? 1.0 : 0.0;
        out[i][j] += hA[i][j];
    }
    // P = hA (power 1)
    double P[STATE][STATE];
    for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) P[i][j] = hA[i][j];
    double fact = 1.0;
    for (int k=2;k<=order;k++) {
        // Pnext = P * hA / hA? We want A^k * h^k: multiply P by hA/h? simpler multiply by hA then divide by k
        double Pnext[STATE][STATE];
        for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) {
            double s=0.0;
            for (int m=0;m<STATE;m++) s += P[i][m] * hA[m][j];
            Pnext[i][j] = s;
        }
        fact *= (double)(k);
        for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) out[i][j] += Pnext[i][j] / fact;
        for (int i=0;i<STATE;i++) for (int j=0;j<STATE;j++) P[i][j] = Pnext[i][j];
    }
}

/* exp_matrix_taylor_B: approximate integral_0^h exp(A s) ds * B
   Python used sum_{i=1..n} h^i / i! * A^{i-1} * B
*/
void exp_matrix_taylor_B(double A[STATE][STATE], double B[STATE][CTL], double h, int order, double out[STATE][CTL]) {
    // initialize: out = h * B  (i=1)
    for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) out[i][j] = h * B[i][j];
    // ApowB = B (A^{0}*B)
    double ApowB[STATE][CTL];
    for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) ApowB[i][j] = B[i][j];
    double hpow = h;
    double fact = 1.0;
    for (int k=2;k<=order+1;k++) {
        // ApowB = A * ApowB  (now A^{k-1} * B)
        double ApowB_next[STATE][CTL];
        for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) {
            double s=0.0;
            for (int m=0;m<STATE;m++) s += A[i][m] * ApowB[m][j];
            ApowB_next[i][j] = s;
        }
        hpow *= h;
        fact *= (double)(k);
        for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) out[i][j] += (hpow / fact) * ApowB_next[i][j];
        for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) ApowB[i][j] = ApowB_next[i][j];
    }
}

void RK5_step(double xk[STATE],
              double uk[CTL],
              double dt,
              double x_next[STATE])
{
    double k1[STATE], k2[STATE], k3[STATE], k4[STATE], k5[STATE];
    double tmp[STATE];

    f_SCVx(xk, uk, k1);

    for (int i=0;i<STATE;i++)
        tmp[i] = xk[i] + (1.0/4.0)*dt*k1[i];
    f_SCVx(tmp, uk, k2);

    for (int i=0;i<STATE;i++)
        tmp[i] = xk[i] + (3.0/8.0)*dt*k2[i];
    f_SCVx(tmp, uk, k3);

    for (int i=0;i<STATE;i++)
        tmp[i] = xk[i] + (12.0/13.0)*dt*k3[i];
    f_SCVx(tmp, uk, k4);

    for (int i=0;i<STATE;i++)
        tmp[i] = xk[i] + dt*k4[i];
    f_SCVx(tmp, uk, k5);

    for (int i=0;i<STATE;i++)
        x_next[i] = xk[i] + dt*(16.0/135.0*k1[i]
                                + 6656.0/12825.0*k3[i]
                                + 28561.0/56430.0*k4[i]
                                - 9.0/50.0*k5[i]);
}

/* f_SCVx: discrete map using SIZE_N_SUB RK5 substeps */
void f_SCVx_map(double xk[STATE], double uk[CTL], double outx[STATE]) {
    double x_aux[STATE];
    for (int i=0;i<STATE;i++) x_aux[i] = xk[i];
    double dt_sub = tau / (double)(SIZE_N_SUB);
    for (int s=0;s<SIZE_N_SUB;s++) {
        double xnext[STATE];
        rk5_step(x_aux, uk, dt_sub, xnext);
        for (int i=0;i<STATE;i++) x_aux[i] = xnext[i];
    }
    for (int i=0;i<STATE;i++) outx[i] = x_aux[i];
}

/* norms used for J/L cost (vector L2 and L1) */
double l2_norm_vec(double *v, int n) {
    double s=0.0;
    for (int i=0;i<n;i++) s += v[i]*v[i];
    return sqrt(s);
}
double l1_norm_vec(double *v, int n) {
    double s=0.0;
    for (int i=0;i<n;i++) s += fabs(v[i]);
    return s;
}

double J_SCVx_c(double x_seq[STATE][TP1], double u_seq[CTL][T]){
    double cost = 0.0;
    /* ---- Cost 1: quadratic control effort ---- */
    for (int k = 0; k < T; k++) {
        double norm2 = 0.0;
        for (int i = 0; i < 3; i++) norm2 += u_seq[i][k] * u_seq[i][k];
        cost += tau * norm2;
    }

    printf("Cost 1: %f\n", cost);

    scaling_end(u_seq, u_qw_scaling);

    for (int k = 0; k < T; k++) {
        double x_next_pred[STATE];
        double xk[STATE];
        double uk[CTL];
        for (int j =0; j<STATE; j++) xk[j] = x_seq[j][k];
        for (int j =0; j<CTL; j++) uk[j] = u_seq[j][k];
        f_SCVx_map(xk, uk, x_next_pred);
        for (int i = 0; i < STATE; i++) {
            double diff = x_seq[i][k+1] - x_next_pred[i];
            cost += tau * fabs(lamb_double * diff);
        }
    }

    /* ---- scaling begin (restore u) ---- */
    scaling_begin(u_seq, u_qw_scaling);

    printf("Cost 2: %f\n", cost);

    return cost;
}


/* L_SCVx: cost L used in SCVx (includes vc vector penalty) */
double L_SCVx_c(double x[STATE][TP1], double u[CTL][T], double vc[STATE][T]){
    double cost = 0.0;

    for (int k = 0; k < T; k++){
        /* tau * ||u(:,k)||_2^2 */
        double u_norm2_sq = 0.0;
        for (int i = 0; i < 3; i++){
            u_norm2_sq += (u[i][k]) * (u[i][k]);
        }
        cost += tau * u_norm2_sq;

        /* tau * || lamb * vc(:,k) ||_1 */
        double vc_norm1 = 0.0;
        for (int i = 0; i < STATE; i++){
            vc_norm1 += fabs(lamb_double * vc[i][k]);
        }
        cost += tau * vc_norm1;
    }

    return cost;
}

/* Flatten helpers (column-major): store x[STATE][TP1] into flat array of length STATE*TP1
   Order: for col=0..TP1-1, for row=0..STATE-1 -> flat[index++]=x[row][col]
*/
void flatten_state(double x[STATE][TP1], double out_flat[STATE * TP1]) {
    int idx=0;
    for (int b=0;b<TP1;b++) {
        for (int a=0;a<STATE;a++) {
            out_flat[idx++] = x[a][b];
        }
    }
}
void flatten_control(double u[CTL][T], double out_flat[CTL * T]) {
    int idx=0;
    for (int b=0;b<T;b++) {
        for (int a=0;a<CTL;a++) {
            out_flat[idx++] = u[a][b];
        }
    }
}

/* Unflatten from solver result arrays (assumed column-major flattened) into x[STATE][TP1] */
void unflatten_state_from_solver(double *flat_in, double x_out[STATE][TP1]) {
    int idx=0;
    for (int b=0;b<TP1;b++) {
        for (int a=0;a<STATE;a++) {
            x_out[a][b] = flat_in[idx++];
        }
    }
}
void unflatten_control_from_solver(double *flat_in, double u_out[CTL][T]) {
    int idx=0;
    for (int b=0;b<T;b++) {
        for (int a=0;a<CTL;a++) {
            u_out[a][b] = flat_in[idx++];
        }
    }
}
void unflatten_virtual_control_from_solver(double *flat_in, double vc_out[STATE][T]) {
    int idx=0;
    for (int b=0;b<T;b++) {
        for (int a=0;a<STATE;a++) {
            vc_out[a][b] = flat_in[idx++];
        }
    }
}
double get_time_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
/* ---------------------- MAIN SCVx loop tied to CPG ---------------------- */
/* ADAPT: the cpg_update_* function names must match your cpg_workspace.h */
void COGU_guidance(double u_global[CTL][T], double x_global[STATE][TP1], double start_quat[4], double end_quat[4]) {
    // initialize tau
    tau = tf / (double)(T);
    // scaling diagonal matrix
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) u_qw_scaling[i][j] = 0.0;
    u_qw_scaling[0][0] = 1.0 / u_max_torq_double;
    u_qw_scaling[1][1] = 1.0 / u_max_torq_double;
    u_qw_scaling[2][2] = 1.0 / u_max_torq_double;

    // Replace with actual Euler->quat or direct quaternions like Python's Rotation.from_euler outputs
    // For now we leave identity; user must set these to desired values before calling generator.

    // build startpos & endpos (7x1): q(4) then ω zeros
    double startpos[STATE], endpos[STATE];
    for (int i=0;i<4;i++) { startpos[i] = start_quat[i]; endpos[i] = end_quat[i]; }
    for (int i=4;i<STATE;i++) { startpos[i] = 0.0; endpos[i] = 0.0; }

    // Initialize ox via SLERP for quaternions and angular velocity approx
    double qseq[4][TP1];

    printf("q_startpos: \n");
    for (int k=0;k<7;k++){
        printf("%.6f, ",startpos[k]);
    }
    printf("\n");

    printf("q_endpos: \n");
    for (int k=0;k<7;k++){
        printf("%.6f, ",endpos[k]);
    }
    printf("\n");

    slerp_quat(startpos, endpos, qseq);

    printf("qseq: \n");
    for (int j=0;j<TP1;j++){
        for (int k=0;k<4;k++){
            printf("%.6f, ",qseq[k][j]);
        }
        printf("\n");
    }

    double wseq[3][TP1];
    compute_angvel_from_quat(qseq, tau, wseq);

    for (int k=0;k<TP1;k++) {
        for (int i=0;i<4;i++) ox[i][k] = qseq[i][k];
        for (int j=0;j<3;j++) ox[4+j][k] = wseq[j][k];
    }
    // initialize controls to zero (scaled)
    for (int k=0;k<T;k++) for (int i=0;i<CTL;i++) ou[i][k] = 0.0;

    // set solver parameters - ADAPT calls to your param update API
    cpg_update_omega_max(omega_max_double); // ADAPT if your API differs
    cpg_update_lamb(lamb_double);
    cpg_update_etta(etta_double);

    // set start / end pos into solver (example API used in your sample)
    for (int i=0;i<STATE;i++) {
        cpg_update_start_pos(i, startpos[i]);  // ADAPT name if necessary
        cpg_update_end_pos(i,   endpos[i]);
    }

    scaling_begin(ou,u_qw_scaling); //scaling begin

    int aux_count = 0;
    for (int k=0;k<TP1;k++) {
        for (int j=0;j<STATE;j++){
        cpg_update_ox_cvxpy(aux_count, (double)(ox[j][k]));
        aux_count++;
        }
    }

    aux_count = 0;
    for (int k=0;k<T;k++) {
        for (int j=0;j<CTL;j++){
        cpg_update_ou_cvxpy(aux_count, (double)(ou[j][k]));
        aux_count++;
        }
    }
    /* Main SCVx iterative loop */
    int iter = 1;
    int no_first_iterations = 0;
    double t_start = 0.0; // optional timing

    while (iter <= KMAX) {
        // Unscale controls (scaling_end) before computing linearization matrices

        scaling_end(ou,u_qw_scaling); //scaling end

        // invert diagonal scaling: inverse is diag(u_max)
        double inv_scal[3][3];
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) inv_scal[i][j] = 0.0;
        inv_scal[0][0] = 1.0 / u_qw_scaling[0][0];
        inv_scal[1][1] = 1.0 / u_qw_scaling[1][1];
        inv_scal[2][2] = 1.0 / u_qw_scaling[2][2];

        // Build discrete linearization matrices column by column (k=0..T-1)
        for (int k = 0; k < T; k++) {
            // extract ox column xk
            double xk[STATE];

            for (int i=0;i<STATE;i++){xk[i] = ox[i][k];};

            double Acont[STATE][STATE];
            double Bcont[STATE][CTL];
            A_qw(xk, Acont);
            B_qw(Bcont);
            /*printf("\n A_qw \n");
            for (int k=0;k<STATE;k++){
                for (int j=0;j<STATE;j++){
                    printf("%.6f, ",Acont[k][j]);
                }
                printf("\n");
            }*/
            
            double Adis[STATE][STATE];
            double Bdis[STATE][CTL];
            exp_matrix_taylor_A(Acont, tau, TAYLOR_ORDER, Adis);

            // compute B discrete scaled by inv(scaling) as python: B_qw * inv(u_qw_scaling)
            double Bscaled[STATE][CTL];
            for (int i=0;i<STATE;i++) for (int j=0;j<CTL;j++) {
                double s = 0.0;
                for (int m=0;m<CTL;m++) s += Bcont[i][m] * inv_scal[m][j];
                Bscaled[i][j] = s;
            }
            exp_matrix_taylor_B(Acont, Bscaled, tau, TAYLOR_ORDER, Bdis);
            
            // compute w_discrete = integral_0^h exp(A s) ds * w_local
            double fq[STATE], Ax[STATE], w_local[STATE];
            double zero_u[CTL] = {0.0,0.0,0.0};
            f_qw(xk, zero_u, fq);
            for (int i=0;i<STATE;i++) {
                double sum = 0.0;
                for (int j=0;j<STATE;j++) sum += Acont[i][j] * xk[j];
                Ax[i] = sum;
                w_local[i] = fq[i] - Ax[i];
            }
            // compute wdis = exp_matrix_taylor_B(Acont, w_local_as_Bcol, tau, order)  (7x1)
            // We'll call a small local adaptation: treat w_local as single-column B and compute same series
            double wdis[STATE];
            // i=1 term
            for (int i=0;i<STATE;i++) wdis[i] = tau * w_local[i];
            // iterative
            double Apow[STATE][1];
            for (int i=0;i<STATE;i++) Apow[i][0] = w_local[i];
            double hpow = tau;
            double fact = 1.0;
            for (int ord=2; ord<=TAYLOR_ORDER; ord++) {
                double Apow_next[STATE][1];
                for (int i=0;i<STATE;i++) {
                    double s=0.0;
                    for (int m=0;m<STATE;m++) s += Acont[i][m] * Apow[m][0];
                    Apow_next[i][0] = s;
                }
                hpow *= tau;
                fact *= (double)(ord);
                for (int i=0;i<STATE;i++) wdis[i] += (hpow / fact) * Apow_next[i][0];
                for (int i=0;i<STATE;i++) Apow[i][0] = Apow_next[i][0];
            }


            // copy Adis (7x7) to aux_A_discrete_qw column block 7*k..7*k+6 using column-major sweep:
            // solver expects aux_A_discrete_qw as 7 x (7*T) matrix; we fill aux_A_discrete_qw[row][7*k + col] = Adis[row][col]
            for (int col = 0; col < STATE; col++) {
                for (int row = 0; row < STATE; row++) {
                    aux_A_discrete_qw[row][STATE * k + col] = Adis[row][col];
                }
            }
            // copy Bdis (7x3) to aux_B_discrete_qw_scaled column block 3*k..3*k+2
            for (int col = 0; col < CTL; col++) {
                for (int row = 0; row < STATE; row++) {
                    aux_B_discrete_qw_scaled[row][CTL * k + col] = Bdis[row][col];
                }
            }
            // copy wdis into aux_w_discrete_qw column k
            for (int row=0; row<STATE; row++) aux_w_discrete_qw[row][k] = wdis[row];
        }

        // Now push these arrays to the solver via cpg_update_... column-major, using the same sweep order
        // ADAPT: we assume cpg_update_A_discrete_qw(aux_idx, value) updates sequentially with aux_idx from 0..(7*7*T-1)
        // Following your rule: sweep columns then rows -> index = row + STATE * colBlock where colBlock iterates 0..(7*T -1)
        // We'll use a single aux_count that increases in same order used in cpg_example.

        // Update A_discrete: there are (STATE)*(STATE*T) elements, but typical API expects count over flattened 7x(7*T)
        {
            int aux_count = 0;
            for (int colBlock = 0; colBlock < STATE * T; colBlock++) {
                for (int row = 0; row < STATE; row++) {
                    double val = aux_A_discrete_qw[row][colBlock];
                    cpg_update_A_discrete_qw(aux_count, val); // ADAPT: your API name
                    aux_count++;
                }
            }
        }
        // Update B_discrete scaled: shape 7 x (3*T)
        {
            int aux_count = 0;
            for (int colBlock = 0; colBlock < CTL * T; colBlock++) {
                for (int row = 0; row < STATE; row++) {
                    double val = aux_B_discrete_qw_scaled[row][colBlock];
                    cpg_update_B_discrete_qw_scaled(aux_count, val); // ADAPT name
                    aux_count++;
                }
            }
        }
        // Update w_discrete: shape 7 x T
        {
            int aux_count = 0;
            for (int colBlock = 0; colBlock < T; colBlock++) {
                for (int row = 0; row < STATE; row++) {
                    double val = aux_w_discrete_qw[row][colBlock];
                    cpg_update_w_discrete_qw(aux_count, val); // ADAPT name
                    aux_count++;
                }
            }
        }

        // Also update initial guess parameters for ox_cvxpy and ou_cvxpy as flat arrays

        // Solve problem instance
        printf("Hello world \n");
        cpg_solve();
        printf("obj = %.12f\n", CPG_Result.info->obj_val);
        printf("Hello world 2 \n");
        
        // Read solver's primal result arrays and unflatten into candidate x_opt and u_opt
        // ADAPT: names used below must match your generated C code from CVXPYgen
        // We assume CPG_Result.prim->new_x_qw is length STATE * TP1
        
        double x_opt[STATE][TP1];
        double u_opt[CTL][T];
        double vc_opt[STATE][T];

        int aux_counter = 0;
        for (int j=0;j<TP1;j++){
            for (int k=0;k<STATE;k++){
            x_opt[k][j]=CPG_Result.prim->nx[aux_counter];
            aux_counter++;
            }
        }

        aux_counter = 0;
        for (int j=0;j<T;j++){
            for (int k=0;k<CTL;k++){
            u_opt[k][j]=CPG_Result.prim->u[aux_counter];
            aux_counter++;
            }
        }

        aux_counter = 0;
        for (int j=0;j<T;j++){
            for (int k=0;k<STATE;k++){
            vc_opt[k][j]=CPG_Result.prim->vc[aux_counter];
            aux_counter++;
            }
        }

        scaling_begin(ou,u_qw_scaling); //scaling begin
        printf("%.12f u00 a\n",ou[1][0]);
        double oJ = J_SCVx_c(ox, ou);
        printf("%.12f u11 a\n",u_opt[1][0]);
        double Jopt = J_SCVx_c(x_opt, u_opt);
        printf("%.12f u11 b\n",u_opt[1][0]);
        double Lopt = L_SCVx_c(x_opt, u_opt, vc_opt);
        printf("%.12f u11 c\n",u_opt[1][0]);
        

        double Delta_J = oJ - Jopt;
        double Delta_L = oJ - Lopt;

        // print iteration summary
        // Norm_x_diff = max norm across columns
        double max_norm_x_diff = 0.0;
        for (int col=0; col<TP1; col++) {
            double diff_col[STATE];
            for (int r=0;r<STATE;r++) diff_col[r] = x_opt[r][col] - ox[r][col];
            double ncol = l2_norm_vec(diff_col, STATE);
            if (ncol > max_norm_x_diff) max_norm_x_diff = ncol;
        }
        printf("Iter %d: solver cost (L) reported, L:%.6f, oJ:%.6f Jopt:%.6f Lopt:%.6f Delta_J:%.6f Delta_L:%.6f Norm_x_diff:%.6f\n",
               iter, Lopt, oJ, Jopt, Lopt, Delta_J, Delta_L, max_norm_x_diff);

        // SCVx acceptance and etta update rules (same logic as Python)
        if ((Delta_L < e_tol * fabs(oJ) || max_norm_x_diff < epsilon_stop_norm) && no_first_iterations) {
            // convergence achieved
            // store global results (unscale the internal ou before returning)
            scaling_end(ou,u_qw_scaling); //Scaling end

            for (int col=0; col<TP1; col++) for (int row=0; row<STATE; row++) x_global[row][col] = ox[row][col];
            for (int col=0; col<T; col++) for (int row=0; row<CTL; row++) u_global[row][col] = ou[row][col];

            printf("SCVx converged at iter %d\n", iter);
            break;
        } else {
            double rho_i = (Delta_L == 0.0) ? 0.0 : Delta_J / Delta_L;
            if (rho_i < rho0) {
                // reject step
                double new_etta = fmax(etta0, etta_double / beta_sh);
                etta_double = new_etta;
                // keep ox/ou unchanged
            } else if (rho_i >= rho0 && rho_i < rho1) {
                etta_double = fmax(etta0, etta_double / beta_sh);
                // accept x_opt/u_opt
                for (int col=0; col<TP1; col++) for (int row=0; row<STATE; row++) ox[row][col] = x_opt[row][col];
                for (int col=0; col<T; col++) for (int row=0; row<CTL; row++) ou[row][col] = u_opt[row][col]; // store unscaled into ou
            } else if (rho_i >= rho1 && rho_i < rho2) {
                // accept
                for (int col=0; col<TP1; col++) for (int row=0; row<STATE; row++) ox[row][col] = x_opt[row][col];
                for (int col=0; col<T; col++) for (int row=0; row<CTL; row++) ou[row][col] = u_opt[row][col];
            } else { // rho_i >= rho2
                etta_double = fmin(etta1, beta_gr * etta_double);
                for (int col=0; col<TP1; col++) for (int row=0; row<STATE; row++) ox[row][col] = x_opt[row][col];
                for (int col=0; col<T; col++) for (int row=0; row<CTL; row++) ou[row][col] = u_opt[row][col];
            }
            cpg_update_etta(etta_double);
            printf("Iter %d: updated etta=%.6f rho=%.6f\n", iter, etta_double, rho_i);

            int aux_count = 0;
            for (int k=0;k<TP1;k++) {
                for (int j=0;j<STATE;j++){
                cpg_update_ox_cvxpy(aux_count, ox[j][k]);
                aux_count++;
                }
            }

            aux_count = 0;
            for (int k=0;k<T;k++) {
                for (int j=0;j<CTL;j++){
                cpg_update_ou_cvxpy(aux_count, ou[j][k]);
                aux_count++;
                }
            }
        }

        // set no_first_iterations after first few iters as in Python
        if (iter == 3) no_first_iterations = 1;
        iter++;
    } // end SCVx loop

    printf("SCVx finished (CPG wrapper). Iterations: %d\n", iter-1);
}

int main()
{
    double start_quat[4] = {0.0, 1.0, 0.0, 0.0};
    double end_quat[4]   = {0.0, 0.0, 0.0, 1.0};
    double u_global[CTL][T];
    double x_global[STATE][TP1];
    double t0 = get_time_sec();
    COGU_guidance(u_global,x_global,start_quat,end_quat);
    double t1 = get_time_sec();
    printf("Total time: %.6f s\n", t1 - t0);
    
    printf("\n u_global \n");
    for (int k=0;k<T;k++) {
        for (int j=0;j<CTL;j++){
        printf("%e, ", u_global[j][k]);
        }
        printf("\n");
    }
    printf("\n x_global \n");
    for (int k=0;k<TP1;k++) {
        for (int j=0;j<STATE;j++){
        printf("%e, ", x_global[j][k]);
        }
        printf("\n");
    }
}