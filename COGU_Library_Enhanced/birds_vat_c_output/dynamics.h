#ifndef DYNAMICS_H
#define DYNAMICS_H

/* nx=7, nu=4, ng=0 */

void f_dynamics(double *hx, double *hu, double t, double *dyn_par, double *out);
void A_jacobian(double *hx, double *hu, double t, double *dyn_par, double *out);
void B_jacobian(double *hx, double *hu, double t, double *dyn_par, double *out);
void y_affine(double *hx, double *hu, double t, double *dyn_par, double *out);
void g_constraints(double *hx, double *hu, double t, double *dyn_par, double *out);
void C_jac(double *hx, double *hu, double t, double *dyn_par, double *out);
void D_jac(double *hx, double *hu, double t, double *dyn_par, double *out);
void z_affine(double *hx, double *hu, double t, double *dyn_par, double *out);

#endif /* DYNAMICS_H */
