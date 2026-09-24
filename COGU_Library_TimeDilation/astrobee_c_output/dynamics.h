#ifndef DYNAMICS_H
#define DYNAMICS_H

/* nx=13, nu=6, ng=3, np=0 */

void f_dynamics(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void A_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void B_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void F_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void y_affine(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void g_constraints(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void C_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void D_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void G_jacobian(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);
void z_affine(const double *hx, const double *hu, const double *hp, double t, const double *dyn_par, double *out);

#endif /* DYNAMICS_H */
