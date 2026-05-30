/*
 * scvx_main.c  —  SCVx loop generado por COGU pipeline
 *
 * Dimensiones: NX={NX}, NU={NU}, NG={NG}, NP={NP}, T={T}
 *
 * Compilar (desde el directorio de salida):
 *   gcc -O0 -std=c99 -D_USE_MATH_DEFINES \
 *       scvx_main.c dynamics.c cpg_compat.h \
 *       solver/c/src/cpg_solve.c solver/c/src/cpg_workspace.c \
 *       solver/c/solver_code/src/*.c \
 *       solver/c/solver_code/external/amd/src/*.c \
 *       solver/c/solver_code/external/ldl/src/ldl.c \
 *       -I solver/c/include \
 *       -I solver/c/solver_code/include \
 *       -I solver/c/solver_code/external/amd/include \
 *       -I solver/c/solver_code/external/ldl/include \
 *       -I solver/c/solver_code/external/SuiteSparse_config \
 *       -lm -o scvx_cogu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "cpg_workspace.h"
#include "cpg_compat.h"
#include "dynamics.h"

/* ─────────────────────────────────────────────────────────────────────────────
 * Dimensiones
 * ───────────────────────────────────────────────────────────────────────────*/
#define STATES_SIZE   {NX}
#define INPUTS_SIZE   {NU}
#define N_OBS         {NG}
#define DYN_PAR_SIZE  {NP}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ─────────────────────────────────────────────────────────────────────────────
 * Algebra lineal basica (row-major internamente)
 * ───────────────────────────────────────────────────────────────────────────*/
static void mat_zero(double *A, int n) { memset(A,0,n*sizeof(double)); }
static void mat_copy(const double *s, double *d, int n) { memcpy(d,s,n*sizeof(double)); }

static void mat_eye(double *A, int n) {
    mat_zero(A,n*n);
    for(int i=0;i<n;i++) A[i*n+i]=1.0;
}
static void mat_mul(const double *A,const double *B,double *C,int ra,int ca,int cb){
    for(int i=0;i<ra;i++) for(int j=0;j<cb;j++){
        double s=0; for(int k=0;k<ca;k++) s+=A[i*ca+k]*B[k*cb+j]; C[i*cb+j]=s;
    }
}
static void mat_sub(const double *A,const double *B,double *C,int n){
    for(int i=0;i<n;i++) C[i]=A[i]-B[i];
}
static double vec_norm1(const double *v,int n){
    double s=0; for(int i=0;i<n;i++) s+=fabs(v[i]); return s;
}
static int mat_inv(const double *A,double *Ai,int n){
    double *aug=(double*)malloc(n*2*n*sizeof(double));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) aug[i*2*n+j]=A[i*n+j];
        for(int j=0;j<n;j++) aug[i*2*n+n+j]=(i==j)?1.0:0.0;
    }
    for(int col=0;col<n;col++){
        int piv=col;
        for(int r=col+1;r<n;r++) if(fabs(aug[r*2*n+col])>fabs(aug[piv*2*n+col])) piv=r;
        if(fabs(aug[piv*2*n+col])<1e-14){free(aug);return -1;}
        for(int j=0;j<2*n;j++){double tmp=aug[col*2*n+j];aug[col*2*n+j]=aug[piv*2*n+j];aug[piv*2*n+j]=tmp;}
        double d=aug[col*2*n+col];
        for(int j=0;j<2*n;j++) aug[col*2*n+j]/=d;
        for(int r=0;r<n;r++){
            if(r==col) continue;
            double f=aug[r*2*n+col];
            for(int j=0;j<2*n;j++) aug[r*2*n+j]-=f*aug[col*2*n+j];
        }
    }
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) Ai[i*n+j]=aug[i*2*n+n+j];
    free(aug); return 0;
}

static void to_colmaj(const double *M_row, double *M_col, int ra, int ca){
    for(int i=0;i<ra;i++)
        for(int j=0;j<ca;j++)
            M_col[j*ra+i] = M_row[i*ca+j];
}

/* ─────────────────────────────────────────────────────────────────────────────
 * Dinamica — delegada a dynamics.h (generado por generate_c_functions)
 * ───────────────────────────────────────────────────────────────────────────*/
static void f_SCVx(const double *hx,const double *hu,double t,const double *dp,double *fo){
    f_dynamics((double*)hx,(double*)hu,t,(double*)dp,fo);
}
static void g_SCVx(const double *hx,const double *hu,double t,const double *dp,double *go){
    g_constraints((double*)hx,(double*)hu,t,(double*)dp,go);
}
static void A_nn(const double *hx,const double *hu,double t,const double *dp,double *Ao){
    A_jacobian((double*)hx,(double*)hu,t,(double*)dp,Ao);
}
static void B_nm(const double *hx,const double *hu,double t,const double *dp,double *Bo){
    B_jacobian((double*)hx,(double*)hu,t,(double*)dp,Bo);
}
static void y_ns(const double *hx,const double *hu,double t,const double *dp,double *yo){
    y_affine((double*)hx,(double*)hu,t,(double*)dp,yo);
}
static void C_on(const double *hx,const double *hu,double t,const double *dp,double *Co){
    C_jac((double*)hx,(double*)hu,t,(double*)dp,Co);
}
static void D_om(const double *hx,const double *hu,double t,const double *dp,double *Do){
    D_jac((double*)hx,(double*)hu,t,(double*)dp,Do);
}
static void z_ns_g(const double *hx,const double *hu,double t,const double *dp,double *zo){
    z_affine((double*)hx,(double*)hu,t,(double*)dp,zo);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * Escalado
 * ───────────────────────────────────────────────────────────────────────────*/
static void sc_x(const double *x,const double *iSx,const double *cx,double *xs){
    double tmp[STATES_SIZE];
    for(int i=0;i<STATES_SIZE;i++) tmp[i]=x[i]-cx[i];
    mat_mul(iSx,tmp,xs,STATES_SIZE,STATES_SIZE,1);
}
static void isc_x(const double *xs,const double *Sx,const double *cx,double *x){
    double tmp[STATES_SIZE];
    mat_mul(Sx,xs,tmp,STATES_SIZE,STATES_SIZE,1);
    for(int i=0;i<STATES_SIZE;i++) x[i]=tmp[i]+cx[i];
}
static void sc_u(const double *u,const double *iSu,const double *cu,double *us){
    double tmp[INPUTS_SIZE];
    for(int i=0;i<INPUTS_SIZE;i++) tmp[i]=u[i]-cu[i];
    mat_mul(iSu,tmp,us,INPUTS_SIZE,INPUTS_SIZE,1);
}
static void isc_u(const double *us,const double *Su,const double *cu,double *u){
    double tmp[INPUTS_SIZE];
    mat_mul(Su,us,tmp,INPUTS_SIZE,INPUTS_SIZE,1);
    for(int i=0;i<INPUTS_SIZE;i++) u[i]=tmp[i]+cu[i];
}

static void As(const double *A_ns,const double *iSx,const double *Sx,double *Ao){
    double t1[STATES_SIZE*STATES_SIZE];
    mat_mul(iSx,A_ns,t1,STATES_SIZE,STATES_SIZE,STATES_SIZE);
    mat_mul(t1,Sx,Ao,STATES_SIZE,STATES_SIZE,STATES_SIZE);
}
static void Bs(const double *B_ns,const double *iSx,const double *Su,double *Bo){
    double t1[STATES_SIZE*INPUTS_SIZE];
    mat_mul(iSx,B_ns,t1,STATES_SIZE,STATES_SIZE,INPUTS_SIZE);
    mat_mul(t1,Su,Bo,STATES_SIZE,INPUTS_SIZE,INPUTS_SIZE);
}
static void ys(const double *A_ns,const double *B_ns,const double *y_n,
               const double *iSx,const double *cx,const double *cu,double *yo){
    double Acx[STATES_SIZE],Bcu[STATES_SIZE],t1[STATES_SIZE],t2[STATES_SIZE],t3[STATES_SIZE];
    mat_mul(A_ns,cx,Acx,STATES_SIZE,STATES_SIZE,1);
    mat_mul(iSx,Acx,t1,STATES_SIZE,STATES_SIZE,1);
    mat_mul(B_ns,cu,Bcu,STATES_SIZE,INPUTS_SIZE,1);
    mat_mul(iSx,Bcu,t2,STATES_SIZE,STATES_SIZE,1);
    mat_mul(iSx,y_n,t3,STATES_SIZE,STATES_SIZE,1);
    for(int i=0;i<STATES_SIZE;i++) yo[i]=t1[i]+t2[i]+t3[i];
}
static void Cs(const double *C_ns,const double *Sx,double *Co){
    mat_mul(C_ns,Sx,Co,N_OBS,STATES_SIZE,STATES_SIZE);
}
static void Ds(const double *D_ns,const double *Su,double *Do){
    mat_mul(D_ns,Su,Do,N_OBS,INPUTS_SIZE,INPUTS_SIZE);
}
static void zs(const double *C_ns,const double *D_ns,const double *z_n,
               const double *cx,const double *cu,double *zo){
    double Ccx[N_OBS],Dcu[N_OBS];
    mat_mul(C_ns,cx,Ccx,N_OBS,STATES_SIZE,1);
    mat_mul(D_ns,cu,Dcu,N_OBS,INPUTS_SIZE,1);
    for(int i=0;i<N_OBS;i++) zo[i]=Ccx[i]+Dcu[i]+z_n[i];
}

/* ─────────────────────────────────────────────────────────────────────────────
 * RK4 variacional → Phi, G, zv
 * ───────────────────────────────────────────────────────────────────────────*/
static void disc_ABY(const double *ohx,const double *ohu,double t0,const double *dp,double tau,
                     const double *Sx,const double *iSx,const double *cx,
                     const double *Su,const double *iSu,const double *cu,
                     int sN, double *Phi,double *G,double *zv)
{
    int ns=STATES_SIZE, mi=INPUTS_SIZE;
    mat_eye(Phi,ns); mat_zero(G,ns*mi); mat_zero(zv,ns);
    double dt=tau/(sN-1), t=t0;
    double An[STATES_SIZE*STATES_SIZE],Bn[STATES_SIZE*INPUTS_SIZE],yn[STATES_SIZE];
    double Asc[STATES_SIZE*STATES_SIZE],Bsc[STATES_SIZE*INPUTS_SIZE],ysc[STATES_SIZE];
    double k1P[STATES_SIZE*STATES_SIZE],k1G[STATES_SIZE*INPUTS_SIZE],k1z[STATES_SIZE];
    double k2P[STATES_SIZE*STATES_SIZE],k2G[STATES_SIZE*INPUTS_SIZE],k2z[STATES_SIZE];
    double k3P[STATES_SIZE*STATES_SIZE],k3G[STATES_SIZE*INPUTS_SIZE],k3z[STATES_SIZE];
    double k4P[STATES_SIZE*STATES_SIZE],k4G[STATES_SIZE*INPUTS_SIZE],k4z[STATES_SIZE];
    double P2[STATES_SIZE*STATES_SIZE],G2[STATES_SIZE*INPUTS_SIZE],z2[STATES_SIZE];
    double P3[STATES_SIZE*STATES_SIZE],G3[STATES_SIZE*INPUTS_SIZE],z3[STATES_SIZE];
    double P4[STATES_SIZE*STATES_SIZE],G4[STATES_SIZE*INPUTS_SIZE],z4[STATES_SIZE];
    double tm2[STATES_SIZE*INPUTS_SIZE],tm3[STATES_SIZE];

    for(int s=0;s<sN-1;s++){
        /* k1 */
        A_nn(ohx,ohu,t,dp,An); B_nm(ohx,ohu,t,dp,Bn); y_ns(ohx,ohu,t,dp,yn);
        As(An,iSx,Sx,Asc); Bs(Bn,iSx,Su,Bsc); ys(An,Bn,yn,iSx,cx,cu,ysc);
        mat_mul(Asc,Phi,k1P,ns,ns,ns);
        mat_mul(Asc,G,tm2,ns,ns,mi); for(int i=0;i<ns*mi;i++) k1G[i]=tm2[i]+Bsc[i];
        mat_mul(Asc,zv,tm3,ns,ns,1); for(int i=0;i<ns;i++) k1z[i]=tm3[i]+ysc[i];
        /* k2 */
        double t2=t+0.5*dt;
        for(int i=0;i<ns*ns;i++) P2[i]=Phi[i]+0.5*dt*k1P[i];
        for(int i=0;i<ns*mi;i++) G2[i]=G[i]+0.5*dt*k1G[i];
        for(int i=0;i<ns;i++) z2[i]=zv[i]+0.5*dt*k1z[i];
        A_nn(ohx,ohu,t2,dp,An); B_nm(ohx,ohu,t2,dp,Bn); y_ns(ohx,ohu,t2,dp,yn);
        As(An,iSx,Sx,Asc); Bs(Bn,iSx,Su,Bsc); ys(An,Bn,yn,iSx,cx,cu,ysc);
        mat_mul(Asc,P2,k2P,ns,ns,ns);
        mat_mul(Asc,G2,tm2,ns,ns,mi); for(int i=0;i<ns*mi;i++) k2G[i]=tm2[i]+Bsc[i];
        mat_mul(Asc,z2,tm3,ns,ns,1); for(int i=0;i<ns;i++) k2z[i]=tm3[i]+ysc[i];
        /* k3 (mismo t2) */
        for(int i=0;i<ns*ns;i++) P3[i]=Phi[i]+0.5*dt*k2P[i];
        for(int i=0;i<ns*mi;i++) G3[i]=G[i]+0.5*dt*k2G[i];
        for(int i=0;i<ns;i++) z3[i]=zv[i]+0.5*dt*k2z[i];
        mat_mul(Asc,P3,k3P,ns,ns,ns);
        mat_mul(Asc,G3,tm2,ns,ns,mi); for(int i=0;i<ns*mi;i++) k3G[i]=tm2[i]+Bsc[i];
        mat_mul(Asc,z3,tm3,ns,ns,1); for(int i=0;i<ns;i++) k3z[i]=tm3[i]+ysc[i];
        /* k4 */
        double t4=t+dt;
        for(int i=0;i<ns*ns;i++) P4[i]=Phi[i]+dt*k3P[i];
        for(int i=0;i<ns*mi;i++) G4[i]=G[i]+dt*k3G[i];
        for(int i=0;i<ns;i++) z4[i]=zv[i]+dt*k3z[i];
        A_nn(ohx,ohu,t4,dp,An); B_nm(ohx,ohu,t4,dp,Bn); y_ns(ohx,ohu,t4,dp,yn);
        As(An,iSx,Sx,Asc); Bs(Bn,iSx,Su,Bsc); ys(An,Bn,yn,iSx,cx,cu,ysc);
        mat_mul(Asc,P4,k4P,ns,ns,ns);
        mat_mul(Asc,G4,tm2,ns,ns,mi); for(int i=0;i<ns*mi;i++) k4G[i]=tm2[i]+Bsc[i];
        mat_mul(Asc,z4,tm3,ns,ns,1); for(int i=0;i<ns;i++) k4z[i]=tm3[i]+ysc[i];
        /* update */
        for(int i=0;i<ns*ns;i++) Phi[i]+=(dt/6.0)*(k1P[i]+2*k2P[i]+2*k3P[i]+k4P[i]);
        for(int i=0;i<ns*mi;i++) G[i]+=(dt/6.0)*(k1G[i]+2*k2G[i]+2*k3G[i]+k4G[i]);
        for(int i=0;i<ns;i++) zv[i]+=(dt/6.0)*(k1z[i]+2*k2z[i]+2*k3z[i]+k4z[i]);
        t+=dt;
    }
}

static void disc_CDZ(const double *ohx,const double *ohu,double t,const double *dp,
                     const double *Sx,const double *cx,const double *Su,const double *cu,
                     double *Co,double *Do,double *zo){
    double Cn[N_OBS*STATES_SIZE],Dn[N_OBS*INPUTS_SIZE],zn[N_OBS];
    C_on(ohx,ohu,t,dp,Cn); D_om(ohx,ohu,t,dp,Dn); z_ns_g(ohx,ohu,t,dp,zn);
    Cs(Cn,Sx,Co); Ds(Dn,Su,Do); zs(Cn,Dn,zn,cx,cu,zo);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * RK4 dinamica no lineal
 * ───────────────────────────────────────────────────────────────────────────*/
static void rk4_step(const double *xk,const double *uk,double t,double dt,const double *dp,
                     const double *Sx,const double *cx,const double *Su,const double *cu,
                     double *xn){
    double xns[STATES_SIZE],uns[INPUTS_SIZE];
    double k1[STATES_SIZE],k2[STATES_SIZE],k3[STATES_SIZE],k4[STATES_SIZE],tmp[STATES_SIZE];
    isc_x(xk,Sx,cx,xns); isc_u(uk,Su,cu,uns);
    f_SCVx(xns,uns,t,dp,k1);
    for(int i=0;i<STATES_SIZE;i++) tmp[i]=xns[i]+0.5*dt*k1[i]; f_SCVx(tmp,uns,t+0.5*dt,dp,k2);
    for(int i=0;i<STATES_SIZE;i++) tmp[i]=xns[i]+0.5*dt*k2[i]; f_SCVx(tmp,uns,t+0.5*dt,dp,k3);
    for(int i=0;i<STATES_SIZE;i++) tmp[i]=xns[i]+dt*k3[i];      f_SCVx(tmp,uns,t+dt,dp,k4);
    for(int i=0;i<STATES_SIZE;i++) xn[i]=xns[i]+(dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * SLERP y velocidad angular (utiles para problemas de actitud)
 * ───────────────────────────────────────────────────────────────────────────*/
static void slerp(const double *q1,const double *q2,int N,double *out){
    double qa[4],qb[4];
    memcpy(qa,q1,32); memcpy(qb,q2,32);
    double dot=0;
    for(int i=0;i<4;i++) dot+=qa[i]*qb[i];
    if(dot<0){for(int i=0;i<4;i++) qb[i]=-qb[i]; dot=-dot;}
    if(dot>1) dot=1;
    double th0=acos(dot);
    for(int s=0;s<N;s++){
        double t=(N>1)?(double)s/(N-1):0.0;
        if(fabs(th0)<1e-6){
            for(int i=0;i<4;i++) out[s*4+i]=(1-t)*qa[i]+t*qb[i];
        } else {
            double st0=sin(th0),th=th0*t;
            double s0=cos(th)-dot*sin(th)/st0, s1=sin(th)/st0;
            for(int i=0;i<4;i++) out[s*4+i]=s0*qa[i]+s1*qb[i];
        }
    }
}

static void ang_vel(const double *q,int N,double dt,double *w){
    for(int i=0;i<3;i++) w[i]=0;
    for(int k=0;k<N-1;k++){
        const double *qc=q+k*4,*qn=q+(k+1)*4;
        double c0=qc[0],c1=-qc[1],c2=-qc[2],c3=-qc[3];
        double d0=c0*qn[0]-c1*qn[1]-c2*qn[2]-c3*qn[3];
        double d1=c0*qn[1]+c1*qn[0]+c2*qn[3]-c3*qn[2];
        double d2=c0*qn[2]-c1*qn[3]+c2*qn[0]+c3*qn[1];
        double d3=c0*qn[3]+c1*qn[2]-c2*qn[1]+c3*qn[0];
        double vn=sqrt(d1*d1+d2*d2+d3*d3),ang=2*atan2(vn,d0);
        double *wk=w+(k+1)*3;
        if(vn>1e-12){wk[0]=(ang/vn)*d1/dt;wk[1]=(ang/vn)*d2/dt;wk[2]=(ang/vn)*d3/dt;}
        else{wk[0]=wk[1]=wk[2]=0;}
    }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * Funcion de costo J
 * ───────────────────────────────────────────────────────────────────────────*/
static double J_cost(const double *x,const double *u,int T,double tau,const double *dp,
                     double lam,const double *Sx,const double *cx,
                     const double *Su,const double *cu,double tf){
    int ns=STATES_SIZE,mi=INPUTS_SIZE;
    double cost=0;
    for(int k=0;k<T;k++){
        const double *uk=u+k*mi; double n1=0,n2=0;
        for(int i=0;i<3;i++) n1+=uk[i]*uk[i];
        for(int i=3;i<mi;i++) n2+=uk[i]*uk[i];
        cost+=tau*(n1+n2);
    }
    for(int k=0;k<T;k++){
        double xn[STATES_SIZE],xk1n[STATES_SIZE],def[STATES_SIZE];
        rk4_step(x+k*ns,u+k*mi,(double)k/T*tf,tau,dp,Sx,cx,Su,cu,xn);
        isc_x(x+(k+1)*ns,Sx,cx,xk1n);
        for(int i=0;i<ns;i++) def[i]=xk1n[i]-xn[i];
        cost+=tau*vec_norm1(def,ns)*lam;
    }
    for(int k=0;k<=T;k++){
        double xkn[STATES_SIZE],ukn[INPUTS_SIZE],gv[N_OBS];
        isc_x(x+k*ns,Sx,cx,xkn);
        if(k==T) isc_u(u+(T-1)*mi,Su,cu,ukn);
        else     isc_u(u+k*mi,Su,cu,ukn);
        g_SCVx(xkn,ukn,(double)k/T*tf,dp,gv);
        for(int i=0;i<N_OBS;i++){
            double v=gv[i]>0?gv[i]:0;
            cost+=tau*fabs(lam*v);
        }
    }
    return cost;
}

static double now_ms(void){
    return 1000.0 * (double)clock() / (double)CLOCKS_PER_SEC;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * main — ADAPTAR al problema concreto (condiciones iniciales, dp[], escalado)
 * ───────────────────────────────────────────────────────────────────────────*/
int main(void)
{
    int T={T}; double tf=200.0, tau_d=tf/T; int sN=20;
    int ns=STATES_SIZE, mi=INPUTS_SIZE;

    /* Condiciones de frontera */
    double sr[3]={0.0,-1.0,1.0},   er[3]={5.0,2.0,1.6};
    double sv[3]={0.001,-0.002,0.0},ev[3]={0,0,0};
    double sq[4]={0,1,0,0},         eq[4]={0,0,1,0};

    double start_pos[STATES_SIZE]={0}, end_pos[STATES_SIZE]={0};
    for(int i=0;i<3;i++){start_pos[i]=sr[i];start_pos[3+i]=sv[i];end_pos[i]=er[i];end_pos[3+i]=ev[i];}
    for(int i=0;i<4;i++){start_pos[6+i]=sq[i];end_pos[6+i]=eq[i];}

    double dp[DYN_PAR_SIZE]={
        0.1083,0.1083,0.1083, 0.01,0.02,0.01,
        1.1,-0.5,1.0, 2.6,0.5,1.1, 4.0,1.6,1.2, 0.8,0.8,0.8
    };

    double acc_max=20.0/7.2*0.001, torq_max=100.0*1e-6;
    double vel_max=5.0, omega_max=5.0*M_PI/180.0;

    /* Matrices de escalado */
    double Su[INPUTS_SIZE*INPUTS_SIZE]; mat_zero(Su,mi*mi);
    Su[0*mi+0]=Su[1*mi+1]=Su[2*mi+2]=2*(10*acc_max);
    Su[3*mi+3]=Su[4*mi+4]=Su[5*mi+5]=2*(0.1*torq_max);
    double cu[INPUTS_SIZE]={-(10*acc_max),-(10*acc_max),-(10*acc_max),
                             -(0.1*torq_max),-(0.1*torq_max),-(0.1*torq_max)};

    double Sx[STATES_SIZE*STATES_SIZE]; mat_zero(Sx,ns*ns);
    Sx[0*ns+0]=Sx[1*ns+1]=Sx[2*ns+2]=2*30;
    Sx[3*ns+3]=Sx[4*ns+4]=Sx[5*ns+5]=2*vel_max;
    Sx[6*ns+6]=Sx[7*ns+7]=Sx[8*ns+8]=Sx[9*ns+9]=2;
    Sx[10*ns+10]=Sx[11*ns+11]=Sx[12*ns+12]=2*(0.1*omega_max);
    double cx[STATES_SIZE]={-30,-30,-30,-vel_max,-vel_max,-vel_max,
                             -1,-1,-1,-1,-(0.1*omega_max),-(0.1*omega_max),-(0.1*omega_max)};

    double iSx[STATES_SIZE*STATES_SIZE], iSu[INPUTS_SIZE*INPUTS_SIZE];
    mat_inv(Sx,iSx,ns); mat_inv(Su,iSu,mi);

    int COGU_max=25;
    double rho0=0.0, rho1=0.1, rho2=0.7;
    double etta0=1e-8, etta1=10.0, beta_sh=2.0, beta_gr=2.0;
    double lam=1000.0, etta=1.0, e_tol=0.005, eps_stop=1e-5;

    /* Trayectoria inicial */
    double *ohx=(double*)calloc(ns*(T+1),sizeof(double));
    double *ohu=(double*)calloc(mi*T,    sizeof(double));
    for(int t=0;t<=T;t++){
        double a=(double)t/T;
        for(int i=0;i<3;i++) ohx[i+t*ns]=(1-a)*sr[i]+a*er[i];
    }
    for(int t=1;t<T;t++)
        for(int i=0;i<3;i++)
            ohx[(3+i)+t*ns]=(ohx[i+t*ns]-ohx[i+(t-1)*ns])/tau_d;
    for(int i=0;i<3;i++){ohx[(3+i)+0*ns]=sv[i]; ohx[(3+i)+T*ns]=ev[i];}
    double *qt=(double*)malloc((T+1)*4*sizeof(double));
    slerp(sq,eq,T+1,qt);
    for(int t=0;t<=T;t++) for(int i=0;i<4;i++) ohx[(6+i)+t*ns]=qt[t*4+i];
    double *av=(double*)calloc((T+1)*3,sizeof(double));
    ang_vel(qt,T+1,tau_d,av);
    for(int t=0;t<=T;t++) for(int i=0;i<3;i++) ohx[(10+i)+t*ns]=av[t*3+i];

    double *ox=(double*)malloc(ns*(T+1)*sizeof(double));
    double *ou=(double*)malloc(mi*T*    sizeof(double));
    for(int t=0;t<=T;t++) sc_x(ohx+t*ns,iSx,cx,ox+t*ns);
    for(int t=0;t<T;t++)  sc_u(ohu+t*mi,iSu,cu,ou+t*mi);

    cpg_update_vel_max(vel_max);
    cpg_update_omega_max(omega_max);
    cpg_update_acc_max(acc_max);
    cpg_update_torq_max(torq_max);

    {
        double Sx_cm[STATES_SIZE*STATES_SIZE];
        double Su_cm[INPUTS_SIZE*INPUTS_SIZE];
        to_colmaj(Sx,Sx_cm,ns,ns);
        to_colmaj(Su,Su_cm,mi,mi);
        for(int i=0;i<ns*ns;i++) cpg_update_S_x_scaling(i,Sx_cm[i]);
        for(int i=0;i<mi*mi;i++) cpg_update_S_u_scaling(i,Su_cm[i]);
    }
    for(int i=0;i<ns;i++) cpg_update_c_x_scaling(i,cx[i]);
    for(int i=0;i<mi;i++) cpg_update_c_u_scaling(i,cu[i]);

    double sps[STATES_SIZE],eps_s[STATES_SIZE];
    sc_x(start_pos,iSx,cx,sps);
    sc_x(end_pos,  iSx,cx,eps_s);
    for(int i=0;i<ns;i++) cpg_update_start_pos(i,sps[i]);
    for(int i=0;i<ns;i++) cpg_update_end_pos(i,eps_s[i]);

    double *x_opt =(double*)calloc(ns*(T+1),sizeof(double));
    double *u_opt =(double*)calloc(mi*T,    sizeof(double));
    double *vc_opt=(double*)calloc(ns*T,    sizeof(double));
    double *vi_opt=(double*)calloc(N_OBS*(T+1),sizeof(double));

    int no_first=0, iter=1;
    double rho_i=0.0;

    /* ── Bucle SCVx ─────────────────────────────────────────────────────────*/
    for(int cogu=0;cogu<COGU_max;cogu++){

        double t_disc_0 = now_ms();

        for(int t=0;t<=T;t++) isc_x(ox+t*ns,Sx,cx,ohx+t*ns);
        for(int t=0;t<T;t++)  isc_u(ou+t*mi,Su,cu,ohu+t*mi);

        for(int k=0;k<T;k++){
            double Phi[STATES_SIZE*STATES_SIZE];
            double G[STATES_SIZE*INPUTS_SIZE];
            double zv[STATES_SIZE];
            double Phi_cm[STATES_SIZE*STATES_SIZE];
            double G_cm[STATES_SIZE*INPUTS_SIZE];

            disc_ABY(ohx+k*ns,ohu+k*mi,(double)k/T*tf,dp,tau_d,
                     Sx,iSx,cx,Su,iSu,cu,sN,Phi,G,zv);

            to_colmaj(Phi,Phi_cm,ns,ns);
            to_colmaj(G,G_cm,ns,mi);

            for(int j=0;j<ns;j++)
                for(int i=0;i<ns;i++)
                    cpg_update_A_discrete(i+(k*ns+j)*ns, Phi_cm[j*ns+i]);

            for(int j=0;j<mi;j++)
                for(int i=0;i<ns;i++)
                    cpg_update_B_discrete(i+(k*mi+j)*ns, G_cm[j*ns+i]);

            for(int i=0;i<ns;i++)
                cpg_update_y_discrete(k*ns+i, zv[i]);
        }

        for(int k=0;k<=T;k++){
            double Co[N_OBS*STATES_SIZE];
            double Do[N_OBS*INPUTS_SIZE];
            double zo[N_OBS];
            double Co_cm[N_OBS*STATES_SIZE];
            double Do_cm[N_OBS*INPUTS_SIZE];
            double *uk_r=(k<T)?ohu+k*mi:ohu+(T-1)*mi;

            disc_CDZ(ohx+k*ns,uk_r,(double)k/T*tf,dp,Sx,cx,Su,cu,Co,Do,zo);

            to_colmaj(Co,Co_cm,N_OBS,ns);
            to_colmaj(Do,Do_cm,N_OBS,mi);

            for(int j=0;j<ns;j++)
                for(int i=0;i<N_OBS;i++)
                    cpg_update_C_discrete(i+(k*ns+j)*N_OBS, Co_cm[j*N_OBS+i]);

            for(int j=0;j<mi;j++)
                for(int i=0;i<N_OBS;i++)
                    cpg_update_D_discrete(i+(k*mi+j)*N_OBS, Do_cm[j*N_OBS+i]);

            for(int i=0;i<N_OBS;i++)
                cpg_update_z_discrete(k*N_OBS+i, zo[i]);
        }

        double t_disc_1 = now_ms();
        printf("Discretization time [ms]:  %.6f\n", t_disc_1 - t_disc_0);

        cpg_update_sqrt_tau(sqrt(tau_d));
        cpg_update_tau_lamb(tau_d*lam);
        cpg_update_etta(etta);

        for(int t=0;t<=T;t++)
            for(int i=0;i<ns;i++)
                cpg_update_ox_cvxpy(i+t*ns, ox[t*ns+i]);

        for(int t=0;t<T;t++)
            for(int i=0;i<mi;i++)
                cpg_update_ou_cvxpy(i+t*mi, ou[t*mi+i]);

        double t_solve_0 = now_ms();
        cpg_solve();
        double t_solve_1 = now_ms();
        printf("Solver time [ms]:  %.6f\n", t_solve_1 - t_solve_0);

        double L_opt = CPG_Result.info->obj_val;

        for(int t=0;t<=T;t++)
            for(int i=0;i<ns;i++)
                x_opt[t*ns+i] = CPG_Result.prim->x[i+t*ns];

        for(int t=0;t<T;t++)
            for(int i=0;i<mi;i++)
                u_opt[t*mi+i] = CPG_Result.prim->u[i+t*mi];

        for(int t=0;t<T;t++)
            for(int i=0;i<ns;i++)
                vc_opt[t*ns+i] = CPG_Result.prim->vc[i+t*ns];

        for(int t=0;t<=T;t++)
            for(int i=0;i<N_OBS;i++)
                vi_opt[t*N_OBS+i] = CPG_Result.prim->vi[i+t*N_OBS];

        double J_opt = J_cost(x_opt,u_opt,T,tau_d,dp,lam,Sx,cx,Su,cu,tf);
        double oJ    = J_cost(ox,   ou,   T,tau_d,dp,lam,Sx,cx,Su,cu,tf);
        double Delta_J=oJ-J_opt, Delta_L=oJ-L_opt;

        double nd=0;
        for(int t=0;t<=T;t++){
            double d[STATES_SIZE];
            for(int i=0;i<ns;i++) d[i]=x_opt[t*ns+i]-ox[t*ns+i];
            double n=vec_norm1(d,ns); if(n>nd) nd=n;
        }
        double vc_norm=0;
        for(int i=0;i<ns*T;i++) vc_norm+=fabs(vc_opt[i]);
        double vi_norm=0;
        for(int i=0;i<N_OBS*(T+1);i++) vi_norm+=fabs(vi_opt[i]);

        printf("Iteration number:  %d  Etta:  %.15g  Rho:  ", iter, etta);
        if(iter == 1) printf("None\n");
        else          printf("%.15g\n", rho_i);

        printf("L_SCVx_opt: %.15f J_SCVx_opt: %.15f oJ_SCVx: %.15f Norm_x_diff: %.15f \n\n",
            L_opt, J_opt, oJ, nd);

        if((Delta_L<e_tol*fabs(oJ)||nd<eps_stop)&&no_first){
            printf("\nCOGU convergio! :)\n"); break;
        }

        rho_i=(Delta_L!=0.0)?Delta_J/Delta_L:0.0;
        if(rho_i<rho0){
            etta=etta>etta0?etta/beta_sh:etta0;
        } else if(rho_i<rho1){
            etta=etta>etta0?etta/beta_sh:etta0;
            mat_copy(x_opt,ox,ns*(T+1)); mat_copy(u_opt,ou,mi*T);
        } else if(rho_i<rho2){
            mat_copy(x_opt,ox,ns*(T+1)); mat_copy(u_opt,ou,mi*T);
        } else {
            etta=etta<etta1?beta_gr*etta:etta1;
            mat_copy(x_opt,ox,ns*(T+1)); mat_copy(u_opt,ou,mi*T);
        }

        if(iter==3) no_first=1;
        iter++;
    }

    /* Trayectoria final sin escalar */
    double *xv=(double*)malloc(ns*(T+1)*sizeof(double));
    double *uv=(double*)malloc(mi*T*    sizeof(double));
    for(int t=0;t<=T;t++) isc_x(ox+t*ns,Sx,cx,xv+t*ns);
    for(int t=0;t<T;t++)  isc_u(ou+t*mi,Su,cu,uv+t*mi);

    printf("\n=== Trayectoria optima ===\n");
    printf("Posicion final: %.4f  %.4f  %.4f\n",xv[0+T*ns],xv[1+T*ns],xv[2+T*ns]);

    free(ohx);free(ohu);free(ox);free(ou);free(qt);free(av);
    free(x_opt);free(u_opt);free(vc_opt);free(vi_opt);
    free(xv);free(uv);
    return 0;
}
