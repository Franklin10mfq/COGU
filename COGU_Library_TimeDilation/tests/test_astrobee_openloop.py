"""
===========================================================================
.QUE HACE ESTA PRUEBA, EN CRISTIANO?
===========================================================================

El solver te entrega un PLAN: una lista de controles que, segun el, llevan al
Astrobee del punto A al punto B.

Esta prueba comprueba si ese plan se sostiene en la realidad. Toma esos
controles, los aplica a las ecuaciones de movimiento de verdad (integrandolas
aqui mismo, a mano, SIN usar la libreria), y mira donde termina el satelite.

Eso es lo que significa "lazo abierto": aplicar el plan a ciegas, sin ir
midiendo ni corrigiendo por el camino. Lo contrario seria "lazo cerrado",
donde un piloto automatico va ajustando segun lo que mide.

---------------------------------------------------------------------------
.POR QUE HACE FALTA, SI YA COMPARAMOS CONTRA EL NOTEBOOK Y CONTRA LA
LIBRERIA VIEJA?
---------------------------------------------------------------------------

Porque esas comparaciones responden a "..calcula lo mismo que antes?", no a
"..el resultado es correcto?". Si las tres implementaciones comparten el mismo
punto ciego, las tres coinciden y ninguna se entera.

Y hay un punto ciego concreto. El metodo SCvx solo obliga a que se cumplan
las ecuaciones en 31 instantes concretos (los "nodos"), y encima usa una
version simplificada de esas ecuaciones. Entre nodo y nodo no vigila nada.
Esta prueba mira precisamente lo que pasa ENTRE los nodos.

---------------------------------------------------------------------------
QUE MIDE, Y QUE SE ESPERA
---------------------------------------------------------------------------

1. .Llega al objetivo?
   El plan dice que llega exacto. Aplicandolo de verdad se llega con un error
   de decimas de milimetro. Se exige menos de 1 mm.

2. .Choca con los obstaculos ENTRE nodos?
   Si, un poco, y es inevitable: la optimizacion pega los nodos a la
   superficie de la esfera (holgura exactamente 0), asi que la curva que los
   une por fuerza corta por dentro. Son ~2.6 mm sobre un radio de 800 mm.
   Esto NO es un fallo de la libreria: pasa igual con el notebook de
   referencia (medido: 1.95 mm con su T=50). Es geometria.
   Se exige que no empeore de 1 cm, para que una regresion de verdad salte.

3. .Se pasa de los limites de velocidad entre nodos?
   No deberia, y no se pasa.

---------------------------------------------------------------------------
COMO CORRERLA
---------------------------------------------------------------------------
   Desde la raiz del paquete:   python tests/test_astrobee_openloop.py
===========================================================================
"""

import os
import sys
import warnings

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'examples'))

import numpy as np

from astrobee_model import (DYN_PARAMS_VAL, END_POS, OBSTACLES, OMEGA_MAX,
                            START_POS, VEL_MAX, astrobee_problem)
from scvx_toolbox import solve_trajectory

# --- umbrales -------------------------------------------------------------
TOL_POSICION_FINAL = 1e-3      # 1 mm
TOL_PENETRACION = 1e-2         # 1 cm; hoy son ~2.6 mm y es inherente
MARGEN_LIMITES = 1.02          # 2% de holgura sobre los limites de velocidad

# --- inercia, de DYN_PARAMS_VAL[0:6] = Ixx Iyy Izz Ixy Ixz Iyz ------------
_Ixx, _Iyy, _Izz, _Ixy, _Ixz, _Iyz = DYN_PARAMS_VAL[0:6]
I_MAT = np.array([[_Ixx, _Ixy, _Ixz],
                  [_Ixy, _Iyy, _Iyz],
                  [_Ixz, _Iyz, _Izz]])
I_INV = np.linalg.inv(I_MAT)


def f_astrobee(x, u):
    """
    Las ecuaciones de movimiento del Astrobee, escritas a mano aqui.

    Deliberadamente NO se usa la libreria para esto: si la libreria tuviera un
    error en las ecuaciones, usarla para comprobarse a si misma no detectaria
    nada.

    x = [posicion(3), velocidad(3), cuaternion(4), velocidad angular(3)]
    u = [aceleracion(3), torque(3)]
    """
    v = x[3:6]
    q0, q1, q2, q3 = x[6], x[7], x[8], x[9]
    w = x[10:13]
    a, tau = u[0:3], u[3:6]

    wx, wy, wz = w
    qdot = 0.5 * np.array([
        -wx * q1 - wy * q2 - wz * q3,
        +wx * q0 + wz * q2 - wy * q3,
        +wy * q0 - wz * q1 + wx * q3,
        +wz * q0 + wy * q1 - wx * q2,
    ])
    wdot = I_INV @ (tau - np.cross(w, I_MAT @ w))
    return np.concatenate([v, a, qdot, wdot])


def _rk4(x, u, dt):
    """Un paso de Runge-Kutta de orden 4: integra las ecuaciones un ratito."""
    k1 = f_astrobee(x, u)
    k2 = f_astrobee(x + 0.5 * dt * k1, u)
    k3 = f_astrobee(x + 0.5 * dt * k2, u)
    k4 = f_astrobee(x + dt * k3, u)
    return x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def propagar(x0, u, tf, T, subpasos=50):
    """
    Aplica el plan a ciegas y devuelve la trayectoria REAL, muestreada fina.

    El control se mantiene constante dentro de cada intervalo (eso es lo que
    significa ZOH), pero se integra con 50 pasitos por intervalo en vez de uno
    solo, para que la curva entre nodos se vea de verdad.
    """
    tau = tf / T
    h = tau / subpasos
    x = np.asarray(x0, dtype=float).reshape(-1).copy()
    traza = [x.copy()]
    for k in range(T):
        for _ in range(subpasos):
            x = _rk4(x, u[:, k], h)
            traza.append(x.copy())
    return np.array(traza).T          # (13, n_puntos)


def main():
    warnings.filterwarnings('ignore')
    T, tf = 30, 200.0

    print("Resolviendo (formulacion legacy: la misma del notebook y la "
          "libreria vieja)...")
    kwargs = astrobee_problem(T=T, tf=tf, time_mode='ZOH',
                              legacy_quaternion=True)
    r = solve_trajectory(size_N=20, max_iter=25, pin_aware_trust_region=False,
                         defect_penalty_space='unscaled', verbose=False,
                         **kwargs)
    x_nodos = r['x_opt_unscaled']
    u = r['u_opt_unscaled']
    print(f"  {r['status']}, {r['iterations']} iteraciones")

    xc = propagar(START_POS[:, 0], u, tf, T)
    print(f"  trayectoria real integrada con paso {tf / T / 50:.4f} s "
          f"(los nodos van cada {tf / T:.2f} s)")

    ok = True

    # --- 1. .llega? -------------------------------------------------------
    err = np.linalg.norm(xc[0:3, -1] - END_POS[0:3, 0])
    err_plan = np.linalg.norm(x_nodos[0:3, -1] - END_POS[0:3, 0])
    print()
    print("1. .LLEGA AL OBJETIVO?")
    print(f"   el plan dice que llega con error   {err_plan:.3e} m")
    print(f"   aplicandolo de verdad, el error es {err:.3e} m "
          f"({err * 1000:.4f} mm)")
    if err > TOL_POSICION_FINAL:
        print(f"   >>> FALLA: pasa del umbral de {TOL_POSICION_FINAL * 1000:.1f} mm")
        ok = False

    # --- 2. .entra en los obstaculos entre nodos? -------------------------
    print()
    print("2. .SE METE EN LOS OBSTACULOS ENTRE NODOS?")
    print(f"   {'obstaculo':<12} {'holgura en nodos':>18} {'holgura real':>15}")
    for i, (centro, radio) in enumerate(OBSTACLES, start=1):
        d_nodos = np.linalg.norm(x_nodos[0:3] - centro.reshape(3, 1), axis=0)
        d_real = np.linalg.norm(xc[0:3] - centro.reshape(3, 1), axis=0)
        cn, cc = d_nodos.min() - radio, d_real.min() - radio
        nota = f"   se mete {-cc * 1000:.2f} mm" if cc < 0 else ""
        print(f"   obstaculo {i:<3} {cn:18.5f} {cc:15.5f}{nota}")
        if cc < -TOL_PENETRACION:
            print(f"   >>> FALLA: pasa del umbral de {TOL_PENETRACION * 1000:.0f} mm")
            ok = False

    # --- 3. .respeta los limites de velocidad? ---------------------------
    print()
    print("3. .RESPETA LOS LIMITES DE VELOCIDAD ENTRE NODOS?")
    v_max = np.max(np.linalg.norm(xc[3:6], axis=0))
    w_max = np.max(np.linalg.norm(xc[10:13], axis=0))
    print(f"   max |v| = {v_max:.6f} m/s     (limite {VEL_MAX})")
    print(f"   max |w| = {w_max:.6f} rad/s   (limite {OMEGA_MAX:.6f})")
    if v_max > VEL_MAX * MARGEN_LIMITES or w_max > OMEGA_MAX * MARGEN_LIMITES:
        print("   >>> FALLA: se pasa de los limites")
        ok = False

    print()
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
