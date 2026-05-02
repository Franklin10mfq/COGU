import os

_TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), 'templates', 'scvx_loop_template.c')


def _make_switch(fname, stride, n_cases):
    """Genera una funcion cpg_update_X con switch de n_cases casos."""
    lines = [
        f'static inline void {fname}(int idx, double val) {{',
        f'    int k = idx / {stride};',
        f'    int inner = idx % {stride};',
        f'    switch(k) {{',
    ]
    for k in range(n_cases):
        lines.append(f'        case {k}: {fname.replace("cpg_update_", "cpg_update_", 1)}'
                     f'_{k}(inner, val); break;')
    lines += ['    }', '}', '']
    return '\n'.join(lines)


def generate_cpg_compat(nx, nu, ng, T, output_path):
    """
    Genera cpg_compat.h para dimensiones arbitrarias.

    Mapea las funciones stacked (indice global) a las per-timestep
    que genera CVXPYgen (cpg_update_A_0 ... cpg_update_A_{T-1}).

    Parametros:
        nx:          numero de estados
        nu:          numero de controles
        ng:          numero de restricciones no convexas
        T:           numero de pasos de tiempo
        output_path: ruta donde escribir el archivo
    """
    header = f"""\
/* cpg_compat.h — Adapta API stacked a per-timestep (generado automaticamente)
 *
 * Dimensiones: NS={nx}, MI={nu}, NG={ng}, T={T}
 */
#ifndef CPG_COMPAT_H
#define CPG_COMPAT_H

#include "cpg_solve.h"

/* Parametros horneados como constantes — no-op */
static inline void cpg_update_sqrt_tau(double val)              {{ (void)val; }}
static inline void cpg_update_tau_lamb(double val)              {{ (void)val; }}
static inline void cpg_update_vel_max(double val)               {{ (void)val; }}
static inline void cpg_update_omega_max(double val)             {{ (void)val; }}
static inline void cpg_update_acc_max(double val)               {{ (void)val; }}
static inline void cpg_update_torq_max(double val)              {{ (void)val; }}
static inline void cpg_update_S_x_scaling(int idx, double val)  {{ (void)idx; (void)val; }}
static inline void cpg_update_c_x_scaling(int idx, double val)  {{ (void)idx; (void)val; }}
static inline void cpg_update_S_u_scaling(int idx, double val)  {{ (void)idx; (void)val; }}
static inline void cpg_update_c_u_scaling(int idx, double val)  {{ (void)idx; (void)val; }}

"""

    def switch_block(outer, inner_prefix, stride, n_cases):
        """outer: nombre de la funcion wrapper; inner_prefix: prefijo de las CPG per-timestep."""
        lines = [f'static inline void {outer}(int idx, double val) {{']
        lines.append(f'    int k = idx / {stride};')
        lines.append(f'    int inner = idx % {stride};')
        lines.append('    switch(k) {')
        for k in range(n_cases):
            lines.append(f'        case {k}: {inner_prefix}_{k}(inner, val); break;')
        lines.append('    }')
        lines.append('}')
        lines.append('')
        return '\n'.join(lines)

    blocks = [header]
    blocks.append(switch_block('cpg_update_A_discrete', 'cpg_update_A', nx * nx,  T))
    blocks.append(switch_block('cpg_update_B_discrete', 'cpg_update_B', nx * nu,  T))
    blocks.append(switch_block('cpg_update_y_discrete', 'cpg_update_y', nx,       T))
    if ng > 0:
        blocks.append(switch_block('cpg_update_C_discrete', 'cpg_update_C', ng * nx,  T + 1))
        blocks.append(switch_block('cpg_update_D_discrete', 'cpg_update_D', ng * nu,  T + 1))
        blocks.append(switch_block('cpg_update_z_discrete', 'cpg_update_z', ng,       T + 1))
    blocks.append(switch_block('cpg_update_ox_cvxpy',   'cpg_update_ox', nx,       T + 1))
    blocks.append(switch_block('cpg_update_ou_cvxpy',   'cpg_update_ou', nu,       T))
    blocks.append('#endif /* CPG_COMPAT_H */\n')

    content = '\n'.join(blocks)
    with open(output_path, 'w', encoding='utf-8') as fp:
        fp.write(content)


def fill_scvx_template(nx, nu, ng, np_, T, output_dir):
    """
    Genera scvx_main.c y cpg_compat.h en output_dir.

    Parametros:
        nx:         numero de estados
        nu:         numero de controles
        ng:         numero de restricciones no convexas
        np_:        numero de parametros dinamicos
        T:          numero de pasos de tiempo
        output_dir: directorio de salida (se crea si no existe)
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(_TEMPLATE_PATH, 'r', encoding='utf-8') as fp:
        content = fp.read()

    content = content.replace('{NX}', str(nx))
    content = content.replace('{NU}', str(nu))
    content = content.replace('{NG}', str(ng))
    content = content.replace('{NP}', str(np_))
    content = content.replace('{T}',  str(T))

    with open(os.path.join(output_dir, 'scvx_main.c'), 'w', encoding='utf-8') as fp:
        fp.write(content)

    generate_cpg_compat(nx, nu, ng, T,
                        os.path.join(output_dir, 'cpg_compat.h'))

    print(f'[template_filler] scvx_main.c generado  (NX={nx}, NU={nu}, NG={ng}, NP={np_}, T={T})')
    print(f'[template_filler] cpg_compat.h generado')
