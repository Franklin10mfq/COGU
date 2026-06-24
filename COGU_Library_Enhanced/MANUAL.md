# Manual de usuario — COGU_Library_Enhanced

Librería para calcular **trayectorias óptimas** de un vehículo (satélite, robot) mediante
**SCVx** (Successive Convexification). El usuario define el problema en SymPy + numpy y la librería
resuelve en Python y, opcionalmente, **genera código C embebible** para correr a bordo.

> Astrobee y BIRDS-TM son **ejemplos** incluidos, no casos especiales: la librería es genérica.

---

## 1. Introducción

COGU resuelve el problema: *"¿cómo debe moverse el vehículo del estado inicial al final, gastando el
mínimo esfuerzo, respetando su física y sus límites?"*. Internamente linealiza la dinámica, resuelve
un problema convexo en cada iteración (con región de confianza) y repite hasta converger.

Flujo de dos etapas:
- **Python** — prototipo y validación en tierra.
- **C** — código autocontenido que corre en la computadora embebida (se genera solo si lo pides).

---

## 2. Requisitos

| Entorno | Paquetes / herramientas |
|---|---|
| Python | `cvxpy`, `cvxpygen`, `ecos`, `sympy`, `numpy`, `scipy` |
| C (solo si `generate_c=True`) | `gcc`, `cmake`, `ninja` (en Windows: msys64/ucrt64) |

---

## 3. Diagrama de bloques

```
            INPUTS                          PROCESO  (solve_trajectory)            OUTPUTS
 ┌──────────────────────────┐     ┌────────────────────────────────────┐  ┌────────────────────┐
 │ Modelo SymPy             │     │ 1. codegen: SymPy → funciones num. │  │ Python:            │
 │  states, controls,       │ ──▶ │ 2. build_problem: arma CVXPY       │─▶│  result (dict):    │
 │  dynamics, dyn_params    │     │ 3. solve_scvx: loop SCVx (ECOS)    │  │  x_opt_unscaled,   │
 │ Frontera (start, end)    │     │ 4. [si generate_c] genera código C │  │  u_opt, converged, │
 │ Costo (cost_terms)       │     │    en subproceso (CVXPYgen)        │  │  iterations,history│
 │ Constraints (bounds)     │     └────────────────────────────────────┘  ├────────────────────┤
 │ Escalado (S_x,c_x,…)     │                                             │ C (opcional):      │
 │ Warm start               │                                             │  scvx_main.c,      │
 │ Params algoritmo (T,tf…) │                                             │  dynamics.c/.h,    │
 └──────────────────────────┘                                             │  cpg_compat.h,     │
                                                                          │  solver/, CMake →  │
                                                                          │  ejecutable        │
                                                                          └────────────────────┘
```

El **único punto de entrada** es `solve_trajectory(...)` (en `sympy_cvx_converter/pipeline.py`).
Resuelve primero en Python; solo si el solve corre sin error **y** `generate_c=True`, genera el C.

---

## 4. Inputs — referencia de `solve_trajectory(...)`

### Obligatorios

| Parámetro | Tipo | Qué es |
|---|---|---|
| `states` | lista de símbolos SymPy | Variables de estado `x` (ej. posición, velocidad, quaternion). `nx = len(states)` |
| `controls` | lista de símbolos SymPy | Variables de control `u` (ej. aceleración, torque). `nu = len(controls)` |
| `dynamics` | lista SymPy | `dx/dt = f(x,u)`, una expresión por estado |
| `start`, `end` | numpy `(nx,1)` | Estado inicial y final |
| `T` | int | Número de pasos de tiempo |
| `tf` | float | Tiempo final (s). `tau = tf/T` |
| `cost_terms` | lista de dicts | **Función de costo** → sección 5. **Obligatorio** (sin default) |

### Opcionales

| Parámetro | Default | Qué es |
|---|---|---|
| `dynamic_parameters_sym` / `_val` | None | Símbolos SymPy de parámetros físicos (inercia, etc.) y sus valores numpy |
| `nonconvex_constraints` | None | Restricciones no convexas `g(x,u) ≤ 0` (ej. obstáculos) |
| `state_bounds`, `control_bounds` | None | Límites sobre estado/control → sección 6 |
| `S_x`, `c_x`, `S_u`, `c_u` | identidad | Matrices de escalado `x_scaled = S_x⁻¹(x − c_x)` |
| `warm_start_x` `(nx,T+1)`, `warm_start_u` `(nu,T)` | interpolación | Trayectoria inicial. Para actitud usa `utils.slerp` / `utils.angular_vel` |
| `size_N` | 20 | Sub-pasos RK4 en la discretización |
| `lamb` | 1000 | Peso de la penalización de defecto (vc) |
| `max_iter` | 25 | Iteraciones SCVx |
| `rho0,rho1,rho2` | 0, 0.1, 0.7 | Umbrales de la región de confianza |
| `etta0,etta1,etta_init` | 1e-8, 10, 1 | Límites y valor inicial del radio de confianza |
| `beta_sh,beta_gr` | 2, 2 | Factores de contracción/expansión del radio |
| `e_tol,epsilon_stop` | 0.005, 1e-5 | Tolerancias de convergencia |
| `solver` | `'ECOS'` | Solver convexo |
| `generate_c` | False | Si True, genera el C tras un solve Python exitoso |
| `c_output_dir` | `"c_output"` | Carpeta destino del C |

---

## 5. `cost_terms` — la función de costo (sección clave)

**El usuario siempre debe declarar el costo.** No hay costo "de fábrica": si omites `cost_terms`,
la librería lanza un error claro. El costo es una **lista de términos**; el costo total es la suma.

### Formato de un término

```python
{
  'kind':    'sumsq' | 'norm1' | 'norm2',   # tipo: cuadrado, norma-1, norma-2
  'var':     'u' | 'x',                       # penaliza control o estado
  'slice':   slice(a, b),                     # qué componentes (default: todos)
  'coeff':   'sqrt_tau' | 'tau' | 'tau_lamb' | float,   # escalar multiplicador (default 'sqrt_tau')
  'weight':  float,                            # peso del término (default 1.0)
  'offset':  numpy (n,1) | None,               # objetivo a restar: penaliza (var − offset)
  'k_range': 'T' | 'T+1',                      # rango temporal (default: 'T' para u, 'T+1' para x)
}
```

- `sumsq` → `weight · Σ‖coeff·(var−offset)‖²`
- `norm1` → `weight · Σ‖coeff·(var−offset)‖₁`
- `norm2` → `weight · Σ‖coeff·(var−offset)‖₂`

(definición única en `sympy_cvx_converter/cost_dsl.py`, usada por las versiones CVXPY, numpy y C).

### Cómo sacarlo del notebook de Jupyter

La función objetivo del notebook se traduce **1-a-1**. Ejemplo real de BIRDS-TM:

| En el notebook | `cost_term` equivalente |
|---|---|
| `cp.sum_squares(sqrt_tau*u[0:3])` | `{'kind':'sumsq','var':'u','slice':slice(0,3),'coeff':'sqrt_tau'}` |
| `1.0*cp.norm(tau*u[3:7],1)` | `{'kind':'norm1','var':'u','slice':slice(3,7),'coeff':'tau','weight':1.0}` |
| `100*cp.sum_squares(sqrt_tau*x[0:10]-sqrt_tau_xq_goal)` | `{'kind':'sumsq','var':'x','slice':slice(0,10),'coeff':'sqrt_tau','weight':100,'offset':x_goal,'k_range':'T+1'}` |

> **No incluyas** las penalizaciones de defecto (`vc`) ni de obstáculos (`vi`): esas son parte
> interna del algoritmo SCVx y la librería las añade automáticamente.

El `offset` (`x_goal`) es el estado objetivo **ya escalado**: `inv(S_x) @ (end − c_x)`.

---

## 6. Constraints — `state_bounds` / `control_bounds`

Cada bound limita un tramo (`slice`) de la variable. Dos formatos conviven:

**Tuple (norma con tope):**
```python
(slice(3,6), 'norm2', vel_max)        # ‖x[3:6]‖₂ ≤ vel_max
```

**Dict (más expresivo):**
```python
{'kind':'norm', 'slice':slice(0,3), 'norm_type':'norminf', 'limit':torq_max}
{'kind':'box',  'slice':slice(3,7), 'lower':0.0, 'upper':Fth_max}     # 0 ≤ u[3:7] ≤ Fth_max
{'kind':'norm', 'slice':slice(4,7), 'norm_type':'norm2', 'limit':omega_max,
                'S':inv_Jb_S, 'c':inv_Jb_c}   # con transformación propia (no el escalado global)
```

`norm_type`: `'norm1'` | `'norm2'` | `'norminf'`. Los campos `S`/`c` opcionales permiten una matriz
de transformación distinta del escalado global (ej. `inv(J)` para pasar de momento a velocidad
angular en BIRDS). Definición en `_build_bound_constraints` (problem.py).

Ejemplo de traducción (BIRDS-TM):

| En el notebook | bound |
|---|---|
| `norm(S_u[0:3]@u[0:3]+c_u[0:3],'inf') ≤ torqW_max` | `(slice(0,3),'norminf',torqW_max)` |
| `0 ≤ S_u[i]@u[i]+c_u[i] ≤ Fth_max` | `{'kind':'box','slice':slice(3,7),'lower':0.0,'upper':Fth_max}` |
| `norm(S_u[3:5]@u[3:5]+c_u[3:5],1) ≤ Fth_max` | `{'kind':'norm','slice':slice(3,5),'norm_type':'norm1','limit':Fth_max}` |
| `norm(inv(J)@S_x[4:7]@x[4:7]+...,2) ≤ omegab_max` | `{'kind':'norm','slice':slice(4,7),'norm_type':'norm2','limit':omegab_max,'S':inv_Jb_S,'c':inv_Jb_c}` |

---

## 7. Outputs

### Python — `result` (dict)

| Clave | Qué es |
|---|---|
| `x_opt_unscaled` | Estados óptimos `(nx, T+1)` en unidades físicas |
| `u_opt_unscaled` | Controles óptimos `(nu, T)` en unidades físicas |
| `x_opt`, `u_opt` | Igual pero escalados |
| `converged` | bool |
| `iterations` | nº de iteraciones realizadas |
| `history` | lista por iteración con `L_SCVx`, `J_SCVx`, etc. (para comparar con la referencia) |

### C (si `generate_c=True`) — en `c_output_dir/`

| Archivo | Qué es |
|---|---|
| `scvx_main.c` | Loop SCVx + datos del problema embebidos + J_cost generado |
| `dynamics.c` / `dynamics.h` | `f`, jacobianos `A`/`B`, etc. (de la dinámica SymPy) |
| `cpg_compat.h` | Adaptador a la API del solver CVXPYgen |
| `solver/` | Solver convexo generado por CVXPYgen |
| `CMakeLists.txt` | Para compilar |

Al compilar y correr, el ejecutable imprime los `L/J` por iteración y el **Estado final**.

---

## 8. Tutorial — añadir un problema nuevo

Usa `examples/birds_tm.py` como plantilla. Pasos:

1. **Modelo SymPy.** Define `states`, `controls`, `dynamics` (una expresión `dx/dt` por estado) y,
   si aplica, `dynamic_parameters_sym/val`. Copia la dinámica **tal cual** del notebook (no la reinventes).
2. **Frontera y tiempo.** `start`, `end` (nx,1); `T`, `tf`.
3. **Escalado.** `S_x`, `c_x`, `S_u`, `c_u` (cópialos del notebook; si no hay, deja identidad).
4. **Warm start.** Construye `warm_start_x/u`. Para actitud usa `utils.slerp` y `utils.angular_vel`.
5. **Costo.** Traduce la función objetivo del notebook → `cost_terms` (sección 5).
6. **Restricciones.** Traduce las constraints → `state_bounds` / `control_bounds` (sección 6).
7. **Resuelve en Python** con `generate_c=False`:
   ```python
   result = solve_trajectory(states, controls, dynamics, start, end, T, tf,
                             cost_terms=cost_terms, state_bounds=..., control_bounds=...,
                             S_x=S_x, c_x=c_x, S_u=S_u, c_u=c_u,
                             warm_start_x=warm_x, warm_start_u=warm_u,
                             lamb=..., max_iter=..., generate_c=False)
   ```
8. **Valida** contra el notebook: compara `result['history']` (L/J por iteración) y el estado final.
9. **Genera el C**: cambia a `generate_c=True` y `c_output_dir="..."`. Compila y corre:
   ```powershell
   python tu_problema.py
   $env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
   cmake -S <dir> -B <dir>/build -G "Ninja" -DCMAKE_C_COMPILER=gcc
   cmake --build <dir>/build
   .\<dir>\build\cogu_scvx.exe
   ```

---

## 9. Alcances y limitaciones

**Alcances**
- Cualquier problema de trayectoria óptima con **dinámica continua** definida en SymPy.
- Costo combinando términos `sumsq` / `norm1` / `norm2`, con coeficientes, pesos y objetivo (`offset`).
- Restricciones de norma (`1`/`2`/`inf`), de **rango/caja**, y con transformación propia.
- Obstáculos no convexos (`nonconvex_constraints`).
- Genera **C embebible** del problema completo. Astrobee y BIRDS-TM son ejemplos, **no** hardcodeo.

**Limitaciones (ser consciente)**
- Algoritmo fijo: **SCVx** con región de confianza; solver **ECOS**.
- El costo se limita a `sumsq`/`norm1`/`norm2` (no cualquier expresión arbitraria).
- El `offset` (objetivo de seguimiento) es **constante**, no se puede cambiar en caliente sin
  reconstruir el problema.
- **Parámetros del algoritmo no ajustables desde el C.** *Toda la física de tu problema* (dinámica,
  costo, restricciones, frontera, escalado, warm start, `lam`, `tf`) **sí** se embebe: el C resuelve
  tu problema, cualquiera que sea. Lo que **no** se embebe son los parámetros del *algoritmo* SCVx
  (`max_iter`, `rho*`, `etta*`, `beta_*`, tolerancias): el C los usa con sus **valores por defecto**,
  aunque los cambies en Python. Solo importa si tu problema necesita valores distintos a los defaults
  para converger; en ese caso, edítalos a mano en el `scvx_main.c` generado (≈4 líneas) o amplía la
  librería para embeberlos (deuda técnica, igual que se hizo con `lam`/`tf`). Astrobee y BIRDS usan
  los defaults, por eso no les afecta.
- El problema debe ser **DPP-compliant** para que CVXPYgen genere el solver C.
- Los escalares del costo se hornean como `float` (no `cp.Parameter`) por un bug de CVXPY 1.8.
- Generar el solver C (CVXPYgen) puede tardar varios minutos.

---

## 10. Referencia rápida

```powershell
# Astrobee (Python + C, validación de referencia)
python COGU_Library_Enhanced/verify_final.py

# BIRDS-TM (Python + C)
python COGU_Library_Enhanced/examples/birds_tm.py

# Compilar y correr un C generado (ej. Astrobee)
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
cmake -S COGU_Library_Enhanced/astrobee_c_output -B COGU_Library_Enhanced/astrobee_c_output/build -G "Ninja" -DCMAKE_C_COMPILER=gcc
cmake --build COGU_Library_Enhanced/astrobee_c_output/build
.\COGU_Library_Enhanced\astrobee_c_output\build\cogu_scvx.exe
```

**Plantillas:** `examples/astrobee_simple.py` (con obstáculos, `ng>0`) y `examples/birds_tm.py`
(actitud pura, `ng=0`, con `offset` y transformación propia). Empieza copiando la más parecida a tu problema.
