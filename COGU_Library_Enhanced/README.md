# COGU Library Enhanced

Versión mejorada de `COGU_Dev_Library`. Misma funcionalidad base — toma un problema de trayectoria definido en SymPy/CVXPY y genera código C autónomo que lo resuelve con SCVx + ECOS — pero con un pipeline más robusto y limpio.

---

## Diferencias respecto a `COGU_Dev_Library`

| | `COGU_Dev_Library` | `COGU_Library_Enhanced` |
|---|---|---|
| **Penalización `vi`** | Loop duplicado (fiel al notebook de referencia) | Loop único (SCVx canónico) |
| **Orden del pipeline** | Genera C en medio del solve Python | **Python primero, C solo si Python corre sin error** |
| **CVXPYgen** | Directo en el proceso principal (puede fallar por RAM fragmentada) | **Subproceso limpio** — el hijo arranca sin SymPy en memoria |
| **`solver_source_dir`** | Requerido como workaround de memoria | **No necesario** — se genera automáticamente |
| **Caché manual** | `regen_solver.py` para regenerar el solver | Sin caché manual — el pipeline lo gestiona solo |

### Por qué Python primero

CVXPYgen necesita un bloque de RAM contiguo de ~200-500 MB durante la canonicalización DPP. En el proceso principal, SymPy (cálculo de Jacobianos) ocupa y fragmenta la RAM antes de que CVXPYgen corra → `MemoryError`.

La solución: resolver en Python primero, y si no hay excepción, lanzar un **subproceso hijo limpio** que reconstruye solo el problema CVXPY y genera el C. El hijo arranca sin residuo de SymPy en memoria.

```
Pipeline Enhanced:
  [A] SymPy → funciones numpy        (solve Python)
  [D] build_problem() → CVXPY
  [E] solve_scvx()                   ← Python resuelve aquí
       └─ si no hubo excepción:
  [F] subprocess limpio → CVXPYgen   ← C se genera aquí
      dynamics.c, scvx_main.c, ...
```

---

## Requisitos

### Python
- Python 3.10+
- `numpy`, `sympy`, `cvxpy`, `cvxpygen`, `scipy`

```bash
pip install numpy sympy cvxpy cvxpygen scipy
```

### C (Windows via MSYS2)
```bash
pacman -S mingw-w64-ucrt-x86_64-gcc
pacman -S mingw-w64-ucrt-x86_64-cmake
pacman -S mingw-w64-ucrt-x86_64-ninja
```

Agregar al PATH antes de compilar:
```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
```

---

## Uso — 3 pasos

Todos los comandos desde la raíz del repo (`COGU/`).

### Paso 1 — Resolver y generar C

```bash
python COGU_Library_Enhanced/verify_final.py
```

Output esperado:
```
COGU pipeline (Enhanced) -- Verificacion final
...
Iter   1 | etta=... | L=... | J=...
...
Iter  25 | ...
[CVXPYgen] Python OK -> generando solver C en subproceso (RAM limpia)...
[CVXPYgen] Solver C generado en COGU_Library_Enhanced/astrobee_c_output/solver/c/
[template_filler] scvx_main.c generado  (NX=13, NU=6, NG=3, NP=18, T=30)
...
End pos: [5.  2.  1.6]
```

> Si el solve Python lanza una excepción (ECOS fatal, NaN, etc.), el pipeline
> se detiene **antes** de generar C — no se produce código de un solve fallido.

Genera en `astrobee_c_output/`:
- `dynamics.c / .h` — dinámica del satélite (SymPy → C)
- `scvx_main.c` — loop SCVx completo
- `cpg_compat.h` — adaptador del solver
- `CMakeLists.txt` — build multiplataforma
- `solver/` — solver ECOS generado por CVXPYgen

### Paso 2 — Compilar

```powershell
cmake -S COGU_Library_Enhanced/astrobee_c_output `
      -B COGU_Library_Enhanced/astrobee_c_output/build `
      -G "Ninja" -DCMAKE_C_COMPILER=gcc

cmake --build COGU_Library_Enhanced/astrobee_c_output/build
```

### Paso 3 — Ejecutar

```powershell
.\COGU_Library_Enhanced\astrobee_c_output\build\cogu_scvx.exe
```

Output esperado:
```
COGU convergio! Posicion final: 5.0000  2.0000  1.6000
```

---

## API

```python
from sympy_cvx_converter import solve_trajectory

result = solve_trajectory(
    states=states,               # lista de símbolos SymPy
    controls=controls,
    dynamics=dynamics,           # lista SymPy dx/dt = f(x, u)
    start=start_pos,             # (nx, 1) numpy array
    end=end_pos,
    T=30,                        # pasos de tiempo
    tf=200.0,                    # horizonte temporal
    nonconvex_constraints=...,   # opcional: g(x,u) <= 0
    dynamic_parameters_sym=...,
    dynamic_parameters_val=...,
    S_x=S_x, c_x=c_x,
    S_u=S_u, c_u=c_u,
    state_bounds=[...],
    control_bounds=[...],
    generate_c=True,             # genera C si Python no lanza excepción
    c_output_dir="mi_output/",   # carpeta destino
    # solver_source_dir ya no es necesario
)

x = result['x_opt_unscaled']    # trayectoria óptima (nx, T+1)
u = result['u_opt_unscaled']    # controles óptimos (nu, T)
```

---

## Estructura

```
COGU_Library_Enhanced/
├── sympy_cvx_converter/
│   ├── pipeline.py           # Entry point: solve_trajectory()
│   ├── _cpg_worker.py        # Worker del subproceso de CVXPYgen
│   ├── problem.py            # CVXPY problem builder (vi loop único)
│   ├── codegen.py            # SymPy → funciones numpy y C
│   ├── template_filler.py    # Genera scvx_main.c, cpg_compat.h
│   ├── solver.py             # Loop SCVx + discretización RK4
│   └── templates/
│       └── scvx_loop_template.c
├── examples/
│   └── astrobee_simple.py
├── references/
│   ├── astrobee_t30/         # Loop C manual (referencia histórica)
│   └── solver_t30/           # Solver ECOS pre-generado (backup opcional)
├── verify_final.py           # Verificación end-to-end
└── plot_trajectory.py        # Graficador de trayectoria
```

---

## Notas

- `astrobee_c_output/` está en `.gitignore` — se regenera corriendo `verify_final.py`.
- El solver en `references/solver_t30/` es un backup. El pipeline ya **no depende de él** — genera el solver automáticamente en cada corrida.
- Para compilación cruzada ARM (COGU embedded): agregar `-DCMAKE_TOOLCHAIN_FILE=arm-none-eabi.cmake` al primer comando CMake.
- `COGU_Dev_Library` se conserva sin cambios como referencia fiel al notebook de Mariel (penalización `vi` doble, workaround de caché manual).
