# COGU Dev Library

Librería Python que toma un problema de trayectoria espacial definido con SymPy/CVXPY y genera código C autónomo que lo resuelve usando el algoritmo SCVx con ECOS. El resultado es un ejecutable C sin dependencia de Python — diseñado para correr en la unidad embebida COGU.

---

## Requisitos

### Python
- Python 3.10+
- `numpy`
- `sympy`
- `cvxpy`
- `cvxpygen`
- `scipy`

Instalar todo de una vez:
```
pip install numpy sympy cvxpy cvxpygen scipy
```

### Herramientas de compilación C (Windows via MSYS2)
- **GCC** — instalar desde [msys2.org](https://www.msys2.org), luego correr `pacman -S mingw-w64-ucrt-x86_64-gcc`
- **CMake** — `pacman -S mingw-w64-ucrt-x86_64-cmake`
- **Ninja** — `pacman -S mingw-w64-ucrt-x86_64-ninja`

Agregar MSYS2 al PATH antes de compilar:
```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
```

En Linux/macOS: instalar `gcc`, `cmake` y `ninja` con el gestor de paquetes del sistema.

---

## Inicio rápido — 3 pasos

Todos los comandos se corren desde la **raíz del repositorio** (`COGU/`).

### Paso 1 — Correr el pipeline Python

Resuelve la trayectoria Astrobee con SCVx y genera todos los archivos C en `COGU_Dev_Library/astrobee_c_output/`:

```
python COGU_Dev_Library/verify_final.py
```

Salida esperada:
```
COGU pipeline -- Verificacion final
...
[Python] Convergio en X iteraciones
[Python] Posicion final: 5.0000  2.0000  1.6000
Archivos C generados en: COGU_Dev_Library/astrobee_c_output/
```

Esto genera:
- `astrobee_c_output/dynamics.c / .h` — dinámica del satélite (SymPy → C)
- `astrobee_c_output/scvx_main.c` — loop SCVx completo
- `astrobee_c_output/cpg_compat.h` — adaptador del solver
- `astrobee_c_output/CMakeLists.txt` — archivo de build multiplataforma
- `astrobee_c_output/solver/` — código fuente de ECOS (copiado desde `references/solver_t30/`)

### Paso 2 — Compilar con CMake

```powershell
cmake -S COGU_Dev_Library/astrobee_c_output `
      -B COGU_Dev_Library/astrobee_c_output/build `
      -G "Ninja" -DCMAKE_C_COMPILER=gcc

cmake --build COGU_Dev_Library/astrobee_c_output/build
```

### Paso 3 — Correr

```
.\COGU_Dev_Library\astrobee_c_output\build\cogu_scvx.exe
```

Salida esperada:
```
COGU convergio! Posicion final: 5.0000  2.0000  1.6000
```

---

## VS Code

Si tienes la extensión CMake Tools instalada, apunta el directorio fuente a `COGU_Dev_Library/astrobee_c_output` y usa:
- `F7` — compilar
- `Shift+F5` — correr

---

## Estructura del proyecto

```
COGU_Dev_Library/
├── sympy_cvx_converter/      # Código fuente de la librería (Python)
│   ├── pipeline.py           # Punto de entrada: solve_trajectory()
│   ├── codegen.py            # Generación de código SymPy → C
│   ├── template_filler.py    # Rellena el template scvx_main.c
│   ├── problem.py            # Constructor del problema CVXPY
│   ├── solver.py             # Loop SCVx + discretización RK4
│   └── templates/
│       └── scvx_loop_template.c
├── examples/
│   └── astrobee_simple.py    # Ejemplo mínimo de uso
├── references/
│   ├── astrobee_t30/         # Loop C manual (referencia pre-pipeline)
│   ├── solver_t30/           # Solver ECOS pre-generado para T=30
│   └── test_cvxpygen_t30.py  # Script para regenerar solver_t30/ si es necesario
└── verify_final.py           # Verificación end-to-end (corre los pasos 1–3)
```

---

## Notas

- `references/solver_t30/` está commiteado en el repo para evitar volver a correr CVXPYgen (que puede dar `MemoryError` con T=30). No borrarlo.
- `astrobee_c_output/` está en el gitignore — se regenera completamente con `verify_final.py`.
- Para compilación cruzada a ARM (target embebido COGU): agregar `-DCMAKE_TOOLCHAIN_FILE=arm-none-eabi.cmake` al primer comando de CMake.
