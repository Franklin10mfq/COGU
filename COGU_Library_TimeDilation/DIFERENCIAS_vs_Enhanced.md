# `scvx_toolbox` frente a `COGU_Library_Enhanced` — análisis de diferencias

Comparación entre la librería original (`sympy_cvx_converter`, en
`COGU_Library_Enhanced/`) y la que escribió Franklin (`scvx_toolbox`, esta).

Todo lo medido aquí sale de la librería de Franklin **tal como la entregó**: los
once módulos que producen números (`pipeline`, `problem`, `solver`,
`discretization`, `symbolic`, `cost_dsl`, `scaling`, `modes`, `quaternion`,
`utils`, `_cpg_worker`) están byte a byte iguales a su entrega. Lo único que se
tocó fue la generación de C (`cgen.py`, `templates/scvx_loop_template.c`) y la
restauración de `scvx_toolbox/__init__.py`, que faltaba.

---

## 1. Estructura

| `sympy_cvx_converter` | `scvx_toolbox` | qué cambió |
| --- | --- | --- |
| `pipeline.py` (265) | `pipeline.py` (311) | + modos de tiempo, escalado de `p` |
| `problem.py` (230) | `problem.py` (330) | + `F`/`G`/`op` como parámetros, cotas de `tf` |
| `solver.py` (487) | `solver.py` (400) + `discretization.py` (241) | el RK4 se separó del bucle SCvx |
| `translator.py` (380) | `symbolic.py` (296) | mismo trabajo, + Jacobianas de dilatación |
| `codegen.py` (473) + `template_filler.py` (173) | `cgen.py` (1173) | generación de C reescrita |
| `cost_dsl.py` (41) | `cost_dsl.py` (97) | + `k_range='U'`, resolución compartida |
| `utils/quaternion.py` (73) | `quaternion.py` (103) | + corrección de convención |
| — | `modes.py` (104) | nuevo: los cuatro modos y sus alias |

## 2. API

`solve_trajectory` es un **superconjunto estricto**: 34 parámetros en común,
3 renombrados (`etta0/etta1/etta_init` → `eta0/eta1/eta_init`, con los nombres
viejos aceptados como alias) y 14 nuevos:

`time_mode`, `tf_init`, `tf_min`, `tf_max`, `alpha_p`, `S_p`, `c_p`,
`dilation_symbol`, `pin_aware_trust_region`, `defect_penalty_space`,
`defect_substeps`, `max_stall`, `solver_options`, `**legacy`.

Es decir: **cualquier llamada válida en Enhanced sigue siendo válida aquí.**

## 3. Diferencias de comportamiento (verificadas en código, no en el README)

### 3.1 Región de confianza y nodos fijados — `pin_aware_trust_region`

`x[:,0]` y `x[:,T]` están fijados por restricciones de igualdad, así que
`‖x_0 − ox_0‖∞` es una **constante**. Enhanced la cobra contra `eta`; si el warm
start discrepa de `start`/`end`, esa constante es un piso duro bajo `eta` y el
subproblema se vuelve infactible cuando la región se encoge por debajo.

- Enhanced: siempre la cobra.
- `scvx_toolbox`: `pin_aware_trust_region=True` por defecto la omite.
- **Efecto en Astrobee: ninguno.** Su warm start ya coincide en ambos extremos.

### 3.2 Test de convergencia — con signo vs valor absoluto

```python
# Enhanced,  solver.py:437
if no_first_iterations and (Delta_L < e_tol * np.abs(oJ_SCVx) or ...):
# scvx_toolbox, solver.py:333
if past_startup and (abs(delta_L) < e_tol * abs(J_ref) or ...):
```

Un paso rechazado tiene `Delta_L` muy negativo, con lo que el test con signo se
cumple trivialmente y declara convergencia justo cuando la linealización es menos
confiable. **Es una corrección real.** En Astrobee no cambia nada porque ahí
`delta_L` se mantiene positivo (~28) todo el tiempo.

### 3.3 Unidades del defecto en `J` — `defect_penalty_space`

`J` mide el defecto en unidades físicas; `L` penaliza el slack `vc` en unidades
escaladas. Con un `S_x` agresivo no están en el mismo pie y `rho` puede salir
negativo en pasos que sí mejoran `J`.

- Enhanced: siempre `'unscaled'`.
- `scvx_toolbox`: `'unscaled'` por defecto (reproduce Enhanced), `'scaled'` opcional.
- **Efecto en Astrobee: sólo cambia las unidades de `J`** (753 → 762). Su `S_x` es
  suave. Esos dos números **no son comparables entre sí**.

### 3.4 Convención de cuaternión en el warm start

Enhanced hace `R.from_quat(quats)` pasando cuaterniones *scalar-first*
`(q0,q1,q2,q3)` a una función de scipy que los espera *scalar-last*. La velocidad
angular del warm start no es la del slerp.

- `scvx_toolbox`: `angular_velocity(..., scalar_first=True)` por defecto lo corrige;
  `scalar_first=False` reproduce el comportamiento viejo.
- **Para paridad con Enhanced o con el notebook hay que pasar `scalar_first=False`.**

### 3.5 Subproblema infactible y estancamiento

- Enhanced: no comprueba `x.value is None`. Un subproblema infactible tira la
  corrida entera.
- `scvx_toolbox`: lo detecta, conserva el último iterado aceptado y reporta
  `status='subproblem_infeasible'`. Añade `max_stall` para cortar cuando el paso
  se sigue rechazando con `eta` ya en el piso.

### 3.6 Modos de tiempo (lo nuevo de verdad)

`FOH`, `ZOH+dilation`, `FOH+dilation`. Con dilatación el problema se replantea en
`s ∈ [0,1]` con `t = p·s`, `p ≡ tf` como variable de decisión, Jacobianas extra
`F = ∂f/∂p` y `G = ∂g/∂p`, y una ODE más en el RK4 (`H' = AH + F`).

---

## 4. Resultados medidos

### 4.1 Astrobee — paridad con las dos referencias

Referencias, ambas verificadas:

1. `COGU_Library_Enhanced/examples/astrobee_simple.py` (T=30, tf=200)
2. `COGU_algorithms/Astrobee_template/SCVx_official_template_cvxpy_cvxpygen_ZOH.ipynb`
   (**T = 51−1 = 50**, tf=200). La copia en `ASTROBEE_SCVx_TEMPLATE_C/` es
   numéricamente idéntica.

| referencia | peor diferencia relativa | test |
| --- | --- | --- |
| Enhanced, T=30, it. 1–12 | 1.0e-06 (it. 1: **1e-12**) | `tests/test_astrobee_vs_enhanced.py` |
| Notebook, T=50, it. 1–8 | 1.1e-05 | `tests/test_astrobee_zoh.py` |
| Python vs C generado, it. 1–4 | 1.1e-07 | `tests/test_astrobee_c.py` |

Para reproducir cualquiera de las dos hacen falta los tres interruptores legacy:
`pin_aware_trust_region=False`, `defect_penalty_space='unscaled'` y el warm start
con la convención vieja de cuaternión.

**Por qué antes parecía que Astrobee fallaba:** el ejemplo viejo
`COGU_Library_Enhanced/examples/astrobee_simple.py` corre `T=30` pero imprime como
referencia los números del notebook, que usa `T=50`. Nunca fue un fallo de ninguna
librería.

### 4.2 Ninguna de las dos converge en Astrobee

Las dos agotan `max_iter=25` con `rho` negativo. `delta_L = J_ref − L` se estanca
en ~28 porque el defecto de discretización ZOH a `tau = 6.67 s` no baja, mientras
el criterio lo compara contra `e_tol·|J_ref| ≈ 3.8`. El **paso** sí converge
(`‖dx‖ = 1e-5` en la it. 25) y el estado terminal se alcanza a 1e-9 m: la
trayectoria está convergida aunque la bandera diga otra cosa.

Ninguna de las correcciones de 3.1–3.4 lo cambia — ver
`examples/astrobee_compare_settings.py`.

### 4.3 Coste de la reescritura, a igualdad de problema

Astrobee T=30 ZOH, formulación idéntica (interruptores legacy puestos), un proceso
por librería y ejecutados en serie:

| librería | wall [s] | it | s/iter | J final |
| --- | --- | --- | --- | --- |
| `sympy_cvx_converter` (Enhanced) | 38.79 | 25 | 1.55 | 753.5529 |
| `scvx_toolbox` | 43.18 | 25 | 1.73 | 753.1169 |

**Mismas iteraciones, +11.3% de tiempo.** La reescritura cuesta algo de overhead
por iteración (el discretizador es más general: lleva las ramas de FOH y de
dilatación aunque estén apagadas). No cambia el resultado.

Nota: no se puede cruzar BIRDS entre las dos librerías porque **son modelos
distintos** — `COGU_Library_Enhanced/examples/birds_tm.py` tiene 13 estados con
ruedas de reacción (`hw1,hw2,hw3`) y torques de rueda, y el de aquí es VAT-only
con 7 estados y sólo thrusters. Astrobee sí es idéntico en ambas, y por eso es el
banco de pruebas correcto para el overhead.

### 4.4 Qué cuesta y qué compra la dilatación del tiempo

BIRDS-VAT, misma librería, mismo modelo, mismo warm start, mismo ajuste de SCvx;
lo único que cambia es `time_mode` (`examples/birds_vat_timing.py`):

| modo | it | tf [s] | wall [s] | s/iter |
| --- | --- | --- | --- | --- |
| ZOH (fijo) | 18 | 480.000 | 27.5 | 1.53 |
| ZOH+dilation | 18 | **420.244** | 28.6 | 1.59 |

**La dilatación no acelera el cálculo.** Mismas iteraciones, 4.1% más lento por
iteración — lo esperable, porque añade una variable, una Jacobiana por nodo, una
ODE en el RK4 y dos restricciones. Lo que compra es el `tf` óptimo: la maniobra
baja de 480 s a 420.2 s, **12.4% menos tiempo de maniobra física**.
