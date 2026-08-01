# Waskiri-Sat — Subsistema COGU
## Documento de Control de Interfaces (ICD)

> **Borrador en español.** Este documento se traducirá al inglés una vez revisado y
> completado. Los campos marcados con **[A COMPLETAR]** requieren datos que no están
> definidos en el repositorio de software (dimensiones físicas, conectores, voltajes,
> corrientes) — corresponden a decisiones de hardware/PCB que no están a cargo del
> autor de este borrador. La estructura de este documento sigue el formato ya usado
> por el equipo en el ICD de OBC-CDH (`Waskiri_Sat_OBC_ICD_260705`), para mantener
> consistencia entre los ICD de los distintos subsistemas del satélite.

---

## Registro de cambios

| Versión | Fecha | Modificado por | Comentarios |
|---|---|---|---|
| 0.1 | 2026-07-04 | [Nombre del usuario] | Borrador inicial en español |
| 0.2 | 2026-07-05 | [Nombre del usuario] | Reestructurado siguiendo el formato del ICD de OBC-CDH: se agregó la tabla de asignación de requisitos (SR/DR/VR), se separaron las especificaciones de interfaz por subsistema con el mismo campo estándar (Interface Name, Connection Type, etc.), se confirmó I2C como protocolo COGU↔ADCS, y se convirtió el diagrama de flujo a formato de flowchart |

---

## Tabla de contenidos

1. Alcance del documento
2. Introducción
3. Tabla de asignación de requisitos
4. Identificación de interfaces
   - 4.1 Interfaces mecánicas
   - 4.2 Interfaces eléctricas
   - 4.3 Interfaces de software
5. Especificaciones de interfaz por subsistema
   - 5.1 COGU <-> ADCS
   - 5.2 COGU <-> EPS

Apéndice A — Cómo funciona COGU por dentro
Apéndice B — DSL de costo y restricciones

---

## 1. Alcance del documento

Este documento define las interfaces técnicas entre el subsistema **COGU** y los demás
subsistemas de Waskiri-Sat.

---

## 2. Introducción

### 2.1 ¿Qué es COGU?

COGU es el subsistema de cómputo encargado de calcular, en tiempo real y a bordo, la
**maniobra de actitud óptima** del satélite. Recibe el estado actual y el estado
objetivo de actitud (provistos por el ADCS) y devuelve la secuencia de control óptima
que debe ejecutar el actuador — en el caso de Waskiri-Sat, el **dipolo magnético** de
los magnetorquers (Waskiri-Sat no usa ruedas de reacción como actuador).

### 2.2 Funcionalidad principal

- Resolver, a bordo, un problema de **optimización convexa** (Sequential Convex
  Programming / SCVx) que minimiza el consumo de energía del actuador y el tiempo de
  asentamiento de la maniobra, sujeto a los límites físicos del satélite (velocidad
  angular máxima, dipolo magnético máximo).
- Incorporar el modelo del campo magnético terrestre a lo largo de la órbita (World
  Magnetic Model, WMM) para calcular el torque disponible en cada instante
  (`torque = dipolo × campo magnético`).
- Entregar la secuencia óptima de control al **ADCS** — el único subsistema con el
  que COGU tiene interfaz directa (ver §4 y §5).
- El algoritmo corre embebido en C, generado a partir de un modelo validado primero en
  Python/CVXPY, para que el código a bordo sea liviano y determinístico en tiempo de
  ejecución.

### 2.3 Diagrama de bloques del subsistema (conceptual)

```
  ┌─────────┐   comando: activar ADCS  ┌─────────┐   estado actual (q, ω/h_B)   ┌─────────┐
  │         │ ────────────────────────▶│         │   vía I2C                    │         │
  │   OBC   │   y solicitar envío      │  ADCS   │ ─────────────────────────────▶│  COGU   │
  │         │   de datos a COGU        │         │   [futuro] warm start (NN)    │         │
  └─────────┘                          └─────────┘ ─────────────────────────────▶└─────────┘
                                                    ◀─────────────────────────────
                                                     secuencia óptima de control m(t)
                                                     (vía I2C — ver §5.1)
```

**Confirmado por el equipo:** COGU se comunica **únicamente con el ADCS, vía I2C**. El
ADCS es el puente de comunicación entre COGU y el OBC (a través de su propio CPLD /
Mission Boss) — **no existe una interfaz directa COGU↔OBC**. El OBC activa al ADCS por
comando; el ADCS, una vez activado, es quien efectivamente envía los datos a COGU y
recibe de vuelta la secuencia óptima de control.

El **modelo del satélite, el costo, las restricciones y el objetivo de la maniobra ya
están embebidos** en el programa COGU — lo único que varía en tiempo real es la
**condición inicial** (actitud y momento angular actuales, medidos por el ADCS) y, a
futuro, el **warm start** (trayectoria inicial sugerida, generada por una red
neuronal). Detalle completo en §4.3.4.

### 2.4 Composición interna

- **Boards:** `[A COMPLETAR]` — a definir con el equipo de hardware/PCB (no es parte de
  este borrador porque esa documentación la llevan mis compañeros).
- **Especificaciones de componentes principales:** `[A COMPLETAR]` (procesador/MCU,
  memoria, etc.)
- **Especificaciones de conectores:** `[A COMPLETAR]` (fabricante, modelo, número de
  pines).

---

## 3. Tabla de asignación de requisitos

Sigue el mismo esquema **SR → DR → VR → Método de verificación** que el ICD de
OBC-CDH, para mantener trazabilidad consistente entre subsistemas de Waskiri-Sat.
Varias filas quedan `TBD` porque dependen de la interfaz física/protocolo con el ADCS
(§4.2/§5.1), aún no definida con el equipo de hardware.

| ID | Requisito de Sistema | ID | Requisito de Diseño | ID | Requisito de Verificación | Método de Verificación |
|---|---|---|---|---|---|---|
| SR1 | Calcular la maniobra de actitud óptima | DR1.1 | Implementar el algoritmo SCVx reproduciendo el modelo del notebook de referencia | VR1.1.1 | Los valores de costo (`L_SCVx`, `J_SCVx`) por iteración coinciden con el notebook de referencia | Correr `solve_trajectory()` en Python y comparar contra el notebook — **ya hecho para Astrobee y BIRDS-TM** |
| | | DR1.2 | Validar visualmente la trayectoria óptima | VR1.2.1 | Las gráficas de la trayectoria (actitud, velocidad angular, control) coinciden con lo esperado físicamente | Graficar y comparar visualmente — **hecho para Astrobee** (comparación gráfica guardada); **pendiente aplicar el mismo nivel a BIRDS-TM y a la librería de Waskiri** |
| | | DR1.3 | Generar código C que reproduzca la solución de Python | VR1.3.1 | El "Estado final" impreso por el ejecutable C coincide con `x_opt_unscaled` de Python dentro de tolerancia numérica | Compilar y correr el C generado, comparar contra el resultado de Python — **ya hecho para Astrobee y BIRDS-TM** |
| SR2 | Intercambiar datos con el ADCS | DR2.1 | Establecer comunicación I2C entre COGU y ADCS | VR2.1.1 | `TBD` | `TBD` — pendiente de hardware/protocolo (ver §4.2, §5.1) |
| | | DR2.2 | Recibir condición inicial (y, a futuro, warm start) desde el ADCS en tiempo real | VR2.2.1 | `TBD` | `TBD` — desarrollo pendiente, ver nota de implementación en §4.3.4 |
| | | DR2.3 | Entregar la secuencia óptima de control al ADCS | VR2.3.1 | `TBD` | `TBD` |
| SR3 | Reportar el estado del solver | DR3.1 | Comunicar convergencia y métricas de iteración hacia el ADCS/OBC | VR3.1.1 | `TBD` | `TBD` — ver manejo de errores en §4.3.5 |

---

## 4. Identificación de interfaces

| Source | Destination | Tipo de interfaz | Notas/Observaciones |
|---|---|---|---|
| ADCS | COGU | Software (datos), vía **I2C** | El ADCS, activado momentáneamente por comando del OBC, envía a COGU la condición inicial actual (cuaternión, velocidad/momento angular) y, a futuro, un warm start generado por red neuronal |
| OBC | ADCS | Software (comando) | El OBC comanda al ADCS para que se active y transmita sus datos a COGU (ver §2.3) |
| COGU | ADCS | Software (datos), vía **I2C** | Entrega la secuencia óptima de control (dipolo magnético `m(t)`) |
| COGU | OBC | **Ninguna directa** | Confirmado por el equipo: COGU no tiene interfaz directa con el OBC. Toda la comunicación pasa por el ADCS, que actúa de puente hacia el OBC (vía su propio CPLD/Mission Boss) |
| COGU | EPS | Eléctrica | Línea de alimentación — voltaje y corriente `[A COMPLETAR]` según el hardware final |
| COGU | Magnetorquers | Ninguna directa | COGU **no** comanda los magnetorquers directamente: entrega la secuencia óptima al ADCS, que es quien la ejecuta sobre el actuador |

### 4.1 Interfaces mecánicas

`[A COMPLETAR]` — pendiente de la definición de hardware/PCB del subsistema:
- Dimensiones físicas de la(s) placa(s) de COGU.
- Orientación del subsistema respecto a los ejes del satélite (si aplica; COGU es
  cómputo, no tiene una orientación funcional propia salvo la que imponga el chasis).
- Interfaz mecánica con el Backplane Board (conector de 50 pines, posición del pin 1).
- Interfaz mecánica con la estructura del satélite.

### 4.2 Interfaces eléctricas

**Confirmado:** el protocolo de datos entre COGU y el ADCS es **I2C**. Niveles
lógicos, resistencias pull-up, velocidad de bus (100kHz/400kHz) y pines específicos
del conector siguen `[A COMPLETAR]`, pendientes de la definición de hardware/PCB:

| Pin no. | Nombre de la conexión | Nivel de voltaje [V] | Corriente [A] | Frecuencia/baud rate | Protocolo | Source | Destination | Observaciones |
|---|---|---|---|---|---|---|---|---|
| `[A COMPLETAR]` | SDA | `[A COMPLETAR]` | `[A COMPLETAR]` | `[A COMPLETAR]` | I2C | COGU | ADCS | Línea de datos |
| `[A COMPLETAR]` | SCL | `[A COMPLETAR]` | `[A COMPLETAR]` | `[A COMPLETAR]` | I2C | COGU | ADCS | Línea de reloj |

### 4.3 Interfaces de software

> **Cómo está organizada esta sección.** El cuerpo principal (4.3.1 a 4.3.6) es el
> **contrato de interfaz** propiamente dicho: trata a COGU como caja negra y documenta,
> en tablas, exactamente qué datos cruzan la frontera, en qué formato y con qué
> unidades — siguiendo la práctica estándar de un ICD (NASA SWE-112 / IDD; ver
> Referencias en 4.3.6). La explicación de **cómo funciona COGU por dentro** (el
> algoritmo, el lenguaje para declarar costo/restricciones) se deja en los
> **Apéndices A y B**, al final del documento, para quien necesite implementar un
> problema nuevo — no forma parte del contrato de interfaz en sí. El detalle de campo
> por campo (Interface Name, Connection Type, Data Protocol, etc.) del enlace físico
> con el ADCS está en §5.1, siguiendo el mismo formato que el resto del proyecto.

#### 4.3.1 Diagrama de contexto

Ver diagrama de bloques en §2.3. Resumen de la interfaz vista desde afuera: COGU
recibe una condición inicial (y, a futuro, un warm start) desde el ADCS por I2C, y
devuelve por el mismo enlace una secuencia óptima de control. El modelo del satélite,
el costo, las restricciones y el objetivo de la maniobra **no** se transmiten en cada
activación — ya están embebidos en el programa de COGU desde la etapa de generación
en tierra (§4.3.2).

#### 4.3.2 Diccionario de datos — generación del problema (etapa de tierra, previa al vuelo)

Antes de que COGU exista como programa embebido, un ingeniero define la maniobra
completa mediante la función `solve_trajectory(...)`
(`COGU_Library_Enhanced/sympy_cvx_converter/pipeline.py`). Esta tabla es el
**diccionario de datos** de esa definición — no es una interfaz en tiempo de vuelo,
pero determina el contenido fijo que después queda embebido en COGU.

**Datos obligatorios:**

| Nombre | Tipo / forma | Unidades | Descripción |
|---|---|---|---|
| `states` | lista de símbolos, tamaño `nx` | — | Variables de estado del satélite (para Waskiri: cuaternión + momento angular) |
| `controls` | lista de símbolos, tamaño `nu` | — | Variables de control (para Waskiri: 3 componentes del dipolo magnético) |
| `dynamics` | lista de expresiones, tamaño `nx` | — | `ẋ = f(x,u,t)` simbólico, una expresión por estado |
| `start` | numpy `(nx,1)` | físicas (rad, kg·m²·rad/s, etc.) | Condición inicial de la maniobra |
| `end` | numpy `(nx,1)` | físicas | Condición objetivo de la maniobra |
| `T` | entero | pasos | Número de pasos de discretización temporal |
| `tf` | float | s | Duración de la maniobra |
| `cost_terms` | lista de dicts | — | Función de costo (sin default — ver Apéndice B) |

**Datos opcionales:**

| Nombre | Tipo / forma | Default | Descripción |
|---|---|---|---|
| `dynamic_parameters_sym` / `_val` | símbolos / numpy | `None` | Parámetros físicos (inercia, momento residual, tabla `B(t)` del WMM) |
| `nonconvex_constraints` | lista de expresiones | `None` | Restricciones no convexas `g(x,u)≤0` (no usado en Waskiri hoy) |
| `state_bounds`, `control_bounds` | lista de dicts/tuplas | `None` | Límites físicos (ver Apéndice B) |
| `S_x`, `c_x`, `S_u`, `c_u` | numpy | identidad | Escalado numérico interno |
| `warm_start_x` `(nx,T+1)`, `warm_start_u` `(nu,T)` | numpy | interpolación lineal | Trayectoria inicial del algoritmo |
| `size_N` | entero | 20 | Sub-pasos de integración (RK4) por paso de discretización |
| `lamb` | float | 1000.0 | Peso de penalización del control virtual interno |
| `max_iter` | entero | 25 | Iteraciones máximas del algoritmo |
| `rho0,rho1,rho2` | float | 0, 0.1, 0.7 | Umbrales de aceptación de paso |
| `etta0,etta1,etta_init` | float | 1e-8, 10.0, 1.0 | Límites y valor inicial de la región de confianza |
| `beta_sh,beta_gr` | float | 2.0, 2.0 | Factores de contracción/expansión de la región de confianza |
| `e_tol,epsilon_stop` | float | 0.005, 1e-5 | Tolerancias de convergencia |
| `solver` | string | `'ECOS'` | Solver convexo interno |
| `generate_c` | bool | `False` | Si genera el programa en C |
| `c_output_dir` | string | `"c_output"` | Carpeta de salida del código C |

(Detalle conceptual de `cost_terms`, `state_bounds`, `control_bounds`, y del
algoritmo en general → **Apéndice A y B**.)

#### 4.3.3 Diccionario de datos — salidas de la etapa de generación

| Nombre | Tipo / forma | Descripción |
|---|---|---|
| `x_opt_unscaled` | numpy `(nx, T+1)` | Trayectoria óptima de estado, unidades físicas |
| `u_opt_unscaled` | numpy `(nu, T)` | Trayectoria óptima de control, unidades físicas |
| `x_opt`, `u_opt` | numpy | Igual, en variables escaladas internas |
| `converged` | bool | Si el algoritmo cumplió el criterio de convergencia |
| `iterations` | entero | Iteraciones realmente ejecutadas |
| `history` | lista de dicts | Un registro por iteración (`L_SCVx`, `J_SCVx`, `rho`, `etta`, `norm_diff`) |

Si `generate_c=True`, además se produce una carpeta con el programa en C
(`scvx_main.c`, `dynamics.c/.h`, `cpg_compat.h`, `CMakeLists.txt`, `solver/`) —
detalle de cada archivo en `MANUAL.md` §7.

#### 4.3.4 Diccionario de datos — interfaz en tiempo de vuelo (ADCS ↔ COGU, por I2C)

**Contexto operacional** (según arquitectura confirmada en §2.3): el OBC activa
momentáneamente al ADCS, que a su vez envía datos a COGU por I2C. De todo el problema
que resuelve COGU, **solo dos cosas cambian en cada activación**; todo lo demás queda
fijo desde la generación en tierra:

| Dato | Origen | ¿Varía en tiempo real? | Tipo / forma | Unidades |
|---|---|---|---|---|
| Cuaternión actual | ADCS | **Sí** | `(4,1)` | adimensional |
| Velocidad angular / momento angular actual | ADCS | **Sí** | `(3,1)` | rad/s ó kg·m²·rad/s |
| Warm start (trayectoria inicial `x`,`u`) | Red neuronal (futuro, ver notas de proyecto internas) | **Sí** (a futuro) | `(nx,T+1)` / `(nu,T)` | según variable |
| Modelo dinámico, costo, restricciones | Definidos en tierra (§4.3.2) | No | — | — |
| Objetivo de la maniobra (`end`) | Definido en tierra (§4.3.2) | No (fijo por tipo de maniobra) | `(nx,1)` | físicas |
| Parámetros físicos (inercia, `m_RMM`, tabla `B(t)` WMM) | Definidos en tierra (§4.3.2) | No | — | — |
| Secuencia óptima de control `m(t)` | COGU → ADCS | (salida) | `(3,T)` | A·m² |

**Estado de implementación (importante):** esta tabla describe el **diseño previsto**
de la interfaz en tiempo de vuelo. Hoy, el programa en C generado **todavía no lee
ninguno de estos datos desde el bus I2C real** — la condición inicial también queda
horneada como constante en tiempo de generación, igual que el resto. Implementar la
lectura real de estos datos desde la interfaz física (§4.2, §5.1) es desarrollo
pendiente.

**Mensajes / comandos** (formato exacto `TBD`, a definir con el equipo de OBC/ADCS
sobre el detalle del protocolo I2C — ver §5.1):

| Mensaje | Dirección | Formato | Estado |
|---|---|---|---|
| Activar ADCS y solicitar envío de datos a COGU | OBC → ADCS | `TBD` | Pendiente de definir |
| Estado actual (cuaternión, velocidad/momento angular) [+ warm start, a futuro] | ADCS → COGU (I2C) | `TBD` | Pendiente de definir |
| Secuencia óptima de control `m(t)` | COGU → ADCS (I2C) | `TBD` | Pendiente de definir |

#### 4.3.5 Manejo de errores / modos de falla

| Condición | Qué pasa hoy (etapa de tierra, Python) | Qué pasa hoy (programa en C, a bordo) | Qué falta definir |
|---|---|---|---|
| El algoritmo no converge dentro de `max_iter` | `result['converged'] = False`; igual devuelve la última trayectoria calculada | No se imprime el mensaje de convergencia; el programa igual imprime el "Estado final" alcanzado (el mejor disponible, no necesariamente válido) | Un código/bandera de error explícito que el ADCS/OBC pueda leer para decidir qué hacer (ej. reintentar, usar un warm start distinto) |
| El solve lanza una excepción (infactible, `NaN`, error del solver) | La excepción se propaga; **no se genera C** de una solución inválida (barrera de seguridad intencional) | No aplica — el C ya fue generado y validado en tierra antes de compilarse | — |
| El solver interno (ECOS embebido) no reporta óptimo en alguna iteración | No aplica (Python usa CVXPY, con su propio manejo) | **Hallazgo:** el programa en C no verifica el estado de salida del solver en cada iteración — asume éxito siempre | Verificar el código de estado del solver en cada iteración y decidir qué hacer si falla (abortar, devolver el último válido, etc.) |
| Datos de entrada fuera de rango (ej. cuaternión no normalizado, `NaN`) | No hay validación explícita | No aplica todavía (no hay entrada en tiempo real) | Validar la entrada una vez que exista la interfaz I2C real ADCS→COGU |

#### 4.3.6 Referencias

*A completar.* Categorías previstas: el/los paper(s) de SCVx / Successive
Convexification en los que se basa el algoritmo; la documentación del modelo WMM
usado para el campo magnético; y las herramientas de software empleadas (SymPy,
CVXPY, CVXPYgen, ECOS). Además, referencia interna a `MANUAL.md` (manual de usuario
de la librería) y a los Apéndices A/B de este documento.

#### 4.3.7 Diagrama de flujo de operación de COGU

Siguiendo el mismo formato de flowchart usado en los escenarios de operación del
ICD de OBC-CDH (§3.3.6 de ese documento):

```
                    ┌───────┐
                    │ START │
                    └───┬───┘
                        ▼
          ┌───────────────────────────┐
          │ Recibir activación del     │
          │ ADCS (I2C)                 │
          └─────────────┬──────────────┘
                         ▼
          ┌───────────────────────────┐
          │ Leer condición inicial:    │
          │ q, ω/h_B [+ warm start,    │
          │ a futuro]                  │
          └─────────────┬──────────────┘
                         ▼
          ┌───────────────────────────┐
          │ Iterar SCVx: discretizar + │◀─────────────┐
          │ resolver subproblema       │              │
          │ convexo                    │              │
          └─────────────┬──────────────┘              │
                         ▼                             │
                  ╱─────────────╲    no                │
                 ╱ ¿Convergió o   ╲───────────────────-┘
                 ╲ max_iter?      ╱
                  ╲──────┬───────╱
                         │ sí
                         ▼
          ┌───────────────────────────┐
          │ Entregar m(t) óptimo al    │
          │ ADCS (I2C)                 │
          └─────────────┬──────────────┘
                         ▼
                     ┌───────┐
                     │  END  │
                     └───────┘
```

---

## 5. Especificaciones de interfaz por subsistema

Sigue el mismo formato de campos fijos usado en el ICD de OBC-CDH (§4 de ese
documento: `Interface Name`, `Connection Type`, `Electrical Characteristics`,
`Connector Type`, `Pinout`, `Data Protocol`, `Message Format`, `Timeout`,
`Retry Logic`), para que todos los ICD de Waskiri-Sat sean comparables entre sí.

### 5.1 COGU <-> ADCS

- **Interface Name:** COGU_ADCS
- **Connection Type:** Hardware I2C (**confirmado**)
- **Electrical Characteristics:** `[A COMPLETAR]` (niveles lógicos, resistencias
  pull-up, corriente máxima — pendiente de hardware/PCB)
- **Connector Type:** `[A COMPLETAR]` (vía backplane, pendiente de definición física)
- **Pinout:** `[A COMPLETAR]` — referirse al mapeo de pines del backplane una vez
  definido (ver §4.2)
- **Data Protocol:** I2C — velocidad de bus (100kHz/400kHz) y rol maestro/esclavo
  `[A COMPLETAR]`
- **Message Format:** `TBD` — el contenido de datos a intercambiar ya está definido
  (§4.3.4: condición inicial, warm start futuro, secuencia óptima de control); el
  formato exacto de bytes queda pendiente
- **Timeout:** `TBD`
- **Retry Logic:** `TBD`

### 5.2 COGU <-> EPS

- **Interface Name:** COGU_EPS
- **Connection Type:** Eléctrica (alimentación, sin protocolo de datos)
- **Electrical Characteristics:** `[A COMPLETAR]` (voltaje, corriente)
- **Connector Type:** `[A COMPLETAR]`
- **Pinout:** `[A COMPLETAR]`

> COGU no tiene interfaz directa con el OBC ni con los magnetorquers (ver tabla de
> identificación de interfaces en §4).

---

## Apéndice A — Cómo funciona COGU por dentro

*(Nota: los demás ICD de Waskiri-Sat, como el de OBC-CDH, no incluyen un apéndice de
este tipo porque documentan protocolos de comunicación estándar. COGU sí lo amerita
por tratarse de un subsistema de optimización con un algoritmo propio que alguien
deberá poder extender — ver Apéndice B. Este material es de referencia para quien
vaya a extender la librería o depurar su comportamiento; no es parte del contrato de
interfaz de §4.3/§5.)*

COGU resuelve, para una maniobra de actitud dada, un problema de control óptimo:

```
  min_{x(t), u(t)}   J(x, u)                     (costo: energía, seguimiento, etc.)
  s.a.               ẋ(t) = f(x(t), u(t), t)     (dinámica no lineal del satélite)
                      x(t) ∈ X,  u(t) ∈ U          (límites físicos)
                      x(0) = x0,  x(tf) = x_goal    (condición inicial y objetivo)
```

Es un problema **no convexo** (la dinámica de actitud lo es). COGU lo resuelve con
**SCVx (Sequential Convex Programming)**: en cada iteración, linealiza la dinámica
alrededor de la mejor trayectoria conocida, resuelve el subproblema **convexo**
resultante con un solver tipo ECOS, y actualiza la trayectoria de referencia. Se
repite hasta que la solución deja de cambiar significativamente o se alcanza
`max_iter`.

Dos mecanismos hacen que el proceso sea numéricamente bien comportado:

- **Control virtual (`vc`):** variable de holgura que absorbe temporalmente el error
  de la linealización, penalizada fuertemente en el costo para que tienda a cero al
  converger. Sin esto, el subproblema linealizado podría no tener solución factible
  en las primeras iteraciones.
- **Región de confianza (`η`):** la nueva trayectoria no puede alejarse más de `η` de
  la anterior (norma infinito). `η` se contrae o expande según un criterio `ρ` que
  mide qué tan bien predijo el modelo linealizado la mejora real del costo:

  | `ρ` | Interpretación | Acción sobre `η` |
  |---|---|---|
  | `ρ < rho0` | El paso fue malo, se rechaza | Se contrae (`η/beta_sh`), no se actualiza la trayectoria |
  | `rho0 ≤ ρ < rho1` | Paso aceptado, pero la predicción fue mala | Se contrae, se actualiza la trayectoria |
  | `rho1 ≤ ρ < rho2` | Paso aceptado, predicción razonable | `η` sin cambios |
  | `ρ ≥ rho2` | Paso muy bueno | Se expande (`η·beta_gr`) |

Esta arquitectura (SCVx + `vc` + región de confianza adaptativa) es lo que implementa
la función `solve_trajectory(...)` en tierra, y es también la lógica que queda
"congelada" dentro del programa en C que corre a bordo.

**Cómo se valida esto en la práctica (ver también §3, VR1.1.1/VR1.2.1):** para
Astrobee, el caso más completo, se comparó no solo el valor del costo por iteración
contra el notebook de referencia, sino también gráficas de la trayectoria óptima
(posición, velocidad, actitud, control) para verificar visualmente que el resultado
tiene sentido físico. Para BIRDS-TM se comparó contra el notebook de referencia
(valores `L_SCVx`/`J_SCVx` por iteración) sin gráficas adicionales. Se espera aplicar
el mismo criterio — idealmente al nivel de Astrobee — a la librería customizada de
Waskiri (`COGU_Library_Waskiri`) una vez que exista.

---

## Apéndice B — DSL de costo y restricciones (para implementar un problema nuevo)

*(Es el mecanismo concreto para traducir la función objetivo y las restricciones de
un modelo matemático — ej. de un notebook de diseño — al formato que entiende la
librería.)*

### B.1 Costo (`cost_terms`)

`cost_terms` es una **lista de términos**; el costo total es la suma de todos. Cada
término es un diccionario (definición única en `sympy_cvx_converter/cost_dsl.py`,
compartida por las tres materializaciones del costo: CVXPY en Python, evaluación
numérica, y C generado — para que nunca se desincronicen):

```python
{
  'kind':    'sumsq' | 'norm1' | 'norm2',   # tipo de penalización
  'var':     'u' | 'x',                      # penaliza el control o el estado
  'slice':   slice(a, b),                    # qué componentes (default: todos)
  'coeff':   'sqrt_tau' | 'tau' | 'tau_lamb' | float,  # escalar (default: 'sqrt_tau')
  'weight':  float,                          # peso relativo (default: 1.0)
  'offset':  numpy (n,1) | None,             # valor objetivo a restar (default: None)
  'k_range': 'T' | 'T+1',                    # rango de pasos (default: 'T' para u, 'T+1' para x)
}
```

Con `v = coeff · (var[slice, k] − offset)`, cada término se traduce a:

| `kind` | Fórmula | Uso típico |
|---|---|---|
| `'sumsq'` | `weight · Σₖ ‖v‖²` | energía / mínimo esfuerzo cuadrático |
| `'norm1'` | `weight · Σₖ ‖v‖₁` | favorece soluciones dispersas ("on-off"), aunque no lo garantiza siempre |
| `'norm2'` | `weight · Σₖ ‖v‖₂` | magnitud, sin elevar al cuadrado |

**Ejemplos:**

```python
# Energía del dipolo magnético (caso típico de Waskiri)
cost_terms = [
    {'kind': 'sumsq', 'var': 'u', 'coeff': 'sqrt_tau'},
]

# + seguimiento de una actitud de referencia (ej. fine guidance / sun-pointing)
cost_terms = [
    {'kind': 'sumsq', 'var': 'u', 'coeff': 'sqrt_tau'},
    {'kind': 'sumsq', 'var': 'x', 'slice': slice(0,4), 'coeff': 'sqrt_tau',
     'weight': 100, 'offset': q_goal_escalado, 'k_range': 'T+1'},
]
```

> No declarar las penalizaciones del control virtual `vc` (ni de `vi`, si hay
> restricciones no convexas) — son mecanismo interno del algoritmo y la librería las
> agrega automáticamente. `cost_terms` describe solo el objetivo de ingeniería.

### B.2 Restricciones (`state_bounds` / `control_bounds`)

Cada restricción limita una porción (`slice`) del estado o del control:

```python
# Formato simple — norma con un tope
(slice(0,3), 'norminf', m_max)     # ‖m‖∞ ≤ m_max → así se limita el dipolo en Waskiri

# Formato explícito — admite caja y transformación propia
{'kind': 'norm', 'slice': slice(4,7), 'norm_type': 'norm2', 'limit': omega_max,
 'S': inv_J, 'c': ceros}
# ‖inv(J)·h_B‖₂ ≤ ω_max — límite de velocidad angular en Waskiri, transformando
# el momento angular (el estado real) a velocidad angular con inv(J)

{'kind': 'box', 'slice': slice(0,3), 'lower': 0.0, 'upper': f_max}
# 0 ≤ u[0:3] ≤ f_max (límite de caja, ej. actuadores unidireccionales)
```

`norm_type` acepta `'norm1'`, `'norm2'` o `'norminf'`. Los campos opcionales `S`/`c`
aplican una transformación lineal propia antes de medir la norma, en vez del escalado
global — es el mecanismo que usa Waskiri para expresar el límite de velocidad angular
en términos del momento angular.

### B.3 Cómo traducir un modelo desde un notebook/paper a este formato

| Expresión matemática / en notebook | `cost_term` o `bound` equivalente |
|---|---|
| `Σ‖√τ·u_k‖²` (energía de control) | `{'kind':'sumsq','var':'u','coeff':'sqrt_tau'}` |
| `w·Σ‖√τ·(x_k − x_goal)‖²` (seguimiento) | `{'kind':'sumsq','var':'x','coeff':'sqrt_tau','weight':w,'offset':x_goal,'k_range':'T+1'}` |
| `‖u‖∞ ≤ u_max` | `(slice(None), 'norminf', u_max)` |
| `0 ≤ u ≤ u_max` | `{'kind':'box','slice':slice(None),'lower':0.0,'upper':u_max}` |
| `‖S·x + c‖₂ ≤ x_max` (con transformación propia) | `{'kind':'norm','norm_type':'norm2','limit':x_max,'S':S,'c':c}` |

### B.4 Ejemplo completo de referencia

`COGU_Library_Enhanced/examples/birds_tm.py` es, hoy, el ejemplo más parecido en
estructura a lo que necesitará Waskiri (actitud pura, sin traslación, con límites
mixtos de norma y de caja) — usarlo como plantilla de partida es más rápido que
escribir desde cero. El modelo específico de Waskiri (magnetorquers + WMM) requiere
además una librería customizada (`COGU_Library_Waskiri`, pendiente) porque usa un
esquema de discretización distinto (First Order Hold) al que soporta
`COGU_Library_Enhanced` (Zero Order Hold) — ver notas de proyecto internas.

---

## Notas para completar este borrador

1. Obtener del equipo de hardware/PCB: dimensiones, conector, voltajes, corrientes,
   y los parámetros exactos del bus I2C (velocidad, pull-ups) — secciones 2.4, 4.1,
   4.2, 5.1.
2. Definir con el equipo de OBC/ADCS el formato exacto de los mensajes marcados `TBD`
   en §4.3.4 y §5.1, una vez definidos los detalles del protocolo I2C.
3. Completar la Tabla de asignación de requisitos (§3) — varias filas de SR2/SR3
   quedan `TBD` hasta que exista la interfaz I2C real.
4. Completar §4.3.6 con las citas exactas (papers de SCVx, WMM, herramientas usadas).
5. Aplicar a `COGU_Library_Waskiri` (cuando exista) el mismo nivel de validación que
   Astrobee (comparación numérica + gráficas) — ver Apéndice A.
6. Una vez completado y revisado en español, traducir a inglés siguiendo el formato
   del template oficial y el ICD de OBC-CDH como referencia de estilo.
