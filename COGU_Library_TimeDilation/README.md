# scvx_toolbox

SCvx trajectory optimization from a SymPy model, with selectable time
discretization. Successor to `sympy_cvx_converter`, extended to free final
time (time dilation) and first-order hold.

Write the dynamics once, pick a time mode:


| `time_mode`      | control hold | final time |
| ---------------- | ------------ | ---------- |
| `'ZOH'`          | zero-order   | fixed      |
| `'FOH'`          | first-order  | fixed      |
| `'ZOH+dilation'` | zero-order   | **free**   |
| `'FOH+dilation'` | first-order  | **free**   |

Aliases are accepted, so `'zoh_time_dilation'`, `'ZOH + dilation'` and
`'zoh-td'` all resolve to the same mode.

---

## Install

```
pip install numpy scipy sympy cvxpy ecos
# only for C generation (also needs cmake and a C compiler):
pip install cvxpygen
```

Put `scvx_toolbox/` on your path and import it. (`scvx_toolbox/__init__.py`
was missing from an earlier drop — without it `from scvx_toolbox import
solve_trajectory` fails, since the directory then resolves as a namespace
package with no exports.)

---

## Quick start

```python
from scvx_toolbox import solve_trajectory

result = solve_trajectory(
    states=states, controls=controls, dynamics=dynamics,   # SymPy
    start=start_pos, end=end_pos, T=30,

    time_mode='ZOH+dilation',
    tf_init=480.0, tf_min=60.0, tf_max=900.0, alpha_p=0.05,

    cost_terms=[{'kind': 'norm1', 'var': 'u', 'coeff': 'tau', 'weight': 0.5}],
    state_bounds=..., control_bounds=...,
)

result['tf']              # optimized final time
result['x_opt_unscaled']  # (nx, T+1)
result['u_opt_unscaled']  # (nu, T)   -- (nu, T+1) under FOH
result['status']          # 'converged' | 'stalled' | 'max_iter' | 'subproblem_...'
```

---

## Time dilation, in one paragraph

The problem is restated on the normalized interval `s ∈ [0,1]` using the
dilation parameter `p ≡ tf`:

```
t = p·s        dx/ds = (dx/dt)(dt/ds) = p·f(x,u) =: f(x,u,p)
```

`p` becomes a decision variable, so `min α_p·tf` is finally something the
solver can act on. Linearization needs one extra Jacobian,

```
F = ∂f/∂p        y = f - A·x̄ - B·ū - F·p̄
G = ∂g/∂p        z = g - C·x̄ - D·ū - G·p̄
```

scaled as `F = S_x⁻¹ F̂ S_p`, `G = Ĝ S_p`, with `p̂ = S_p p + c_p` chosen so
scaled `p` lives in `[0,1]` across `[tf_min, tf_max]`. The discrete update
becomes

```
ZOH:  x_{k+1} = A_k x_k + Bm_k u_k                + F_k p + y_k + vc_k
FOH:  x_{k+1} = A_k x_k + Bm_k u_k + Bp_k u_{k+1} + F_k p + y_k + vc_k
```

obtained by RK4 on `Φ' = AΦ`, `Gm' = A·Gm + λ⁻B`, `Gp' = A·Gp + λ⁺B`,
`H' = AH + F`, `z' = Az + y`, where `λ⁻ = (t_{k+1}-t)/τ`, `λ⁺ = (t-t_k)/τ`
under FOH and `λ⁻ ≡ 1`, `λ⁺ ≡ 0` under ZOH. The trust region picks up
`‖p - p̄‖∞`, and `tf_min ≤ S_p p + c_p ≤ tf_max` is added.

**You write `dx/dt = f(x,u)` in absolute time.** The library forms `p·f`
itself, which is what keeps one model definition valid in all four modes. (If
you prefer to pre-multiply it yourself, pass `dilation_symbol=` and the
library will detect it and not double up.)

The running cost stays in normalized time and deliberately omits the dilation
factor — `p·‖u‖₁` would be bilinear, so neither convex nor DPP. The final
time enters only through the convex term `α_p·tf`.

---

## `α_p` is the knob

- `α_p = 0` — `tf` is free to be whatever the rest of the cost prefers.
- `α_p` large — drives `tf` towards `tf_min` (minimum-time).

`tf_min`/`tf_max` default to `0.1·tf_init` and `3·tf_init` if you omit them.

---

## Cost DSL

```python
{'kind':    'sumsq' | 'norm1' | 'norm2',
 'var':     'u' | 'x',
 'slice':   slice(...),                        # default slice(None)
 'coeff':   'sqrt_tau' | 'tau' | 'tau_lamb' | float,
 'weight':  float,
 'offset':  (n,1) array,                       # tracking: coeff*(v - offset)
 'k_range': 'T' | 'T+1' | 'U'}
```

`'U'` means "every control node" — `T` under ZOH, `T+1` under FOH — and is the
default for `var='u'`. That is what makes a cost definition portable across
holds.

`'tau'` resolves to the step of the formulation being solved: `tf/T` without
dilation, the normalized `1/T` with it.

The DSL has three backends kept in sync through `cost_dsl.resolve_term`:
CVXPY (`problem.py`), numpy (`solver.py`), and C (`cgen.py`).

---

## Bounds DSL

```python
(slice(3,6), 'norm2', 5.0)                                    # shorthand
{'kind': 'norm', 'slice': ..., 'norm_type': ..., 'limit': ...,
 'S': optional, 'c': optional}                                # own transform
{'kind': 'box', 'slice': ..., 'lower': ..., 'upper': ...}
{'kind': 'box_affine', 'c_affine': (1,n), 'lower': ..., 'upper': ...}
```

`'S'`/`'c'` replace the global scaling sub-block, which is how you bound a
derived quantity — e.g. bounding angular *velocity* `inv(J)·h` when the state
is angular *momentum* `h`.

---

## Migrating from `sympy_cvx_converter`

Mostly mechanical:


| old                                                        | new                                                       |
| ---------------------------------------------------------- | --------------------------------------------------------- |
| `from sympy_cvx_converter import solve_trajectory`         | `from scvx_toolbox import solve_trajectory`               |
| `from sympy_cvx_converter.utils import slerp, angular_vel` | `from scvx_toolbox.utils import slerp, angular_velocity`  |
| `etta0`, `etta1`, `etta_init`                              | `eta0`, `eta1`, `eta_init` (old spellings still accepted) |
| —                                                         | add`time_mode=` (defaults to `'ZOH'`)                     |

Two behavioural differences to know about, both deliberate:

1. `angular_velocity` fixes a quaternion-convention bug. The old helper fed
   scalar-first `(q0,q1,q2,q3)` straight into `scipy`'s scalar-last
   `from_quat`, so the warm-start angular velocity was wrong. Pass
   `scalar_first=False` for the old behaviour.
2. `pin_aware_trust_region=True` by default (see below). Set it `False` to
   reproduce the legacy formulation bit for bit.

---

## Three bugs found in the original formulation

These are present in both the notebook and `sympy_cvx_converter`.

### 1. The trust region charged for boundary-pinned states

`x[:,0]` and `x[:,T]` are fixed by equality constraints, so
`‖x₀ − ox₀‖∞` is a **constant**. Charging it against `eta` does nothing
useful and is actively harmful: whenever the warm start disagrees with
`start`/`end` — which is the norm, since any warm start whose angular
velocity comes from finite differences starts at zero — that constant becomes
a hard floor on `eta`, and the subproblem goes **infeasible** the moment the
trust region shrinks below it.

In the BIRDS VAT problem the floor is `0.0667`; the solve died at `eta = 0.0625`. Fixed by `pin_aware_trust_region=True` (default), which omits the
state term at the two pinned nodes.

### 2. Convergence could be declared on a rejected step

The test was `ΔL < e_tol·|J_ref|`. When a step is rejected, `ΔL` is a large
**negative** number and the test passes trivially — declaring success exactly
when the linearization is least trustworthy. Observed: a dilation run
"converging" at iteration 4 with `tf` never having moved off its initial
guess. Now uses `|ΔL|`.

### 3. `J` and `L` penalize constraint violation in different units

`J` measures the defect in **unscaled** (physical) units; `L` penalizes the
slack `vc` in **scaled** units. With an aggressive `S_x` these are far apart —
in the BIRDS problem `S_x` compresses the momentum states by `~1.7e-3`, so
`J`'s defect penalty is much weaker than `L`'s. The result is `rho < 0` on
steps that clearly improve `J`, every step rejected, and `tf` frozen.

Default is `defect_penalty_space='unscaled'` so existing results reproduce.
**Set `'scaled'` if `rho` is negative on steps that improve `J`.** On the
BIRDS VAT dilation problem this is the difference between a useless answer
and a good one:


|                                | `'unscaled'`           | `'scaled'`   |
| ------------------------------ | ---------------------- | ------------ |
| optimal`tf`                    | 480.00 s (never moved) | **420.24 s** |
| open-loop`‖q_f − q_target‖` | 1.15                   | **2.4e-3**   |

Also added: stall detection (`max_stall`), and an infeasible subproblem now
ends the loop with a diagnostic instead of throwing away the run.

---

## Validation

Run from the package root.

**`tests/test_astrobee_vs_enhanced.py`** — fixed-time ZOH parity with
`sympy_cvx_converter` (COGU_Library_Enhanced), on the Astrobee example that
library shipped: T=30, tf=200, three keep-out spheres. Iteration 1 gives
`L=4919.371870 J=14199.891376` from both, agreeing to 1e-12; the first twelve
iterations hold to 1e-6. Enhanced's generated C binary gives the same numbers,
so one reference covers Python and C.

Neither library converges on this problem — both reach `max_iter` with `rho`
negative. That is a property of the example, not a regression:
`delta_L = J_ref - L` plateaus near 28 because the ZOH defect at `tau = 6.67 s`
does not shrink, while the test compares it against `e_tol*|J_ref| = 3.8`. The
*step* does shrink (`||dx||` is 1e-5 by iteration 25) and the terminal state is
hit to 1e-9 m, so the trajectory is converged even though the flag says
otherwise. See `examples/astrobee_compare_settings.py`.

**`tests/test_astrobee_zoh.py`** — the same model against the reference
notebook, which runs `T = 51-1 = 50`. An earlier Astrobee example in
`sympy_cvx_converter` compared a T=30 run against the notebook's T=50 numbers,
which is why the library used to look wrong when it was not. At T=50 the match
is exact:

```
  it |        L (lib)      L (nb) |        J (lib)      J (nb)
   1 |     3267.22617  3267.22617 |     8899.07051  8899.07040
   5 |      723.57423   723.57440 |      733.68120   733.68164
   8 |      723.53931   723.53931 |      729.46697   729.46696
```

Known notebook bug that does not affect this comparison: cell 36 adds the `vi`
penalty to `L` twice — a guarded block followed by an unguarded copy. The
obstacle constraints are satisfied at the solution, so `vi = 0` and the
duplicate vanishes.

**`tests/test_astrobee_c.py`** — generate, compile with cmake, run, and compare
the C loop against the Python solve. Astrobee is the case that exercises what
BIRDS-TM cannot: `ng > 0` (the C/D/z Jacobians, the `vi` slack, the `g()` term
of `J`) and `nu > 1` (indexed `cpg_update_ou_k`).

**`tests/test_dilation_vs_notebook.py`** — reproduces
`BIRDS_TM_Th_cvxpy_cvxpygen_l2_time_dilation.ipynb` exactly, with
`pin_aware_trust_region=False`:

```
  it |     L (mine)   L(nb)+0.5 |     J (mine)      J (nb) |    tf (mine)    tf (nb)
   1 |      35.4272     35.4272 |     302.1056    302.1056 |     298.6339    298.6339
   2 |      49.8521     49.8521 |     251.4258    251.4258 |     585.6905    585.6905
   3 |      34.5228     34.5228 |     155.0715    155.0715 |     410.7513    410.7513
   4 |      35.8089     35.8089 |      43.7126     43.7126 |     445.8738    445.8738
final tf = 443.9383 s   (notebook: ~443.9 s)
```

The uniform `+0.5` on `L` is intentional: the notebook's `L` drops the
constant `α_p·c_p = 0.5` that its `J` includes, which biases `rho`. This
library adds it back so `L` and `J` are comparable.

_Not in this drop:_ `tests/test_all_modes.py`,
`examples/birds_vat_compare_modes.py` and `tests/test_c_generation.py`
are described below but were not included with the sources.

---

## Examples

**`examples/astrobee_simple.py`** — the Astrobee example from
COGU_Library_Enhanced, ported. `astrobee_model.py` holds the model and data
(shared with the tests), `astrobee_compare_settings.py` measures what the
formulation corrections do on it, and `plot_astrobee.py` draws the six
state/control figures plus the 3D view, with `--reference` to overlay another
run:

```
python examples/plot_astrobee.py --legacy --reference ref_enhanced.npz --save out/
```

**`examples/birds_vat_only_dilation.py`** — `birds_vat_only.py` with ZOH +
time dilation. Includes open-loop RK4 validation against the original
nonlinear dynamics on the *optimized* grid, plus plots.

```
Optimal tf:       420.2440 s     (initial guess 480.0, bounds [60, 900])
Converged:        True  (18 iterations)
Final quaternion OL: [-0.00047974  0.99999709 -0.00228366 -0.00062022]
Target quaternion:   [0 1 0 0]
||q_final - q_target||: 2.41e-03
||omega_final||:        2.58e-05 rad/s
Max ||omega|| OL:  0.01328935 rad/s   (limit 0.02617994)
```

**`examples/birds_vat_compare_modes.py`** — the same model in all four modes:

```
mode            it      status    tf [s]         J   OL q err  OL |w| err    max|w|
ZOH             18   converged    480.00   701.132   2.73e-03    2.20e-05  0.008242
FOH             19   converged    480.00   560.332   2.95e-04    3.17e-06  0.008314
ZOH+dilation    18   converged    420.24    22.880   1.49e-03    1.36e-05  0.013288
FOH+dilation    12   converged    420.36    22.346   4.69e-04    7.89e-06  0.013107
```

FOH buys about an order of magnitude in open-loop accuracy. Both dilation
modes independently find `tf ≈ 420 s`, a 12% reduction on the fixed 480 s.

`J` is **not** comparable between the fixed and dilation rows: `tau` is
`tf/T = 16` in one and `1/T = 0.033` in the other.

---

## C generation (`generate_c=True`)

Fully working, verified end to end: generated, compiled with `cmake`/`gcc`,
executed, and checked against the Python result.

Emits into `c_output_dir`:


| file                | contents                                                    |
| ------------------- | ----------------------------------------------------------- |
| `scvx_main.c`       | the SCvx loop, all four modes behind`#if` guards            |
| `dynamics.c` / `.h` | model and every Jacobian, incl.`F = ∂f/∂p`, `G = ∂g/∂p` |
| `solver/`           | convex solver, via CVXPYgen in a clean subprocess           |
| `cpg_compat.h`      | stacked update API over CVXPYgen's per-step setters         |
| `embedded_data.h`   | problem data and SCvx tuning as C constants                 |
| `CMakeLists.txt`    | builds everything                                           |

```cd
cd c_output_dir; cmake -B build -DCMAKE_C_FLAGS="-O3"; cmake --build build --config Release; .\build\Release\scvx_app.exe
```

### C reproduces Python

`examples/birds_vat_only_dilation.py` with `generate_c=True`, then running the
compiled binary:

```
Iter   1 | rho=   None | J=    128.140040 | tf= 342.9272s
Iter  18 | rho= 0.0102 | J=     22.880151 | tf= 420.2448s
SCvx converged.
optimal tf: 420.243953 s
```

Same 18 iterations, same `rho` sequence, same `J` to six digits, same
`tf = 420.244 s` as the Python solve. Discretization runs in 1.7 ms and the
solve in ~15 ms per iteration.

`tests/test_c_generation.py` builds and runs all four modes on a problem that
also has a nonconvex constraint, so every preprocessor branch is exercised:

```
mode            files  cmake  build   run   tf (py)    tf (C)
ZOH              True   True   True  True   10.0000   10.0000
FOH              True   True   True  True   10.0000   10.0000
ZOH+dilation     True   True   True  True    5.5293    5.5293
FOH+dilation     True   True   True  True    5.2397    5.2397
```

### Two things the template has to get right

**CVXPYgen prunes unused parameters.** With `pin_aware_trust_region=True`,
`ox_0` and `ox_T` appear in no constraint, so CVXPYgen emits no
`cpg_update_ox_0`. `cpg_compat.h` turns those switch cases into no-ops; a
naive generator produces a header that will not link.

**Size-1 parameters get a different signature.** CVXPYgen emits
`cpg_update_NAME(cpg_float val)` — no index — whenever a parameter's total
size is 1. So `ou_k` is indexed when `nu > 1` but scalar when `nu == 1`, and
`op` (always 1×1) is always scalar. The compat generator switches on the
stride to match.

### The SCvx tuning travels with the code

`size_N`, `max_iter`, `rho*`, `eta*`, `beta*`, `e_tol`, `epsilon_stop`,
`max_stall` and `defect_penalty_space` are emitted as `#define`s, so the
generated program runs the same algorithm as the Python call that produced
it. No hand-editing of `main()`.

### If you want your own loop template

Drop it at `scvx_toolbox/templates/scvx_loop_template.c`. Placeholders
(`cgen.TEMPLATE_PLACEHOLDERS`):

```
{NX} {NU} {NG} {NP} {T} {U_NODES} {DILATION} {FOH} {TIME_MODE}
{USER_COST} {EMBEDDED_DATA}
```

`{DILATION}` and `{FOH}` are `0`/`1` and `{U_NODES}` is `T` or `T+1`, so one
template serves all four modes. `{NP}` is the number of *dynamic parameters*
(the size of `dp[]`), not the number of dilation parameters.

### Notes

CVXPYgen runs in a subprocess because DPP canonicalization needs a few hundred
MB of *contiguous* memory; in the parent, SymPy and the SCvx solve have
already fragmented the heap. The worker also stubs out `pdaqp`, whose import
chain pulls in `juliacall` and tries to download a Julia toolchain — unused
here, and fatal on a machine without network access.

## Module map

```
scvx_toolbox/
  modes.py           mode names, aliases, control-node count
  symbolic.py        SymPy -> A, B, F, y, C, D, G, z; numpy + C emitters
  discretization.py  RK4 of Phi/Gm/Gp/H/z; nonlinear flow map
  problem.py         CVXPY subproblem (per-step params, DPP-sparse)
  solver.py          SCvx loop, trust region, rho logic
  pipeline.py        solve_trajectory()
  cost_dsl.py        cost term resolution, shared by all three backends
  scaling.py         x, u, p scaling
  quaternion.py      slerp, angular_velocity  (re-exported by utils.py)
  cgen.py            C generation
  _cpg_worker.py     CVXPYgen subprocess
  templates/
    scvx_loop_template.c   SCvx loop in C, all four modes
```

---

## Tuning notes

- **`rho < 0` on steps that improve `J`** → `defect_penalty_space='scaled'`.
- **Subproblem goes infeasible as `eta` shrinks** → something constant is in
  the trust region. Keep `pin_aware_trust_region=True` and check the warm
  start against `start`/`end`.
- **`L` and `J` plateau with a gap `e_tol` can't clear** → that gap is
  discretization error, not a solver failure. Stop on step size
  (`epsilon_stop`) or raise `T`/`size_N`. The BIRDS example uses
  `epsilon_stop=1e-3` for this reason.
- **`tf` jumps wildly early on** → `eta_init=1` with `p` scaled to `[0,1]`
  lets `tf` cross its whole range in one step. Lower `eta_init`, or narrow
  `[tf_min, tf_max]`.
