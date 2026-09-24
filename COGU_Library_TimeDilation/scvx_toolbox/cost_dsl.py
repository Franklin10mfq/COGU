"""
Cost DSL -- the single source of truth for the semantics of `cost_terms`.

The user cost is materialized by three backends:

    CVXPY  (problem.py: build_convex_cost)   -- builds the convex subproblem
    numpy  (solver.py:  eval_user_cost)      -- evaluates J in the SCvx loop
    C      (cgen.py:    generate_c_user_cost)-- emits the J_cost body in C

All three resolve fields and defaults through `resolve_term`, so the backends
cannot drift apart.

Format of a cost term
---------------------
    {'kind':    'sumsq' | 'norm1' | 'norm2',
     'var':     'u' | 'x',
     'slice':   slice | None,                        # default: slice(None)
     'coeff':   'sqrt_tau' | 'tau' | 'tau_lamb' | float,   # default 'sqrt_tau'
     'weight':  float,                               # default: 1.0
     'offset':  numpy (n,1) | None,                  # default: None
     'k_range': 'T' | 'T+1' | 'U'}                   # default: 'U' (u), 'T+1' (x)

`k_range`
    'T'    -> nodes 0 .. T-1
    'T+1'  -> nodes 0 .. T
    'U'    -> every control node: T columns under ZOH, T+1 under FOH.
              This is the default for `var='u'` and is what makes a cost
              definition portable across holds.

Note on `coeff` under time dilation
-----------------------------------
`tau` is the discretization step of the formulation actually being solved:
tf/T without dilation, and the normalized 1/T with dilation.  The running
cost is deliberately evaluated in normalized time and the dilation factor is
omitted from it, because p*||u||_1 would be bilinear (neither convex nor DPP).
The final time enters the objective through the separate convex term
alpha_p * tf.
"""

from collections import namedtuple

ResolvedTerm = namedtuple('ResolvedTerm',
                          'var slc coeff weight offset kind k_range')

VALID_KINDS = ('sumsq', 'norm1', 'norm2')
VALID_RANGES = ('T', 'T+1', 'U')


def resolve_term(term):
    """Resolve the fields and defaults of one cost term."""
    if 'var' not in term:
        raise ValueError(f"cost term is missing 'var': {term!r}")
    if 'kind' not in term:
        raise ValueError(f"cost term is missing 'kind': {term!r}")

    var = term['var']
    if var not in ('x', 'u'):
        raise ValueError(f"cost term 'var' must be 'x' or 'u', got {var!r}")

    slc = term.get('slice', slice(None))
    if slc is None:
        slc = slice(None)

    coeff = term.get('coeff', 'sqrt_tau')
    weight = float(term.get('weight', 1.0))
    offset = term.get('offset')

    kind = term['kind']
    if kind not in VALID_KINDS:
        raise ValueError(f"unknown cost term kind: {kind!r}")

    k_range = term.get('k_range') or ('T+1' if var == 'x' else 'U')
    if k_range not in VALID_RANGES:
        raise ValueError(f"unknown cost term k_range: {k_range!r}")

    return ResolvedTerm(var, slc, coeff, weight, offset, kind, k_range)


def node_count(k_range, T, u_nodes):
    """Number of nodes a term is summed over."""
    if k_range == 'T':
        return T
    if k_range == 'T+1':
        return T + 1
    return u_nodes


def coefficient_map(sqrt_tau, tau, tau_lamb):
    """Numeric values of the symbolic coefficient names."""
    return {'sqrt_tau': float(sqrt_tau),
            'tau': float(tau),
            'tau_lamb': float(tau_lamb)}


def resolve_coefficient(rt, coeff_map):
    """Numeric coefficient of a resolved term."""
    return coeff_map[rt.coeff] if isinstance(rt.coeff, str) else float(rt.coeff)
