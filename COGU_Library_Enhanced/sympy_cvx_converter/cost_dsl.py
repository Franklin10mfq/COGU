"""
DSL de costo — única fuente de verdad de la semántica de cost_terms.

El costo del usuario se materializa en tres backends:
    - CVXPY  (problem.py: _build_cost_from_terms)   — construye el problema
    - numpy  (solver.py:  eval_user_cost)           — evalua J en el loop SCVx
    - C      (codegen.py: generate_c_user_cost)     — genera J_cost en C

Los tres comparten la resolución de campos y defaults de cada término vía
`resolve_term`. Así no se duplican los defaults ni la inferencia de k_range
(evita que los backends se desincronicen).

Formato de un cost_term:
    {'kind':   'sumsq'|'norm1'|'norm2',
     'var':    'u'|'x',
     'slice':  slice | None,                 # default: slice(None)
     'coeff':  'sqrt_tau'|'tau'|'tau_lamb'|float,   # default: 'sqrt_tau'
     'weight': float,                          # default: 1.0
     'offset': numpy (n,1) | None,             # default: None (resta x_goal)
     'k_range':'T'|'T+1'}                       # default: 'T' (u), 'T+1' (x)
"""

from collections import namedtuple

ResolvedTerm = namedtuple('ResolvedTerm', 'var slc coeff weight offset kind k_range')


def resolve_term(term):
    """Resuelve campos + defaults de un cost_term. Lanza ValueError si kind es inválido."""
    var = term['var']
    slc = term.get('slice', slice(None))
    if slc is None:
        slc = slice(None)
    coeff = term.get('coeff', 'sqrt_tau')
    weight = float(term.get('weight', 1.0))
    offset = term.get('offset')
    kind = term['kind']
    if kind not in ('sumsq', 'norm1', 'norm2'):
        raise ValueError(f"cost_term kind desconocido: {kind!r}")
    k_range = term.get('k_range') or ('T+1' if var == 'x' else 'T')
    return ResolvedTerm(var, slc, coeff, weight, offset, kind, k_range)
