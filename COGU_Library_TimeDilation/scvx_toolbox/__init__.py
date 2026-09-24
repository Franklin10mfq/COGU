"""
scvx_toolbox -- SCvx trajectory optimization from a SymPy model.

Define the dynamics symbolically, pick a time mode, and call
`solve_trajectory`.  Four time modes are supported:

    'ZOH'            zero-order hold,  fixed final time
    'FOH'            first-order hold, fixed final time
    'ZOH+dilation'   zero-order hold,  free final time
    'FOH+dilation'   first-order hold, free final time

Under a dilation mode the problem is restated on the normalized interval
s in [0,1] with t = p*s and p = tf a decision variable, so the final time can
be optimized through the convex term alpha_p*tf.
"""

from .modes import (ALL_MODES, FOH, FOH_DILATION, ZOH, ZOH_DILATION,
                    parse_time_mode)
from .pipeline import solve_trajectory
from .problem import build_problem
from .solver import solve_scvx
from .symbolic import SymbolicModel

__version__ = '1.0.0'

__all__ = [
    'solve_trajectory',
    'build_problem',
    'solve_scvx',
    'SymbolicModel',
    'parse_time_mode',
    'ZOH', 'FOH', 'ZOH_DILATION', 'FOH_DILATION', 'ALL_MODES',
]
