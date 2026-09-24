"""
Time-discretization modes.

The library supports four combinations of a control hold and a time
parameterization:

    ZOH             zero-order hold, fixed final time
    FOH             first-order hold, fixed final time
    ZOH+dilation    zero-order hold, free final time (time dilation)
    FOH+dilation    first-order hold, free final time (time dilation)

Control hold
------------
ZOH: u(t) = u_k for t in [t_k, t_k+1).  The control is a variable at the
     T interval starts, so `u` has shape (nu, T).

FOH: u(t) = lam_minus(t) u_k + lam_plus(t) u_k+1, with
     lam_minus = (t_k+1 - t)/tau and lam_plus = (t - t_k)/tau.  The control
     is a variable at every node, so `u` has shape (nu, T+1).

Time dilation (free final time)
-------------------------------
Without dilation the problem is stated in absolute time t in [0, tf] with a
fixed step tau = tf/T.

With dilation the problem is restated on the normalized interval s in [0, 1]
by means of the dilation parameter p = tf:

    t = p * s,    dx/ds = (dx/dt)(dt/ds) = p * f(x, u) =: f(x, u, p)

so p becomes an optimization variable and the step is the dimensionless
tau = 1/T.  The final time then enters the objective through the convex
linear term alpha_p * tf.
"""

ZOH = 'ZOH'
FOH = 'FOH'
ZOH_DILATION = 'ZOH+dilation'
FOH_DILATION = 'FOH+dilation'

ALL_MODES = (ZOH, FOH, ZOH_DILATION, FOH_DILATION)

# Accepted spellings -> canonical name.  Comparison is done on a normalized
# key (lowercase, separators stripped) so that 'zoh_time_dilation',
# 'ZOH + dilation', 'zoh-td', ... all resolve to the same mode.
_ALIASES = {
    'zoh': ZOH,
    'zeroorderhold': ZOH,
    'foh': FOH,
    'firstorderhold': FOH,

    'zohdilation': ZOH_DILATION,
    'zohtimedilation': ZOH_DILATION,
    'zohdil': ZOH_DILATION,
    'zohtd': ZOH_DILATION,
    'zohfreetime': ZOH_DILATION,
    'zohfreefinaltime': ZOH_DILATION,

    'fohdilation': FOH_DILATION,
    'fohtimedilation': FOH_DILATION,
    'fohdil': FOH_DILATION,
    'fohtd': FOH_DILATION,
    'fohfreetime': FOH_DILATION,
    'fohfreefinaltime': FOH_DILATION,
}


def _normalize(name):
    key = str(name).lower()
    for ch in (' ', '_', '-', '+', '.'):
        key = key.replace(ch, '')
    return key


def parse_time_mode(mode):
    """
    Resolve a user-supplied time mode to (canonical_name, hold, dilation).

    Returns
    -------
    name : str
        One of ZOH, FOH, 'ZOH+dilation', 'FOH+dilation'.
    hold : str
        'ZOH' or 'FOH'.
    dilation : bool
        True if the final time is a free variable.
    """
    if mode is None:
        mode = ZOH
    key = _normalize(mode)
    if key not in _ALIASES:
        raise ValueError(
            f"unknown time_mode {mode!r}. Valid options are: "
            f"{', '.join(ALL_MODES)}"
        )
    name = _ALIASES[key]
    hold = FOH if name in (FOH, FOH_DILATION) else ZOH
    dilation = name in (ZOH_DILATION, FOH_DILATION)
    return name, hold, dilation


def control_node_count(hold, T):
    """Number of columns of the control variable: T for ZOH, T+1 for FOH."""
    return T + 1 if hold == FOH else T
