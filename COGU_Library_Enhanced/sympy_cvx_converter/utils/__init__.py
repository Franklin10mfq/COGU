"""
sympy_cvx_converter.utils
Helpers reutilizables para construir warm starts y manipular
representaciones de estado especificas (cuaterniones, etc.).

Acceso recomendado:
    from sympy_cvx_converter.utils import slerp, angular_vel
Acceso explicito (auto-documenta el dominio):
    from sympy_cvx_converter.utils.quaternion import slerp, angular_vel
"""

from .quaternion import slerp, angular_vel

__all__ = ["slerp", "angular_vel"]
