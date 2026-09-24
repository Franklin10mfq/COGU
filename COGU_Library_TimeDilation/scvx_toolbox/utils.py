"""Convenience re-exports of the warm-start helpers."""

from .quaternion import (angular_vel, angular_velocity, quat_from_axis_angle,
                         slerp)

__all__ = ['slerp', 'angular_velocity', 'angular_vel', 'quat_from_axis_angle']
