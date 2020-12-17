# -*- coding: utf-8 -*-
import numpy as np

from ..fun.constants import rd, g


################################################################################
# Convert pressure levels to geom. height or the other way around
#
# IN: x   pressure or height vector
#     p0  surface pressure    (default: 1013.25 hPa as in International Std.)
#     t0  surface temperature (default: 288.15 K )
#     invert  Option Flag  (true: height to p, false: p to height)
################################################################################
def pressure_to_height(x,
                       p0=1013.25,
                       t0=288.15):
    if isinstance(x, np.ndarray):
        func = np.vectorize(pressure_to_height)
        return func(x, p0=p0, t0=t0)

    dt = 0.0065  # K/m
    a = g / (rd * dt)

    # p to h
    h = (t0 / dt) * (1 - (x / p0) ** (1 / a))
    return np.where(x < 226.32, (11000 + np.log(x / 226.32) * (rd * 216.65) / -g), h)


def height_to_pressure(x,
                       p0=1013.25,
                       t0=288.15):
    if isinstance(x, np.ndarray):
        func = np.vectorize(height_to_pressure)
        return func(x, p0=p0, t0=t0)

    dt = 0.0065  # K/m
    a = g / (rd * dt)
    # h to p
    p = p0 * (1 - (dt * x / t0)) ** a
    return np.where(x > 11000, (226.32 * np.exp((-g * (x - 11000)) / (rd * 216.65))), p)
