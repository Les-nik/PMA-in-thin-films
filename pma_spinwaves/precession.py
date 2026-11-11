"""Precession orbits and their sense of rotation (Larmor or anti-Larmor).

A spin-wave mode gives complex amplitudes m_x and m_z at every depth z. In
time, the magnetization traces

    m_x(t) = Re[m_x exp(i w t)],   m_z(t) = Re[m_z exp(i w t)],

an ellipse. Its sense of rotation follows from the sign of

    s = Im(conj(m_x) m_z),

which is the signed area swept per period divided by pi. Precession in the
same sense as the free (Larmor) precession around the effective field
corresponds to one sign of s, anti-Larmor precession to the other. Which sign
is which depends on the axis convention, so ``handedness`` returns the raw sign
and the caller fixes the convention from a reference point (k = 0, where the
precession is always Larmor).
"""

from __future__ import annotations

import numpy as np


def orbit(mx: complex, mz: complex, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Points of the precession ellipse over one period."""
    t = np.linspace(0, 1, n)
    phase = np.exp(2j * np.pi * t)
    return np.real(mx * phase), np.real(mz * phase)


def handedness(mx: np.ndarray, mz: np.ndarray) -> np.ndarray:
    """Sign of the rotation sense at every depth (+1, -1, or 0 for linear motion)."""
    return np.sign(np.imag(np.conj(mx) * mz))


def ellipticity(mx: np.ndarray, mz: np.ndarray) -> np.ndarray:
    """|m_z| / |m_x| at every depth."""
    return np.abs(mz) / np.abs(mx)
