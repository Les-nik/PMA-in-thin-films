"""Isofrequency maps from MuMax3 single-frequency fields.

The ``Mf0`` arrays in the data release are complex maps m(x, y) of one
frequency component, taken from a MuMax3 run with a continuous-wave source.
Their 2D Fourier transform shows which wavevectors are present at that
frequency, i.e. the isofrequency contour seen by the simulation.

Sign convention: ``Mf0`` is a forward FFT in time (exp(-i w t)), and NumPy's
``fft2`` uses exp(-i k r) in space. For a wave cos(k.r - w t) the forward time
FFT keeps the exp(-i k.r) part, so the wave would appear at -k.
``isofrequency_map`` returns |F(-k)| so that the bright points sit at the
wavevectors along which the waves actually travel. (The original notebooks
did the same by reversing one axis and plotting with ``origin="upper"``.)
"""

from __future__ import annotations

import numpy as np


def isofrequency_map(field: np.ndarray, dx: float, dy: float):
    """|FFT| of a complex (ny, nx) field on ascending (ky, kx) grids in rad/m.

    Use with ``imshow(amp, origin="lower", extent=extent)``.
    """
    ny, nx = field.shape
    amp = np.fft.fftshift(np.abs(np.fft.fft2(field)))
    amp = amp[::-1, ::-1]  # k -> -k, see module docstring
    kx = np.fft.fftshift(np.fft.fftfreq(nx, dx)) * 2 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(ny, dy)) * 2 * np.pi
    # After reversing, pixel i sits at -k[i]; flip the grids to match.
    kx, ky = -kx[::-1], -ky[::-1]
    extent = [kx[0], kx[-1], ky[0], ky[-1]]
    return kx, ky, amp, extent
